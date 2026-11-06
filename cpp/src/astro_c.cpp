#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include "piazzi/version.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <sstream>

using namespace astro;

#include "session_internal.hpp"

static double parse_angle(const std::string &s, bool hours) {
    double a = 0, b = 0, c = 0;
    int n = std::sscanf(s.c_str(), "%lf %lf %lf", &a, &b, &c);
    if (n <= 1)
        return a;
    double sg = s.find('-') != std::string::npos ? -1 : 1;
    return sg * (std::fabs(a) + b / 60 + c / 3600) * (hours ? 15 : 1);
}

extern "C" {
AstroSession *astro_new(void) {
    try {
        return new AstroSession();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return nullptr;
    } catch (...) {
        return nullptr;
    }
}
void astro_free(AstroSession *s) {
    try {
        delete s;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}
void astro_load_cfg(AstroSession *s, const char *p) {
    try {
        if (!s)
            return;
        s->s.load(p);
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}
void astro_set(AstroSession *s, const char *k, const char *v) {
    try {
        if (!s)
            return;
        s->s.set(k, v);
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}
const char *astro_get(AstroSession *s, const char *k) {
    try {
        if (!s)
            return "";
        s->getv = s->s.str(k);
        return s->getv.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}
int astro_add_image(AstroSession *s, const char *p) {
    try {
        if (!s)
            return 0;
        s->images.push_back(p);
        return int(s->images.size());
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
void astro_clear_images(AstroSession *s) {
    try {
        if (!s)
            return;
        s->images.clear();
        s->ok = false;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

const char *astro_catalog_url(AstroSession *s, const char *ra, const char *dec) {
    try {
        if (!s)
            return "";
        // box = image width + 2 × Pointing (pointing margin); without an image, 30′
        double box = 30.0;
        if (!s->images.empty()) {
            try {
                Image im = load_image(s->images[0], s->s.num("flip_horizontal") != 0, s->s.num("flip_vertical") != 0);
                double sx = s->s.num("pixel_width_mm") / s->s.num("focal_length_mm") * 206264.8;
                if (im.w > 0 && sx > 0)
                    box = im.w * sx / 60 + 2 * s->s.num("pointing_tolerance");
            } catch (...) {
            }
        }
        s->buf = gaia_box_url(parse_angle(ra, true), parse_angle(dec, false), box, s->s.num("catalog_faint_mag"));
        return s->buf.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}

static std::string esc(const std::string &x) {
    return session_util::json_esc(x);
}

} // extern "C"

const char *astro_version(void) {
    return PIAZZI_VERSION;
}

int astro_reduce(AstroSession *s, const char *catalog, const char *ra, const char *dec, const char *object) {
    try {
        if (!s)
            return -1;
        s->ok = false;
        s->manual.clear();
        s->notes.clear();
        s->excluded.clear();
        if (s->images.empty()) {
            s->json = "{\"error\":\"No images — add at least two exposures of the same field\",\"code\":\"no_images\"}";
            return -1;
        }
        std::vector<Image> ims;
        for (size_t k = 0; k < s->images.size(); ++k) {
            try {
                ims.push_back(session_load_image(s, k, true));
            } catch (...) {
                ims.push_back(Image{});
            }
            if (ims.back().w == 0) {
                s->json = "{\"error\":\"Image unreadable: " + esc(s->images[k]) +
                          " — use FITS or SBIG files\",\"code\":\"image_unreadable\",\"detail\":\"" + esc(s->images[k]) + "\"}";
                return -1;
            }
        }
        double t = ims[0].jd_mid;
        for (auto &im : ims)
            t = std::min(t, im.jd_mid);
        Catalog cat = load_vizier_tsv(catalog, 2000.0 + (t - 2451545.0) / 365.25, s->s.num("catalog_faint_mag"),
                                      s->s.str("plate_ruwe").empty() ? 1.4 : s->s.num("plate_ruwe"));
        if (cat.size() < 3) {
            s->json = "{\"error\":\"Catalog empty or unreadable — check the internet connection or load a Gaia DR3 TSV\","
                      "\"code\":\"catalog_empty\"}";
            return -1;
        }
        s->R = reduce(ims, cat, parse_angle(ra, true), parse_angle(dec, false), s->s, true);
        s->cat = Catalog{}; // as inside reduce: only stars fainter than catalog_bright_mag (Plate::cat_idx indices)
        for (size_t i = 0; i < cat.size(); ++i)
            if (cat.mag[i] >= s->s.num("catalog_bright_mag"))
                s->cat.push_from(cat, i);
        // with no object given the designation is empty: ADES carries only trkSub and MPC uses each track's
        // temporary designation (the UI generates the display label)
        s->object = object ? object : "";
        for (auto &T : s->R.movers)
            T.desig = s->object;
        s->ok = true;
        session_rebuild_json(s);
        return int(s->R.movers.size());
    } catch (const std::exception &e) { // C++ exceptions must not cross the extern "C" boundary
        if (s) {
            s->ok = false;
            s->json = "{\"error\":\"Internal error: " + esc(e.what()) + "\",\"code\":\"internal\",\"detail\":\"" + esc(e.what()) + "\"}";
        }
        return -1;
    } catch (...) {
        return -1;
    }
}

namespace {
// JSON number: NaN/Inf do not exist in JSON — emitted as null
struct Num {
    double v;
};
std::ostream &operator<<(std::ostream &o, Num n) {
    if (std::isfinite(n.v))
        o << n.v;
    else
        o << "null";
    return o;
}
} // namespace

void session_rebuild_json(AstroSession *s) {
    std::ostringstream j;
    j.precision(10);
    j << "{\"catalog\":" << s->cat.size() << ",\"frames\":[";
    for (size_t k = 0; k < s->R.frames.size(); ++k) {
        const Frame &f = s->R.frames[k];
        const Plate &p = f.plate;
        if (k)
            j << ",";
        j << "{\"name\":\"" << esc(f.im.path.substr(f.im.path.find_last_of("/\\") + 1)) << "\",\"time\":\"" << jd_to_isot(f.im.jd_mid)
          << "\",\"jd\":" << Num{f.im.jd_mid} << ",\"ndet\":" << f.det.size() << ",\"nref\":" << p.nref << ",\"nmatch\":" << p.nmatch
          << ",\"dra\":" << Num{p.dra_mean} << ",\"dde\":" << Num{p.dde_mean} << ",\"zp\":" << Num{p.zp} << ",\"nphot\":" << p.nphot
          << ",\"seeing\":" << Num{f.seeing} << ",\"w\":" << f.im.w << ",\"h\":" << f.im.h
          << ",\"stacked\":" << (f.stacked ? "true" : "false") << ",\"dets\":[";
        for (size_t i = 0; i < f.det.size(); ++i)
            j << (i ? "," : "") << "[" << Num{f.det[i].x} << "," << Num{f.det[i].y} << "]";
        j << "],\"refs\":[";
        for (size_t i = 0; i < p.ref_idx.size(); ++i)
            j << (i ? "," : "") << "[" << Num{f.det[p.ref_idx[i]].x} << "," << Num{f.det[p.ref_idx[i]].y} << "," << int(p.ref_ok[i]) << "]";
        j << "]}";
    }
    j << "],\"movers\":[";
    size_t nm = s->R.movers.size();
    // tracks as they enter the report (no excluded measurements, no empty tracks): trkSub and σ match the ADES
    std::vector<Track> rep = session_report_tracks(s);
    size_t ri = 0; // position of the current track in `rep`
    for (size_t m = 0; m < nm + s->manual.size(); ++m) {
        const Track &T = m < nm ? s->R.movers[m] : s->manual[m - nm];
        bool any = false;
        for (const auto &q : T.meas) {
            auto key = std::make_pair(int(m), q.frame);
            if (!(s->excluded.count(key) && s->excluded[key]))
                any = true;
        }
        const Track *RT = any && ri < rep.size() ? &rep[ri] : nullptr;
        std::string trk = RT ? temp_id(*RT, ri) : "";
        if (RT)
            ++ri;
        if (m)
            j << ",";
        j << "{\"desig\":\"" << esc(T.desig) << "\",\"trksub\":\"" << esc(trk) << "\",\"manual\":" << (m < nm ? "false" : "true")
          << ",\"rate\":" << Num{T.rate} << ",\"pa\":" << Num{T.pa} << ",\"meas\":[";
        for (size_t i = 0; i < T.meas.size(); ++i) {
            const auto &q = T.meas[i];
            auto key = std::make_pair(int(m), q.frame);
            std::string note = s->notes.count(key) ? s->notes[key] : q.note;
            bool ex = s->excluded.count(key) && s->excluded[key];
            double era = NAN, ede = NAN, tl = NAN;
            if (q.frame >= 0 && q.frame < int(s->R.frames.size())) {
                const Plate &p = s->R.frames[q.frame].plate;
                tl = q.trail_len_px * 0.5 * (p.scale_x() + p.scale_y());
                measurement_sigma(s->R, RT ? *RT : T, q, s->s, era, ede);
            }
            j << (i ? "," : "") << "{\"frame\":" << q.frame << ",\"x\":" << Num{q.x} << ",\"y\":" << Num{q.y} << ",\"ra\":\""
              << fmt_ra(q.ra) << "\",\"dec\":\"" << fmt_dec(q.dec) << "\",\"mag\":" << (std::isfinite(q.mag) ? q.mag : 0.0)
              << ",\"rad\":" << Num{q.ra} << ",\"decd\":" << Num{q.dec} << ",\"jd\":" << Num{q.jd} << ",\"time\":\"" << jd_to_isot(q.jd)
              << "\",\"note\":\"" << esc(note) << "\",\"excluded\":" << (ex ? "true" : "false") << ",\"trksub\":\"" << esc(trk)
              << "\",\"sig_ra\":" << Num{era} << ",\"sig_dec\":" << Num{ede} << ",\"trail_len\":" << Num{tl} << ",\"remarks\":\""
              << esc(q.remark) << "\"}";
        }
        j << "]}";
    }
    j << "]}";
    s->json = j.str();
}

extern "C" {
const char *astro_result_json(AstroSession *s) {
    try {
        if (!s)
            return "";
        return s->json.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}
const char *astro_report(AstroSession *s, int which) {
    try {
        if (!s || !s->ok)
            return "";
        auto tr = session_report_tracks(s);
        s->report = which == 0 ? ades_report(s->R, tr, s->s) : which == 1 ? phot_report(s->R, tr, s->s) : mpc_report(s->R, tr, s->s);
        return s->report.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}
int astro_frame_count(AstroSession *s) {
    try {
        if (!s)
            return 0;
        return s->ok ? int(s->R.frames.size()) : 0;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
int astro_frame_size(AstroSession *s, int k, int *w, int *h) {
    try {
        if (!s || !w || !h || !s->ok || k < 0 || k >= int(s->R.frames.size()))
            return 0;
        *w = s->R.frames[k].im.w;
        *h = s->R.frames[k].im.h;
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
void astro_render(AstroSession *s, int k, unsigned char *rgba) {
    try {
        if (!s || !rgba || !s->ok || k < 0 || k >= int(s->R.frames.size()))
            return;
        const Image &im = s->R.frames[k].im;
        if (im.pix.empty())
            return;
        std::vector<int> v(im.pix.begin(), im.pix.end()), t = v;
        std::nth_element(t.begin(), t.begin() + t.size() / 2, t.end());
        double med = t[t.size() / 2];
        std::nth_element(t.begin(), t.begin() + size_t(t.size() * 0.995), t.end());
        double hi = t[size_t(t.size() * 0.995)];
        double lo = med - (hi - med) * 0.08;
        if (hi <= lo)
            hi = lo + 1;
        for (size_t i = 0; i < v.size(); ++i) {
            double q = std::clamp((v[i] - lo) / (hi - lo), 0.0, 1.0);
            unsigned char g = (unsigned char)(255 * std::sqrt(q));
            rgba[i * 4] = rgba[i * 4 + 1] = rgba[i * 4 + 2] = g;
            rgba[i * 4 + 3] = 255;
        }
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

int astro_sky_to_pix(AstroSession *s, int k, double ra, double dec, double *x, double *y) {
    try {
        if (!s || !s->ok || k < 0 || k >= int(s->R.frames.size()) || !s->R.frames[k].solved)
            return 0;
        const Plate &p = s->R.frames[k].plate;
        double xi, eta;
        radec_to_std(ra, dec, p.ra0, p.dec0, xi, eta);
        // plate inverse: linear part + Newton (valid for any order)
        double a = p.cx[1], b = p.cx[2], c = p.cy[1], d = p.cy[2], det = a * d - b * c;
        if (det == 0 || !std::isfinite(det)) // degenerate plate: no inverse
            return 0;
        double px = (d * (xi - p.cx[0]) - b * (eta - p.cy[0])) / det, py = (-c * (xi - p.cx[0]) + a * (eta - p.cy[0])) / det;
        for (int it = 0; it < 10; ++it) {
            auto t = p.terms(px, py);
            double fx = 0, fy = 0;
            for (size_t i = 0; i < t.size(); ++i) {
                fx += t[i] * p.cx[i];
                fy += t[i] * p.cy[i];
            }
            double ex = xi - fx, ey = eta - fy;
            px += (d * ex - b * ey) / det;
            py += (-c * ex + a * ey) / det;
            if (std::fabs(ex) + std::fabs(ey) < 1e-6)
                break;
        }
        *x = px;
        *y = py;
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
int astro_frame_center(AstroSession *s, int k, double *ra, double *dec) {
    try {
        if (!s || !s->ok || k < 0 || k >= int(s->R.frames.size()) || !s->R.frames[k].solved)
            return 0;
        const Frame &f = s->R.frames[k];
        f.plate.pix_to_radec(f.im.w / 2.0, f.im.h / 2.0, *ra, *dec);
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

void astro_set_mover_desig(AstroSession *s, int i, const char *desig) {
    try {
        if (!s)
            return;
        int nm = int(s->R.movers.size());
        if (s->ok && i >= 0 && i < nm)
            s->R.movers[i].desig = desig;
        else if (s->ok && i >= nm && i < nm + int(s->manual.size()))
            s->manual[i - nm].desig = desig;
        if (s->ok)
            session_rebuild_json(s);
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

int astro_pix_to_sky(AstroSession *s, int k, double x, double y, double *ra, double *dec) {
    try {
        if (!s || !s->ok || k < 0 || k >= int(s->R.frames.size()) || !s->R.frames[k].solved)
            return 0;
        s->R.frames[k].plate.pix_to_radec(x, y, *ra, *dec);
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
}
