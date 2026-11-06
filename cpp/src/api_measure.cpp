// Full parity — calibration (dark/flat), time/header, manual measurement and notes/exclusions (Track & Stack in api_stack.cpp).
// Espelha piazzi/extras.py (load_calibrated, measure, track_and_stack) e session.py (measure, set_time, stack, accept).
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include "session_internal.hpp"
#include <exception>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>

using namespace astro;

// ---------------------------------------------------------------- utilities shared with api_stack.cpp
namespace session_util {

// RA/Dec -> pixel through the plate (linear part + Newton), like astro_sky_to_pix / extras.radec_to_pix.
void sky2pix(const Plate &p, double ra, double dec, double &px, double &py) {
    double xi, eta;
    radec_to_std(ra, dec, p.ra0, p.dec0, xi, eta);
    double a = p.cx[1], b = p.cx[2], c = p.cy[1], d = p.cy[2], det = a * d - b * c;
    px = (d * (xi - p.cx[0]) - b * (eta - p.cy[0])) / det;
    py = (-c * (xi - p.cx[0]) + a * (eta - p.cy[0])) / det;
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
}

double plate_scale(const Plate &p) {
    return std::sqrt(std::fabs(p.cx[1] * p.cy[2] - p.cx[2] * p.cy[1]));
}

double nanmedian(std::vector<float> &v) {
    v.erase(std::remove_if(v.begin(), v.end(), [](float q) { return std::isnan(q); }), v.end());
    if (v.empty())
        return std::numeric_limits<double>::quiet_NaN();
    size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n / 2, v.end());
    double hi = v[n / 2];
    if (n % 2)
        return hi;
    double lo = *std::max_element(v.begin(), v.begin() + n / 2);
    return 0.5 * (lo + hi);
}

} // namespace session_util
using namespace session_util;

namespace {

double numd(const Settings &s, const char *k, double def) {
    return s.str(k).empty() ? def : s.num(k);
}

bool is_fits(const std::string &p) {
    std::string l;
    for (char c : p)
        l += char(std::tolower((unsigned char)c));
    auto ends = [&](const char *e) {
        size_t n = std::strlen(e);
        return l.size() >= n && l.compare(l.size() - n, n, e) == 0;
    };
    return ends(".fit") || ends(".fits") || ends(".fts");
}

std::string rtrim(std::string x) {
    while (!x.empty() && (x.back() == ' ' || x.back() == '\r' || x.back() == '\n' || x.back() == '\t'))
        x.pop_back();
    return x;
}

// Text header: FITS = 80-column cards up to END; SBIG = lines of the 2048-byte block. kv receives key -> value.
std::string read_header(const std::string &path, std::map<std::string, std::string> *kv) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return "";
    std::string out;
    if (is_fits(path)) {
        char card[81];
        card[80] = 0;
        while (f.read(card, 80)) {
            std::string c = rtrim(std::string(card, 80));
            out += c + "\n";
            std::string key = rtrim(c.substr(0, std::min<size_t>(8, c.size())));
            if (key == "END")
                break;
            if (kv && c.size() > 10 && c[8] == '=') {
                std::string v = c.substr(10);
                size_t a = v.find_first_not_of(' ');
                if (a != std::string::npos && v[a] == '\'') {
                    size_t b = v.find('\'', a + 1);
                    v = rtrim(v.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1));
                } else {
                    v = rtrim(v.substr(0, v.find('/')));
                    a = v.find_first_not_of(' ');
                    v = a == std::string::npos ? "" : v.substr(a);
                }
                if (!kv->count(key))
                    (*kv)[key] = v;
            }
        }
    } else {
        std::vector<char> raw(2048, 0);
        f.read(raw.data(), 2048);
        size_t n = 0;
        while (n < raw.size() && raw[n])
            ++n;
        std::string cur;
        for (size_t i = 0; i <= n; ++i) {
            char c = i < n ? raw[i] : '\n';
            if (c == '\r' || c == '\n') {
                std::string l = rtrim(cur);
                size_t a = l.find_first_not_of(" \t");
                l = a == std::string::npos ? "" : l.substr(a);
                if (!l.empty()) {
                    out += l + "\n";
                    size_t e = l.find('=');
                    if (kv && e != std::string::npos)
                        (*kv)[rtrim(l.substr(0, e))] = [&] {
                            std::string v = l.substr(e + 1);
                            size_t b = v.find_first_not_of(' ');
                            return b == std::string::npos ? std::string() : v.substr(b);
                        }();
                }
                cur.clear();
            } else
                cur += c;
        }
    }
    return out;
}

bool frame_ok(AstroSession *s, int k) {
    return s && s->ok && k >= 0 && k < int(s->R.frames.size()) && s->R.frames[k].solved;
}

// Measurement on frame k (like detection): fills Detection and Measurement; false on failure.
bool do_measure(AstroSession *s, int k, double x, double y, Detection &d, Measurement &m) {
    if (!frame_ok(s, k))
        return false;
    const Frame &f = s->R.frames[k];
    if (!measure_at(f.im, x, y, detect_params(s->s), d))
        return false;
    m = Measurement{};
    m.frame = k;
    m.x = d.x;
    m.y = d.y;
    f.plate.pix_to_radec(d.x, d.y, m.ra, m.dec);
    m.flux = d.flux;
    m.mag = d.flux > 0 ? f.plate.zp - 2.5 * std::log10(d.flux) : std::numeric_limits<double>::quiet_NaN();
    m.snr = d.snr;
    m.snr_tot = d.snr_tot;
    m.fwhm_as = d.fwhm * plate_scale(f.plate);
    m.jd = f.im.jd_mid;
    m.note = f.stacked ? "K" : "";
    return true;
}

Track *track_at(AstroSession *s, int t) {
    int nm = int(s->R.movers.size());
    if (t >= 0 && t < nm)
        return &s->R.movers[t];
    if (t >= nm && t < nm + int(s->manual.size()))
        return &s->manual[t - nm];
    return nullptr;
}

void update_motion(AstroSession *s, Track &T) {
    if (T.meas.size() < 2) {
        T.rate = T.pa = 0;
        return;
    }
    const auto &a = T.meas.front();
    const auto &b = T.meas.back();
    double dt = (b.jd - a.jd) * 1440.0;
    if (dt == 0)
        return;
    double xi, eta;
    radec_to_std(b.ra, b.dec, a.ra, a.dec, xi, eta);
    T.rate = std::hypot(xi, eta) / dt;
    T.pa = std::fmod(std::atan2(xi, eta) * 180.0 / PI + 360.0, 360.0);
    (void)s;
}

} // namespace

// ---------------------------------------------------------------- loading with time + calibration (extras.load_calibrated)
Image session_load_image(AstroSession *s, size_t idx, bool calibrate) {
    const Settings &S = s->s;
    const std::string &path = s->images.at(idx);
    int tmode = int(S.num("time_reference"));
    Image im = load_image(path, S.num("flip_horizontal") != 0, S.num("flip_vertical") != 0, tmode == 0);
    bool fits = is_fits(path);
    int eu = int(numd(S, "ExpUnit", 1));
    if (fits && eu != 1 && (eu == 0 || eu == 2)) { // "Exposure from FITS": min / s / ms
        double e_old = im.exptime;
        im.exptime = e_old * (eu == 0 ? 60.0 : 0.001);
        if (tmode == 0)
            im.jd_mid += (im.exptime - e_old) / 2 / 86400.0;
    }
    if (tmode == 2)
        im.jd_mid -= im.exptime / 2 / 86400.0; // timestamp = end
    im.jd_mid += -numd(S, "TimeOffsetH", 0) / 24.0 + numd(S, "TimeDelay", 0) / 86400.0;
    if (fits && numd(S, "AutoRotatePier", 0) != 0) {
        std::map<std::string, std::string> kv;
        read_header(path, &kv);
        std::string ps = kv["PIERSIDE"];
        size_t a = ps.find_first_not_of(' ');
        if (a != std::string::npos && std::toupper((unsigned char)ps[a]) == 'W')
            std::reverse(im.pix.begin(), im.pix.end());
    }
    im.jd_mid += numd(S, "time_offset_s", 0) / 86400.0;
    auto ov = s->jd_override.find(path);
    if (ov != s->jd_override.end())
        im.jd_mid = ov->second;
    if (calibrate && (!s->dark_path.empty() || !s->flat_path.empty())) {
        std::vector<float> v(im.pix.begin(), im.pix.end());
        if (!s->dark_path.empty()) {
            Image d = load_image(s->dark_path, S.num("flip_horizontal") != 0, S.num("flip_vertical") != 0);
            if (d.w == im.w && d.h == im.h)
                for (size_t i = 0; i < v.size(); ++i)
                    v[i] = float(v[i] - float(d.pix[i]));
        }
        if (!s->flat_path.empty()) {
            Image fl = load_image(s->flat_path, S.num("flip_horizontal") != 0, S.num("flip_vertical") != 0);
            if (fl.w == im.w && fl.h == im.h) {
                std::vector<float> t(fl.pix.begin(), fl.pix.end());
                double med = nanmedian(t);
                if (med > 0)
                    for (size_t i = 0; i < v.size(); ++i)
                        if (fl.pix[i] > 0)
                            v[i] = float(v[i] / float(fl.pix[i] / med));
            }
        }
        for (size_t i = 0; i < v.size(); ++i) {
            double r = std::nearbyint(double(v[i]));
            im.pix[i] = uint16_t(r > 0 ? std::min(r, 65535.0) : 0.0);
        }
    }
    return im;
}

std::vector<Track> session_report_tracks(AstroSession *s) {
    std::vector<Track> out;
    int nm = int(s->R.movers.size());
    for (int t = 0; t < nm + int(s->manual.size()); ++t) {
        Track T = *track_at(s, t);
        T.meas.clear();
        for (auto m : track_at(s, t)->meas) {
            auto key = std::make_pair(t, m.frame);
            auto ex = s->excluded.find(key);
            if (ex != s->excluded.end() && ex->second)
                continue;
            auto no = s->notes.find(key);
            if (no != s->notes.end())
                m.note = no->second;
            T.meas.push_back(m);
        }
        if (!T.meas.empty())
            out.push_back(T);
    }
    return out;
}

extern "C" {

int astro_set_dark(AstroSession *s, const char *path) {
    try {
        if (!s)
            return 0;
        if (!path || !*path) {
            s->dark_path.clear();
            return 1;
        }
        try {
            Image d = load_image(path, false, false);
            if (d.w <= 0)
                return 0;
        } catch (...) {
            return 0;
        }
        s->dark_path = path;
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
int astro_set_flat(AstroSession *s, const char *path) {
    try {
        if (!s)
            return 0;
        if (!path || !*path) {
            s->flat_path.clear();
            return 1;
        }
        try {
            Image d = load_image(path, false, false);
            if (d.w <= 0)
                return 0;
        } catch (...) {
            return 0;
        }
        s->flat_path = path;
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

int astro_image_count(AstroSession *s) {
    try {
        if (!s)
            return 0;
        return int(s->images.size());
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
const char *astro_image_header(AstroSession *s, int i) {
    try {
        if (!s)
            return "";
        s->hdr = (i >= 0 && i < int(s->images.size())) ? read_header(s->images[i], nullptr) : std::string();
        return s->hdr.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}
double astro_image_jd(AstroSession *s, int i) {
    try {
        if (!s)
            return 0;
        if (i < 0 || i >= int(s->images.size()))
            return 0;
        try {
            return session_load_image(s, size_t(i), false).jd_mid;
        } catch (...) {
            return 0;
        }
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}
void astro_set_image_jd(AstroSession *s, int i, double jd) {
    try {
        if (!s)
            return;
        if (i < 0 || i >= int(s->images.size()))
            return;
        if (jd == 0)
            s->jd_override.erase(s->images[i]);
        else
            s->jd_override[s->images[i]] = jd;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

int astro_measure(AstroSession *s, int k, double x, double y, double *out) {
    try {
        if (!s)
            return 0;
        Detection d;
        Measurement m;
        if (!do_measure(s, k, x, y, d, m))
            return 0;
        const Plate &p = s->R.frames[k].plate;
        const Frame &f = s->R.frames[k];
        double dra = p.rms_ra, dde = p.rms_dec;
        for (size_t i = 0; i < p.ref_idx.size(); ++i) { // is it a reference star? its O−C residual
            const Detection &r = f.det[p.ref_idx[i]];
            int c = i < p.cat_idx.size() ? p.cat_idx[i] : -1;
            if (std::hypot(r.x - d.x, r.y - d.y) < 0.5 && c >= 0 && c < int(s->cat.size())) {
                double xi, eta, cxi, ceta;
                auto t = p.terms(d.x, d.y);
                xi = eta = 0;
                for (size_t q = 0; q < t.size(); ++q) {
                    xi += t[q] * p.cx[q];
                    eta += t[q] * p.cy[q];
                }
                radec_to_std(s->cat.ra[c], s->cat.dec[c], p.ra0, p.dec0, cxi, ceta);
                dra = xi - cxi;
                dde = eta - ceta;
                break;
            }
        }
        double v[12] = {m.ra, m.dec, m.mag, d.snr, m.fwhm_as, d.x, d.y, d.fit, d.flux, d.snr_tot, dra, dde};
        std::copy(v, v + 12, out);
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

int astro_add_manual(AstroSession *s, int track, int frame, double x, double y) {
    try {
        if (!s)
            return -1;
        Detection d;
        Measurement m;
        if (!do_measure(s, frame, x, y, d, m))
            return -1;
        int nm = int(s->R.movers.size());
        Track *T = nullptr;
        int idx = track;
        if (track < 0) {
            Track n;
            n.desig = s->object;
            s->manual.push_back(n);
            idx = nm + int(s->manual.size()) - 1;
        }
        T = track_at(s, idx);
        if (!T)
            return -1;
        bool done = false;
        for (auto &q : T->meas)
            if (q.frame == frame) {
                q = m;
                done = true;
            }
        if (!done) {
            T->meas.push_back(m);
            std::stable_sort(T->meas.begin(), T->meas.end(), [](const Measurement &a, const Measurement &b) { return a.jd < b.jd; });
        }
        if (idx >= nm)
            update_motion(s, *T);
        session_rebuild_json(s);
        return idx;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return -1;
    } catch (...) {
        return -1;
    }
}
void astro_set_note(AstroSession *s, int track, int frame, const char *note) {
    try {
        if (!s)
            return;
        if (!track_at(s, track))
            return;
        s->notes[{track, frame}] = note ? std::string(note).substr(0, 1) : std::string();
        session_rebuild_json(s);
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}
void astro_set_excluded(AstroSession *s, int track, int frame, int excluded) {
    try {
        if (!s)
            return;
        if (!track_at(s, track))
            return;
        s->excluded[{track, frame}] = excluded != 0;
        session_rebuild_json(s);
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

} // extern "C"
