// C API: alternative catalogs (URL), display with levels/histogram, and 16-bit FITS writing with WCS.
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

using namespace astro;
#include "session_internal.hpp"

static double parse_ang(const std::string &s, bool hours) {
    double a = 0, b = 0, c = 0;
    int n = std::sscanf(s.c_str(), "%lf %lf %lf", &a, &b, &c);
    if (n <= 1)
        return a;
    double sg = s.find('-') != std::string::npos ? -1 : 1;
    return sg * (std::fabs(a) + b / 60 + c / 3600) * (hours ? 15 : 1);
}

static bool frame_ok(AstroSession *s, int k) {
    return s && s->ok && k >= 0 && k < int(s->R.frames.size());
}

static std::vector<uint16_t> sorted_pix(const Image &im) {
    std::vector<uint16_t> t(im.pix);
    std::sort(t.begin(), t.end());
    return t;
}
static double pct(const std::vector<uint16_t> &t, double p) {
    if (t.empty())
        return 0;
    double q = std::clamp(p, 0.0, 100.0) / 100.0 * (t.size() - 1);
    size_t i = size_t(q);
    double f = q - i;
    return i + 1 < t.size() ? t[i] * (1 - f) + t[i + 1] * f : t[i];
}

// TAN (+SIP if order > 1) exactly equivalent to the plate: xi, eta (″) = poly in 0-based (x, y), tangent at ra0/dec0.
struct Wcs {
    double crpix1, crpix2, cd[2][2];
    int sip = 0;
    std::vector<std::vector<double>> A, B;
};
static Wcs plate_wcs(const Plate &p) {
    // exponents of each plate term
    std::vector<std::pair<int, int>> ex{{0, 0}, {1, 0}, {0, 1}};
    if (p.order >= 2)
        for (auto e : {std::pair{2, 0}, {1, 1}, {0, 2}})
            ex.push_back(e);
    if (p.order >= 3)
        for (auto e : {std::pair{3, 0}, {2, 1}, {1, 2}, {0, 3}})
            ex.push_back(e);
    if (p.order >= 4)
        for (auto e : {std::pair{4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4}})
            ex.push_back(e);
    auto ev = [&](const std::vector<double> &c, double x, double y) {
        double s = 0;
        for (size_t i = 0; i < c.size(); ++i)
            s += c[i] * std::pow(x, ex[i].first) * std::pow(y, ex[i].second);
        return s;
    };
    // point where xi = eta = 0 (Newton on the linear part)
    double a = p.cx[1], b = p.cx[2], c = p.cy[1], d = p.cy[2], det = a * d - b * c;
    double x0 = (-d * p.cx[0] + b * p.cy[0]) / det, y0 = (c * p.cx[0] - a * p.cy[0]) / det;
    for (int it = 0; it < 20; ++it) {
        double ex_ = -ev(p.cx, x0, y0), ey = -ev(p.cy, x0, y0);
        x0 += (d * ex_ - b * ey) / det;
        y0 += (-c * ex_ + a * ey) / det;
    }
    // re-expand the polynomial in u = x - x0, v = y - y0
    int N = std::max(p.order, 1);
    auto binom = [](int n, int k) {
        double r = 1;
        for (int i = 1; i <= k; ++i)
            r = r * (n - k + i) / i;
        return r;
    };
    auto shift = [&](const std::vector<double> &cf) {
        std::vector<std::vector<double>> q(N + 1, std::vector<double>(N + 1, 0.0));
        for (size_t t = 0; t < cf.size(); ++t) {
            int i = ex[t].first, j = ex[t].second;
            for (int u = 0; u <= i; ++u)
                for (int v = 0; v <= j; ++v)
                    q[u][v] += cf[t] * binom(i, u) * binom(j, v) * std::pow(x0, i - u) * std::pow(y0, j - v);
        }
        return q;
    };
    auto qx = shift(p.cx), qy = shift(p.cy);
    Wcs w;
    w.crpix1 = x0 + 1;
    w.crpix2 = y0 + 1;
    w.cd[0][0] = qx[1][0] / 3600;
    w.cd[0][1] = qx[0][1] / 3600;
    w.cd[1][0] = qy[1][0] / 3600;
    w.cd[1][1] = qy[0][1] / 3600;
    if (p.order >= 2) {
        w.sip = p.order;
        w.A.assign(N + 1, std::vector<double>(N + 1, 0.0));
        w.B = w.A;
        double ia = qy[0][1], ib = -qx[0][1], ic = -qy[1][0], id = qx[1][0], dt = qx[1][0] * qy[0][1] - qx[0][1] * qy[1][0];
        for (int u = 0; u <= N; ++u)
            for (int v = 0; u + v <= N; ++v)
                if (u + v >= 2) {
                    w.A[u][v] = (ia * qx[u][v] + ib * qy[u][v]) / dt;
                    w.B[u][v] = (ic * qx[u][v] + id * qy[u][v]) / dt;
                }
    }
    return w;
}

static std::string iso_ms(double jd) {
    double z = jd + 0.5;
    long Z = long(std::floor(z));
    double F = z - Z;
    long A = Z;
    if (Z >= 2299161) {
        long al = long((Z - 1867216.25) / 36524.25);
        A = Z + 1 + al - al / 4;
    }
    long B = A + 1524, C = long((B - 122.1) / 365.25), D = long(365.25 * C), E = long((B - D) / 30.6001);
    long day = B - D - long(30.6001 * E), mon = E < 14 ? E - 1 : E - 13, yr = mon > 2 ? C - 4716 : C - 4715;
    long ms = std::lround(F * 86400000.0);
    if (ms >= 86400000)
        ms = 86399999;
    char b[40];
    std::snprintf(b, sizeof b, "%04ld-%02ld-%02ldT%02ld:%02ld:%02ld.%03ld", yr, mon, day, ms / 3600000, ms / 60000 % 60, ms / 1000 % 60,
                  ms % 1000);
    return b;
}

extern "C" {

const char *astro_catalog_url_for(AstroSession *s, int catalog, const char *ra, const char *dec) {
    try {
        if (!s)
            return "";
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
        s->buf =
            catalog_box_url(catalog, parse_ang(ra ? ra : "", true), parse_ang(dec ? dec : "", false), box, s->s.num("catalog_faint_mag"));
        return s->buf.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "";
    } catch (...) {
        return "";
    }
}

void astro_render_levels(AstroSession *s, int frame, double lo_pct, double hi_pct, int mode, int invert, unsigned char *rgba) {
    try {
        if (!frame_ok(s, frame) || !rgba)
            return;
        const Image &im = s->R.frames[frame].im;
        auto t = sorted_pix(im);
        double lo = pct(t, lo_pct), hi = pct(t, hi_pct);
        if (hi <= lo)
            hi = lo + 1;
        for (size_t i = 0; i < im.pix.size(); ++i) {
            double q = std::clamp((im.pix[i] - lo) / (hi - lo), 0.0, 1.0);
            if (mode == 1)
                q = std::sqrt(q);
            else if (mode == 2)
                q = std::log10(1 + 1000 * q) / 3;
            if (invert)
                q = 1 - q;
            unsigned char g = (unsigned char)std::lround(255 * q);
            rgba[i * 4] = rgba[i * 4 + 1] = rgba[i * 4 + 2] = g;
            rgba[i * 4 + 3] = 255;
        }
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return;
    } catch (...) {
        return;
    }
}

int astro_histogram(AstroSession *s, int frame, int bins, double *counts, double *minmax2) {
    try {
        if (!frame_ok(s, frame) || bins <= 0 || !counts)
            return 0;
        const Image &im = s->R.frames[frame].im;
        if (im.pix.empty())
            return 0;
        auto t = sorted_pix(im);
        double lo = pct(t, 0.1), hi = pct(t, 99.9);
        if (hi <= lo)
            hi = lo + 1; // robust min/max
        std::fill(counts, counts + bins, 0.0);
        for (uint16_t v : im.pix) {
            if (v < lo || v > hi)
                continue;
            int b = int((v - lo) / (hi - lo) * bins);
            counts[std::min(b, bins - 1)] += 1;
        }
        if (minmax2) {
            minmax2[0] = lo;
            minmax2[1] = hi;
        }
        return 1;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

int astro_save_fits(AstroSession *s, int frame, const char *path) {
    try {
        if (!frame_ok(s, frame) || !path)
            return 0;
        const Frame &f = s->R.frames[frame];
        const Image &im = f.im;
        std::vector<std::string> cards;
        auto card = [&](const std::string &k, const std::string &v, const std::string &c = "") {
            char b[81];
            std::snprintf(b, sizeof b, "%-8s= %20s%s%s", k.c_str(), v.c_str(), c.empty() ? "" : " / ", c.c_str());
            std::string x = b;
            x.resize(80, ' ');
            cards.push_back(x);
        };
        auto str = [](const std::string &v) {
            std::string x = "'" + v + "'";
            while (x.size() < 10)
                x.insert(x.size() - 1, " ");
            return x;
        };
        auto num = [](double v) {
            char b[32];
            std::snprintf(b, sizeof b, "%.15G", v);
            std::string x = b;
            if (x.find_first_of(".E") == std::string::npos)
                x += ".0";
            return x;
        };
        auto strc = [&](const std::string &k, const std::string &v, const std::string &c = "") {
            std::string x = k;
            x.resize(8, ' ');
            x += "= " + str(v);
            if (!c.empty())
                x += " / " + c;
            x.resize(80, ' ');
            cards.push_back(x);
        };
        card("SIMPLE", "T");
        card("BITPIX", "16");
        card("NAXIS", "2");
        card("NAXIS1", std::to_string(im.w));
        card("NAXIS2", std::to_string(im.h));
        card("BZERO", "32768");
        card("BSCALE", "1");
        strc("DATE-OBS", iso_ms(im.jd_mid), "UT mid-exposure");
        card("EXPTIME", num(im.exptime), "s");
        std::string base = im.path.substr(im.path.find_last_of("/\\") + 1);
        strc("ORIGFILE", base.substr(0, 60));
        if (f.solved && f.plate.cx.size() >= 3) {
            const Plate &p = f.plate;
            Wcs w = plate_wcs(p);
            strc("CTYPE1", w.sip ? "RA---TAN-SIP" : "RA---TAN");
            strc("CTYPE2", w.sip ? "DEC--TAN-SIP" : "DEC--TAN");
            card("CRVAL1", num(p.ra0));
            card("CRVAL2", num(p.dec0));
            card("CRPIX1", num(w.crpix1));
            card("CRPIX2", num(w.crpix2));
            card("CD1_1", num(w.cd[0][0]));
            card("CD1_2", num(w.cd[0][1]));
            card("CD2_1", num(w.cd[1][0]));
            card("CD2_2", num(w.cd[1][1]));
            strc("CUNIT1", "deg");
            strc("CUNIT2", "deg");
            card("EQUINOX", "2000.0");
            strc("RADESYS", "ICRS");
            if (w.sip) {
                card("A_ORDER", std::to_string(w.sip));
                card("B_ORDER", std::to_string(w.sip));
                for (int u = 0; u <= w.sip; ++u)
                    for (int v = 0; u + v <= w.sip; ++v)
                        if (u + v >= 2) {
                            card("A_" + std::to_string(u) + "_" + std::to_string(v), num(w.A[u][v]));
                            card("B_" + std::to_string(u) + "_" + std::to_string(v), num(w.B[u][v]));
                        }
            }
            double rc, dc;
            p.pix_to_radec(im.w / 2.0, im.h / 2.0, rc, dc);
            std::string r = fmt_ra(rc), d = fmt_dec(dc);
            strc("OBJCTRA", r);
            strc("OBJCTDEC", d);
        }
        {
            std::string x = "HISTORY Piazzi";
            x.resize(80, ' ');
            cards.push_back(x);
        }
        {
            std::string x = "END";
            x.resize(80, ' ');
            cards.push_back(x);
        }
        std::string hdr;
        for (auto &c : cards)
            hdr += c;
        hdr.resize((hdr.size() + 2879) / 2880 * 2880, ' ');
        std::ofstream o(path, std::ios::binary);
        if (!o)
            return 0;
        o.write(hdr.data(), hdr.size());
        std::vector<unsigned char> data(im.pix.size() * 2);
        for (size_t i = 0; i < im.pix.size(); ++i) {
            uint16_t v = uint16_t(im.pix[i] ^ 0x8000);
            data[2 * i] = v >> 8;
            data[2 * i + 1] = v & 0xff;
        }
        data.resize((data.size() + 2879) / 2880 * 2880, 0);
        o.write((const char *)data.data(), data.size());
        return o.good() ? 1 : 0;
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

} // extern "C"
