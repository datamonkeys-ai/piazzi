// Offline identification of solar-system bodies (MPCORB.DAT / COMET.DAT), like piazzi/ident.py:
// unperturbed two-body orbit, light time, topocentric position from [ObservingSite].
// Earth comes from a Chebyshev fit to erfa.epv00 (the same model as astropy 'builtin' used by Python),
// valid from 1950 to 2100 (error < 1e-7 AU).
#include "session_internal.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "earth_cheb.inc"

namespace {
using Orbit = AstroSession::Orbit;
constexpr double PI = 3.14159265358979323846, D2R = PI / 180, R2D = 180 / PI;
constexpr double KG = 0.01720209895;      // constante gaussiana
constexpr double EPS0 = 23.4392911 * D2R; // J2000 obliquity (same as Python)
constexpr double LT_DAY = 0.0057755183;   // days per AU
constexpr double MU_MOON = 1 / (1 + 81.30056);
constexpr double J2000 = 2451545.0;

double cal_jd(int y, int m, double d) { // gregoriano -> JD
    if (m <= 2) {
        y -= 1;
        m += 12;
    }
    int A = y / 100, B = 2 - A + A / 4;
    return std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1)) + d + B - 1524.5;
}

// TT - UTC (s): leap seconds since 1972; approximation before that
double tt_minus_utc(double jd) {
    static const double tab[][2] = {{2441317.5, 10}, {2441499.5, 11}, {2441683.5, 12}, {2442048.5, 13}, {2442413.5, 14}, {2442778.5, 15},
                                    {2443144.5, 16}, {2443509.5, 17}, {2443874.5, 18}, {2444239.5, 19}, {2444786.5, 20}, {2445151.5, 21},
                                    {2445516.5, 22}, {2446247.5, 23}, {2447161.5, 24}, {2447892.5, 25}, {2448257.5, 26}, {2448804.5, 27},
                                    {2449169.5, 28}, {2449534.5, 29}, {2450083.5, 30}, {2450630.5, 31}, {2451179.5, 32}, {2453736.5, 33},
                                    {2454832.5, 34}, {2456109.5, 35}, {2457204.5, 36}, {2457754.5, 37}};
    if (jd < tab[0][0])
        return 42.184;
    double ls = 10;
    for (auto &r : tab)
        if (jd >= r[0])
            ls = r[1];
    return 32.184 + ls;
}

void moon_simple(double t, double v[3]) { // identical to the table generator's (only needs to match)
    double L = (218.316 + 13.176396 * t) * D2R, M = (134.963 + 13.064993 * t) * D2R, F = (93.272 + 13.229350 * t) * D2R,
           D = (297.850 + 12.190749 * t) * D2R;
    double lon = L + 6.289 * D2R * std::sin(M) + 1.274 * D2R * std::sin(2 * D - M) + 0.658 * D2R * std::sin(2 * D);
    double lat = 5.128 * D2R * std::sin(F), r = (385001 - 20905 * std::cos(M)) / 149597870.7;
    double x = r * std::cos(lat) * std::cos(lon), y = r * std::cos(lat) * std::sin(lon), z = r * std::sin(lat);
    v[0] = x;
    v[1] = y * std::cos(EPS0) - z * std::sin(EPS0);
    v[2] = y * std::sin(EPS0) + z * std::cos(EPS0);
}

// heliocentric Earth, equatorial J2000 (AU), jd in TDB(=TT)
void earth_helio(double jd, double e[3]) {
    double t = jd - EC_START;
    int k = int(std::floor(t / EC_SEG));
    k = std::max(0, std::min(EC_NSEG - 1, k));
    double x = 2 * (t - k * EC_SEG) / EC_SEG - 1;
    const long long *c = EC_COEF + size_t(k) * 3 * (EC_DEG + 1);
    for (int i = 0; i < 3; ++i) {
        double b0 = 0, b1 = 0, b2; // Clenshaw
        for (int j = EC_DEG; j >= 1; --j) {
            b2 = b1;
            b1 = b0;
            b0 = 2 * x * b1 - b2 + c[i * (EC_DEG + 1) + j] * 1e-11;
        }
        e[i] = x * b0 - b1 + c[i * (EC_DEG + 1)] * 1e-11;
    }
    double m[3];
    moon_simple(jd - J2000, m);
    for (int i = 0; i < 3; ++i)
        e[i] -= MU_MOON * m[i];
}

// geocentric observer position (AU, equatorial J2000): WGS84 + GMST + IAU 1976 precession (no nutation: < 1 km)
void observer_gcrs(double jd_utc, double jd_tt, double lon_e, double lat, double h_m, double o[3]) {
    const double a = 6378.137, f = 1 / 298.257223563, e2 = f * (2 - f);
    double sl = std::sin(lat * D2R), cl = std::cos(lat * D2R), N = a / std::sqrt(1 - e2 * sl * sl), hk = h_m / 1000;
    double rxy = (N + hk) * cl, z = (N * (1 - e2) + hk) * sl;
    double Tu = (jd_utc - J2000) / 36525;
    double gmst = std::fmod(280.46061837 + 360.98564736629 * (jd_utc - J2000) + 0.000387933 * Tu * Tu - Tu * Tu * Tu / 38710000, 360);
    double th = (gmst + lon_e) * D2R, v[3] = {rxy * std::cos(th), rxy * std::sin(th), z};
    double T = (jd_tt - J2000) / 36525, as = D2R / 3600;
    double zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) * as,
           zz = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) * as,
           theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) * as;
    double cz = std::cos(zeta), sz = std::sin(zeta), cZ = std::cos(zz), sZ = std::sin(zz), ct = std::cos(theta), st = std::sin(theta);
    double P[3][3] = {{cz * ct * cZ - sz * sZ, -sz * ct * cZ - cz * sZ, -st * cZ}, // J2000 -> data
                      {cz * ct * sZ + sz * cZ, -sz * ct * sZ + cz * cZ, -st * sZ},
                      {cz * st, -sz * st, ct}};
    for (int i = 0; i < 3; ++i)
        o[i] = (P[0][i] * v[0] + P[1][i] * v[1] + P[2][i] * v[2]) / 149597870.7; // P^T
}

// heliocentric equatorial J2000 position (AU)
void helio(const Orbit &el, double jd, double out[3]) {
    double xv, yv;
    if (!el.comet) {
        double M = (el.M + el.n * (jd - el.epoch)) * D2R, e = el.e, E = M;
        for (int i = 0; i < 30; ++i) {
            double d = (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
            E -= d;
            if (std::fabs(d) < 1e-13)
                break;
        }
        xv = el.a * (std::cos(E) - e);
        yv = el.a * std::sqrt(1 - e * e) * std::sin(E);
    } else {
        double q = el.q, e = el.e, dt = jd - el.Tp, nu;
        if (std::fabs(e - 1) < 1e-4) {
            double W = 3 * KG * dt / std::sqrt(2 * q * q * q), Y = std::cbrt(W / 2 + std::sqrt(W * W / 4 + 1));
            nu = 2 * std::atan(Y - 1 / Y);
        } else if (e < 1) {
            double a = q / (1 - e), M = KG * dt / std::pow(a, 1.5), E = M;
            for (int i = 0; i < 50; ++i) {
                double d = (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
                E -= d;
                if (std::fabs(d) < 1e-13)
                    break;
            }
            nu = 2 * std::atan2(std::sqrt(1 + e) * std::sin(E / 2), std::sqrt(1 - e) * std::cos(E / 2));
        } else {
            double a = q / (e - 1), M = KG * dt / std::pow(a, 1.5), H = std::asinh(M / e);
            for (int i = 0; i < 50; ++i) {
                double d = (e * std::sinh(H) - H - M) / (e * std::cosh(H) - 1);
                H -= d;
                if (std::fabs(d) < 1e-13)
                    break;
            }
            nu = 2 * std::atan(std::sqrt((e + 1) / (e - 1)) * std::tanh(H / 2));
        }
        double r = q * (1 + e) / (1 + e * std::cos(nu));
        xv = r * std::cos(nu);
        yv = r * std::sin(nu);
    }
    double w = el.w * D2R, O = el.node * D2R, i = el.inc * D2R;
    double cO = std::cos(O), sO = std::sin(O), cw = std::cos(w), sw = std::sin(w), ci = std::cos(i), si = std::sin(i);
    double xh = xv * (cO * cw - sO * sw * ci) - yv * (cO * sw + sO * cw * ci);
    double yh = xv * (sO * cw + cO * sw * ci) - yv * (sO * sw - cO * cw * ci);
    double zh = xv * sw * si + yv * cw * si;
    out[0] = xh;
    out[1] = yh * std::cos(EPS0) - zh * std::sin(EPS0);
    out[2] = yh * std::sin(EPS0) + zh * std::cos(EPS0);
}

struct Obs {
    double jd_tt, E[3];
}; // heliocentric observer at an epoch
Obs make_obs(const AstroSession *s, double jd_utc) {
    Obs o;
    o.jd_tt = jd_utc + tt_minus_utc(jd_utc) / 86400;
    earth_helio(o.jd_tt, o.E);
    // in the .cfg longitude is positive West (negative = East)
    double lon = -s->s.num("longitude"), lat = s->s.num("latitude"), h = s->s.num("altitude_m");
    if (lon != 0 || lat != 0) {
        double g[3];
        observer_gcrs(jd_utc, o.jd_tt, lon, lat, h, g);
        for (int i = 0; i < 3; ++i)
            o.E[i] += g[i];
    }
    return o;
}

struct Eph {
    double ra, dec, V, delta, rh;
};
Eph ephem(const Orbit &el, const Obs &o) {
    double r[3], g[3], delta = 0;
    for (int it = 0; it < 3; ++it) {
        helio(el, o.jd_tt - delta * LT_DAY, r);
        for (int i = 0; i < 3; ++i)
            g[i] = r[i] - o.E[i];
        delta = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    }
    Eph p;
    p.delta = delta;
    p.ra = std::fmod(std::atan2(g[1], g[0]) * R2D + 360, 360);
    p.dec = std::asin(g[2] / delta) * R2D;
    double rh = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    p.rh = rh;
    if (el.comet) {
        p.V = el.M1 + 5 * std::log10(delta) + el.K1 * std::log10(rh);
        return p;
    }
    double cph = (r[0] * g[0] + r[1] * g[1] + r[2] * g[2]) / (rh * delta), ph = std::acos(std::max(-1.0, std::min(1.0, cph)));
    double t = std::tan(ph / 2), p1 = std::exp(-3.33 * std::pow(t, 0.63)), p2 = std::exp(-1.87 * std::pow(t, 1.22));
    p.V = el.H + 5 * std::log10(rh * delta) - 2.5 * std::log10((1 - el.G) * p1 + el.G * p2);
    return p;
}

double sep_deg(double r1, double d1, double r2, double d2) {
    double c = std::sin(d1 * D2R) * std::sin(d2 * D2R) + std::cos(d1 * D2R) * std::cos(d2 * D2R) * std::cos((r1 - r2) * D2R);
    return std::acos(std::max(-1.0, std::min(1.0, c))) * R2D;
}

// ---------------------------------------------------------------- leitura
int c2n(char ch) {
    return std::isdigit((unsigned char)ch) ? ch - '0' : (ch >= 'a' ? ch - 'a' + 36 : ch - 'A' + 10);
}
bool num_field(const std::string &l, size_t a, size_t n, double &v, bool allow_empty = false, double dflt = 0) {
    if (a >= l.size()) {
        v = dflt;
        return allow_empty;
    }
    std::string f = l.substr(a, n);
    size_t b = f.find_first_not_of(' '), e = f.find_last_not_of(" \r\n");
    if (b == std::string::npos) {
        v = dflt;
        return allow_empty;
    }
    f = f.substr(b, e - b + 1);
    char *end;
    v = std::strtod(f.c_str(), &end);
    return *end == 0;
}
std::string trim(const std::string &x) {
    size_t b = x.find_first_not_of(" \r\n\t"), e = x.find_last_not_of(" \r\n\t");
    return b == std::string::npos ? "" : x.substr(b, e - b + 1);
}
std::string unpack_desig(const std::string &p) {
    if (p.size() == 7 && (p[0] == 'I' || p[0] == 'J' || p[0] == 'K') && std::isdigit((unsigned char)p[1]) &&
        std::isdigit((unsigned char)p[5])) {
        int y = (p[0] == 'I' ? 18 : p[0] == 'J' ? 19 : 20) * 100 + std::atoi(p.substr(1, 2).c_str());
        int n = c2n(p[4]) * 10 + (p[5] - '0');
        return std::to_string(y) + " " + p[3] + p[6] + (n ? std::to_string(n) : "");
    }
    return p;
}
long long orbit_number(const Orbit &o) { // numbered: "(N) Name" or numeric packed designation
    if (o.comet)
        return 0;
    if (!o.name.empty() && o.name[0] == '(')
        return std::atoll(o.name.c_str() + 1);
    const std::string &p = o.packed;
    if (p.size() == 5 && std::isalnum((unsigned char)p[0])) {
        for (size_t i = 1; i < 5; ++i)
            if (!std::isdigit((unsigned char)p[i]))
                return 0;
        return c2n(p[0]) * 10000LL + std::atoi(p.c_str() + 1);
    }
    return 0;
}

int read_mpcorb(const std::string &path, std::vector<Orbit> &out) {
    std::ifstream fh(path, std::ios::binary);
    if (!fh)
        return 0;
    int n = 0;
    std::string l;
    while (std::getline(fh, l)) {
        if (l.size() < 103)
            continue;
        std::string ep = l.substr(20, 5);
        if (!(ep[0] == 'I' || ep[0] == 'J' || ep[0] == 'K') || !std::isdigit((unsigned char)ep[1]) || !std::isdigit((unsigned char)ep[2]))
            continue;
        Orbit o;
        double dummy;
        if (!num_field(l, 8, 5, o.H, true, 99) || !num_field(l, 14, 5, o.G, true, 0.15) || !num_field(l, 26, 9, o.M) ||
            !num_field(l, 37, 9, o.w) || !num_field(l, 48, 9, o.node) || !num_field(l, 59, 9, o.inc) || !num_field(l, 70, 9, o.e) ||
            !num_field(l, 80, 11, o.n) || !num_field(l, 92, 11, o.a))
            continue;
        (void)dummy;
        if (o.e >= 1 || o.a <= 0)
            continue;
        o.epoch = cal_jd((ep[0] == 'I' ? 1800 : ep[0] == 'J' ? 1900 : 2000) + std::atoi(ep.substr(1, 2).c_str()), c2n(ep[3]), c2n(ep[4]));
        o.packed = trim(l.substr(0, 7));
        o.name = l.size() > 166 ? trim(l.substr(166, 28)) : "";
        if (o.name.empty())
            o.name = unpack_desig(o.packed);
        o.q = o.a * (1 - o.e);
        out.push_back(std::move(o));
        ++n;
    }
    return n;
}

int read_comets(const std::string &path, std::vector<Orbit> &out) {
    std::ifstream fh(path, std::ios::binary);
    if (!fh)
        return 0;
    int n = 0;
    std::string l;
    while (std::getline(fh, l)) {
        if (l.size() < 100)
            continue;
        Orbit o;
        o.comet = true;
        double y, m, d;
        if (!num_field(l, 14, 4, y) || !num_field(l, 19, 2, m) || !num_field(l, 22, 7, d) || !num_field(l, 30, 9, o.q) ||
            !num_field(l, 41, 8, o.e) || !num_field(l, 51, 8, o.w) || !num_field(l, 61, 8, o.node) || !num_field(l, 71, 8, o.inc) ||
            !num_field(l, 91, 4, o.M1, true, 99) || !num_field(l, 96, 4, o.K1, true, 10))
            continue;
        o.Tp = cal_jd(int(y), int(m), d);
        o.H = o.M1;
        o.G = 0.15;
        o.packed = trim(l.substr(0, 12));
        o.name = l.size() > 102 ? trim(l.substr(102, 56)) : "";
        if (o.name.empty())
            o.name = o.packed;
        if (o.e < 1)
            o.a = o.q / (1 - o.e);
        out.push_back(std::move(o));
        ++n;
    }
    return n;
}

// ---------------------------------------------------------------- JSON
std::string jesc(const std::string &x) {
    return session_util::json_esc(x);
}
struct Hit {
    const Orbit *o;
    Eph p;
    double sep, rate, pa, x = 0, y = 0;
};
void motion(const Orbit &el, const AstroSession *s, double jd, const Eph &p, double &rate, double &pa) {
    Obs o2 = make_obs(s, jd + 10.0 / 1440);
    Eph q = ephem(el, o2);
    double dra = q.ra - p.ra;
    if (dra > 180)
        dra -= 360;
    if (dra < -180)
        dra += 360;
    double vx = dra * std::cos(p.dec * D2R) * 3600 / 10, vy = (q.dec - p.dec) * 3600 / 10;
    rate = std::hypot(vx, vy);
    pa = std::fmod(std::atan2(vx, vy) * R2D + 360, 360);
}
std::string to_json(const std::vector<Hit> &v, bool xy) {
    std::string out = "[";
    char b[512];
    for (size_t i = 0; i < v.size(); ++i) {
        const Hit &h = v[i];
        long long num = orbit_number(*h.o);
        out += i ? ",{" : "{";
        out += "\"name\":\"" + jesc(h.o->name) + "\",\"packed\":\"" + jesc(h.o->packed) + "\",";
        if (num) {
            std::snprintf(b, sizeof b, "\"number\":%lld,", num);
            out += b;
        }
        std::snprintf(b, sizeof b,
                      "\"kind\":\"%s\",\"ra\":%.7f,\"dec\":%.7f,\"sep_arcsec\":%.2f,\"vmag\":%.2f,\"rate_arcsec_min\":%.4f,\"pa_deg\":%.2f,"
                      "\"delta_au\":%.6f",
                      h.o->comet ? "comet" : "asteroid", h.p.ra, h.p.dec, h.sep * 3600, h.p.V, h.rate, h.pa, h.p.delta);
        out += b;
        if (xy) {
            std::snprintf(b, sizeof b, ",\"x\":%.2f,\"y\":%.2f", h.x, h.y);
            out += b;
        }
        out += "}";
    }
    return out + "]";
}
std::string g_buf;

// candidates within rad_deg of (ra, dec): geometric pre-filter without light time, then full ephemeris
std::vector<Hit> cone(AstroSession *s, double jd, double ra, double dec, double rad_deg, double vmax) {
    std::vector<Hit> out;
    Obs o = make_obs(s, jd);
    double u[3] = {std::cos(dec * D2R) * std::cos(ra * D2R), std::cos(dec * D2R) * std::sin(ra * D2R), std::sin(dec * D2R)};
    double cmin = std::cos(std::min(180.0, rad_deg + 3.0) * D2R);
    for (const Orbit &el : s->orbits) {
        double r[3];
        helio(el, o.jd_tt, r);
        double g[3] = {r[0] - o.E[0], r[1] - o.E[1], r[2] - o.E[2]}, d = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
        if (!(d > 0))
            continue;
        // margin: 3° + motion during light time (very close objects)
        if ((g[0] * u[0] + g[1] * u[1] + g[2] * u[2]) / d < cmin && d > 0.02)
            continue;
        Eph p = ephem(el, o);
        if (!std::isfinite(p.ra) || !std::isfinite(p.dec))
            continue;
        double sp = sep_deg(p.ra, p.dec, ra, dec);
        if (sp > rad_deg || p.V > vmax)
            continue;
        Hit h{&el, p, sp, 0, 0};
        motion(el, s, jd, p, h.rate, h.pa);
        out.push_back(h);
    }
    return out;
}
} // namespace

extern "C" {
int astro_load_orbits(AstroSession *s, const char *mpcorb_path, const char *comet_path) {
    try {
        if (!s)
            return 0;
        s->orbits.clear();
        if (mpcorb_path && *mpcorb_path)
            read_mpcorb(mpcorb_path, s->orbits);
        if (comet_path && *comet_path)
            read_comets(comet_path, s->orbits);
        return int(s->orbits.size());
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return 0;
    } catch (...) {
        return 0;
    }
}

const char *astro_identify_offline(AstroSession *s, double jd, double ra, double dec, double radius_deg) {
    try {
        if (!s) {
            g_buf = "[]";
            return g_buf.c_str();
        }
        std::string &B = s->buf3;
        B = "[]";
        auto v = cone(s, jd, ra, dec, radius_deg, 1e9);
        std::sort(v.begin(), v.end(), [](const Hit &a, const Hit &b) { return a.sep < b.sep; });
        B = to_json(v, false);
        return B.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "[]";
    } catch (...) {
        return "[]";
    }
}

const char *astro_known_offline(AstroSession *s, int frame, double vmax) {
    try {
        if (!s) {
            g_buf = "[]";
            return g_buf.c_str();
        }
        std::string &B = s->buf2;
        B = "[]";
        if (!s->ok || frame < 0 || frame >= int(s->R.frames.size()) || !s->R.frames[frame].solved)
            return B.c_str();
        const astro::Frame &f = s->R.frames[frame];
        const astro::Plate &p = f.plate;
        int w = f.im.w, h = f.im.h;
        double cra, cde;
        p.pix_to_radec(w / 2.0, h / 2.0, cra, cde);
        double rad = 0;
        for (int cx = 0; cx <= 1; ++cx)
            for (int cy = 0; cy <= 1; ++cy) {
                double r, d;
                p.pix_to_radec(cx * w, cy * h, r, d);
                rad = std::max(rad, sep_deg(r, d, cra, cde));
            }
        auto v = cone(s, f.im.jd_mid, cra, cde, rad + 0.05, vmax);
        std::vector<Hit> out;
        for (Hit &hh : v) {
            double xi, eta;
            astro::radec_to_std(hh.p.ra, hh.p.dec, p.ra0, p.dec0, xi, eta);
            // plate inverse (like astro_sky_to_pix): linear part + Newton
            double a = p.cx[1], b = p.cx[2], c = p.cy[1], d = p.cy[2], det = a * d - b * c;
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
            if (px < -5 || py < -5 || px > w + 5 || py > h + 5)
                continue;
            hh.x = px;
            hh.y = py;
            out.push_back(hh);
        }
        std::sort(out.begin(), out.end(), [](const Hit &a, const Hit &b) { return a.p.V < b.p.V; });
        B = to_json(out, true);
        return B.c_str();
    } catch (const std::exception &) { // C++ exceptions must not cross the extern "C" boundary
        return "[]";
    } catch (...) {
        return "[]";
    }
}
}
