// ADES (PSV), photometric and MPC 80-column reports, per the public MPC specifications
// (ADES: https://minorplanetcenter.net/iau/info/ADES.html; 80-column format: https://minorplanetcenter.net/iau/info/OpticalObs.html).
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace astro {

static std::string f(const char *fmt, double v) {
    char b[64];
    std::snprintf(b, sizeof b, fmt, v);
    return b;
}

static void jd_to_cal(double jd, int &y, int &mo, int &d, double &dayfrac) {
    double z = std::floor(jd + 0.5), F = jd + 0.5 - z;
    long Z = long(z), A = Z;
    if (Z >= 2299161) {
        long al = long((Z - 1867216.25) / 36524.25);
        A = Z + 1 + al - al / 4;
    }
    long B = A + 1524, C = long((B - 122.1) / 365.25), D = long(365.25 * C), E = long((B - D) / 30.6001);
    d = int(B - D - long(30.6001 * E));
    mo = int(E < 14 ? E - 1 : E - 13);
    y = int(mo > 2 ? C - 4716 : C - 4715);
    dayfrac = F;
}

std::string jd_to_isot(double jd) {
    int y, mo, d;
    double fr;
    jd_to_cal(jd, y, mo, d, fr);
    long ms = std::lround(fr * 86400000.0);
    if (ms >= 86400000) {
        jd_to_cal(jd + 0.5 / 86400.0, y, mo, d, fr);
        ms = 0;
    }
    long s = ms / 1000;
    char b[32];
    std::snprintf(b, sizeof b, "%04d-%02d-%02dT%02ld:%02ld:%02ld", y, mo, d, s / 3600, (s / 60) % 60, s % 60);
    return b;
}

// "YYYY-MM-DDTHH:MM:SS.sss" (rounded milliseconds) for the ADES obsTime: truncating to the whole second shifts
// the epoch by 0.5 s on average, which for an object at v ″/min is a v/120 ″ bias along the motion.
static std::string jd_to_isot_ms(double jd) {
    int y, mo, d;
    double fr;
    jd_to_cal(jd, y, mo, d, fr);
    long ms = std::lround(fr * 86400000.0);
    if (ms >= 86400000) {
        jd_to_cal(jd + 0.5 / 86400.0, y, mo, d, fr);
        ms = 0;
    }
    long s = ms / 1000;
    char b[40];
    std::snprintf(b, sizeof b, "%04d-%02d-%02dT%02ld:%02ld:%02ld.%03ld", y, mo, d, s / 3600, (s / 60) % 60, s % 60, ms % 1000);
    return b;
}

std::string fmt_ra(double ra, bool extra) {
    int nd = extra ? 3 : 2;
    double h = ra / 15;
    int hh = int(h);
    double m = (h - hh) * 60;
    int mm = int(m);
    double s = (m - mm) * 60;
    double p = std::pow(10.0, nd);
    if (std::nearbyint(s * p) / p >= 60) {
        s = std::max(s - 60, 0.0); // 59.97″ -> 0, not "-0.0" (spurious sign that breaks the MPC format columns)
        mm += 1;
    }
    if (mm >= 60) {
        mm -= 60;
        hh += 1;
    }
    char b[32];
    std::snprintf(b, sizeof b, "%02d %02d %0*.*f", hh % 24, mm, nd + 3, nd, s);
    return b;
}
std::string fmt_dec(double de, bool extra) {
    int nd = extra ? 2 : 1;
    char sg = de < 0 ? '-' : '+';
    de = std::fabs(de);
    int dd = int(de);
    double m = (de - dd) * 60;
    int mm = int(m);
    double s = (m - mm) * 60;
    double p = std::pow(10.0, nd);
    if (std::nearbyint(s * p) / p >= 60) {
        s = std::max(s - 60, 0.0); // 59.97″ -> 0, not "-0.0" (spurious sign that breaks the MPC format columns)
        mm += 1;
    }
    if (mm >= 60) {
        mm -= 60;
        dd += 1;
    }
    char b[32];
    std::snprintf(b, sizeof b, "%c%02d %02d %0*.*f", sg, dd, mm, nd + 3, nd, s);
    return b;
}

static bool isnum(const std::string &p) {
    for (char c : p)
        if (c < '0' || c > '9')
            return false;
    return !p.empty();
}

// designation without control characters (line breaks etc. would break the report lines)
static std::string clean_desig(const std::string &x) {
    std::string o;
    for (char c : x)
        o += (unsigned char)c < 0x20 || c == 0x7f || c == '|' ? ' ' : c;
    return o;
}

std::string packed_from(const std::string &desig) {
    std::string d;
    for (char c : clean_desig(desig))
        if (c != '(' && c != ')')
            d += c;
    while (!d.empty() && d.front() == ' ')
        d.erase(0, 1);
    while (!d.empty() && d.back() == ' ')
        d.pop_back();
    if (isnum(d)) {
        // number with unbounded digits: only up to 619999 is packable (A–Z = 10–35, a–z = 36–61); above that, or with
        // too many digits for a long, returns the designation unchanged (no std::stol, which would throw)
        size_t z = d.find_first_not_of('0');
        std::string sig = z == std::string::npos ? "0" : d.substr(z);
        if (sig.size() > 6)
            return d;
        long n = std::strtol(sig.c_str(), nullptr, 10);
        if (n > 619999)
            return d;
        char b[16];
        if (n < 100000)
            std::snprintf(b, sizeof b, "%05ld", n);
        else {
            long v = n / 10000;
            std::snprintf(b, sizeof b, "%c%04ld", v < 36 ? char('A' + v - 10) : char('a' + v - 36), n % 10000);
        }
        return b;
    }
    // "2001 SP75" -> K01S75P
    std::string s;
    for (char c : d)
        if (c != ' ')
            s += char(std::toupper(c));
    if (s.size() >= 6 && isnum(s.substr(0, 4)) && std::isalpha(s[4]) && std::isalpha(s[5])) {
        std::string cent = s.substr(0, 2) == "18" ? "I" : s.substr(0, 2) == "19" ? "J" : "K";
        // cycle number: digits only, at most 3 (up to 619 fits the format); otherwise do not pack
        std::string tail = s.substr(6);
        if (!tail.empty() && (!isnum(tail) || tail.size() > 3))
            return d.substr(0, 7);
        int n = tail.empty() ? 0 : std::atoi(tail.c_str());
        if (n > 619)
            return d.substr(0, 7);
        char cyc[4];
        if (n < 100)
            std::snprintf(cyc, sizeof cyc, "%02d", n);
        else
            std::snprintf(cyc, sizeof cyc, "%c%d", n / 10 < 36 ? char('A' + n / 10 - 10) : char('a' + n / 10 - 36), n % 10);
        return cent + s.substr(2, 2) + s[4] + cyc + s[5];
    }
    return d.substr(0, 7);
}

// Temporary designation/trkSub for an unidentified track (ADES trkSub ≤ 8 characters; MPC 80-col., columns 6–12
// ≤ 7): 'P' + 4 base-36 digits of the UTC minute of the first measurement (~3.2-year cycle) + 2 for the track index.
std::string temp_id(const Track &T, size_t idx) {
    static const char *B36 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    double jd = T.meas.empty() ? 2451544.5 : T.meas.front().jd;
    for (auto &m : T.meas)
        jd = std::min(jd, m.jd);
    long long mn = (long long)std::floor((jd - 2451544.5) * 1440.0);
    mn = ((mn % 1679616) + 1679616) % 1679616; // 36⁴
    std::string o = "P";
    for (long long d = 36 * 36 * 36; d >= 1; d /= 36)
        o += B36[(mn / d) % 36];
    o += B36[(idx / 36) % 36];
    o += B36[idx % 36];
    return o;
}

static std::string unpack_desig(const std::string &p) {
    if (p.size() == 7 && (p[0] == 'I' || p[0] == 'J' || p[0] == 'K') && std::isdigit((unsigned char)p[1]) &&
        std::isdigit((unsigned char)p[2]) && std::isalnum((unsigned char)p[4]) && std::isdigit((unsigned char)p[5])) {
        int y = (p[0] == 'I' ? 18 : p[0] == 'J' ? 19 : 20) * 100 + std::stoi(p.substr(1, 2));
        int n = (std::isdigit((unsigned char)p[4]) ? (p[4] - '0') * 10 : (std::isupper((unsigned char)p[4]) ? p[4] - 55 : p[4] - 61) * 10) +
                (p[5] - '0');
        return std::to_string(y) + " " + p[3] + p[6] + (n ? std::to_string(n) : "");
    }
    return p;
}

struct Col {
    const char *name;
    int w;
    char al;
};
static const Col ADES[] = {{"permID", 7, 'l'},   {"provID", 11, 'l'}, {"trkSub", 8, 'l'}, {"mode", 4, 'r'},   {"stn", 4, 'l'},
                           {"obsTime", 24, 'l'}, {"ra", 11, 'l'},     {"dec", 11, 'l'},   {"rmsRA", 6, 'l'},  {"rmsDec", 6, 'l'},
                           {"astCat", 8, 'r'},   {"mag", 5, 'l'},     {"rmsMag", 6, 'l'}, {"band", 4, 'r'},   {"fltr", 4, 'l'},
                           {"photCat", 8, 'r'},  {"photAp", 6, 'l'},  {"logSNR", 6, 'l'}, {"seeing", 6, 'l'}, {"exp", 4, 'r'},
                           {"rmsFit", 6, 'l'},   {"nStars", 6, 'r'},  {"notes", 5, 'l'},  {"remarks", 0, 'l'}};

static std::string pad(const std::string &v, int w, char al, bool header) {
    if (w == 0 || int(v.size()) >= w)
        return v;
    return (al == 'r' && !header) ? std::string(w - v.size(), ' ') + v : v + std::string(w - v.size(), ' ');
}

// Uncertainties (ADES rmsRA/rmsDec), per axis: σ² = σ_c² + σ_sys² + σ_pred², with σ_c the centroid Cramér–Rao bound
// for a Gaussian PSF (King 1983, PASP 95, 163; Mighell 2005, MNRAS 361, 861), σ_sys the systematic plate component
// estimated on the references (reduced χ² = 1 after removing their centroid noise) and σ_pred the plate's own
// uncertainty at the point (coefficient covariance). Without that model (manual measurement / old plate):
// σ_c ≈ FWHM/(2.355·SNR) ⊕ out-of-sample plate rms. cra/cde (optional): centroid term only.
void measurement_sigma(const Reduction &R, const Track &T, const Measurement &m, const Settings &s, double &era, double &ede, double *cra,
                       double *cde) {
    const Frame &fr = R.frames[m.frame];
    const Plate &p = fr.plate;
    double st = m.snr_tot > 0 ? m.snr_tot : m.snr, sn = std::max(st, 1e-3);
    int di = m.frame < int(T.idx.size()) ? T.idx[m.frame] : -1;
    double scr, scd;
    // anisotropic centroid covariance (trails): C_sky = J·C_px·Jᵀ, with J the local plate Jacobian (″/px).
    // No fit (cxx missing) but a detected trail: analytic Cramér–Rao model for segment ⊗ PSF — along the trail
    // σ² = (σ_PSF² + L²/12)/N·(1 + 2·A·B/N), across σ_PSF²/N·(1 + 2·A·B/N), A = 4πσ_PSF²·(1 + L/(2√π σ_PSF)) the
    // effective area (4πσ² for a point source, as in centroid_sigma) and B the sky per pixel (e⁻)
    double C[3] = {m.cxx, m.cxy, m.cyy};
    const bool trailcov = s.str("trail_cov").empty() || s.num("trail_cov") != 0;
    if (trailcov && !std::isfinite(C[0]) && di >= 0 && di < int(fr.det.size()) && fr.det[di].trail_len > 0 && fr.seeing > 0) {
        const Detection &d = fr.det[di];
        const double g = R.gain > 0 ? R.gain : 1.0, scd = 0.5 * (p.scale_x() + p.scale_y());
        const double sp = fr.seeing / scd / 2.3548, Lp = d.trail_len, Ne = std::max(d.flux, 1e-9) * g, B = d.sigma * d.sigma * g * g;
        const double A = 4 * PI * sp * sp * (1 + Lp / (2 * std::sqrt(PI) * sp)), k = (1 + 2 * A * B / Ne) / Ne;
        const double va = (sp * sp + Lp * Lp / 12) * k, vc = sp * sp * k, ux = std::cos(d.trail_ang), uy = std::sin(d.trail_ang);
        C[0] = va * ux * ux + vc * uy * uy;
        C[1] = (va - vc) * ux * uy;
        C[2] = va * uy * uy + vc * ux * ux;
    }
    if (trailcov && std::isfinite(C[0]) && std::isfinite(C[1]) && std::isfinite(C[2]) && !p.cov_ra.empty()) {
        double J[4];
        plate_jacobian(p, m.x, m.y, J);
        const double sxx = J[0] * J[0] * C[0] + 2 * J[0] * J[1] * C[1] + J[1] * J[1] * C[2];
        const double syy = J[2] * J[2] * C[0] + 2 * J[2] * J[3] * C[1] + J[3] * J[3] * C[2];
        double pr, pd;
        p.pred_sigma(m.x, m.y, pr, pd);
        if (cra)
            *cra = std::sqrt(std::max(sxx, 0.0));
        if (cde)
            *cde = std::sqrt(std::max(syy, 0.0));
        era = std::sqrt(std::max(sxx, 0.0) + p.sys_ra * p.sys_ra + pr * pr);
        ede = std::sqrt(std::max(syy, 0.0) + p.sys_dec * p.sys_dec + pd * pd);
        return;
    }
    if (di >= 0 && di < int(fr.det.size()) && !p.cov_ra.empty()) {
        double scale = 0.5 * (p.scale_x() + p.scale_y());
        double sc = internal::centroid_sigma(fr.det[di], scale, R.gain > 0 ? R.gain : (s.str("gain").empty() ? 1.0 : s.num("gain")));
        double pr, pd;
        p.pred_sigma(m.x, m.y, pr, pd);
        scr = scd = sc;
        era = std::sqrt(sc * sc + p.sys_ra * p.sys_ra + pr * pr);
        ede = std::sqrt(sc * sc + p.sys_dec * p.sys_dec + pd * pd);
    } else {
        double sc = m.fwhm_as / (2.355 * sn);
        double pra = p.rms_cv_ra > 0 ? p.rms_cv_ra : p.rms_ra, pde = p.rms_cv_dec > 0 ? p.rms_cv_dec : p.rms_dec;
        scr = scd = sc;
        era = std::hypot(pra, sc);
        ede = std::hypot(pde, sc);
    }
    if (cra)
        *cra = scr;
    if (cde)
        *cde = scd;
}

std::string ades_report(const Reduction &R, const std::vector<Track> &tracks, const Settings &s) {
    std::vector<std::string> L{"# version=2022", "# observatory", "! mpcCode " + s.str("mpc_code"), "# submitter"};
    std::string sub = !s.str("submitter").empty()  ? s.str("submitter")
                      : !s.str("measurer").empty() ? s.str("measurer")
                      : !s.str("observer").empty() ? s.str("observer")
                                                   : "Unknown";
    L.push_back("! name " + sub);
    if (!s.str("observer").empty()) {
        L.push_back("# observers");
        L.push_back("! name " + s.str("observer"));
    }
    L.push_back("# measurers");
    L.push_back("! name " + (!s.str("measurer").empty()   ? s.str("measurer")
                             : !s.str("observer").empty() ? s.str("observer")
                                                          : std::string("Unknown")));
    // # telescope: design and aperture are mandatory in ADES; without telescope_design/telescope_aperture_m they are derived from the
    // free-text "telescope" (e.g. "1.2-m f/2.4 Schmidt + CCD"), and when that is missing "Unknown"/omitted aperture
    // is not allowed — the aperture is computed from focal length and focal ratio (f/N) when present
    std::string tel = s.str("telescope"), tl, design = s.str("telescope_design");
    for (char c : tel)
        tl += char(std::tolower((unsigned char)c));
    if (design.empty())
        design = tl.find("schmidt") != std::string::npos                                                 ? "Schmidt"
                 : tl.find("refractor") != std::string::npos                                             ? "Refractor"
                 : tl.find("ritchey") != std::string::npos || tl.find("cassegrain") != std::string::npos ? "Cassegrain"
                 : tl.find("reflector") != std::string::npos || tl.find("newton") != std::string::npos   ? "Reflector"
                                                                                                         : "Unknown";
    double ap = s.num("telescope_aperture_m");
    if (!(ap > 0)) {
        double v = 0, N = 0;
        for (size_t i = 0; i < tl.size() && !(v > 0); ++i)
            if ((std::isdigit((unsigned char)tl[i]) || tl[i] == '.') && (i == 0 || !std::isalnum((unsigned char)tl[i - 1]))) {
                char *e = nullptr;
                double x = std::strtod(tl.c_str() + i, &e);
                size_t j = size_t(e - tl.c_str());
                while (j < tl.size() && (tl[j] == ' ' || tl[j] == '-'))
                    ++j;
                if (j < tl.size() && tl[j] == 'm' && (j + 1 == tl.size() || !std::isalpha((unsigned char)tl[j + 1])))
                    v = x;
                else if (tl.compare(j, 2, "cm") == 0)
                    v = x / 100;
                i = std::max(i, j > 0 ? j - 1 : 0);
            }
        size_t q = tl.find("f/");
        if (q != std::string::npos)
            N = std::strtod(tl.c_str() + q + 2, nullptr);
        ap = v > 0 ? v : (N > 0 && s.num("focal_length_mm") > 0 ? s.num("focal_length_mm") / N / 1000 : 0);
    }
    L.push_back("# telescope");
    L.push_back("! design " + design);
    L.push_back("! aperture " + f(ap > 0 && ap < 0.95 ? "%.2f" : "%.1f", ap));
    L.push_back("! detector " + (s.str("telescope_detector").empty() ? std::string("CCD") : s.str("telescope_detector")));
    L.push_back("# software");
    L.push_back("! astrometry Piazzi");
    L.push_back("! photometry Piazzi");
    std::string h;
    for (size_t i = 0; i < 24; ++i) {
        if (i)
            h += '|';
        h += pad(ADES[i].name, ADES[i].w, ADES[i].al, true);
    }
    L.push_back(h);
    bool e6 = s.num("ades_six_digit_positions") != 0, m2 = s.num("ades_two_digit_mags") != 0;
    for (size_t ti = 0; ti < tracks.size(); ++ti) {
        const Track &T = tracks[ti];
        std::string packed = packed_from(T.desig), trk = temp_id(T, ti);
        bool num = isnum(packed) || (packed.size() == 5 && isnum(packed.substr(1)));
        for (auto &m : T.meas) {
            const Frame &fr = R.frames[m.frame];
            const Plate &p = fr.plate;
            double st = m.snr_tot > 0 ? m.snr_tot : m.snr, sn = std::max(st, 1e-3);
            double era, ede;
            measurement_sigma(R, T, m, s, era, ede);
            double emag = std::hypot(p.zp_rms, 1.0857 / sn);
            // position precision: 5 decimal places in degrees (step 0.036″, rms rounding error 0.036/√12 = 0.010″
            // per axis) only when that is negligible relative to σ; with σ < 0.1″ (sharp images, high SNR) 6
            // places are used, as in ades_six_digit_positions — otherwise rounding becomes a relevant part of the reported error
            const bool fine = std::min(era, ede) < 0.1;
            std::string v[24] = {num ? packed : "", num ? "" : (packed.size() == 7 ? unpack_desig(packed) : clean_desig(T.desig)), trk,
                                 "CCD", s.str("mpc_code"), jd_to_isot_ms(m.jd) + "Z", f(e6 || fine ? "%.6f" : "%.5f", m.ra),
                                 f(e6 || fine ? "%+.6f" : "%+.5f", m.dec),
                                 // 3 decimals: with σ of 0.01–0.03″ (SDSS, sharp images) 2 decimals round σ by up to ±35%
                                 f("%.3f", era), f("%.3f", ede), "Gaia3", f(m2 ? "%.2f" : "%.1f", m.mag), f("%.2f", emag), "Gr", "",
                                 "Gaia3", f("%.1f", s.num("aperture_radius") * 0.5 * (p.scale_x() + p.scale_y())), // aperture radius in ″
                                 f("%.2f", std::log10(std::max(st, 1.0))), f("%.1f", fr.seeing), f("%.0f", fr.im.exptime),
                                 f("%.2f", std::hypot(p.rms_ra, p.rms_dec)), std::to_string(p.nref), m.note, m.remark};
            if (num)
                while (!v[0].empty() && v[0][0] == '0')
                    v[0].erase(0, 1);
            std::string line;
            for (size_t i = 0; i < 24; ++i) {
                if (i)
                    line += '|';
                line += pad(v[i], ADES[i].w, ADES[i].al, false);
            }
            L.push_back(line);
        }
    }
    std::string out;
    for (auto &l : L)
        out += l + "\r\n";
    return out;
}

std::string phot_report(const Reduction &R, const std::vector<Track> &tracks, const Settings &s) {
    std::string sep(63, '-');
    std::vector<std::string> L{"OBSERVER:    " + s.str("observer"),
                               "OBS CODE:    " + s.str("mpc_code"),
                               "TELESCOPE:   " + f("%.1f", s.num("telescope_aperture_m")) + "m " + s.str("telescope_design"),
                               "EXPOSURE JD: Mid-exposure, not corrected for light time",
                               sep,
                               "      JD         mag     Flt   SNR    ZeroPt      Cat   Design.",
                               sep};
    for (auto &T : tracks)
        for (auto &m : T.meas) {
            char b[200];
            std::snprintf(b, sizeof b, "%.5f   %6.3f %2s C  %6.2f   %6.3f   %6s   %s", m.jd, m.mag, "Gr", m.snr_tot > 0 ? m.snr_tot : m.snr,
                          R.frames[m.frame].plate.zp, "Gaia3", clean_desig(T.desig).c_str());
            L.push_back(b);
        }
    L.push_back("----- end -----");
    std::string out;
    for (auto &l : L)
        out += l + "\r\n";
    return out;
}

std::string mpc_report(const Reduction &, const std::vector<Track> &tracks, const Settings &s) {
    std::ostringstream o;
    o << "COD " << s.str("mpc_code") << "\nOBS " << s.str("observer") << "\nMEA " << s.str("measurer") << "\nTEL " << s.str("telescope")
      << "\nNET Gaia-DR3\nACK Piazzi\n";
    bool extra = s.num("extra_time_digit") != 0;
    for (size_t ti = 0; ti < tracks.size(); ++ti) {
        const Track &T = tracks[ti];
        std::string packed = packed_from(T.desig);
        if (packed.empty())
            packed = temp_id(T, ti); // unidentified object: observer's temporary designation
        bool num = isnum(packed) || (packed.size() == 5 && isnum(packed.substr(1)));
        char p12[16];
        if (num)
            std::snprintf(p12, sizeof p12, "%5s       ", packed.c_str());
        else
            std::snprintf(p12, sizeof p12, "     %-7s", packed.c_str());
        for (auto &m : T.meas) {
            int y, mo, d;
            double fr;
            jd_to_cal(m.jd, y, mo, d, fr);
            char date[32];
            std::snprintf(date, sizeof date, "%04d %02d %08.5f", y, mo, d + fr);
            char mag[16];
            if (s.num("report_mags") != 0 && std::isfinite(m.mag))
                std::snprintf(mag, sizeof mag, "%4.1f %c", m.mag, s.str("report_band").empty() ? 'R' : s.str("report_band")[0]);
            else
                std::snprintf(mag, sizeof mag, "      ");
            char line[128];
            std::snprintf(line, sizeof line, "%-12.12s %cC%s %-12s%-12s         %s      %s", p12, m.note.empty() ? ' ' : m.note[0], date,
                          fmt_ra(m.ra, extra).c_str(), fmt_dec(m.dec, extra).c_str(), mag, s.str("mpc_code").c_str());
            o << line << "\n";
        }
    }
    return o.str();
}

} // namespace astro
