// End-to-end accuracy test on the bundled example (cpp/gui/assets/example: 3 ZTF frames of (5157) Hindemith)
// against independent references:
//   - catalog: the example Gaia DR3 (cpp/gui/assets/example/gaia_dr3.tsv);
//   - ephemerides: JPL Horizons at the mid-exposure instants (cpp/tests/data/horizons_hindemith.txt).
// Requires: one moving object identified as (5157) Hindemith (offline orbits from cpp/tests/data/orbits_mpcorb.dat),
// positions close to Horizons in each image (limits ~2–3× the measured error), rate ≈ Horizons, plates with ≥ 50
// references and small rms, mag 16–18.5 with small amplitude across images; repeats with automatic plate order.
// Uso: test_accuracy [raiz-do-repo]
#include "piazzi/astro_c.h"
#include "piazzi/reduce.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace astro;

// example field center (07 11 43.534 +22 24 17.78) and Horizons motion (dRA*cosD -21.02"/h, dDec +2.45"/h)
static const double RA0 = 107.931392, DE0 = 22.404939, HZ_RATE = 0.3527, HZ_PA = 276.65;
// limits ~2–3× measured (max O−C 0.021", rms 0.016"; plate rms ≤ 0.053"; amplitude 0.09 mag): max O−C/rms, plate rms,
// magnitude amplitude
static const double DMAX = 0.06, WMAX = 0.06, RMSMAX = 0.045, PRMS = 0.15, AMP = 0.25;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

int main(int argc, char **argv) {
    std::string root = argc > 1 ? std::string(argv[1]) + "/" : "";
    // Horizons ephemerides: "date time m? RA DEC dRA dDEC ..." (lines without '#')
    std::vector<std::pair<double, double>> hz;
    {
        std::ifstream f(root + "cpp/tests/data/horizons_hindemith.txt");
        std::string l;
        while (std::getline(f, l)) {
            if (l.empty() || l[0] == '#')
                continue;
            std::istringstream ss(l);
            std::string d, t, tok;
            ss >> d >> t;
            std::vector<double> v;
            while (ss >> tok) {
                try {
                    size_t used = 0;
                    double x = std::stod(tok, &used);
                    if (used == tok.size())
                        v.push_back(x);
                } catch (...) {
                }
            }
            if (v.size() >= 2)
                hz.push_back({v[0], v[1]});
        }
    }
    check(hz.size() == 3, "Horizons ephemerides (3 instants)");
    Settings s = Settings::defaults();
    const std::string ex = root + "cpp/gui/assets/example/";
    s.load(ex + "example.cfg");
    std::vector<Image> ims;
    for (int k = 0; k < 3; ++k)
        ims.push_back(load_image(ex + "ztf_hindemith_" + std::to_string(k) + ".fits", s.num("flip_horizontal") != 0,
                                 s.num("flip_vertical") != 0, s.num("time_reference") == 0));
    Catalog cat = load_vizier_tsv(ex + "gaia_dr3.tsv", 2000.0 + (ims[0].jd_mid - 2451545.0) / 365.25, s.num("catalog_faint_mag"));
    Reduction R = reduce(ims, cat, RA0, DE0, s, true);
    char b[256];
    for (auto &l : R.log)
        if (l.rfind("trail", 0) == 0)
            std::printf("     %s\n", l.c_str());
    for (size_t k = 0; k < R.frames.size(); ++k) {
        const Plate &p = R.frames[k].plate;
        std::snprintf(b, sizeof b,
                      "image %zu: %d refs (of %d), order %d, rms RA %.3f\" Dec %.3f\" (out-of-sample %.3f\"/%.3f\"), ZP %.3f ± %.3f, "
                      "seeing %.1f\"",
                      k + 1, p.nref, p.nmatch, p.order, p.rms_ra, p.rms_dec, p.rms_cv_ra, p.rms_cv_dec, p.zp, p.zp_rms, R.frames[k].seeing);
        check(R.frames[k].solved && p.nref >= 50 && p.rms_ra < PRMS && p.rms_dec < PRMS, b);
    }
    std::snprintf(b, sizeof b, "one moving object (%zu)", R.movers.size());
    check(R.movers.size() == 1, b);
    if (R.movers.empty())
        return 1;
    const Track &T = R.movers[0];
    std::snprintf(b, sizeof b, "velocidade %.3f\"/min (Horizons %.3f ± 0.02), AP %.1f° (Horizons %.1f°)", T.rate, HZ_RATE, T.pa, HZ_PA);
    check(std::fabs(T.rate - HZ_RATE) <= 0.02 && std::fabs(T.pa - HZ_PA) <= 3, b);
    std::snprintf(b, sizeof b, "measurements in %zu images", T.meas.size());
    check(T.meas.size() == 3, b);
    double worst = 0, ss = 0;
    for (auto &m : T.meas) {
        if (m.frame < 0 || m.frame >= int(hz.size()))
            continue;
        double dra = (m.ra - hz[m.frame].first) * 3600 * std::cos(m.dec * PI / 180), dde = (m.dec - hz[m.frame].second) * 3600;
        double d = std::hypot(dra, dde);
        worst = std::max(worst, d);
        ss += d * d;
        std::snprintf(b, sizeof b, "image %d: O−C (Horizons) dRA* %+.3f\" dDec %+.3f\" |%.3f\"|, mag %.2f, SNR %.1f", m.frame + 1, dra, dde,
                      d, m.mag, m.snr_tot);
        check(d <= DMAX && m.mag >= 16 && m.mag <= 18.5, b);
    }
    double rms = std::sqrt(ss / T.meas.size()), mlo = 99, mhi = -99;
    for (auto &m : T.meas) {
        mlo = std::min(mlo, m.mag);
        mhi = std::max(mhi, m.mag);
    }
    std::snprintf(b, sizeof b, "max residual %.3f\" (≤ %.2f), rms %.3f\" (≤ %.2f)", worst, WMAX, rms, RMSMAX);
    check(worst <= WMAX && rms <= RMSMAX, b);
    // magnitudes consistent across frames (filters r, i, g, all calibrated to Gaia RP)
    std::snprintf(b, sizeof b, "magnitude amplitude %.2f (≤ %.2f)", mhi - mlo, AMP);
    check(mhi - mlo <= AMP, b);
    // automatic plate order (plate_order = 0, k-fold cross-validation)
    {
        Settings s0 = s;
        s0.set("plate_order", "0");
        Reduction R0 = reduce(ims, cat, RA0, DE0, s0, true);
        double w0 = 0, q0 = 0;
        int n0 = 0;
        std::string ords;
        for (auto &f : R0.frames)
            ords += std::to_string(f.plate.order);
        for (auto &m : R0.movers.empty() ? std::vector<Measurement>{} : R0.movers[0].meas) {
            double dra = (m.ra - hz[m.frame].first) * 3600 * std::cos(m.dec * PI / 180), dde = (m.dec - hz[m.frame].second) * 3600;
            w0 = std::max(w0, std::hypot(dra, dde));
            q0 += dra * dra + dde * dde;
            ++n0;
        }
        std::snprintf(b, sizeof b, "automatic plate_order: orders %s, max residual %.3f\", rms %.3f\"", ords.c_str(), w0,
                      n0 ? std::sqrt(q0 / n0) : NAN);
        check(n0 == 3 && w0 <= WMAX, b);
    }
    // identification: offline orbits (versioned extract) at the position measured in the 1st image, at Palomar (I41)
    AstroSession *a = astro_new();
    int n = astro_load_orbits(a, (root + "cpp/tests/data/orbits_mpcorb.dat").c_str(), (root + "cpp/tests/data/orbits_comets.dat").c_str());
    astro_set(a, "longitude", "116.8598");
    astro_set(a, "latitude", "33.35728");
    astro_set(a, "altitude_m", "1665.5");
    std::string js = n > 0 ? astro_identify_offline(a, T.meas[0].jd, T.meas[0].ra, T.meas[0].dec, 0.05) : "";
    size_t p = js.find("\"packed\":\"05157\"");
    double sep = p == std::string::npos ? NAN : std::atof(js.c_str() + js.find("\"sep_arcsec\":", p) + 13);
    std::snprintf(b, sizeof b, "identified as (5157) Hindemith at %.1f\" from the offline orbit", sep);
    check(p != std::string::npos && sep < 2, b);
    astro_free(a);
    std::printf(fails ? "%d failures\n" : "OK\n", fails);
    return fails ? 1 : 0;
}
