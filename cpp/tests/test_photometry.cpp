// Tests of color photometry and of the ADES report:
//   - fit_zero_point with a color term: references m_RP = zp + k·(BP − RP − 0.82) − 2.5·log F recover zp (solar color)
//     and k; without colors (or with few references) k = 0 and zp = median;
//   - Plate::star_mag returns the predicted in-frame magnitude of a star of known color;
//   - load_vizier_tsv reads BP − RP when BPmag and RPmag are present;
//   - flux_coherent: a constant-flux track passes; a "red" noise track (bright detection in a low-sensitivity frame)
//     is rejected; without color terms the tolerance is larger;
//   - ADES: obsTime with milliseconds and positions with 6 decimals when σ < 0.1″, rmsRA/rmsDec with 3 decimals.
// Uso: test_photometry
#include "../src/reduce_internal.hpp"
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace astro;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

// synthetic frame with n references of uniform color in [0.4, 2.4] and flux given by the transformation equation
static Frame color_frame(int n, double zp, double k, bool with_color, Catalog &cat, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uc(0.4, 2.4), um(14, 19);
    std::normal_distribution<double> noise(0, 0.01);
    Frame f;
    cat = Catalog{};
    for (int i = 0; i < n; ++i) {
        double c = uc(rng), m = um(rng);
        double inst = m - zp - k * (c - SOLAR_BP_RP) + noise(rng); // −2,5·log F
        Detection d;
        d.flux = std::pow(10.0, -0.4 * inst);
        d.peak = 1000;
        f.det.push_back(d);
        cat.ra.push_back(0);
        cat.dec.push_back(0);
        cat.mag.push_back(m);
        if (with_color)
            cat.color.push_back(c);
        f.plate.ref_idx.push_back(i);
        f.plate.cat_idx.push_back(i);
        f.plate.ref_ok.push_back(1);
    }
    return f;
}

int main() {
    // ---- zero-point color term
    {
        Catalog cat;
        Frame f = color_frame(80, 28.0, 0.6, true, cat, 1);
        internal::fit_zero_point(f, cat, 60000);
        check(std::fabs(f.plate.zp - 28.0) < 0.01, "zero point referred to solar color recovered (k = 0.6)");
        check(std::fabs(f.plate.zp_color - 0.6) < 0.02, "color term recovered");
        check(std::fabs(f.plate.star_mag(17.0, 1.82) - (17.0 - 0.6)) < 0.03, "star_mag: star 1 mag redder than the Sun");
        check(f.plate.star_mag(17.0, NAN) == 17.0, "star_mag: unknown color returns the catalog magnitude");
    }
    {
        Catalog cat;
        Frame f = color_frame(80, 28.0, 0.6, false, cat, 2);
        internal::fit_zero_point(f, cat, 60000);
        check(f.plate.zp_color == 0, "no colors in catalog: no color term");
        check(std::fabs(f.plate.zp - 28.0) > 0.05, "no colors: zp = median (shifted by the mean reference color)");
    }
    {
        Catalog cat;
        Frame f = color_frame(6, 28.0, 0.6, true, cat, 3);
        internal::fit_zero_point(f, cat, 60000);
        check(f.plate.zp_color == 0, "few references with color (< 10): no color term");
    }
    // ---- catalog: BP − RP color
    {
        const char *path = "test_photometry_cat.tsv";
        std::ofstream o(path);
        o << "RAJ2000\tDEJ2000\tpmRA\tpmDE\tGmag\tBPmag\tRPmag\ndeg\tdeg\tmas/yr\tmas/"
             "yr\tmag\tmag\tmag\n------\t------\t------\t------\t------\t------\t------\n";
        o << "10.0\t+5.0\t0\t0\t15.0\t15.5\t14.5\n10.01\t+5.01\t0\t0\t16.0\t\t15.6\n";
        o.close();
        Catalog c = load_vizier_tsv(path, 2000.0, 20);
        std::remove(path);
        check(c.has_color() && std::fabs(c.color[0] - 1.0) < 1e-9, "load_vizier_tsv: BP − RP lido");
        check(c.size() == 2 && std::isnan(c.color[1]), "load_vizier_tsv: missing BPmag -> NaN color");
    }
    // ---- fluxo coerente
    {
        check(internal::flux_coherent({20.0, 20.1, 19.9, 20.05}, {0.1, 0.1, 0.2, 0.1}, true), "fluxo constante passa");
        // typical SDSS noise track (r, i, z, g): the least sensitive frame (z) has the brightest detection
        check(!internal::flux_coherent({21.8, 20.9, 20.5, 22.1}, {0.15, 0.17, 0.24, 0.16}, true), "'red' noise track rejected");
        check(internal::flux_coherent({18.0, 18.5, 18.9}, {0.03, 0.03, 0.03}, false) &&
                  !internal::flux_coherent({18.0, 18.5, 18.9}, {0.03, 0.03, 0.03}, true),
              "without color terms the intrinsic tolerance is larger");
        check(internal::flux_coherent({19.0}, {0.1}, true), "single measurement: nothing to test");
    }
    // ---- ADES: precision of time, positions and σ
    {
        Reduction R;
        Frame f;
        f.solved = true;
        f.plate.ra0 = 150;
        f.plate.dec0 = 10;
        f.plate.cx = {0, 0.4, 0};
        f.plate.cy = {0, 0, 0.4};
        f.plate.rms_ra = f.plate.rms_dec = 0.02;
        f.im.jd_mid = 2452910.5 + 0.123456789;
        f.im.exptime = 54;
        R.frames.push_back(f);
        Track T;
        T.idx = {-1};
        Measurement m{};
        m.frame = 0;
        m.ra = 150.0000012;
        m.dec = 10.0000034;
        m.snr = m.snr_tot = 300;
        m.fwhm_as = 1.2;
        m.jd = f.im.jd_mid;
        m.mag = 17;
        T.meas.push_back(m);
        T.desig = "Piazzi";
        Settings s;
        s.set("mpc_code", "645");
        std::string a = ades_report(R, {T}, s);
        check(a.find("2003-09-28T02:57:46.667Z") != std::string::npos, "ADES obsTime with milliseconds");
        check(a.find("150.000001") != std::string::npos && a.find("+10.000003") != std::string::npos,
              "ADES positions with 6 decimals when σ < 0.1″");
        check(a.find("|0.020 |") != std::string::npos || a.find("|0.021 |") != std::string::npos, "ADES rmsRA with 3 decimals");
    }
    std::printf("%s\n", fails ? "FAILED" : "all ok");
    return fails ? 1 : 0;
}
