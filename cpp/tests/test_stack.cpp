// Synthetic tests of the shift-and-stack search and of 2-detection tracks with forced confirmation:
//   - object below the single-frame threshold (SNR ≈ 3 per frame in the matched filter, 5 frames) recovered by
//     stacking and not by detection linking;
//   - no false alarms in noise + stars only images (3 and 5 frames, several seeds);
//   - object detected in 2 of 3 frames (the middle one 1.8 mag fainter, below the detection threshold): track
//     confirmed by forced photometry in the middle frame.
// Uso: test_stack
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

static const double RA0 = 150.0, DEC0 = 10.0, SKY = 300, NOISE = 8, FWHM = 2.6;
static const int W = 420, H = 420;

struct Src {
    double xi, eta, mag; // ″ around (RA0, DEC0); magnitude (ZP 25)
};

// simulated image: sky + Gaussian noise + Gaussian sources (integrated over the pixel by 3×3 subsampling)
static Image render(const std::vector<Src> &S, double jd, unsigned seed) {
    Image im;
    im.w = W;
    im.h = H;
    im.jd_mid = jd;
    im.exptime = 30;
    std::vector<double> f(size_t(W) * H, SKY);
    const double sg = FWHM / 2.3548;
    for (auto &s : S) {
        double x0 = W / 2.0 + s.xi, y0 = H / 2.0 + s.eta, F = std::pow(10.0, -0.4 * (s.mag - 25));
        int R = int(5 * sg) + 2;
        for (int y = int(y0) - R; y <= int(y0) + R; ++y)
            for (int x = int(x0) - R; x <= int(x0) + R; ++x) {
                if (x < 0 || y < 0 || x >= W || y >= H)
                    continue;
                double a = 0;
                for (int j = 0; j < 3; ++j)
                    for (int i = 0; i < 3; ++i) {
                        double dx = x - 1.0 / 3 + i / 3.0 - x0, dy = y - 1.0 / 3 + j / 3.0 - y0;
                        a += std::exp(-(dx * dx + dy * dy) / (2 * sg * sg));
                    }
                f[size_t(y) * W + x] += F * a / 9 / (2 * 3.141592653589793 * sg * sg);
            }
    }
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0, NOISE);
    im.pix.resize(f.size());
    for (size_t i = 0; i < f.size(); ++i)
        im.pix[i] = uint16_t(std::clamp(std::lround(f[i] + nd(rng)), 0L, 65535L));
    return im;
}

static Settings settings() {
    Settings s = Settings::defaults();
    // 1″/px: 0,015 mm / 3093,7 mm · 206264,8
    s.set("focal_length_mm", "3093.7");
    s.set("pixel_width_mm", "0.015");
    s.set("pixel_height_mm", "0.015");
    s.set("flip_horizontal", "0");
    s.set("catalog_faint_mag", "19");
    s.set("saturation", "60000");
    return s;
}

struct Scene {
    std::vector<Src> stars;
    Catalog cat;
};

static Scene make_stars(unsigned seed) {
    Scene sc;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-195, 195), um(14, 19.5);
    for (int k = 0; k < 160; ++k) {
        Src s{u(rng), u(rng), um(rng)};
        sc.stars.push_back(s);
        if (s.mag < 18.5) { // stars fainter than this are left out of the catalog (only cross-matching flags them)
            double ra, dec;
            std_to_radec(s.xi, s.eta, RA0, DEC0, ra, dec);
            sc.cat.ra.push_back(ra);
            sc.cat.dec.push_back(dec);
            sc.cat.mag.push_back(s.mag);
        }
    }
    return sc;
}

struct Mover {
    double xi0, eta0, vx, vy;
    std::vector<double> mag; // per frame
};

static Reduction run(const Scene &sc, const std::vector<double> &tmin, const std::vector<Mover> &mv, unsigned seed, Settings s) {
    std::vector<Image> ims;
    for (size_t k = 0; k < tmin.size(); ++k) {
        std::vector<Src> S = sc.stars;
        for (auto &m : mv)
            S.push_back({m.xi0 + m.vx * tmin[k], m.eta0 + m.vy * tmin[k], m.mag[k]});
        ims.push_back(render(S, 2460000.5 + tmin[k] / 1440, seed + unsigned(k)));
    }
    Reduction R = reduce(ims, sc.cat, RA0, DEC0, s, true);
    if (std::getenv("STACK_VERBOSE"))
        for (auto &l : R.log)
            if (l.find("stack") != std::string::npos || l.find("detection") != std::string::npos)
                std::printf("    %s\n", l.c_str());
    return R;
}

// track T with ≥ nmin measurements, all within < tol ″ of object m
static bool matches(const Track &T, const Mover &m, double tol, size_t nmin) {
    if (T.meas.size() < nmin)
        return false;
    for (auto &me : T.meas) {
        double t = (me.jd - 2460000.5) * 1440, xi, eta;
        radec_to_std(me.ra, me.dec, RA0, DEC0, xi, eta);
        if (std::hypot(xi - (m.xi0 + m.vx * t), eta - (m.eta0 + m.vy * t)) > tol)
            return false;
    }
    return true;
}

int main() {
    char b[256];
    Scene sc = make_stars(42);
    // matched-filter SNR of a point source: F / (σ_sky·√(4π)·σ_PSF) (σ_PSF = FWHM/2.355)
    const double snr1 = std::pow(10.0, -0.4 * (19.8 - 25)) / (NOISE * std::sqrt(4 * 3.141592653589793) * FWHM / 2.3548);
    // ---- 1. faint object, 5 frames, SNR ≈ 3 per frame
    {
        std::vector<double> t{0, 15, 30, 45, 60};
        Mover m{-70, 50, 0.45, -0.35, std::vector<double>(5, 19.8)};
        Settings s = settings();
        s.set("stack_search", "0");
        Reduction R0 = run(sc, t, {m}, 500, s);
        int ok0 = 0;
        for (auto &T : R0.movers)
            ok0 += matches(T, m, 1.5, 2);
        s.set("stack_search", "1");
        {
            double ra, dec;
            std_to_radec(m.xi0, m.eta0, RA0, DEC0, ra, dec);
            std::snprintf(b, sizeof b, "%.8f %.8f 2460000.5 %.4f %.4f", ra, dec, m.vx * 60, m.vy * 60);
            s.set("stack_probe", b);
        }
        Reduction R1 = run(sc, t, {m}, 500, s);
        int ok1 = 0;
        size_t nm = 0;
        for (auto &T : R1.movers)
            if (matches(T, m, 1.5, 3)) {
                ++ok1;
                nm = T.meas.size();
            }
        std::snprintf(b, sizeof b,
                      "object SNR %.1f/frame in 5 frames: linking %d, stacking %d track(s) (%zu measurements), %zu track(s) in total", snr1,
                      ok0, ok1, nm, R1.movers.size());
        check(ok0 == 0 && ok1 == 1 && R1.movers.size() == 1, b);
        std::string ades = ades_report(R1, R1.movers, settings());
        check(ades.find("synthetic tracking") != std::string::npos, "ADES: remarks flags the forced stacking measurements");
    }
    // ---- 2. noise + stars only: no false alarms
    {
        int nfa = 0, nrun = 0;
        for (unsigned seed : {600u, 700u, 800u}) {
            for (int nf : {3, 5}) {
                std::vector<double> t;
                for (int k = 0; k < nf; ++k)
                    t.push_back(60.0 * k / (nf - 1));
                Settings s = settings();
                s.set("stack_search", "1");
                Reduction R = run(sc, t, {}, seed, s);
                nfa += int(R.movers.size());
                ++nrun;
            }
        }
        std::snprintf(b, sizeof b, "noise + stars (%d sequences of 3 and 5 frames): %d false track(s)", nrun, nfa);
        check(nfa == 0, b);
    }
    // ---- 3. 2 detections + forced confirmation: middle frame 1.8 mag fainter (forced SNR ≈ 3)
    {
        std::vector<double> t{0, 30, 60};
        Mover m{40, -60, -0.5, 0.4, {17.8, 19.6, 17.8}};
        Settings s = settings();
        s.set("stack_search", "0");
        s.set("pair_tracks", "0");
        Reduction R0 = run(sc, t, {m}, 900, s);
        int ok0 = 0;
        for (auto &T : R0.movers)
            ok0 += matches(T, m, 1.5, 2);
        s.set("pair_tracks", "1");
        if (std::getenv("STACK_VERBOSE"))
            s.set("stack_debug", "1");
        Reduction R1 = run(sc, t, {m}, 900, s);
        int ok1 = 0;
        size_t nm = 0;
        for (auto &T : R1.movers)
            if (matches(T, m, 1.5, 2)) {
                ++ok1;
                nm = T.meas.size();
            }
        std::snprintf(b, sizeof b, "2 detections + forced: without pairs %d, with pairs %d track(s) (%zu measurements), %zu in total", ok0,
                      ok1, nm, R1.movers.size());
        check(ok0 == 0 && ok1 == 1 && R1.movers.size() == 1, b);
    }
    std::printf(fails ? "\n%d failure(s)\n" : "\nOK\n", fails);
    return fails ? 1 : 0;
}
