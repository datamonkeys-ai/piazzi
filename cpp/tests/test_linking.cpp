// Synthetic tests of detection and moving-object search (full reduction on simulated images):
//   - faint source detected with the matched filter (matched_filter=1) and not with the 3×3 binomial;
//   - normal object (0.5″/min, 3 frames) found at the correct position;
//   - slow object (0.25″/min) with two nearly simultaneous frames at the end of the arc is not taken for a star;
//   - object missing from one of 4 frames: track with 3 detections;
//   - field with many random transient sources (no real object): no chance track;
//   - short arc (4.8 min) with sharp images (FWHM 1.5″): slow object (0.5″/min, 2.4″ over the arc ≈ 1.6 FWHM) found;
//   - same short arc, stars with frame-to-frame centroid jitter ~ σ: no track;
//   - ADES with a unique trkSub per track and no provID when the object was not identified.
// Uso: test_linking
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace astro;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

static const double RA0 = 150.0, DEC0 = 10.0, SKY = 300, NOISE = 8;
static double FWHM = 2.6, JITTER = 0; // FWHM (px = ″) and per-axis Gaussian jitter of the stars in each frame (″)
static const int W = 420, H = 420;

struct Src {
    double xi, eta, mag; // ″ around (RA0, DEC0); magnitude (ZP 25)
    double color = NAN;  // BP − RP (NaN = unknown / solar-colored object)
};
// per-axis chromatic shift (″/mag in BP − RP − 0.82) applied to the stars when rendering (differential chromatic
// refraction / lateral color): the solar-colored target does not shift
static double COLSHIFT_X = 0, COLSHIFT_Y = 0;

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
        const double cc = std::isfinite(s.color) ? s.color - SOLAR_BP_RP : 0;
        double x0 = W / 2.0 + s.xi + COLSHIFT_X * cc, y0 = H / 2.0 + s.eta + COLSHIFT_Y * cc, F = std::pow(10.0, -0.4 * (s.mag - 25));
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

static Scene make_stars(unsigned seed, bool with_color = false) {
    Scene sc;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-195, 195), um(14, 19.5);
    std::uniform_real_distribution<double> uc(0.3, 2.5);
    for (int k = 0; k < 160; ++k) {
        Src s{u(rng), u(rng), um(rng), with_color ? uc(rng) : NAN};
        sc.stars.push_back(s);
        if (s.mag < 18.5) { // stars fainter than this are left out of the catalog (only cross-matching flags them)
            double ra, dec;
            std_to_radec(s.xi, s.eta, RA0, DEC0, ra, dec);
            sc.cat.ra.push_back(ra);
            sc.cat.dec.push_back(dec);
            sc.cat.mag.push_back(s.mag);
            if (with_color)
                sc.cat.color.push_back(s.color);
        }
    }
    return sc;
}

// moving object: position at t (min) from (xi0, eta0) with rate (vx, vy) ″/min
struct Mover {
    double xi0, eta0, vx, vy, mag;
};

static Reduction run(const Scene &sc, const std::vector<double> &tmin, const std::vector<Mover> &mv, const std::vector<int> &skip,
                     const std::vector<std::vector<Src>> &extra, unsigned seed) {
    std::vector<Image> ims;
    for (size_t k = 0; k < tmin.size(); ++k) {
        std::vector<Src> S = sc.stars;
        if (JITTER > 0) {
            std::mt19937 rj(seed * 7 + unsigned(k));
            std::normal_distribution<double> nj(0, JITTER);
            for (auto &q : S) {
                q.xi += nj(rj);
                q.eta += nj(rj);
            }
        }
        for (size_t j = 0; j < mv.size(); ++j)
            if (!(j < skip.size() && skip[j] == int(k)))
                S.push_back({mv[j].xi0 + mv[j].vx * tmin[k], mv[j].eta0 + mv[j].vy * tmin[k], mv[j].mag});
        if (k < extra.size())
            S.insert(S.end(), extra[k].begin(), extra[k].end());
        ims.push_back(render(S, 2460000.5 + tmin[k] / 1440, seed + unsigned(k)));
    }
    return reduce(ims, sc.cat, RA0, DEC0, settings(), true);
}

// does track T pass within < tol ″ of object m in all its measurements?
static bool matches(const Track &T, const Mover &m, const std::vector<double> &tmin, double tol) {
    if (T.meas.size() < 3)
        return false;
    for (auto &me : T.meas) {
        double t = (me.jd - 2460000.5) * 1440, xi, eta;
        radec_to_std(me.ra, me.dec, RA0, DEC0, xi, eta);
        if (std::hypot(xi - (m.xi0 + m.vx * t), eta - (m.eta0 + m.vy * t)) > tol)
            return false;
    }
    (void)tmin;
    return true;
}

int main() {
    char b[256];
    // ---- 1. faint source: matched filter × binomial
    {
        // FWHM 5 px (poor seeing): the 3×3 kernel loses SNR; the PSF-matched filter recovers the source
        std::vector<Src> S;
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> u(-190, 190);
        for (int k = 0; k < 40; ++k)
            S.push_back({u(rng), u(rng), 15.0 + 0.05 * k});
        Image im;
        {
            // render with FWHM 5 by swapping the constant locally
            im.w = W;
            im.h = H;
            std::vector<double> f(size_t(W) * H, SKY);
            const double sg = 5.0 / 2.3548;
            S.push_back({20.5, -30.5, 18.7}); // faint target: F ≈ 330 ADU, optimal SNR ≈ 5.5
            for (auto &s : S) {
                double x0 = W / 2.0 + s.xi, y0 = H / 2.0 + s.eta, F = std::pow(10.0, -0.4 * (s.mag - 25));
                for (int y = int(y0) - 14; y <= int(y0) + 14; ++y)
                    for (int x = int(x0) - 14; x <= int(x0) + 14; ++x)
                        if (x >= 0 && y >= 0 && x < W && y < H)
                            f[size_t(y) * W + x] += F * std::exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (2 * sg * sg)) /
                                                    (2 * 3.141592653589793 * sg * sg);
            }
            std::mt19937 r2(11);
            std::normal_distribution<double> nd(0, NOISE);
            im.pix.resize(f.size());
            for (size_t i = 0; i < f.size(); ++i)
                im.pix[i] = uint16_t(std::clamp(std::lround(f[i] + nd(r2)), 0L, 65535L));
        }
        auto found = [&](const std::vector<Detection> &D) {
            for (auto &d : D)
                if (std::hypot(d.x - (W / 2.0 + 20.5), d.y - (H / 2.0 - 30.5)) < 2)
                    return true;
            return false;
        };
        DetectParams p;
        p.matched = false;
        bool bin = found(detect_sources(im, p));
        p.matched = true;
        bool mf = found(detect_sources(im, p));
        std::snprintf(b, sizeof b, "faint source (optimal SNR ≈ 5.5, FWHM 5 px): binomial %s, matched filter %s", bin ? "yes" : "no",
                      mf ? "yes" : "no");
        check(mf && !bin, b);
    }
    Scene sc = make_stars(42);
    // ---- 2. normal object, 3 frames
    {
        std::vector<double> t{0, 20, 40};
        Mover m{-60, 40, 0.5, -0.3, 17.5};
        Reduction R = run(sc, t, {m}, {}, {}, 100);
        int ok = 0;
        for (auto &T : R.movers)
            ok += matches(T, m, t, 1.0);
        std::snprintf(b, sizeof b, "object 0.58″/min in 3 frames: %zu track(s), %d correct", R.movers.size(), ok);
        check(ok == 1 && R.movers.size() == 1, b);
        std::string ades = ades_report(R, R.movers, settings());
        check(ades.find("|P") != std::string::npos && ades.find("Piazzi     |") == std::string::npos,
              "ADES: trkSub filled and provID empty (object not identified)");
    }
    // ---- 3. slow object with a nearly simultaneous pair at the end of the arc
    {
        std::vector<double> t{0, 30, 57, 60};
        Mover m{50, -70, 0.25, 0.0, 17.5};
        Reduction R = run(sc, t, {m}, {}, {}, 200);
        int ok = 0;
        for (auto &T : R.movers)
            ok += matches(T, m, t, 1.0);
        std::snprintf(b, sizeof b, "slow object 0.25″/min (3 min between the last 2 frames): %d correct track(s)", ok);
        check(ok == 1, b);
    }
    // ---- 4. missing point: object absent from the 2nd of 4 frames
    {
        std::vector<double> t{0, 15, 30, 45};
        Mover m{-100, -120, -0.4, 0.6, 17.0};
        Reduction R = run(sc, t, {m}, {1}, {}, 300);
        int ok = 0;
        size_t nm = 0;
        for (auto &T : R.movers)
            if (matches(T, m, t, 1.0)) {
                ++ok;
                nm = T.meas.size();
            }
        std::snprintf(b, sizeof b, "object missing from 1 of 4 frames: %d correct track(s) with %zu measurements", ok, nm);
        check(ok == 1 && nm == 3, b);
    }
    // ---- 5. random transients (single-frame sources, high density): no chance track
    {
        std::vector<double> t{0, 20, 40};
        std::vector<std::vector<Src>> extra(3);
        std::mt19937 rng(99);
        std::uniform_real_distribution<double> u(-190, 190), um(17.5, 19.0);
        for (int k = 0; k < 3; ++k)
            for (int j = 0; j < 150; ++j)
                extra[k].push_back({u(rng), u(rng), um(rng)});
        Reduction R = run(sc, t, {}, {}, extra, 400);
        std::snprintf(b, sizeof b, "150 random transients per frame, no real object: %zu chance track(s)", R.movers.size());
        check(R.movers.empty(), b);
    }
    // ---- 6. short arc with sharp images (SDSS-like: 4 frames in 4.8 min, FWHM 1.5″): the minimum displacement comes
    // from each centroid's uncertainty (≈ 0.01–0.05″), not from a fixed floor in ″
    FWHM = 1.5;
    {
        std::vector<double> t{0, 1.2, 3.6, 4.8};
        Mover m{30, 60, 0.4, -0.3, 16.5};
        Reduction R = run(sc, t, {m}, {}, {}, 500);
        int ok = 0;
        for (auto &T : R.movers)
            ok += matches(T, m, t, 0.5);
        std::snprintf(b, sizeof b, "short arc 4.8 min, FWHM 1.5″, object 0.5″/min (2.4″ over the arc): %zu track(s), %d correct",
                      R.movers.size(), ok);
        check(ok == 1 && R.movers.size() == 1, b);
    }
    // ---- 7. same short arc, no object, stars with 0.05″ per-axis centroid jitter (chromatic refraction,
    // turbulence) added to the noise: no track
    {
        std::vector<double> t{0, 1.2, 3.6, 4.8};
        int ntr = 0;
        for (unsigned rep = 0; rep < 3; ++rep) {
            JITTER = 0.05;
            Scene s2 = make_stars(1000 + rep);
            Reduction R = run(s2, t, {}, {}, {}, 600 + 10 * rep);
            ntr += int(R.movers.size());
        }
        JITTER = 0;
        std::snprintf(b, sizeof b, "short arc, stars with 0.05″ jitter (3 fields): %d chance track(s)", ntr);
        check(ntr == 0, b);
    }
    // ---- 8. plate color term (Pier et al. 2003): stars shifted by 0.03″/mag in x and −0.02″/mag in y according to
    // color; with plate_color the fit recovers the coefficients and the target (solar color) is unbiased; without the
    // term, the target inherits the shift of the mean reference color (≈ 1.4 − 0.82 = 0.6 mag → ≈ 0.02″ in x)
    {
        std::vector<double> t{0, 1.2, 3.6, 4.8};
        Mover m{-40, 20, 0.5, 0.2, 15.5};
        const Scene scc = make_stars(42, true);
        double bias[2][2] = {};
        double col[2] = {};
        for (int mode = 0; mode < 2; ++mode) {
            COLSHIFT_X = 0.03;
            COLSHIFT_Y = -0.02;
            std::vector<Image> ims;
            for (size_t k = 0; k < t.size(); ++k) {
                std::vector<Src> S = scc.stars;
                S.push_back({m.xi0 + m.vx * t[k], m.eta0 + m.vy * t[k], m.mag});
                ims.push_back(render(S, 2460000.5 + t[k] / 1440, 700 + unsigned(k)));
            }
            COLSHIFT_X = COLSHIFT_Y = 0;
            Settings st = settings();
            st.set("plate_color", mode ? "1" : "0");
            Reduction R = reduce(ims, scc.cat, RA0, DEC0, st, true);
            int n = 0;
            for (auto &T : R.movers)
                if (matches(T, m, t, 0.5))
                    for (auto &me : T.meas) {
                        double tt = (me.jd - 2460000.5) * 1440, xi, eta;
                        radec_to_std(me.ra, me.dec, RA0, DEC0, xi, eta);
                        bias[mode][0] += xi - (m.xi0 + m.vx * tt);
                        bias[mode][1] += eta - (m.eta0 + m.vy * tt);
                        ++n;
                    }
            bias[mode][0] /= std::max(n, 1);
            bias[mode][1] /= std::max(n, 1);
            if (mode) {
                for (auto &f : R.frames)
                    if (f.solved) {
                        col[0] += f.plate.col_ra / 4;
                        col[1] += f.plate.col_dec / 4;
                    }
            }
        }
        // convention: ξ_cat = P(x) + col·(c − c☉), so a star shifted by +Δ on the detector gives col = −Δ
        std::snprintf(b, sizeof b,
                      "color term: coefficients %+.4f/%+.4f ″/mag (expected −0.0300/+0.0200), target bias without term %+.4f/%+.4f″, with "
                      "term %+.4f/%+.4f″",
                      col[0], col[1], bias[0][0], bias[0][1], bias[1][0], bias[1][1]);
        check(std::fabs(col[0] + 0.03) < 0.005 && std::fabs(col[1] - 0.02) < 0.005 && std::fabs(bias[1][0]) < 0.006 &&
                  std::fabs(bias[1][1]) < 0.006 && std::fabs(bias[0][0]) > 2 * std::fabs(bias[1][0]),
              b);
    }
    FWHM = 2.6;
    std::printf(fails ? "\n%d failure(s)\n" : "\nOK\n", fails);
    return fails ? 1 : 0;
}
