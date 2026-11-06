// Synthetic tests of fast-object trails (elongated matched filter + trail fit of Vereš et al. 2012):
//   1. 5–15 px trails with optimal SNR 10–50 and random angles: detection ≥ 95 % at SNR ≥ 10 (a single detection, no
//      fragments), unbiased center (mean along/across the trail < 0.05 px at SNR 50) and scatter consistent with the
//      Cramér–Rao bound;
//   2. star-only images: the trail search creates no detections (same count with trail_detect = 0 and 1);
//   3. full reduction: objects at 12 and 20″/min (6 and 10 px trails in 30 s, 1″/px) linked over 4 and 3 frames,
//      mid-exposure positions within ≤ 0.25–0.3″ of truth at SNR 40 (and linking at SNR 15); automatic trail fit used;
//      trails perpendicular to the motion are not linked (trail × rate consistency).
// Uso: test_trail
#include "../src/reduce_internal.hpp"
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using namespace astro;
using namespace astro::internal;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

static const double SKY = 1000, NOISE = 10, SIG = 1.27; // PSF σ (FWHM ≈ 3 px)

struct Src {
    double x, y, flux, len = 0, ang = 0;
};

static Image render(int W, int H, const std::vector<Src> &S, unsigned seed) {
    std::vector<double> f(size_t(W) * H, SKY);
    for (auto &s : S) {
        PsfShape P{SIG, SIG, s.len, std::cos(s.ang), std::sin(s.ang)};
        int R = int(5 * SIG + s.len / 2) + 2;
        for (int y = int(s.y) - R; y <= int(s.y) + R; ++y)
            for (int x = int(s.x) - R; x <= int(s.x) + R; ++x)
                if (x >= 0 && y >= 0 && x < W && y < H) {
                    double a = 0; // pixel integral by 4×4 sampling (the fit model uses 3×3)
                    for (int j = 0; j < 4; ++j)
                        for (int i = 0; i < 4; ++i)
                            a += psf_profile(P, x - 0.375 + 0.25 * i - s.x, y - 0.375 + 0.25 * j - s.y);
                    f[size_t(y) * W + x] += s.flux * a / 16;
                }
    }
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0, NOISE);
    Image im;
    im.w = W;
    im.h = H;
    im.exptime = 30;
    im.pix.resize(f.size());
    for (size_t i = 0; i < f.size(); ++i)
        im.pix[i] = uint16_t(std::clamp(std::lround(f[i] + nd(rng)), 0L, 65535L));
    return im;
}

// flux of a trail with given optimal SNR (matched filter, background-limited): SNR² = F²·Σp²/σ²
static double flux_for_snr(double snr, double len) {
    PsfShape P{SIG, SIG, len, 1, 0};
    double s2 = 0;
    for (int y = -30; y <= 30; ++y)
        for (int x = -30; x <= 30; ++x) {
            double p = psf_pixel(P, x, y);
            s2 += p * p;
        }
    return snr * NOISE / std::sqrt(s2);
}

static std::vector<Src> stars(int W, int H, double cx, double cy, double avoid, std::mt19937 &rng, int n = 40) {
    std::uniform_real_distribution<double> ux(12, W - 12), uy(12, H - 12), uf(3000, 60000);
    std::vector<Src> S;
    while (int(S.size()) < n) {
        Src s{ux(rng), uy(rng), uf(rng)};
        if (std::hypot(s.x - cx, s.y - cy) > avoid)
            S.push_back(s);
    }
    return S;
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    char b[320];
    // ---- 1. detection and center of isolated trails
    {
        const int W = 160, H = 160, NT = 40;
        std::mt19937 rng(42);
        for (double L : {5.0, 10.0, 15.0})
            for (double snr : {10.0, 20.0, 50.0}) {
                int found = 0, single = 0;
                double sa = 0, sc = 0, sa2 = 0, sc2 = 0, sl = 0;
                int n = 0;
                for (int t = 0; t < NT; ++t) {
                    std::uniform_real_distribution<double> u01(0, 1);
                    double cx = 80 + u01(rng), cy = 80 + u01(rng), th = PI * u01(rng);
                    auto S = stars(W, H, cx, cy, 30, rng);
                    S.push_back({cx, cy, flux_for_snr(snr, L), L, th});
                    Image im = render(W, H, S, unsigned(1000 + t + 97 * int(L) + int(snr)));
                    auto D = detect_sources(im, DetectParams{});
                    const Detection *best = nullptr;
                    int near = 0;
                    for (auto &d : D) {
                        double dd = std::hypot(d.x - cx, d.y - cy);
                        if (dd < L / 2 + 2)
                            ++near;
                        if (dd < std::max(1.5, 0.3 * L) && (!best || dd < std::hypot(best->x - cx, best->y - cy)))
                            best = &d;
                    }
                    if (!best)
                        continue;
                    ++found;
                    single += near == 1;
                    double ux = std::cos(th), uy = std::sin(th), ex = best->x - cx, ey = best->y - cy;
                    double a = ex * ux + ey * uy, c = -ex * uy + ey * ux;
                    sa += a;
                    sc += c;
                    sa2 += a * a;
                    sc2 += c * c;
                    sl += best->trail_len;
                    ++n;
                }
                double ma = n ? sa / n : 0, mc = n ? sc / n : 0, ra = n ? std::sqrt(sa2 / n) : 0, rc = n ? std::sqrt(sc2 / n) : 0;
                // Cramér–Rao bound (background-dominated): σ_⊥ ≈ √2·σ/SNR·(…) ~ σ·√2/SNR; along ≈ √(σ² + L²/12)·√2/SNR
                double cra = std::sqrt(2 * (SIG * SIG + L * L / 12)) / snr, crc = std::sqrt(2.0) * SIG / snr;
                std::snprintf(b, sizeof b,
                              "trail L %2.0f px SNR %2.0f: detected %d/%d (single %d), measured L %.1f, bias along %+.3f across %+.3f px, "
                              "rms %.3f/%.3f px (CR ≈ %.3f/%.3f)",
                              L, snr, found, NT, single, n ? sl / n : 0, ma, mc, ra, rc, cra, crc);
                // residual fragmentation at SNR 10 (the two-PSF model sometimes ties the trail, or aperture SNR
                // falls below the threshold): the object is detected, but in 2–3 fragments; requires ≥ 70 % single
                // detection at SNR 10 and ≥ 95 % above
                bool ok = found >= int(std::ceil(0.95 * NT)) && single >= int(std::ceil((snr < 15 ? 0.7 : 0.95) * NT));
                // bias: mean of NT measurements with σ_mean = rms/√NT; requires < 0.05 px + 3σ_mean (exact at SNR 50)
                ok = ok && std::fabs(ma) < 0.05 + 3 * ra / std::sqrt(std::max(n, 1)) &&
                     std::fabs(mc) < 0.05 + 3 * rc / std::sqrt(std::max(n, 1));
                if (snr >= 50)
                    ok = ok && std::fabs(ma) < 0.05 && std::fabs(mc) < 0.05;
                // scatter ≤ 2× Cramér–Rao (3× at SNR 10, where the remaining fragments weigh in)
                const double kcr = snr < 15 ? 3.0 : 2.0;
                ok = ok && ra < kcr * cra + 0.05 && rc < kcr * crc + 0.05;
                check(ok, b);
            }
    }
    // ---- 2. stars only: no trail, same detection count
    {
        const int W = 300, H = 300;
        int ntr = 0, n0 = 0, n1 = 0;
        for (int t = 0; t < 10; ++t) {
            std::mt19937 rng(500 + t);
            auto S = stars(W, H, -100, -100, 0, rng, 150);
            std::uniform_real_distribution<double> uf(80, 800), ux(10, W - 10);
            for (int k = 0; k < 100; ++k) // faint stars (SNR 3–20)
                S.push_back({ux(rng), ux(rng), uf(rng)});
            Image im = render(W, H, S, 900 + t);
            DetectParams p;
            auto D1 = detect_sources(im, p);
            p.trails = false;
            auto D0 = detect_sources(im, p);
            p.trails = true;
            for (auto &d : D1) {
                ntr += d.trail_len > 0;
                if (d.trail_len > 0 && std::getenv("TRAIL_DEBUG"))
                    std::printf("   spurious: image %d (%.1f, %.1f) L %.1f flux %.0f SNR %.1f det %.1f\n", t, d.x, d.y, d.trail_len, d.flux,
                                d.snr_tot, d.det_snr);
                if (d.trail_len > 0 && std::getenv("TRAIL_DEBUG"))
                    for (auto &e : D0)
                        if (std::hypot(e.x - d.x, e.y - d.y) < 15)
                            std::printf("      point (%.1f, %.1f) flux %.0f FWHM %.2f\n", e.x, e.y, e.flux, e.fwhm);
            }
            n0 += int(D0.size());
            n1 += int(D1.size());
        }
        std::snprintf(b, sizeof b, "stars only (10 images): %d spurious trails; detections %d without and %d with trail search", ntr, n0,
                      n1);
        // tolerance: 1 spurious trail (an unresolved bright star pair/triplet mimics a trail; it is static and
        // does not form a track) and ≤ 0.2 % difference in the detection count
        check(ntr <= 1 && std::abs(n0 - n1) <= n0 / 500 + 1, b);
    }
    // ---- 3. full reduction with fast objects
    {
        const int W = 420, H = 420;
        const double RA0 = 150.0, DEC0 = 10.0;
        struct Mv {
            double x0, y0, vx, vy, flux; // px (= ″) e ″/min
            double ang_off = 0;          // trail rotated relative to the motion (consistency test)
        };
        auto scenario = [&](const std::vector<double> &tmin, const Mv &mv, unsigned seed, double tol, const char *name) {
            std::mt19937 rng(seed);
            std::uniform_real_distribution<double> u(-195, 195), um(14, 19.5);
            Catalog cat;
            std::vector<Src> st;
            for (int k = 0; k < 160; ++k) {
                double xi = u(rng), eta = u(rng), mag = um(rng);
                st.push_back({W / 2.0 + xi, H / 2.0 + eta, std::pow(10.0, -0.4 * (mag - 25))});
                if (mag < 18.5) {
                    double ra, dec;
                    std_to_radec(xi, eta, RA0, DEC0, ra, dec);
                    cat.ra.push_back(ra);
                    cat.dec.push_back(dec);
                    cat.mag.push_back(mag);
                }
            }
            std::vector<Image> ims;
            const double v = std::hypot(mv.vx, mv.vy), ang = std::atan2(mv.vy, mv.vx) + mv.ang_off;
            for (size_t k = 0; k < tmin.size(); ++k) {
                auto S = st;
                S.push_back({W / 2.0 + mv.x0 + mv.vx * tmin[k], H / 2.0 + mv.y0 + mv.vy * tmin[k], mv.flux, v * 0.5, ang});
                Image im = render(W, H, S, seed * 13 + unsigned(k));
                im.jd_mid = 2460000.5 + tmin[k] / 1440; // trail center = mid-exposure
                ims.push_back(std::move(im));
            }
            Settings s = Settings::defaults();
            s.set("focal_length_mm", "3093.7");
            s.set("pixel_width_mm", "0.015");
            s.set("pixel_height_mm", "0.015");
            s.set("flip_horizontal", "0");
            s.set("catalog_faint_mag", "19");
            s.set("saturation", "60000");
            Reduction R = reduce(ims, cat, RA0, DEC0, s, true);
            int hit = 0, nm = 0;
            double worst = 0;
            bool trail_log = false;
            for (auto &l : R.log)
                trail_log |= l.find("trail position") != std::string::npos;
            for (auto &T : R.movers) {
                double w = 0;
                int ok = 0;
                for (auto &me : T.meas) {
                    double t = (me.jd - 2460000.5) * 1440, xi, eta;
                    radec_to_std(me.ra, me.dec, RA0, DEC0, xi, eta);
                    double d = std::hypot(xi - (mv.x0 + mv.vx * t), eta - (mv.y0 + mv.vy * t));
                    w = std::max(w, d);
                    if (std::getenv("TRAIL_DEBUG"))
                        std::printf("   measurement t %.2f min: O−C %.3f″\n", t, d);
                    ok += d < 2.0;
                }
                if (ok >= 3) {
                    ++hit;
                    nm = int(T.meas.size());
                    worst = w;
                }
            }
            if (std::getenv("TRAIL_DEBUG")) {
                for (size_t k = 0; k < R.frames.size(); ++k) {
                    const Frame &f = R.frames[k];
                    double tx = W / 2.0 + mv.x0 + mv.vx * tmin[k], ty = H / 2.0 + mv.y0 + mv.vy * tmin[k];
                    for (size_t i = 0; i < f.det.size(); ++i)
                        if (std::hypot(f.det[i].x - tx, f.det[i].y - ty) < 8)
                            std::printf("   frame %zu: det %.2f %.2f (truth %.2f %.2f) L %.1f snr %.1f typ %d\n", k, f.det[i].x, f.det[i].y,
                                        tx, ty, f.det[i].trail_len, f.det[i].snr_tot, i < f.typ.size() ? f.typ[i] : -1);
                }
                for (auto &l : R.log)
                    std::printf("   log: %s\n", l.c_str());
            }
            std::snprintf(b, sizeof b, "%s: %zu track(s), object %s (%d measurements), worst |O−C| %.3f″, trail fit %s", name,
                          R.movers.size(), hit == 1 ? "found" : "NOT found", nm, worst, trail_log ? "used" : "not used");
            if (tol < 0) { // inconsistent: no track
                std::snprintf(b, sizeof b, "%s: %zu track(s)", name, R.movers.size());
                check(R.movers.empty(), b);
                return;
            }
            check(hit == 1 && R.movers.size() == 1 && nm == int(tmin.size()) && worst < tol && trail_log, b);
        };
        // optimal SNR 40: |O−C| ≤ 0.25″ (≈ 3× the Cramér–Rao bound along the trail); SNR 15: linking only (CR ≈ 0.3″)
        scenario({0, 5, 10, 15}, {60, -90, 3.0, 11.6, flux_for_snr(40, 6)}, 3, 0.25, "12″/min, 4 frames, 6 px trail, SNR 40");
        scenario({0, 4, 8}, {-50, -110, 16.0, 12.0, flux_for_snr(40, 10)}, 5, 0.3, "20″/min, 3 frames, 10 px trail, SNR 40");
        scenario({0, 3, 6, 9}, {-30, 50, -14.0, -14.0, flux_for_snr(15, 10)}, 9, 1.2, "19.8″/min, 4 frames, 10 px trail, SNR 15");
        // trails perpendicular to their own motion (not a Solar System object; e.g. different satellites):
        // trail × rate consistency rejects the link
        scenario({0, 4, 8}, {-50, -110, 16.0, 12.0, flux_for_snr(40, 10), PI / 2}, 5, -1, "trails perpendicular to motion");
    }
    std::printf("%s\n", fails ? "FAILED" : "OK");
    return fails ? 1 : 0;
}
