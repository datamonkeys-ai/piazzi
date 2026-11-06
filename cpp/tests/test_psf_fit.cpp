// Synthetic tests of the trail fit (Vereš et al. 2012) and of the joint fit with a neighboring star (Stetson 1987):
//   1. profile normalization (sum ≈ 1);
//   2. known trail + noise at high SNR: center recovered with error < 0.1 px and flux within < 2 %;
//   3. Monte Carlo at moderate SNR: center scatter with the trail model vs. star PSF (L = 0);
//   4. neighboring star at ~1 FWHM: the joint fit recovers the target flux and center; the target-only fit is biased.
#include "../src/reduce_internal.hpp"
#include "piazzi/core.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace astro;
using namespace astro::internal;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

// deterministic generator (64-bit LCG) + Box–Muller
struct Rng {
    uint64_t s;
    double uni() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return ((s >> 11) + 0.5) / 9007199254740992.0;
    }
    double gauss() { return std::sqrt(-2 * std::log(uni())) * std::cos(2 * PI * uni()); }
};

struct Src {
    double x, y, flux;
    PsfShape P;
};

static Image make(const std::vector<Src> &ss, double bkg, double noise, Rng &r, int W = 64) {
    Image im;
    im.w = im.h = W;
    im.pix.resize(size_t(W) * W);
    for (int y = 0; y < W; ++y)
        for (int x = 0; x < W; ++x) {
            double v = bkg;
            // true profile integrated over the pixel by 5×5 sampling (the model uses 3×3)
            for (auto &s : ss)
                for (int j = 0; j < 5; ++j)
                    for (int i = 0; i < 5; ++i)
                        v += s.flux * psf_profile(s.P, x - 0.4 + 0.2 * i - s.x, y - 0.4 + 0.2 * j - s.y) / 25;
            v += noise * r.gauss();
            im.pix[size_t(y) * W + x] = uint16_t(std::lround(std::clamp(v, 0.0, 65535.0)));
        }
    return im;
}

int main() {
    char b[256];
    // 1. normalization
    {
        PsfShape P{1.1, 0.8, 5.0, std::cos(0.5), std::sin(0.5)};
        double s = 0;
        for (int y = -20; y <= 20; ++y)
            for (int x = -20; x <= 20; ++x)
                s += psf_pixel(P, x, y);
        std::snprintf(b, sizeof b, "trail profile normalized (sum %.4f)", s);
        check(std::fabs(s - 1) < 2e-3, b);
    }
    const double ang = 30 * PI / 180;
    const PsfShape star{1.2, 1.0, 0, 1, 0};
    PsfShape trail = star;
    trail.len = 6;
    trail.dx = std::cos(ang);
    trail.dy = std::sin(ang);
    // 2. SNR alto
    {
        Rng r{12345};
        const double x0 = 31.37, y0 = 32.71, F = 2e5;
        Image im = make({{x0, y0, F, trail}}, 1000, 5, r);
        std::vector<SceneSource> src{{x0 + 0.6, y0 - 0.5, F * 0.7, true, true, false}};
        double bk = 990;
        bool ok = fit_scene(im, x0 + 0.6, y0 - 0.5, 8, trail, star, src, bk);
        double e = std::hypot(src[0].x - x0, src[0].y - y0);
        std::snprintf(b, sizeof b, "trail L=6 px, high SNR: center error %.4f px, flux %.4f of true", e, src[0].flux / F);
        check(ok && e < 0.1 && std::fabs(src[0].flux / F - 1) < 0.02, b);
    }
    // 3. Monte Carlo, SNR ~ 15: trail vs. star PSF scatter
    {
        Rng r{777};
        double s_tr = 0, s_st = 0, bias_st = 0;
        const int N = 60;
        for (int k = 0; k < N; ++k) {
            double x0 = 30 + r.uni(), y0 = 30 + r.uni(), F = 6000;
            Image im = make({{x0, y0, F, trail}}, 1000, 10, r);
            for (int mdl = 0; mdl < 2; ++mdl) {
                std::vector<SceneSource> src{{x0 + 0.3, y0 + 0.2, F, true, true, false}};
                double bk = 1000;
                fit_scene(im, x0, y0, 8, mdl == 0 ? trail : star, star, src, bk);
                double e2 = std::pow(src[0].x - x0, 2) + std::pow(src[0].y - y0, 2);
                (mdl == 0 ? s_tr : s_st) += e2;
                if (mdl == 1)
                    bias_st += src[0].flux / F;
            }
        }
        s_tr = std::sqrt(s_tr / N);
        s_st = std::sqrt(s_st / N);
        std::snprintf(b, sizeof b, "Monte Carlo (%d, SNR≈15, L=6 px): center rms trail %.3f px vs. star PSF %.3f px", N, s_tr, s_st);
        check(s_tr < 0.25 && s_tr <= s_st * 1.05, b);
    }
    // 4. neighboring star at ~1 FWHM
    {
        Rng r{4242};
        const double x0 = 32.2, y0 = 31.6, F = 40000, xs = x0 + 2.4, ys = y0 + 1.3, Fs = 80000;
        Image im = make({{x0, y0, F, trail}, {xs, ys, Fs, star}}, 1000, 5, r);
        std::vector<SceneSource> one{{x0, y0, F, true, true, false}};
        double bk = 1000;
        fit_scene(im, x0, y0, 7, trail, star, one, bk);
        std::vector<SceneSource> two{{x0, y0, F, true, true, false}, {xs, ys, Fs * 0.5, false, false, false}};
        bk = 1000;
        bool ok = fit_scene(im, x0, y0, 7, trail, star, two, bk);
        double e1 = std::hypot(one[0].x - x0, one[0].y - y0), e2 = std::hypot(two[0].x - x0, two[0].y - y0);
        std::snprintf(b, sizeof b,
                      "neighbor at %.1f px: target-only error %.3f px flux %.3f; joint error %.3f px flux %.3f (neighbor %.3f)",
                      std::hypot(xs - x0, ys - y0), e1, one[0].flux / F, e2, two[0].flux / F, two[1].flux / Fs);
        check(ok && e2 < 0.1 && std::fabs(two[0].flux / F - 1) < 0.03 && std::fabs(two[1].flux / Fs - 1) < 0.03 && e1 > 3 * e2, b);
    }
    std::printf(fails ? "%d failures\n" : "OK\n", fails);
    return fails ? 1 : 0;
}
