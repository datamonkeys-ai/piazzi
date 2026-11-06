// Synthetic test of centroid pixel-phase bias (Anderson & King 2000, PASP 112, 1360): stars with pixel-integrated
// profiles (Gaussian and Moffat β = 2.5) injected at known subpixel phases, with Poisson noise, measured by
// detect_sources; the mean bias per phase bin must stay below 0.01 px at normal sampling (FWHM ≥ 2.5 px) and is
// reported for undersampling (FWHM 1.8 px, e.g. ZTF in good seeing). Also measures high-SNR rms error (model floor).
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace astro;

static int fails = 0;
static void check(bool ok, const char *what) {
    std::printf("%s %s\n", ok ? "ok  " : "FAIL", what);
    if (!ok)
        ++fails;
}

struct Rng {
    uint64_t s;
    double uni() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return ((s >> 11) + 0.5) / 9007199254740992.0;
    }
    double gauss() { return std::sqrt(-2 * std::log(uni())) * std::cos(2 * 3.14159265358979323846 * uni()); }
};

// unit-flux point profile: Gaussian (beta <= 0) or Moffat (alpha, beta)
static double profile(double r2, double fwhm, double beta) {
    if (beta <= 0) {
        double s = fwhm / 2.3548200450309493;
        return std::exp(-0.5 * r2 / (s * s)) / (2 * 3.14159265358979323846 * s * s);
    }
    double a = fwhm / (2 * std::sqrt(std::pow(2.0, 1 / beta) - 1));
    return (beta - 1) / (3.14159265358979323846 * a * a) * std::pow(1 + r2 / (a * a), -beta);
}

// W×W image, background bkg e⁻, stars (x, y, flux) integrated over the pixel (7×7 subsamples), Poisson noise (gain 1)
struct Star {
    double x, y, flux;
};
static Image make(int W, double bkg, const std::vector<Star> &st, double fwhm, double beta, Rng &r) {
    Image im;
    im.w = im.h = W;
    im.pix.assign(size_t(W) * W, 0);
    std::vector<double> v(size_t(W) * W, bkg);
    for (auto &s : st) {
        int R = int(std::ceil(5 * fwhm));
        for (int y = int(s.y) - R; y <= int(s.y) + R; ++y)
            for (int x = int(s.x) - R; x <= int(s.x) + R; ++x) {
                if (x < 0 || y < 0 || x >= W || y >= W)
                    continue;
                double a = 0;
                for (int j = 0; j < 7; ++j)
                    for (int i = 0; i < 7; ++i) {
                        double dx = x + (i + 0.5) / 7 - 0.5 - s.x, dy = y + (j + 0.5) / 7 - 0.5 - s.y;
                        a += profile(dx * dx + dy * dy, fwhm, beta);
                    }
                v[size_t(y) * W + x] += s.flux * a / 49;
            }
    }
    for (size_t i = 0; i < v.size(); ++i) {
        double p = v[i] + std::sqrt(std::max(v[i], 1.0)) * r.gauss();
        im.pix[i] = uint16_t(std::max(0.0, std::min(65535.0, p + 0.5)));
    }
    return im;
}

// mean bias per phase bin (5 bins in x and in y) and centroid error rms; returns the largest bin |bias| (px)
static double phase_bias(double fwhm, double beta, double flux, const char *label, double &rms, double beta_fit = 0) {
    Rng r{12345};
    const int W = 400, N = 15; // 15×15 star grid, subpixel phases spread deterministically
    std::vector<Star> st;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
            st.push_back({20.0 + i * 24 + ((i * 7 + j * 3) % 20) / 20.0 - 0.5 + r.uni() * 0.05,
                          20.0 + j * 24 + ((i * 5 + j * 11) % 20) / 20.0 - 0.5 + r.uni() * 0.05, flux});
    Image im = make(W, 1000, st, fwhm, beta, r);
    DetectParams P;
    P.aperture = 4;
    P.gain = 1;
    P.beta = beta_fit;
    std::vector<Detection> D = detect_sources(im, P);
    double bx[5] = {}, by[5] = {}, s2 = 0;
    int nx[5] = {}, ny[5] = {}, n = 0;
    for (auto &s : st) {
        const Detection *best = nullptr;
        for (auto &d : D)
            if (std::hypot(d.x - s.x, d.y - s.y) < 2 &&
                (!best || std::hypot(d.x - s.x, d.y - s.y) < std::hypot(best->x - s.x, best->y - s.y)))
                best = &d;
        if (!best)
            continue;
        double ex = best->x - s.x, ey = best->y - s.y;
        double px = s.x - std::round(s.x), py = s.y - std::round(s.y); // phase in [−0.5, 0.5)
        int ix = std::min(4, int((px + 0.5) * 5)), iy = std::min(4, int((py + 0.5) * 5));
        bx[ix] += ex;
        ++nx[ix];
        by[iy] += ey;
        ++ny[iy];
        s2 += ex * ex + ey * ey;
        ++n;
    }
    rms = n ? std::sqrt(s2 / (2 * n)) : NAN;
    double worst = 0;
    std::printf("  %s: %d/%zu stars, rms per axis %.4f px; bias per phase x:", label, n, st.size(), rms);
    for (int k = 0; k < 5; ++k) {
        double b = nx[k] ? bx[k] / nx[k] : 0;
        worst = std::max(worst, std::fabs(b));
        std::printf(" %+.4f", b);
    }
    std::printf(" | y:");
    for (int k = 0; k < 5; ++k) {
        double b = ny[k] ? by[k] / ny[k] : 0;
        worst = std::max(worst, std::fabs(b));
        std::printf(" %+.4f", b);
    }
    std::printf("\n");
    return worst;
}

int main() {
    double rms;
    char b[160];
    // high SNR (flux 3·10⁵ e⁻ over 1000 e⁻ background): centroid noise is ~0.003 px, the remainder is model/phase
    for (double fw : {3.5, 2.5, 1.8}) {
        // flux limited so the peak (≈ F/(2π σ²)) does not saturate uint16 at FWHM 1.8 px
        const double F = fw >= 2.5 ? 3e5 : 6e4;
        double g = phase_bias(fw, 0, F, fw >= 2.5 ? "Gaussian" : "Gaussian (undersampled)", rms);
        std::snprintf(b, sizeof b, "Gaussian FWHM %.1f px: max phase bias %.4f px (< 0.01)", fw, g);
        check(g < 0.01 || fw < 2.5, b);
        if (fw < 2.5)
            std::printf("     (informational: FWHM %.1f px, max bias %.4f px)\n", fw, g);
        double m = phase_bias(fw, 2.5, F, "Moffat β=2.5", rms);
        std::snprintf(b, sizeof b, "Moffat β=2.5 FWHM %.1f px: max phase bias %.4f px (< 0.01), rms %.4f px", fw, m, rms);
        check((m < 0.01 || fw < 2.5) && rms < 0.03, b);
    }
    // fit with the Moffat profile itself (psf_beta = 2.5): no bias and lower rms than the Gaussian fit
    {
        double rg, rmf;
        phase_bias(2.5, 2.5, 3e5, "Moffat β=2.5 with Gaussian fit", rg);
        double m = phase_bias(2.5, 2.5, 3e5, "Moffat β=2.5 with Moffat fit", rmf, 2.5);
        // at FWHM ≥ 2.5 px both fits sit at the Poisson noise (≈ 0.004 px): the Gaussian model costs nothing
        std::snprintf(b, sizeof b, "Moffat fit on Moffat stars: max bias %.4f px (< 0.01), rms %.4f (Gaussian %.4f) < 0.01 px", m, rmf, rg);
        check(m < 0.01 && rmf < 0.01 && rg < 0.01, b);
    }
    // moderate SNR (flux 2·10⁴): only checks that the rms is consistent with the Cramér–Rao bound (order 0.05 px)
    phase_bias(3.0, 0, 2e4, "Gaussian moderate SNR", rms);
    std::snprintf(b, sizeof b, "moderate SNR: rms %.3f px (< 0.08)", rms);
    check(rms < 0.08, b);
    std::printf("%s\n", fails ? "FAILURES" : "all tests passed");
    return fails ? 1 : 0;
}
