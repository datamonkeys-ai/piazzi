// Source detection, centroiding, PSF fitting and aperture photometry.
//
// References (independent implementation):
//   - Bertin & Arnouts (1996), A&AS 117, 393 — SExtractor: mesh background with sigma-clipped statistics, median filter
//     on the mesh, interpolation; k·σ threshold detection on the filtered image, connected components with min area.
//   - Stetson (1987), PASP 99, 191 — DAOPHOT: aperture photometry with a sky annulus.
//   - Howell (2006), "Handbook of CCD Astronomy", 2nd ed., §4.4 — CCD equation for the SNR:
//       SNR = S / sqrt(S/g + n_pix·(1 + n_pix/n_sky)·σ_sky²).
//   - Levenberg (1944), Marquardt (1963) — nonlinear least squares (circular 2D Gaussian fit).
//   - Centroid by weighted first moments (e.g. Auer & van Altena 1978, AJ 83, 531).
//   - Trails: elongated matched filter (segment ⊗ PSF) over a grid of lengths and angles — Dawson et al. (2016),
//     AJ 152, 112 (matched-filter streak detection); Nir, Zackay & Ofek (2018), AJ 156, 229 ("Optimal and efficient
//     streak detection": a streak's SNR with the point kernel drops as √(FWHM/L)); trail model of Vereš et al.
//     (2012), PASP 124, 1197 (PSF convolved with a segment; the segment center is mid-exposure).
#include "piazzi/core.hpp"
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <vector>

namespace astro {

using internal::clipped_median;

namespace {

constexpr double FWHM_SIGMA = 2.3548200450309493; // 2·sqrt(2·ln 2)

// Mesh background and noise: cells of `cell` px, sigma-clipped (3σ) median in each cell, 3×3 median filter on the
// mesh and bilinear interpolation between cell centers.
struct BackgroundMap {
    int w = 0, h = 0;
    std::vector<float> bkg, rms;
    float b(int x, int y) const { return bkg[size_t(y) * w + x]; }
    float s(int x, int y) const { return rms[size_t(y) * w + x]; }
};

static void median3x3(std::vector<double> &g, int nx, int ny) {
    std::vector<double> o(g.size());
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            std::vector<double> v;
            for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di) {
                    int a = i + di, c = j + dj;
                    if (a >= 0 && c >= 0 && a < nx && c < ny && std::isfinite(g[size_t(c) * nx + a]))
                        v.push_back(g[size_t(c) * nx + a]);
                }
            o[size_t(j) * nx + i] = v.empty() ? NAN : internal::median(v);
        }
    g.swap(o);
}

static BackgroundMap background(const Image &im, bool mesh) {
    BackgroundMap M;
    M.w = im.w;
    M.h = im.h;
    size_t N = size_t(im.w) * im.h;
    M.bkg.assign(N, 0.f);
    M.rms.assign(N, 1.f);
    if (N == 0)
        return M;
    int cell = std::max(im.w, im.h) > 1500 ? 64 : 32;
    if (!mesh)
        cell = std::max(im.w, im.h);
    int nx = std::max(1, (im.w + cell - 1) / cell), ny = std::max(1, (im.h + cell - 1) / cell);
    std::vector<double> gb(size_t(nx) * ny, NAN), gs(size_t(nx) * ny, NAN);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            std::vector<double> v;
            for (int y = j * cell; y < std::min(im.h, (j + 1) * cell); ++y)
                for (int x = i * cell; x < std::min(im.w, (i + 1) * cell); ++x)
                    v.push_back(im.pix[size_t(y) * im.w + x]);
            if (v.size() < 16)
                continue;
            double sg = 0;
            int n = 0;
            double m = clipped_median(v, 3.0, 10, &sg, &n);
            if (n >= 8 && sg > 0) {
                gb[size_t(j) * nx + i] = m;
                gs[size_t(j) * nx + i] = sg;
            }
        }
    // cells without statistics: global median
    std::vector<double> okb, oks;
    for (size_t i = 0; i < gb.size(); ++i)
        if (std::isfinite(gb[i])) {
            okb.push_back(gb[i]);
            oks.push_back(gs[i]);
        }
    double fb = okb.empty() ? 0 : internal::median(okb), fs = oks.empty() ? 1 : internal::median(oks);
    if (nx * ny >= 9) {
        median3x3(gb, nx, ny);
        median3x3(gs, nx, ny);
    }
    for (size_t i = 0; i < gb.size(); ++i) {
        if (!std::isfinite(gb[i]))
            gb[i] = fb;
        if (!std::isfinite(gs[i]) || gs[i] <= 0)
            gs[i] = fs;
    }
    for (int y = 0; y < im.h; ++y) {
        double gy = std::clamp((y + 0.5) / cell - 0.5, 0.0, double(ny - 1));
        int j0 = std::min(int(gy), ny - 1), j1 = std::min(j0 + 1, ny - 1);
        double fy = gy - j0;
        for (int x = 0; x < im.w; ++x) {
            double gx = std::clamp((x + 0.5) / cell - 0.5, 0.0, double(nx - 1));
            int i0 = std::min(int(gx), nx - 1), i1 = std::min(i0 + 1, nx - 1);
            double fx = gx - i0;
            auto bil = [&](const std::vector<double> &g) {
                return (1 - fx) * (1 - fy) * g[size_t(j0) * nx + i0] + fx * (1 - fy) * g[size_t(j0) * nx + i1] +
                       (1 - fx) * fy * g[size_t(j1) * nx + i0] + fx * fy * g[size_t(j1) * nx + i1];
            };
            M.bkg[size_t(y) * im.w + x] = float(bil(gb));
            M.rms[size_t(y) * im.w + x] = float(bil(gs));
        }
    }
    return M;
}

// Circular profile fit by Levenberg–Marquardt, p = {A, x0, y0, w, B}: Gaussian I = B + A·exp(−r²/(2w²)) (beta = 0)
// or Moffat I = B + A·(1 + r²/w²)^−β (beta > 0; Moffat 1969, A&A 3, 455; Trujillo et al. 2001, MNRAS 328, 977:
// β ≈ 2.5–4 for atmospheric seeing), FWHM = 2.355·w or 2·w·√(2^(1/β) − 1).
static double fwhm_of(double w, double beta) {
    return beta > 0 ? 2 * w * std::sqrt(std::pow(2.0, 1 / beta) - 1) : FWHM_SIGMA * w;
}
static bool fit_gauss(const Image &im, int cx, int cy, int R, double p[5], double &rms, double beta = 0) {
    std::vector<double> X, Y, Z;
    for (int y = cy - R; y <= cy + R; ++y)
        for (int x = cx - R; x <= cx + R; ++x)
            if (x >= 0 && y >= 0 && x < im.w && y < im.h) {
                X.push_back(x);
                Y.push_back(y);
                Z.push_back(im.pix[size_t(y) * im.w + x]);
            }
    const int n = int(Z.size());
    if (n < 9)
        return false;
    auto prof = [&](double r2, double w) { return beta > 0 ? std::pow(1 + r2 / (w * w), -beta) : std::exp(-r2 / (2 * w * w)); };
    auto chi2 = [&](const double *q) {
        double c = 0;
        for (int i = 0; i < n; ++i) {
            double dx = X[i] - q[1], dy = Y[i] - q[2];
            double r = Z[i] - (q[4] + q[0] * prof(dx * dx + dy * dy, q[3]));
            c += r * r;
        }
        return c;
    };
    double lambda = 1e-3, c0 = chi2(p);
    bool conv = false;
    for (int it = 0; it < 60; ++it) {
        double JtJ[5][5] = {}, Jtr[5] = {};
        for (int i = 0; i < n; ++i) {
            double dx = X[i] - p[1], dy = Y[i] - p[2], s2 = p[3] * p[3], r2 = dx * dx + dy * dy;
            double e = prof(r2, p[3]);
            double r = Z[i] - (p[4] + p[0] * e);
            // profile derivative with respect to r²/w²: Gaussian −e/2, Moffat −β·e/(1 + r²/w²)
            double g = beta > 0 ? beta * e / (1 + r2 / s2) : 0.5 * e;
            double J[5] = {e, p[0] * g * 2 * dx / s2, p[0] * g * 2 * dy / s2, p[0] * g * 2 * r2 / (s2 * p[3]), 1.0};
            for (int a = 0; a < 5; ++a) {
                Jtr[a] += J[a] * r;
                for (int b = 0; b < 5; ++b)
                    JtJ[a][b] += J[a] * J[b];
            }
        }
        bool improved = false;
        for (int tries = 0; tries < 10 && !improved; ++tries) {
            // (JᵀJ + λ·diag(JᵀJ))·δ = Jᵀr by Gaussian elimination with partial pivoting
            double M[5][6];
            for (int a = 0; a < 5; ++a) {
                for (int b = 0; b < 5; ++b)
                    M[a][b] = JtJ[a][b] + (a == b ? lambda * JtJ[a][a] : 0);
                M[a][5] = Jtr[a];
            }
            bool sing = false;
            for (int c = 0; c < 5; ++c) {
                int piv = c;
                for (int r = c + 1; r < 5; ++r)
                    if (std::fabs(M[r][c]) > std::fabs(M[piv][c]))
                        piv = r;
                if (std::fabs(M[piv][c]) < 1e-300) {
                    sing = true;
                    break;
                }
                for (int k = 0; k < 6; ++k)
                    std::swap(M[c][k], M[piv][k]);
                for (int r = 0; r < 5; ++r)
                    if (r != c) {
                        double f = M[r][c] / M[c][c];
                        for (int k = c; k < 6; ++k)
                            M[r][k] -= f * M[c][k];
                    }
            }
            if (sing)
                return false;
            double q[5];
            for (int a = 0; a < 5; ++a)
                q[a] = p[a] + M[a][5] / M[a][a];
            q[3] = std::fabs(q[3]);
            double c1 = (q[3] > 0.2 && q[3] < 3.0 * R) ? chi2(q) : INFINITY;
            if (c1 < c0) {
                double step = std::fabs(q[1] - p[1]) + std::fabs(q[2] - p[2]);
                bool small = (c0 - c1) < 1e-8 * c0 || step < 1e-5;
                for (int a = 0; a < 5; ++a)
                    p[a] = q[a];
                c0 = c1;
                lambda = std::max(lambda * 0.3, 1e-9);
                improved = true;
                if (small)
                    conv = true;
            } else
                lambda *= 10;
        }
        if (!improved) { // λ exhausted: local minimum
            conv = true;
            break;
        }
        if (conv)
            break;
    }
    rms = std::sqrt(c0 / n);
    return conv && p[0] > 0 && std::fabs(p[1] - cx) < R && std::fabs(p[2] - cy) < R;
}

// Measures a source whose brightest pixel is (px, py): local background in the annulus, moment centroid, Gaussian fit,
// circular aperture and SNR from the CCD equation.
static bool measure_source(const Image &im, int px, int py, const DetectParams &P, Detection &d, double rwin = 1e9) {
    const double rap = std::max(1, P.aperture);
    const double rcen = std::min(rap, rwin); // centroid/fit window
    const double rin = rap + 3, rout = rap + 8;
    // sky annulus
    std::vector<double> sky;
    int Ro = int(std::ceil(rout));
    for (int y = py - Ro; y <= py + Ro; ++y)
        for (int x = px - Ro; x <= px + Ro; ++x) {
            if (x < 0 || y < 0 || x >= im.w || y >= im.h)
                continue;
            double r2 = double(x - px) * (x - px) + double(y - py) * (y - py);
            if (r2 >= rin * rin && r2 <= rout * rout)
                sky.push_back(im.pix[size_t(y) * im.w + x]);
        }
    if (sky.size() < 10)
        return false;
    double sb = 0;
    int nsky = 0;
    double b = clipped_median(sky, 3.0, 10, &sb, &nsky);
    if (!(sb > 0))
        sb = 1;
    // moment centroid (circular window of radius rap, two recentered iterations)
    double xc = px, yc = py;
    double mxx = 0, myy = 0;
    for (int it = 0; it < 3; ++it) {
        double sw = 0, sx = 0, sy = 0, sxx = 0, syy = 0;
        int R = int(std::ceil(rcen));
        for (int y = int(std::lround(yc)) - R; y <= int(std::lround(yc)) + R; ++y)
            for (int x = int(std::lround(xc)) - R; x <= int(std::lround(xc)) + R; ++x) {
                if (x < 0 || y < 0 || x >= im.w || y >= im.h)
                    continue;
                double dx = x - xc, dy = y - yc;
                if (dx * dx + dy * dy > rcen * rcen)
                    continue;
                double v = im.pix[size_t(y) * im.w + x] - b;
                if (v <= 0)
                    continue;
                sw += v;
                sx += v * x;
                sy += v * y;
                sxx += v * dx * dx;
                syy += v * dy * dy;
            }
        if (sw <= 0)
            return false;
        double nx = sx / sw, ny = sy / sw;
        mxx = sxx / sw;
        myy = syy / sw;
        if (std::hypot(nx - px, ny - py) > rcen)
            break;
        xc = nx;
        yc = ny;
    }
    double peak = im.at(px, py);
    double amp = peak - b;
    if (amp <= 0)
        return false;
    // Gaussian fit (optional: if it does not converge, the moments are kept)
    double sig_m = std::sqrt(std::max(0.25, 0.5 * (mxx + myy)));
    // initial width: σ of the moments; for Moffat, the w with the same FWHM
    double w0 = std::min(sig_m, rap);
    if (P.beta > 0)
        w0 = FWHM_SIGMA * w0 / (2 * std::sqrt(std::pow(2.0, 1 / P.beta) - 1));
    double q[5] = {amp, xc, yc, w0, b};
    double frms = 0;
    bool conv = fit_gauss(im, px, py, int(std::max(2.0, std::min(rap + 1, rcen + 1))), q, frms, P.beta);
    d = Detection{};
    if (conv && q[3] > 0.25 && std::hypot(q[1] - xc, q[2] - yc) < 1.5) {
        d.x = q[1];
        d.y = q[2];
        d.fwhm = fwhm_of(q[3], P.beta);
        d.fit = frms / std::max(q[0], 1e-9);
        amp = q[0];
        d.converged = true;
    } else {
        d.x = xc;
        d.y = yc;
        d.fwhm = FWHM_SIGMA * sig_m;
        d.fit = 1;
    }
    // circular aperture with area weights (5×5 subsampling per pixel)
    double flux = 0, npix = 0;
    int R = int(std::ceil(rap)) + 1;
    for (int y = int(d.y) - R; y <= int(d.y) + R + 1; ++y)
        for (int x = int(d.x) - R; x <= int(d.x) + R + 1; ++x) {
            if (x < 0 || y < 0 || x >= im.w || y >= im.h)
                continue;
            int in = 0;
            for (int sy = 0; sy < 5; ++sy)
                for (int sx = 0; sx < 5; ++sx) {
                    double dx = x - 0.4 + 0.2 * sx - d.x, dy = y - 0.4 + 0.2 * sy - d.y;
                    in += dx * dx + dy * dy <= rap * rap;
                }
            if (!in)
                continue;
            double fr = in / 25.0;
            flux += fr * (im.pix[size_t(y) * im.w + x] - b);
            npix += fr;
        }
    double g = P.gain > 0 ? P.gain : 1.0;
    double var = std::max(flux, 0.0) / g + npix * (1 + npix / std::max(nsky, 1)) * sb * sb;
    d.flux = flux;
    d.snr_tot = flux > 0 ? flux / std::sqrt(var) : 0;
    d.snr = amp / sb;
    d.peak = peak;
    d.bkg = b;
    d.sigma = sb;
    d.ix = px;
    d.iy = py;
    return true;
}

} // namespace

DetectParams detect_params(const Settings &s) {
    DetectParams p;
    auto has = [&](const char *k) { return !s.str(k).empty(); };
    if (has("aperture_radius"))
        p.aperture = std::max(1, int(s.num("aperture_radius")));
    if (has("detect_sigma") && s.num("detect_sigma") > 0)
        p.limit = s.num("detect_sigma");
    if (has("saturation") && s.num("saturation") > 0)
        p.saturation = s.num("saturation");
    if (has("min_fwhm"))
        p.min_fwhm = s.num("min_fwhm");
    if (has("plate_max_rms"))
        p.fit_rms = s.num("plate_max_rms");
    if (has("background_mode"))
        p.bkground = s.num("background_mode") != 0;
    if (has("gain") && s.num("gain") > 0)
        p.gain = s.num("gain");
    if (has("matched_filter"))
        p.matched = s.num("matched_filter") != 0;
    if (has("psf_beta") && s.num("psf_beta") > 1)
        p.beta = s.num("psf_beta");
    if (has("trail_detect"))
        p.trails = s.num("trail_detect") != 0;
    if (has("trail_max_len") && s.num("trail_max_len") > 0)
        p.trail_max = s.num("trail_max_len");
    if (has("trail_prescreen"))
        p.trail_prescreen = s.num("trail_prescreen") != 0;
    if (has("trail_fit_iter") && s.num("trail_fit_iter") > 0)
        p.trail_iter = int(s.num("trail_fit_iter"));
    if (has("trail_gate"))
        p.trail_gate = s.num("trail_gate");
    return p;
}

// Detection pass on the image filtered by a kernel `ksig` (Gaussian σ in px; 0 = 3×3 binomial).
// If Fout is non-null, returns the filtered (background-subtracted) image and the kernel noise factor (σ_F = gain·σ).
static std::vector<Detection> detect_pass(const Image &im, const DetectParams &P, const BackgroundMap &M, double ksig,
                                          std::vector<float> *Fout = nullptr, double *gain_out = nullptr) {
    std::vector<Detection> out;
    const int W = im.w, H = im.h;
    std::vector<float> F(size_t(W) * H, 0.f);
    // separable 1D kernel: binomial (1,2,1) or sampled Gaussian truncated at 3σ (filter matched to the PSF: for white
    // noise, convolving with the PSF's own profile maximizes the SNR of a point source — Irwin 1985, MNRAS 214, 575;
    // Bertin & Arnouts 1996, §6.2). Noise gain of the normalized kernel: sqrt(Σw²)/Σw (per axis; product in 2D).
    std::vector<double> k1;
    if (ksig <= 0)
        k1 = {1, 2, 1};
    else {
        int r = std::max(1, int(std::ceil(3 * ksig)));
        for (int i = -r; i <= r; ++i)
            k1.push_back(std::exp(-0.5 * i * i / (ksig * ksig)));
    }
    const int kr = int(k1.size()) / 2;
    double sw = 0, sw2 = 0;
    for (double w : k1) {
        sw += w;
        sw2 += w * w;
    }
    const double noise_gain = sw2 / (sw * sw); // (sqrt(Σw²)/Σw)² = 2D gain of the separable kernel
    {
        std::vector<float> T(size_t(W) * H, 0.f);
        for (int y = 0; y < H; ++y)
            for (int x = 0; x < W; ++x) {
                double a = 0, ws = 0;
                for (int d = -kr; d <= kr; ++d) {
                    int u = x + d;
                    if (u < 0 || u >= W)
                        continue;
                    double w = k1[d + kr];
                    a += w * (im.pix[size_t(y) * W + u] - M.b(u, y));
                    ws += w;
                }
                T[size_t(y) * W + x] = float(a / ws);
            }
        for (int y = 0; y < H; ++y)
            for (int x = 0; x < W; ++x) {
                double a = 0, ws = 0;
                for (int d = -kr; d <= kr; ++d) {
                    int v = y + d;
                    if (v < 0 || v >= H)
                        continue;
                    double w = k1[d + kr];
                    a += w * T[size_t(v) * W + x];
                    ws += w;
                }
                F[size_t(y) * W + x] = float(a / ws);
            }
    }
    std::vector<char> mask(F.size(), 0);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            mask[size_t(y) * W + x] = F[size_t(y) * W + x] > P.limit * noise_gain * M.s(x, y);
    const int minarea = 3;
    const int edge = std::max(2, P.aperture + 1);
    std::vector<char> seen(F.size(), 0);
    std::deque<int> q;
    struct Cand {
        int px, py, n;
        double det;
    };
    std::vector<Cand> cands;
    for (int y0 = 0; y0 < H; ++y0)
        for (int x0 = 0; x0 < W; ++x0) {
            size_t i0 = size_t(y0) * W + x0;
            if (!mask[i0] || seen[i0])
                continue;
            // 8-connected component
            int n = 0, bx = x0, by = y0;
            double bv = -1e30;
            std::vector<int> pix;
            q.clear();
            q.push_back(int(i0));
            seen[i0] = 1;
            while (!q.empty()) {
                int i = q.front();
                q.pop_front();
                int x = i % W, y = i / W;
                ++n;
                pix.push_back(i);
                double v = im.pix[size_t(i)] - M.b(x, y);
                if (v > bv) {
                    bv = v;
                    bx = x;
                    by = y;
                }
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        int a = x + dx, c = y + dy;
                        if (a < 0 || c < 0 || a >= W || c >= H)
                            continue;
                        size_t j = size_t(c) * W + a;
                        if (mask[j] && !seen[j]) {
                            seen[j] = 1;
                            q.push_back(int(j));
                        }
                    }
            }
            if (n < minarea)
                continue;
            // overlapping-source separation (simplified deblending): each local maximum of the filtered image in a
            // 5×5 neighborhood, above threshold, becomes a candidate; the component peak is always a candidate.
            cands.push_back({bx, by, n, F[size_t(by) * W + bx] / (noise_gain * M.s(bx, by))});
            for (int i : pix) {
                int x = i % W, y = i / W;
                if (x == bx && y == by)
                    continue;
                float v = F[size_t(i)];
                bool mx = true;
                for (int dy = -2; dy <= 2 && mx; ++dy)
                    for (int dx = -2; dx <= 2; ++dx) {
                        int a = x + dx, c = y + dy;
                        if ((dx || dy) && a >= 0 && c >= 0 && a < W && c < H && F[size_t(c) * W + a] >= v) {
                            mx = false;
                            break;
                        }
                    }
                double sig = v / (noise_gain * M.s(x, y));
                if (mx && sig > P.limit)
                    cands.push_back({x, y, n, sig});
            }
        }
    // distance to the nearest neighboring candidate (limits the centroid/fit window for overlapping sources)
    std::vector<double> dnear(cands.size(), 1e9);
    for (size_t a = 0; a < cands.size(); ++a)
        for (size_t b = 0; b < cands.size(); ++b)
            if (a != b)
                dnear[a] = std::min(dnear[a], std::hypot(double(cands[a].px - cands[b].px), double(cands[a].py - cands[b].py)));
    for (size_t ci = 0; ci < cands.size(); ++ci) {
        const Cand &c = cands[ci];
        if (c.px < edge + P.margins[0] || c.py < edge + P.margins[2] || c.px >= W - edge - P.margins[1] || c.py >= H - edge - P.margins[3])
            continue; // border
        if (im.at(c.px, c.py) >= P.saturation)
            continue; // saturated
        Detection d;
        if (!measure_source(im, c.px, c.py, P, d, std::max(2.0, 0.5 * dnear[ci] + 0.5)))
            continue;
        if (d.flux <= 0 || d.fwhm < P.min_fwhm || d.fwhm > 6.0 * (P.aperture + 1))
            continue;
        if (d.x < 0 || d.y < 0 || d.x > W - 1 || d.y > H - 1)
            continue;
        d.det_snr = c.det;
        out.push_back(d);
    }
    // sort by flux (desc.) and remove duplicates (two components that converged to the same source)
    std::stable_sort(out.begin(), out.end(), [](const Detection &a, const Detection &b) { return a.flux > b.flux; });
    std::vector<Detection> uniq;
    for (auto &d : out) {
        bool dup = false;
        for (auto &u : uniq)
            if (std::hypot(u.x - d.x, u.y - d.y) < 1.5) {
                dup = true;
                break;
            }
        if (!dup)
            uniq.push_back(d);
    }
    if (Fout)
        Fout->swap(F);
    if (gain_out)
        *gain_out = noise_gain;
    return uniq;
}

// Trail search. S = F/(g·σ) is the SNR map of the point filter (Gaussian kernel with σ = ks, the PSF); the mean of S
// along a segment of length L and angle θ is the matched filter for a trail (segment ⊗ PSF; Dawson et al. 2016,
// Nir et al. 2018), normalized by the empirical σ (MAD) of each kernel. Grid: L = 2·FWHM·2^j up to trail_max_len; angular
// step such that the segment tip moves ≤ σ_PSF. A component of the max map above threshold is a trail if the
// best elongated kernel beats the point kernel (likelihood ratio: for a point source the point kernel is
// optimal); the trail is then fitted (Vereš et al. 2012) with free length and direction, and the point detections that
// were fragments of it are replaced by a single detection at the segment center.
static void trail_search(const Image &im, const DetectParams &P, const BackgroundMap &M, double ks, const std::vector<float> &F, double ng,
                         std::vector<Detection> &D, DetectStats *st, const std::vector<double> *mask = nullptr,
                         std::vector<char> *frag_out = nullptr) {
    using clk = std::chrono::steady_clock;
    auto sec = [](clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); };
    DetectStats dummy;
    if (!st)
        st = &dummy;
    auto t0 = clk::now();
    const int W = im.w, H = im.h;
    const size_t NP = size_t(W) * H;
    const double fw = FWHM_SIGMA * ks;
    std::vector<float> S(NP);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            S[size_t(y) * W + x] = float(F[size_t(y) * W + x] / (ng * std::max(double(M.s(x, y)), 1e-6)));
    std::vector<float> best(NP, 0.f), bestu(NP, 0.f); // clipped (detection) and unclipped (peak, test) response
    std::vector<short> kbest(NP, -1);
    struct Kern {
        double L, th, norm, cap = INFINITY;
        std::vector<int> off; // linear offsets (y·W + x) of the samples
        int r;                // half extent in px
    };
    std::vector<Kern> K;
    for (double L = std::max(4.0, 2 * fw); L <= P.trail_max + 1e-9; L *= 2) {
        int na = std::max(4, int(std::ceil(PI * L / (2 * ks))));
        for (int a = 0; a < na; ++a) {
            Kern k{L, PI * a / na, 1, INFINITY, {}, 0};
            int ns = int(std::ceil(L)) + 1;
            for (int i = 0; i < ns; ++i) {
                double t = -L / 2 + L * i / (ns - 1);
                int dx = int(std::lround(t * std::cos(k.th))), dy = int(std::lround(t * std::sin(k.th)));
                k.off.push_back(dy * W + dx);
                k.r = std::max(k.r, std::max(std::abs(dx), std::abs(dy)));
            }
            K.push_back(std::move(k));
        }
    }
    st->nkern = int(K.size());
    if (K.empty() || K.size() > 30000)
        return;
    const int edge = std::max(2, P.aperture + 1);
    const float gate = float(P.trail_gate);
    for (size_t kk = 0; kk < K.size(); ++kk) {
        Kern &k = K[kk];
        const int r = k.r + edge;
        if (2 * r >= W || 2 * r >= H)
            continue;
        const double inv = 1.0 / k.off.size();
        double cap = INFINITY;
        double gu = 0; // unclipped sum from the last call
        auto G = [&](size_t i) {
            double a = 0, u = 0;
            for (int o : k.off) {
                double v = S[size_t(long(i) + o)];
                a += std::min(v, cap);
                u += v;
            }
            gu = u * inv;
            return a * inv;
        };
        // σ of the mean along the kernel: S is unit white noise filtered by the Gaussian σ = ks, autocorrelation
        // ρ(d) = exp(−d²/4ks²); Var = Σ_ij ρ(d_ij)/n² (analytic: an empirical σ would be inflated by stars)
        double var = 0;
        for (int o1 : k.off)
            for (int o2 : k.off) {
                int dx1 = ((o1 % W) + W + W / 2) % W - W / 2, dx2 = ((o2 % W) + W + W / 2) % W - W / 2;
                double ddx = dx1 - dx2, ddy = double(o1 - dx1) / W - double(o2 - dx2) / W;
                var += std::exp(-(ddx * ddx + ddy * ddy) / (4 * ks * ks));
            }
        double sg = std::sqrt(var) * inv;
        // each sample clipped at 1.5·threshold·σ_G: a trail at threshold has samples ~ threshold·σ_G, stays above it,
        // but a bright star crossed by the segment (few samples) does not create a "trail" (robust filter)
        cap = 1.5 * (P.limit + 1.0) * sg;
        // non-white noise (structured background, resampled images): empirical σ (MAD) of response on a sparse grid;
        // the larger of the two is used
        {
            std::vector<double> smp;
            // only source-free segments (no sample with S > 4): stars do not inflate σ
            for (int y = r; y < H - r; y += 5)
                for (int x = r; x < W - r; x += 5) {
                    size_t i = size_t(y) * W + x;
                    bool clean = true;
                    for (int o : k.off)
                        clean = clean && S[size_t(long(i) + o)] <= 4;
                    if (clean)
                        smp.push_back(G(i));
                }
            double se = smp.size() > 50 ? internal::mad_sigma(smp) : 0;
            if (se > sg) {
                sg = se;
                cap = 1.5 * (P.limit + 1.0) * sg;
            }
        }
        k.norm = 1 / sg;
        k.cap = cap;
        for (int y = r; y < H - r; ++y)
            for (int x = r; x < W - r; ++x) {
                size_t i = size_t(y) * W + x;
                if (S[i] < gate)
                    continue;
                float t = float(G(i) * k.norm), tu = float(gu * k.norm);
                // tie (within 2 %) in the clipped response, which saturates on bright trails: the unclipped one decides
                if (t > 1.02f * best[i] || (t >= 0.98f * best[i] && tu > bestu[i])) {
                    best[i] = std::max(t, best[i]);
                    bestu[i] = tu;
                    kbest[i] = short(kk);
                }
            }
    }
    st->t_kernels = sec(t0);
    t0 = clk::now();
    const double thr = P.limit + 1.0; // hundreds of kernels tested per pixel: threshold above the point one
    // mask: above threshold, elongated kernel better than the point one at the pixel itself, and uniformly high samples
    // along the segment (median ≥ half the mean: a segment that only crosses a star has a high mean and a noise-level
    // median); this keeps the trail from merging with a neighboring star's component
    std::vector<char> mk(NP, 0);
    {
        std::vector<double> v;
        for (size_t i = 0; i < NP; ++i) {
            if (best[i] <= thr || bestu[i] < 1.1f * S[i])
                continue;
            const Kern &k = K[size_t(kbest[i])];
            v.clear();
            double vmax = 0;
            for (int o : k.off) {
                v.push_back(S[size_t(long(i) + o)]);
                vmax = std::max(vmax, v.back());
            }
            // the elongated kernel beats the point one along the whole segment, within 3σ noise at the max sample (else
            // a point source on it explains the response)
            mk[i] = internal::median(v) * k.norm >= 0.5 * bestu[i] && bestu[i] >= 1.1 * vmax - 3;
        }
    }
    std::vector<char> seen(NP, 0);
    std::deque<int> q;
    std::vector<Detection> T;
    struct TCand {
        int ib;
        double ratio, cx, cy;
    };
    std::vector<TCand> cands;
    for (size_t i0 = 0; i0 < NP; ++i0) {
        if (!mk[i0] || seen[i0])
            continue;
        std::vector<int> pix;
        q.clear();
        q.push_back(int(i0));
        seen[i0] = 1;
        int ib = int(i0);
        double smax = 0, sw = 0, sx = 0, sy = 0;
        while (!q.empty()) {
            int i = q.front();
            q.pop_front();
            pix.push_back(i);
            int x = i % W, y = i / W;
            if (bestu[size_t(i)] > bestu[size_t(ib)])
                ib = i;
            smax = std::max(smax, double(S[size_t(i)]));
            double w = best[size_t(i)] - thr;
            sw += w;
            sx += w * x;
            sy += w * y;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    int a = x + dx, c = y + dy;
                    if (a < 0 || c < 0 || a >= W || c >= H)
                        continue;
                    size_t j = size_t(c) * W + a;
                    if (mk[j] && !seen[j]) {
                        seen[j] = 1;
                        q.push_back(int(j));
                    }
                }
        }
        ++st->ncomp;
        if (pix.size() < 6 || kbest[size_t(ib)] < 0)
            continue;
        // a real trail leaves a ridge in the mask along almost its whole length (≥ half the kernel, which
        // is at most twice the trail): few-pixel components with a long kernel are noise peaks
        if (double(pix.size()) < 0.25 * K[size_t(kbest[size_t(ib)])].L)
            continue;
        // the elongated kernel must beat the point one (otherwise it is a point source: keeps the normal detection)
        const Kern &k = K[size_t(kbest[size_t(ib)])];
        const double tu = bestu[size_t(ib)];
        if (tu < 1.1 * smax)
            continue;
        double cx = sx / sw, cy = sy / sw;
        // known static source (deferred search): reference, catalog star, or detection with a counterpart in another
        // frame at the component center — galaxies and halos, not trails (a trail is in one frame only)
        if (mask) {
            bool hitm = false;
            for (size_t q = 0; q + 2 < mask->size() && !hitm; q += 3)
                hitm = std::hypot((*mask)[q] - cx, (*mask)[q + 1] - cy) <= (*mask)[q + 2];
            if (hitm)
                continue;
        }
        // cheap prescreen (the fit costs ~50 ms): a star-shaped point detection at the component center
        // (FWHM < 1.15·PSF and a good Gaussian fit) is a star; fragments of a resolved trail are wider;
        // bright-star halo: the real (non-Gaussian) PSF wings give arc-shaped elongated responses a few FWHM
        // from a star; the component is discarded if within max(3, 2·log10 SNR*)·FWHM + L/2 of a point source with
        // SNR* > 3× the component response
        {
            bool star = false;
            for (auto &e : D) {
                double dd = std::hypot(e.x - cx, e.y - cy);
                star = star || (dd < 1.0 && e.converged && e.fwhm < 1.15 * fw && e.fit < P.fit_rms);
                star = star || (e.snr_tot > 3 * tu && dd < std::max(3.0, 2 * std::log10(e.snr_tot)) * fw + 0.5 * k.L);
            }
            if (star)
                continue;
        }
        cands.push_back({ib, tu / std::max(smax, 1e-3), cx, cy});
    }
    // fits in order of elongation (elongated / point response), at most 40 per image: deep fields (SDSS)
    // have hundreds of elongated galaxies and each fit costs ~10–50 ms
    std::stable_sort(cands.begin(), cands.end(), [](const TCand &a, const TCand &b) { return a.ratio > b.ratio; });
    if (cands.size() > 40)
        cands.resize(40);
    st->ncand = int(cands.size());
    st->t_mask = sec(t0);
    t0 = clk::now();
    for (const TCand &cd : cands) {
        const int ib = cd.ib;
        const double cx = cd.cx, cy = cd.cy;
        const Kern &k = K[size_t(kbest[size_t(ib)])];
        internal::PsfShape tr{ks, ks, k.L, std::cos(k.th), std::sin(k.th), ks >= 1.2};
        const int R = std::min(60, int(std::ceil(3 * ks + 0.75 * k.L)) + 2);
        if (cx < R || cy < R || cx >= W - 1 - R || cy >= H - 1 - R)
            continue;
        // initial flux: sum of (I − background) within ≤ 2σ of the segment
        double f0 = 0, b0 = M.b(int(cx), int(cy)), pk = 0;
        for (int y = int(cy) - R; y <= int(cy) + R; ++y)
            for (int x = int(cx) - R; x <= int(cx) + R; ++x) {
                double u = x - cx, v = y - cy, t = std::clamp(u * tr.dx + v * tr.dy, -k.L / 2, k.L / 2);
                if (std::hypot(u - t * tr.dx, v - t * tr.dy) <= 2 * ks) {
                    f0 += im.pix[size_t(y) * W + x] - M.b(x, y);
                    pk = std::max(pk, double(im.pix[size_t(y) * W + x]));
                }
            }
        if (pk >= P.saturation || !(f0 > 0))
            continue;
        // prescreen (trail_prescreen, default 1): aperture SNR in a capsule of radius = aperture around the kernel
        // segment, with sky in an annulus of the capsule; below threshold + 0.5 the fit would be rejected by the same
        // criterion later (a faint trail keeps the point detections) and would cost ~50 ms
        if (P.trail_prescreen) {
            const double rap0 = std::max(1, P.aperture), h0 = k.L / 2;
            auto dk = [&](double u, double v) {
                double t = std::clamp(u * tr.dx + v * tr.dy, -h0, h0);
                return std::hypot(u - t * tr.dx, v - t * tr.dy);
            };
            const int Ro = int(std::ceil(h0 + rap0 + 8)) + 1;
            if (cx < Ro || cy < Ro || cx >= W - 1 - Ro || cy >= H - 1 - Ro)
                continue;
            std::vector<double> sky;
            double fl = 0, np = 0;
            for (int y = int(cy) - Ro; y <= int(cy) + Ro; ++y)
                for (int x = int(cx) - Ro; x <= int(cx) + Ro; ++x) {
                    double dd = dk(x - cx, y - cy);
                    if (dd >= rap0 + 3 && dd <= rap0 + 8)
                        sky.push_back(im.pix[size_t(y) * W + x]);
                }
            double sb = 0;
            int nsk = 0;
            double bb = clipped_median(sky, 3.0, 10, &sb, &nsk);
            for (int y = int(cy) - Ro; y <= int(cy) + Ro; ++y)
                for (int x = int(cx) - Ro; x <= int(cx) + Ro; ++x)
                    if (dk(x - cx, y - cy) <= rap0) {
                        fl += im.pix[size_t(y) * W + x] - bb;
                        np += 1;
                    }
            const double g = P.gain > 0 ? P.gain : 1.0;
            const double sn = fl > 0 && sb > 0 ? fl / std::sqrt(fl / g + np * (1 + np / std::max(nsk, 1)) * sb * sb) : 0;
            if (sn < P.limit + 0.5) {
                ++st->nskip;
                continue;
            }
        }
        std::vector<internal::SceneSource> src(1);
        src[0].x = cx;
        src[0].y = cy;
        src[0].flux = f0;
        src[0].trail = true;
        src[0].free_len = true;
        src[0].len_sigma = k.L; // weak prior: the kernel only locates the length within a factor of 2
        src[0].pix_noise = std::max(double(M.s(int(cx), int(cy))), 1e-6);
        double bk = b0;
        double rms_t = 0;
        ++st->nfit;
        if (!internal::fit_scene(im, cx, cy, R, tr, tr, src, bk, &rms_t, P.trail_iter, 3 * ks + 3) || !(src[0].flux > 0))
            continue;
        const double L = src[0].fit_len;
        if (L < 1.5 * fw || std::hypot(src[0].x - cx, src[0].y - cy) > 0.5 * k.L + 2 * ks)
            continue;
        // competing model: the point sources already detected on the trail (≥ 1), free positions/fluxes, same window.
        // The trail is kept if it has the lower information criterion χ²/σ² + 4·(no. of params) (penalty between
        // Akaike's 1974, IEEE TAC 19, 716, and Schwarz's 1978, Ann. Stat. 6, 461; chosen from synthetic tests): a star
        // pair 1–2 FWHM apart mimics a short trail, but the two-PSF model describes it better
        {
            internal::PsfShape pt{ks, ks, 0, 1, 0, ks >= 1.2};
            std::vector<internal::SceneSource> ps;
            const double ux = src[0].fit_dx, uy = src[0].fit_dy;
            for (auto &e : D) {
                double u = e.x - src[0].x, v = e.y - src[0].y, t = std::clamp(u * ux + v * uy, -L / 2, L / 2);
                if (std::hypot(u - t * ux, v - t * uy) <= fw && std::fabs(u) < R - 1 && std::fabs(v) < R - 1) {
                    internal::SceneSource q;
                    q.x = e.x;
                    q.y = e.y;
                    q.flux = std::max(e.flux, 1.0);
                    ps.push_back(q);
                }
            }
            if (ps.size() == 1) {
                // a single point detection (pair not separated by detection): the competitor is two PSFs along the axis
                internal::SceneSource a = ps[0], c = ps[0];
                a.x -= 0.25 * L * ux;
                a.y -= 0.25 * L * uy;
                c.x += 0.25 * L * ux;
                c.y += 0.25 * L * uy;
                a.flux = c.flux = 0.5 * ps[0].flux;
                ps = {a, c};
            }
            if (!ps.empty()) {
                double bp = bk, rms_p = 0;
                if (internal::fit_scene(im, cx, cy, R, tr, pt, ps, bp, &rms_p, P.trail_iter, 3 * ks + 3) || rms_p > 0) {
                    // number of capsule pixels (the same window for both fits)
                    double n = 0;
                    for (int y = int(std::lround(cy)) - R; y <= int(std::lround(cy)) + R; ++y)
                        for (int x = int(std::lround(cx)) - R; x <= int(std::lround(cx)) + R; ++x) {
                            const double h = tr.len / 2 + 3 * ks + 3;
                            double u = x - cx, v = y - cy, t = std::clamp(u * tr.dx + v * tr.dy, -h, h);
                            n += std::hypot(u - t * tr.dx, v - t * tr.dy) <= 3 * ks + 3;
                        }
                    const double s2n = std::pow(src[0].pix_noise, 2);
                    const double pen = 4; // between Akaike (2) and Schwarz (ln n ≈ 6)
                    double aic_t = n * rms_t * rms_t / s2n + pen * 6, aic_p = n * rms_p * rms_p / s2n + pen * (1 + 3 * double(ps.size()));
                    if (aic_p <= aic_t)
                        continue;
                }
            }
        }
        // profile along the trail: in a trail the S-map ridge is flat (unit noise only); two neighboring point sources
        // (star pair) have a valley between the peaks. Reject if the minimum over the central 80 % falls below
        // half the mean minus 2.5σ.
        {
            auto Sb = [&](double x, double y) {
                int x0 = int(std::floor(x)), y0 = int(std::floor(y));
                double fx = x - x0, fy = y - y0;
                auto at = [&](int a, int c) { return double(S[size_t(c) * W + a]); };
                return (1 - fy) * ((1 - fx) * at(x0, y0) + fx * at(x0 + 1, y0)) +
                       fy * ((1 - fx) * at(x0, y0 + 1) + fx * at(x0 + 1, y0 + 1));
            };
            double smin = 1e30, ssum = 0;
            int ns = 0;
            for (double t = -0.4 * L; t <= 0.4 * L + 1e-9; t += 0.5) {
                double v = Sb(src[0].x + t * src[0].fit_dx, src[0].y + t * src[0].fit_dy);
                smin = std::min(smin, v);
                ssum += v;
                ++ns;
            }
            if (ns == 0 || smin < 0.5 * ssum / ns - 2.5)
                continue;
        }
        // aperture photometry in a capsule of radius rap around the fitted segment; sky between rap+3 and rap+8
        const double rap = std::max(1, P.aperture), ux = src[0].fit_dx, uy = src[0].fit_dy;
        auto dseg = [&](double u, double v) {
            double t = std::clamp(u * ux + v * uy, -L / 2, L / 2);
            return std::hypot(u - t * ux, v - t * uy);
        };
        const double x0 = src[0].x, y0 = src[0].y;
        const int Ro = int(std::ceil(L / 2 + rap + 8)) + 1;
        if (x0 < Ro || y0 < Ro || x0 >= W - 1 - Ro || y0 >= H - 1 - Ro)
            continue;
        std::vector<double> sky;
        for (int y = int(y0) - Ro; y <= int(y0) + Ro; ++y)
            for (int x = int(x0) - Ro; x <= int(x0) + Ro; ++x) {
                double dd = dseg(x - x0, y - y0);
                if (dd >= rap + 3 && dd <= rap + 8)
                    sky.push_back(im.pix[size_t(y) * W + x]);
            }
        double sb = 0;
        int nsky = 0;
        double b = clipped_median(sky, 3.0, 10, &sb, &nsky);
        double flux = 0, npix = 0, peak = 0;
        for (int y = int(y0) - Ro; y <= int(y0) + Ro; ++y)
            for (int x = int(x0) - Ro; x <= int(x0) + Ro; ++x) {
                int in = 0;
                for (int j = 0; j < 5; ++j)
                    for (int i = 0; i < 5; ++i)
                        in += dseg(x - 0.4 + 0.2 * i - x0, y - 0.4 + 0.2 * j - y0) <= rap;
                if (!in)
                    continue;
                flux += in / 25.0 * (im.pix[size_t(y) * W + x] - b);
                npix += in / 25.0;
                peak = std::max(peak, double(im.pix[size_t(y) * W + x]));
            }
        if (!(sb > 0) || !(flux > 0))
            continue;
        Detection d;
        d.snr_tot = flux / std::sqrt(flux / (P.gain > 0 ? P.gain : 1.0) + npix * (1 + npix / std::max(nsky, 1)) * sb * sb);
        if (d.snr_tot < P.limit + 0.5) // faint trail: keep the point detections (faint source elongated by noise)
            continue;
        d.x = x0;
        d.y = y0;
        d.flux = flux;
        d.fwhm = FWHM_SIGMA * std::sqrt(ks * ks + L * L / 12); // effective width (centroid σ along the trail)
        double g = P.gain > 0 ? P.gain : 1.0;
        d.snr_tot = flux / std::sqrt(flux / g + npix * (1 + npix / std::max(nsky, 1)) * sb * sb);
        d.snr = (peak - b) / sb;
        d.fit = 1;
        d.peak = peak;
        d.bkg = b;
        d.sigma = sb;
        d.det_snr = best[size_t(ib)];
        d.ix = int(std::lround(x0));
        d.iy = int(std::lround(y0));
        d.trail_len = L;
        d.trail_ang = std::atan2(uy, ux);
        // a dominant point source inside the capsule (star on the trail or star + noise): keep the point source
        bool star_dom = false;
        for (auto &e : D)
            if (dseg(e.x - x0, e.y - y0) <= std::max(rap, fw) && e.flux > 0.8 * flux && e.fwhm < 1.3 * fw)
                star_dom = true;
        if (!star_dom)
            T.push_back(d);
    }
    st->t_fits = sec(t0);
    t0 = clk::now();
    if (frag_out)
        frag_out->assign(D.size(), 0);
    if (T.empty())
        return;
    // overlapping trails (same trail found by two components): keep the most significant
    std::stable_sort(T.begin(), T.end(), [](const Detection &a, const Detection &b) { return a.det_snr > b.det_snr; });
    std::vector<Detection> keep;
    for (auto &t : T) {
        bool dup = false;
        for (auto &u : keep)
            dup |= std::hypot(u.x - t.x, u.y - t.y) < 0.5 * std::max(u.trail_len, t.trail_len) + fw;
        if (!dup)
            keep.push_back(t);
    }
    // remove the point fragments of each trail
    const double rap = std::max(1, P.aperture);
    if (frag_out) {
        // deferred mode: mark the fragments and return only the trails (D ← trails); detection order is unchanged
        for (size_t i = 0; i < D.size(); ++i)
            for (auto &t : keep) {
                double ux = std::cos(t.trail_ang), uy = std::sin(t.trail_ang), u = D[i].x - t.x, v = D[i].y - t.y;
                double sp = std::clamp(u * ux + v * uy, -t.trail_len / 2, t.trail_len / 2);
                if (std::hypot(u - sp * ux, v - sp * uy) <= std::max(rap, fw))
                    (*frag_out)[i] = 1;
            }
        D.swap(keep);
        st->ntrail = int(D.size());
        st->t_replace = sec(t0);
        return;
    }
    std::vector<Detection> out;
    for (auto &e : D) {
        bool frag = false;
        for (auto &t : keep) {
            double ux = std::cos(t.trail_ang), uy = std::sin(t.trail_ang), u = e.x - t.x, v = e.y - t.y;
            double s = std::clamp(u * ux + v * uy, -t.trail_len / 2, t.trail_len / 2);
            frag |= std::hypot(u - s * ux, v - s * uy) <= std::max(rap, fw);
        }
        if (!frag)
            out.push_back(e);
    }
    for (auto &t : keep)
        out.push_back(t);
    std::stable_sort(out.begin(), out.end(), [](const Detection &a, const Detection &b) { return a.flux > b.flux; });
    D.swap(out);
    st->ntrail = int(keep.size());
    st->t_replace = sec(t0);
}

std::vector<Detection> detect_sources(const Image &im, const DetectParams &P, DetectStats *st, DetectScratch *scr) {
    const auto tt = std::chrono::steady_clock::now();
    if (im.w < 8 || im.h < 8 || im.pix.size() != size_t(im.w) * im.h)
        return {};
    BackgroundMap M = background(im, P.bkground);
    std::vector<Detection> D = detect_pass(im, P, M, 0);
    if (!P.matched)
        return D;
    // matched filter: PSF σ from the median FWHM of well-measured sources (unsaturated, SNR > 20) from the 1st pass
    std::vector<double> fw;
    for (auto &d : D)
        if (d.converged && d.snr_tot > 20 && d.peak < P.saturation && d.fit < P.fit_rms)
            fw.push_back(d.fwhm);
    if (fw.size() < 5)
        return D;
    double ks = internal::median(fw) / FWHM_SIGMA;
    if (!(ks > 0.6) || ks > 8)
        return D;
    std::vector<float> F;
    double ng = 1;
    D = detect_pass(im, P, M, ks, &F, &ng);
    if (P.trails && scr) {
        scr->bkg.swap(M.bkg);
        scr->rms.swap(M.rms);
        scr->F.swap(F);
        scr->ks = ks;
        scr->ng = ng;
        scr->ok = true;
    } else if (P.trails)
        trail_search(im, P, M, ks, F, ng, D, st);
    if (st)
        st->t_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - tt).count();
    return D;
}

std::vector<Detection> detect_trails(const Image &im, const DetectParams &P, const DetectScratch &scr, const std::vector<Detection> &D,
                                     const std::vector<double> &mask, std::vector<char> &frag, DetectStats *st) {
    frag.assign(D.size(), 0);
    if (!scr.ok || !P.trails || scr.F.size() != size_t(im.w) * im.h)
        return {};
    const auto tt = std::chrono::steady_clock::now();
    BackgroundMap M;
    M.w = im.w;
    M.h = im.h;
    M.bkg = scr.bkg;
    M.rms = scr.rms;
    std::vector<Detection> T = D, out;
    trail_search(im, P, M, scr.ks, scr.F, scr.ng, T, st, &mask, &frag);
    for (auto &t : T)
        if (t.trail_len > 0)
            out.push_back(t);
    if (out.empty())
        frag.assign(D.size(), 0);
    T.swap(out);
    if (st)
        st->t_total = std::chrono::duration<double>(std::chrono::steady_clock::now() - tt).count();
    return T;
}

bool measure_at(const Image &im, double x, double y, const DetectParams &p, Detection &out) {
    if (im.w <= 0 || im.h <= 0 || !std::isfinite(x) || !std::isfinite(y))
        return false;
    // smoothed image (3×3 binomial) in a cutout around the click; pick the local maximum closest to the click
    // within the aperture radius (so a faint object next to a star is measured, not the star)
    const int r = std::max(1, p.aperture), R = r + 3;
    const int x0 = int(std::lround(x)), y0 = int(std::lround(y)), n = 2 * R + 1;
    std::vector<double> S(size_t(n) * n, NAN);
    std::vector<double> vals;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            int xx = x0 - R + i, yy = y0 - R + j;
            if (xx < 1 || yy < 1 || xx >= im.w - 1 || yy >= im.h - 1)
                continue;
            double s = 0;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                    s += (2 - std::abs(dx)) * (2 - std::abs(dy)) * im.pix[size_t(yy + dy) * im.w + xx + dx];
            S[size_t(j) * n + i] = s / 16;
            vals.push_back(s / 16);
        }
    if (vals.empty())
        return false;
    const double med = internal::median(vals);
    struct Mx {
        int x, y;
        double v;
    };
    std::vector<Mx> mx;
    for (int j = 1; j < n - 1; ++j)
        for (int i = 1; i < n - 1; ++i) {
            double v = S[size_t(j) * n + i];
            if (!std::isfinite(v))
                continue;
            bool m = true;
            for (int dy = -1; dy <= 1 && m; ++dy)
                for (int dx = -1; dx <= 1; ++dx)
                    if ((dx || dy) && !(S[size_t(j + dy) * n + i + dx] < v)) {
                        m = false;
                        break;
                    }
            if (m)
                mx.push_back({x0 - R + i, y0 - R + j, v});
        }
    // local maximum within the radius with the highest significance relative to the brightest, favoring proximity
    int sel = -1;
    double vmax = 0;
    for (auto &m : mx)
        vmax = std::max(vmax, m.v - med);
    double bd = 1e9;
    for (size_t k = 0; k < mx.size(); ++k) {
        double d = std::hypot(mx[k].x - x, mx[k].y - y);
        if (d > r || mx[k].v - med < 0.1 * vmax)
            continue;
        if (d < bd) {
            bd = d;
            sel = int(k);
        }
    }
    if (sel < 0)
        return false;
    // window limited by the distance to the neighboring significant maximum
    double dn = 1e9;
    for (size_t k = 0; k < mx.size(); ++k)
        if (int(k) != sel && mx[k].v - med > 0.25 * (mx[sel].v - med))
            dn = std::min(dn, std::hypot(double(mx[k].x - mx[sel].x), double(mx[k].y - mx[sel].y)));
    // brightest (raw) pixel next to the smoothed maximum
    int bx = mx[sel].x, by = mx[sel].y, best = -1;
    for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
            int v = im.at(mx[sel].x + dx, mx[sel].y + dy);
            if (v > best) {
                best = v;
                bx = mx[sel].x + dx;
                by = mx[sel].y + dy;
            }
        }
    Detection d;
    if (!measure_source(im, bx, by, p, d, std::max(2.0, 0.5 * dn + 0.5)))
        return false;
    if (d.flux <= 0 || d.snr < 2)
        return false;
    d.det_snr = d.snr;
    out = d;
    return true;
}

} // namespace astro
