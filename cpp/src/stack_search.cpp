// Faint-object search by stacking along velocities ("shift-and-stack", "synthetic/digital tracking") and
// confirmation of 2-detection tracks by forced photometry.
//
// References (independent implementation):
//   - Shao et al. (2014), ApJ 782, 1 — "synthetic tracking": sum of short frames shifted over a velocity grid;
//     the grid step is chosen so the shift error over the arc stays below ~1 FWHM.
//   - Zhai et al. (2014), ApJ 792, 60 — NEO detection by synthetic tracking; threshold on the stacked SNR.
//   - Heinze, Metchev & Trollo (2015), AJ 150, 125 — "digital tracking": number of independent hypotheses
//     (positions × velocities) and false-alarm threshold corrected by that number.
//   - Whidden et al. (2019), AJ 157, 119 — KBMOD: stacked likelihood with images ψ = Σ K·D/σ² and
//     φ = Σ K²/σ² (per-frame PSF matched filter), SNR = Σψ/√Σφ; static-source mask; light-curve filtering
//     (here: per-frame SNR clipping) against single-frame contributions.
//   - Kubica et al. (2007), Icarus 189, 151; Denneau et al. (2013), PASP 125, 357 — track confirmation by linear
//     prediction; chance-alignment control.
//   - Mesh background with sigma-clipped median and bilinear interpolation (Bertin & Arnouts 1996, A&AS 117, 393).
//   - Per-frame forced measurement with the frame PSF model (Stetson 1987, PASP 99, 191).
//
// False-alarm control: the threshold is the larger of (a) the Gaussian quantile for FA = stack_far per sequence given
// the effective number of hypotheses and (b) the maximum of the same stack with permuted frame times (the resulting
// "velocities" do not correspond to rectilinear motion: a real object does not add up, noise and subtraction
// residuals do), which measures the actual tail of the null distribution in the data itself.
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>

namespace astro {

bool plate_radec_to_pix(const Plate &p, double ra, double dec, double &x, double &y) {
    for (int it = 0; it < 10; ++it) {
        double r0, d0, r1, d1, r2, d2, u0, v0, u1, v1, u2, v2;
        p.pix_to_radec(x, y, r0, d0);
        p.pix_to_radec(x + 1, y, r1, d1);
        p.pix_to_radec(x, y + 1, r2, d2);
        radec_to_std(r0, d0, ra, dec, u0, v0);
        radec_to_std(r1, d1, ra, dec, u1, v1);
        radec_to_std(r2, d2, ra, dec, u2, v2);
        double a = u1 - u0, b = u2 - u0, c = v1 - v0, d = v2 - v0, det = a * d - b * c;
        if (!(std::fabs(det) > 0))
            return false;
        double dx = (d * u0 - b * v0) / det, dy = (-c * u0 + a * v0) / det;
        x -= dx;
        y -= dy;
        if (std::hypot(dx, dy) < 1e-4)
            break;
    }
    return std::isfinite(x) && std::isfinite(y);
}

namespace internal {

namespace {

// upper Gaussian tail and its inverse (bisection)
double gauss_sf(double k) {
    return 0.5 * std::erfc(k / std::sqrt(2.0));
}
double gauss_isf(double p) {
    double lo = 0, hi = 40;
    for (int i = 0; i < 100; ++i) {
        double m = 0.5 * (lo + hi);
        (gauss_sf(m) > p ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
}

double bilinear(const std::vector<float> &a, int w, int h, double x, double y) {
    if (!(x >= 0 && y >= 0 && x <= w - 1 && y <= h - 1))
        return 0;
    int x0 = std::min(int(x), w - 2), y0 = std::min(int(y), h - 2);
    double fx = x - x0, fy = y - y0;
    const float *p = &a[size_t(y0) * w + x0];
    return (1 - fy) * ((1 - fx) * p[0] + fx * p[1]) + fy * ((1 - fx) * p[w] + fx * p[w + 1]);
}

// separable convolution (1D kernels kx, ky with half-widths rx, ry), border = 0
void conv_sep(const std::vector<float> &in, int w, int h, const std::vector<double> &kx, const std::vector<double> &ky,
              std::vector<float> &out) {
    const int rx = int(kx.size() / 2), ry = int(ky.size() / 2);
    std::vector<float> t(in.size(), 0.f);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            double s = 0;
            for (int k = -rx; k <= rx; ++k) {
                int xx = x + k;
                if (xx >= 0 && xx < w)
                    s += kx[k + rx] * in[size_t(y) * w + xx];
            }
            t[size_t(y) * w + x] = float(s);
        }
    out.assign(in.size(), 0.f);
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            double s = 0;
            for (int k = -ry; k <= ry; ++k) {
                int yy = y + k;
                if (yy >= 0 && yy < h)
                    s += ky[k + ry] * t[size_t(yy) * w + x];
            }
            out[size_t(y) * w + x] = float(s);
        }
}

std::vector<double> gauss1d(double s, double power) {
    int r = std::max(2, int(std::ceil(3.5 * s)));
    std::vector<double> k(size_t(2 * r + 1));
    double sum = 0;
    for (int i = -r; i <= r; ++i)
        sum += k[size_t(i + r)] = std::exp(-0.5 * i * i / (s * s));
    for (auto &v : k)
        v = std::pow(v / sum, power);
    return k;
}

} // namespace

// Background (sigma-clipped median per 32 px cell, 3×3 median filter of the mesh, bilinear) and sky σ (MAD).
void background_mesh(const Image &im, std::vector<float> &B, std::vector<float> &S) {
    const int W = im.w, H = im.h, cell = 32;
    const int nx = std::max(1, (W + cell - 1) / cell), ny = std::max(1, (H + cell - 1) / cell);
    std::vector<double> mb(size_t(nx) * ny), ms(size_t(nx) * ny);
    std::vector<double> v;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            v.clear();
            for (int y = j * cell; y < std::min(H, (j + 1) * cell); ++y)
                for (int x = i * cell; x < std::min(W, (i + 1) * cell); ++x)
                    v.push_back(im.pix[size_t(y) * W + x]);
            double sg = 0;
            mb[size_t(j) * nx + i] = clipped_median(v, 3.0, 5, &sg, nullptr);
            ms[size_t(j) * nx + i] = std::max(sg, 0.5);
        }
    auto mfilt = [&](std::vector<double> &a) {
        std::vector<double> o(a.size()), t;
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                t.clear();
                for (int dj = -1; dj <= 1; ++dj)
                    for (int di = -1; di <= 1; ++di)
                        if (i + di >= 0 && i + di < nx && j + dj >= 0 && j + dj < ny)
                            t.push_back(a[size_t(j + dj) * nx + i + di]);
                o[size_t(j) * nx + i] = median(t);
            }
        a = o;
    };
    mfilt(mb);
    mfilt(ms);
    B.assign(size_t(W) * H, 0.f);
    S.assign(size_t(W) * H, 0.f);
    for (int y = 0; y < H; ++y) {
        double gy = std::clamp((y + 0.5) / cell - 0.5, 0.0, double(ny - 1));
        int j0 = std::min(int(gy), std::max(ny - 2, 0)), j1 = std::min(j0 + 1, ny - 1);
        double fy = gy - j0;
        for (int x = 0; x < W; ++x) {
            double gx = std::clamp((x + 0.5) / cell - 0.5, 0.0, double(nx - 1));
            int i0 = std::min(int(gx), std::max(nx - 2, 0)), i1 = std::min(i0 + 1, nx - 1);
            double fx = gx - i0;
            auto L = [&](const std::vector<double> &a) {
                return (1 - fy) * ((1 - fx) * a[size_t(j0) * nx + i0] + fx * a[size_t(j0) * nx + i1]) +
                       fy * ((1 - fx) * a[size_t(j1) * nx + i0] + fx * a[size_t(j1) * nx + i1]);
            };
            B[size_t(y) * W + x] = float(L(mb));
            S[size_t(y) * W + x] = float(L(ms));
        }
    }
}

bool build_mf(const Frame &f, const Catalog &cat, const std::vector<Pt> &extra, double saturation, int aperture, MFMap &M) {
    M = MFMap();
    const int W = f.im.w, H = f.im.h;
    if (W < 16 || H < 16)
        return false;
    M.w = W;
    M.h = H;
    double sc = 0.5 * (f.plate.scale_x() + f.plate.scale_y());
    if (!frame_psf(f, saturation, aperture, M.psf, M.apc)) {
        double s0 = std::max(0.6, (f.seeing > 0 ? f.seeing / sc : 2.5) / 2.3548);
        M.psf = PsfShape{s0, s0, 0, 1, 0};
        M.apc = 1;
    }
    M.fwhm = 2.3548 * 0.5 * (M.psf.sx + M.psf.sy);
    background_mesh(f.im, M.bkg, M.sig);
    std::vector<double> sv(M.sig.begin(), M.sig.end());
    M.noise = median(sv);
    // mask: static sources (fixed detections, catalog stars, saturated pixels) and `extra` (detections of already
    // linked tracks). Radius where the source's Gaussian wing falls to 0.3σ of sky: r = σ_PSF·√(2 ln(A/0.3σ)) + 1 px.
    std::vector<unsigned char> mask(size_t(W) * H, 0);
    const double sp = std::max(M.psf.sx, M.psf.sy);
    auto paint = [&](double x, double y, double amp) {
        double r = sp * std::sqrt(2 * std::log(std::max(amp / (0.3 * M.noise), 1.5))) + 1.0;
        r = std::max(r, 1.5 * M.fwhm);
        int R = int(std::ceil(r));
        for (int v = int(y) - R; v <= int(y) + R + 1; ++v)
            for (int u = int(x) - R; u <= int(x) + R + 1; ++u)
                if (u >= 0 && v >= 0 && u < W && v < H && (u - x) * (u - x) + (v - y) * (v - y) <= r * r)
                    mask[size_t(v) * W + u] = 1;
    };
    for (size_t i = 0; i < f.det.size(); ++i)
        if (i < f.typ.size() && f.typ[i] != 0)
            paint(f.det[i].x, f.det[i].y, std::max(f.det[i].peak - f.det[i].bkg, 0.0));
    for (auto &p : extra)
        paint(p.x, p.y, 1e4 * M.noise);
    const double norm = 1.0 / (2 * PI * M.psf.sx * M.psf.sy);
    for (size_t j = 0; j < cat.size(); ++j) {
        double x = (W - 1) / 2.0, y = (H - 1) / 2.0;
        double xi, eta;
        radec_to_std(cat.ra[j], cat.dec[j], f.plate.ra0, f.plate.dec0, xi, eta);
        if (std::fabs(xi) > 2 * W * sc || std::fabs(eta) > 2 * H * sc)
            continue;
        if (!plate_radec_to_pix(f.plate, cat.ra[j], cat.dec[j], x, y) || x < -20 || y < -20 || x > W + 20 || y > H + 20)
            continue;
        double F = std::isfinite(cat.mag[j])
                       ? std::pow(10.0, 0.4 * (f.plate.zp - f.plate.star_mag(cat.mag[j], cat.has_color() ? cat.color[j] : NAN)))
                       : 1e3 * M.noise;
        paint(x, y, F * norm);
    }
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            if (f.im.pix[size_t(y) * W + x] >= saturation)
                paint(x, y, 1e4 * M.noise);
    // ψ = K ⊗ (D·m/σ²), φ = K² ⊗ (m/σ²), u = ψ/√φ (Whidden et al. 2019, eq. 2–4)
    std::vector<float> a(size_t(W) * H), b(size_t(W) * H);
    for (size_t i = 0; i < a.size(); ++i) {
        double s2 = double(M.sig[i]) * M.sig[i];
        bool ok = !mask[i];
        a[i] = ok ? float((f.im.pix[i] - M.bkg[i]) / s2) : 0.f;
        b[i] = ok ? float(1.0 / s2) : 0.f;
    }
    std::vector<float> psi, phi;
    conv_sep(a, W, H, gauss1d(M.psf.sx, 1), gauss1d(M.psf.sy, 1), psi);
    conv_sep(b, W, H, gauss1d(M.psf.sx, 2), gauss1d(M.psf.sy, 2), phi);
    M.u.assign(a.size(), 0.f);
    std::vector<double> ph;
    for (size_t i = 0; i < a.size(); ++i)
        if (phi[i] > 0) {
            M.u[i] = float(psi[i] / std::sqrt(phi[i]));
            ph.push_back(std::sqrt(phi[i]));
        }
    // typical φ (unmasked part): √φ = frame sensitivity (SNR per unit flux)
    M.sqrt_phi = ph.empty() ? 0 : median(ph);
    M.ok = M.sqrt_phi > 0;
    return M.ok;
}

double MFMap::at(double x, double y) const {
    return bilinear(u, w, h, x, y);
}
double MFMap::at_s(double x, double y) const {
    return bilinear(us.empty() ? u : us, w, h, x, y);
}

bool forced_measure(const Frame &f, const MFMap &M, double x, double y, double gain, Detection &d, double &snr, int fixed_fallback) {
    snr = M.at(x, y);
    const int W = f.im.w, H = f.im.h;
    if (!(x >= 3 && y >= 3 && x < W - 3 && y < H - 3))
        return false;
    size_t i = size_t(int(std::lround(y))) * W + size_t(int(std::lround(x)));
    double bk = M.bkg[i], sg = M.sig[i];
    std::vector<SceneSource> src(1);
    double F0 = std::max(snr, 1.0) * sg * std::sqrt(4 * PI * M.psf.sx * M.psf.sy);
    src[0] = SceneSource{x, y, F0, false, true, false};
    int R = std::max(3, int(std::ceil(3 * std::max(M.psf.sx, M.psf.sy))) + 1);
    bool fixed = false;
    if (fixed_fallback == 2 || !fit_scene(f.im, x, y, R, M.psf, M.psf, src, bk, nullptr) || !(src[0].flux > 0) ||
        std::hypot(src[0].x - x, src[0].y - y) > std::max(1.0, 0.5 * M.fwhm)) {
        // free position did not converge (neighbor, noise): forced measurement with position FIXED at the track prediction —
        // only flux and background fitted (no Eddington bias in position: the position is that of the velocity hypothesis)
        if (fixed_fallback == 0)
            return false;
        bk = M.bkg[i];
        src[0] = SceneSource{x, y, F0, false, false, false};
        if (!fit_scene(f.im, x, y, R, M.psf, M.psf, src, bk, nullptr) || !(src[0].flux > 0))
            return false;
        fixed = true;
    }
    d = Detection();
    d.x = src[0].x;
    d.y = src[0].y;
    d.flux = src[0].flux * M.apc;
    d.fwhm = M.fwhm;
    d.sigma = sg;
    d.bkg = bk;
    d.peak = bk + src[0].flux / (2 * PI * M.psf.sx * M.psf.sy);
    d.snr = (d.peak - bk) / sg;
    d.snr_tot = std::max(M.at(d.x, d.y), snr);
    d.det_snr = snr;
    d.ix = int(std::lround(d.x));
    d.iy = int(std::lround(d.y));
    d.converged = !fixed;
    (void)gain;
    return true;
}

namespace {

// working frame: SNR map, affine pixel(ref) → pixel(frame), time and weight a_n = g_n·√φ_n
struct SF {
    int fi;
    double t, a;
    MFMap M;
    Affine A; // ref → frame
};

struct Cand {
    double x, y, vx, vy, s; // position at t_ref (ref px), velocity (ref px/min), SNR
};

} // namespace

std::vector<Track> stack_search(Reduction &R, const Catalog &cat, const Settings &s) {
    std::vector<Track> out;
    auto get = [&](const char *k, double def) { return s.str(k).empty() ? def : s.num(k); };
    const bool do_stack = get("stack_search", 1) != 0, do_pairs = get("pair_tracks", 0) != 0;
    if (!do_stack && !do_pairs)
        return out;
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<int> F;
    for (size_t k = 0; k < R.frames.size(); ++k)
        if (R.frames[k].solved && !R.frames[k].stacked)
            F.push_back(int(k));
    std::stable_sort(F.begin(), F.end(), [&](int a, int b) { return R.frames[a].im.jd_mid < R.frames[b].im.jd_mid; });
    const int N = int(F.size());
    if (N < 3)
        return out;
    const DetectParams dp = detect_params(s);
    const Frame &fr = R.frames[F[0]];
    const double tA = R.frames[F[0]].im.jd_mid, tB = R.frames[F[N - 1]].im.jd_mid;
    const double span = (tB - tA) * 1440, tref = 0.5 * span;
    if (!(span > 0))
        return out;
    // pixels of already linked tracks: masked (must not reappear or inflate the control)
    std::vector<std::vector<Pt>> trk(R.frames.size());
    for (auto &T : R.movers)
        for (auto &m : T.meas)
            if (m.frame >= 0 && m.frame < int(trk.size()))
                trk[m.frame].push_back({m.x, m.y});
    std::vector<SF> S(N);
    double zref = fr.plate.zp;
    for (int n = 0; n < N; ++n) {
        const Frame &f = R.frames[F[n]];
        S[n].fi = F[n];
        S[n].t = (f.im.jd_mid - tA) * 1440;
        if (!build_mf(f, cat, trk[F[n]], dp.saturation, dp.aperture, S[n].M))
            return out;
        S[n].a = std::pow(10.0, 0.4 * (f.plate.zp - zref)) * S[n].M.sqrt_phi;
        // affine ref → frame from a 5×5 grid of points mapped through the plate
        std::vector<Pt> P, Q;
        for (int j = 0; j < 5; ++j)
            for (int i = 0; i < 5; ++i) {
                double x = (fr.im.w - 1) * i / 4.0, y = (fr.im.h - 1) * j / 4.0, ra, dec;
                fr.plate.pix_to_radec(x, y, ra, dec);
                double u = x - fr.offx + f.offx, v = y - fr.offy + f.offy;
                if (plate_radec_to_pix(f.plate, ra, dec, u, v)) {
                    P.push_back({x, y});
                    Q.push_back({u, v});
                }
            }
        if (P.size() < 6 || !fit_affine(P, Q, S[n].A))
            return out;
    }
    // additional mask for stacking only: every detection with a positional counterpart (any brightness) in another
    // frame — a static source the compatible-brightness criterion missed (color across g/r/i filters, galaxy, variable
    // star) or an object already detectable on its own (handled by linking/pairs). Radius: PSF wing down to 0.3σ + 1 FWHM
    // (extent of the already convolved map).
    {
        auto inv = [](const Affine &A, double u, double v, double &x, double &y) {
            double det = A.a[1] * A.b[2] - A.a[2] * A.b[1], ex = u - A.a[0], ey = v - A.b[0];
            x = (A.b[2] * ex - A.a[2] * ey) / det;
            y = (-A.b[1] * ex + A.a[1] * ey) / det;
        };
        double seeing = 0;
        for (int n = 0; n < N; ++n)
            seeing = std::max(seeing, R.frames[F[n]].seeing);
        const double scr = 0.5 * (fr.plate.scale_x() + fr.plate.scale_y()), rc = std::max(3.0, seeing) / scr;
        std::vector<std::vector<Pt>> Q(N);
        for (int n = 0; n < N; ++n)
            for (auto &d : R.frames[F[n]].det) {
                Pt p;
                inv(S[n].A, d.x, d.y, p.x, p.y);
                Q[n].push_back(p);
            }
        for (int n = 0; n < N; ++n) {
            MFMap &M = S[n].M;
            M.us = M.u;
            const Frame &f = R.frames[F[n]];
            const double sp = std::max(M.psf.sx, M.psf.sy);
            for (size_t i = 0; i < Q[n].size(); ++i) {
                bool intrack = false;
                for (auto &q : trk[F[n]])
                    intrack |= std::hypot(q.x - f.det[i].x, q.y - f.det[i].y) < 1.0;
                (void)intrack;
                bool hit = false;
                for (int m = 0; m < N && !hit; ++m)
                    if (m != n)
                        for (auto &q : Q[m])
                            if (std::hypot(q.x - Q[n][i].x, q.y - Q[n][i].y) < rc) {
                                hit = true;
                                break;
                            }
                if (!hit)
                    continue;
                const Detection &d = f.det[i];
                double amp = std::max(d.peak - d.bkg, 0.0);
                double r = sp * std::sqrt(2 * std::log(std::max(amp / (0.3 * M.noise), 1.5))) + 1.0 + M.fwhm;
                int Rr = int(std::ceil(r));
                for (int v = int(d.y) - Rr; v <= int(d.y) + Rr + 1; ++v)
                    for (int u = int(d.x) - Rr; u <= int(d.x) + Rr + 1; ++u)
                        if (u >= 0 && v >= 0 && u < M.w && v < M.h && (u - d.x) * (u - d.x) + (v - d.y) * (v - d.y) <= r * r)
                            M.us[size_t(v) * M.w + u] = 0.f;
            }
        }
    }
    // static sky seen in the maps themselves (sources detection discarded: galaxies, blends, non-stellar sources):
    // pixels with combined static SNR Σa·u/‖a‖ > 4 and u > 1 in every frame (present from start to end of the arc;
    // a moving object is not in the same place in the first and last frame) are zeroed, with 1 FWHM dilation
    {
        double a2 = 0;
        for (auto &q : S)
            a2 += q.a * q.a;
        const double an = std::sqrt(a2);
        std::vector<std::vector<unsigned char>> sm(N);
        for (int n = 0; n < N; ++n) {
            const MFMap &M = S[n].M;
            const Affine &A = S[n].A;
            double det = A.a[1] * A.b[2] - A.a[2] * A.b[1];
            sm[n].assign(size_t(M.w) * M.h, 0);
            for (int y = 0; y < M.h; ++y)
                for (int x = 0; x < M.w; ++x) {
                    if (M.u[size_t(y) * M.w + x] <= 1.f)
                        continue;
                    double ex = x - A.a[0], ey = y - A.b[0];
                    double rx = (A.b[2] * ex - A.a[2] * ey) / det, ry = (-A.b[1] * ex + A.a[1] * ey) / det;
                    double z = 0;
                    bool all = true;
                    for (int m = 0; m < N && all; ++m) {
                        double u, v;
                        S[m].A.apply(rx, ry, u, v);
                        double um = m == n ? M.u[size_t(y) * M.w + x] : S[m].M.at(u, v);
                        all = um > 1.0;
                        z += S[m].a * um;
                    }
                    if (all && z / an > 4)
                        sm[n][size_t(y) * M.w + x] = 1;
                }
        }
        for (int n = 0; n < N; ++n) {
            MFMap &M = S[n].M;
            int Rr = int(std::ceil(M.fwhm));
            for (int y = 0; y < M.h; ++y)
                for (int x = 0; x < M.w; ++x)
                    if (sm[n][size_t(y) * M.w + x])
                        for (int dy = -Rr; dy <= Rr; ++dy)
                            for (int dx = -Rr; dx <= Rr; ++dx)
                                if (x + dx >= 0 && y + dy >= 0 && x + dx < M.w && y + dy < M.h && dx * dx + dy * dy <= Rr * Rr)
                                    M.us[size_t(y + dy) * M.w + x + dx] = 0.f;
        }
    }
    // extended regions (bright-star halo, galaxy, background residual): connected components of u > 3 with area
    // more than twice that expected for a point source of the same peak, A(u₀) = 4πσ²·ln(u₀/3) (the filtered map has
    // σ = √2·σ_PSF), are zeroed in the stacking map
    for (int n = 0; n < N; ++n) {
        MFMap &M = S[n].M;
        const int W = M.w, H = M.h;
        const double sp = 0.5 * (M.psf.sx + M.psf.sy);
        std::vector<int> lab(size_t(W) * H, 0), stk, comp;
        int nl = 0;
        for (int y0 = 0; y0 < H; ++y0)
            for (int x0 = 0; x0 < W; ++x0) {
                size_t i0 = size_t(y0) * W + x0;
                if (lab[i0] || !(M.us[i0] > 3.f))
                    continue;
                ++nl;
                comp.clear();
                stk.assign(1, int(i0));
                lab[i0] = nl;
                float pk = 0;
                while (!stk.empty()) {
                    int i = stk.back();
                    stk.pop_back();
                    comp.push_back(i);
                    pk = std::max(pk, M.us[size_t(i)]);
                    int x = i % W, y = i / W;
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            int xx = x + dx, yy = y + dy;
                            if (xx < 0 || yy < 0 || xx >= W || yy >= H)
                                continue;
                            size_t j = size_t(yy) * W + xx;
                            if (!lab[j] && M.us[j] > 3.f) {
                                lab[j] = nl;
                                stk.push_back(int(j));
                            }
                        }
                }
                double Aexp = 4 * PI * sp * sp * std::log(std::max(double(pk) / 3.0, 1.0));
                if (double(comp.size()) > 2 * Aexp + 10)
                    for (int i : comp) {
                        int x = i % W, y = i / W, Rr = int(std::ceil(M.fwhm));
                        for (int dy = -Rr; dy <= Rr; ++dy)
                            for (int dx = -Rr; dx <= Rr; ++dx)
                                if (x + dx >= 0 && y + dy >= 0 && x + dx < W && y + dy < H)
                                    M.us[size_t(y + dy) * W + x + dx] = 0.f;
                    }
            }
    }
    double asum = 0;
    for (auto &q : S)
        asum += q.a * q.a;
    const double anorm = std::sqrt(asum);
    const double clip = get("stack_clip", 6.0);
    const double sc = 0.5 * (fr.plate.scale_x() + fr.plate.scale_y());
    double fwhm_px = 0;
    for (auto &q : S)
        fwhm_px = std::max(fwhm_px, q.M.fwhm);
    // velocity limits (″/min): stack_min_rate/stack_max_rate; physical minimum: shift ≥ 1.5 FWHM over the arc
    // (below that an unmasked star adds up like the target)
    double vmin = std::max(get("stack_min_rate", get("min_rate", 0)), 1.5 * fwhm_px * sc / span);
    double vmax = std::min(get("stack_max_rate", 1.0), get("max_rate", 100.0));
    const double alpha = get("stack_far", 0.01), min_snr = get("stack_frame_snr", 3.0);
    // minimum number of measured frames (individual SNR ≥ stack_frame_snr) to emit a stacking track
    const int min_frames = std::max(2, int(get("stack_min_frames", 3)));
    // stack_fixed_pos (default 0; with 1, test_linking with 150 transients per frame yields 12 chance tracks): forced
    // fixed-position measurement at the prediction when the free-position fit fails
    const int fixed_fb = get("stack_fixed_pos", 0) != 0 ? 1 : 0;
    // stack_frame_snr_extra (default 0 = off): minimum SNR of extra fixed-position measurements after acceptance
    const double snr_extra = get("stack_frame_snr_extra", 0.0);
    char b[300];

    // exact stacked SNR for one hypothesis (x, y at t_ref in ref; v in ref px/min); tperm: permuted times (control)
    auto stat = [&](double x, double y, double vx, double vy, const std::vector<double> &tt) {
        double sum = 0;
        for (int n = 0; n < N; ++n) {
            double px = x + vx * (tt[n] - tref), py = y + vy * (tt[n] - tref), u, v;
            S[n].A.apply(px, py, u, v);
            sum += S[n].a * std::min(S[n].M.at_s(u, v), clip);
        }
        return sum / anorm;
    };
    std::vector<double> treal(N);
    for (int n = 0; n < N; ++n)
        treal[n] = S[n].t;
    // pattern-search refinement in (x, y, vx, vy)
    auto refine = [&](Cand c, const std::vector<double> &tt, double dv) {
        double st[4] = {0.5, 0.5, 0.5 * dv, 0.5 * dv};
        c.s = stat(c.x, c.y, c.vx, c.vy, tt);
        for (int h = 0; h < 6; ++h) {
            bool imp = true;
            for (int it = 0; it < 10 && imp; ++it) {
                imp = false;
                for (int k = 0; k < 4; ++k)
                    for (int sg = -1; sg <= 1; sg += 2) {
                        Cand d = c;
                        (k == 0 ? d.x : k == 1 ? d.y : k == 2 ? d.vx : d.vy) += sg * st[k];
                        d.s = stat(d.x, d.y, d.vx, d.vy, tt);
                        if (d.s > c.s + 1e-6) {
                            c = d;
                            imp = true;
                        }
                    }
            }
            for (double &q : st)
                q *= 0.5;
        }
        return c;
    };
    auto to_sky_v = [&](double vx, double vy, double &rate, double &pa) {
        double r0, d0, r1, d1, x0, y0, x1, y1;
        double xc = (fr.im.w - 1) / 2.0, yc = (fr.im.h - 1) / 2.0;
        fr.plate.pix_to_radec(xc, yc, r0, d0);
        fr.plate.pix_to_radec(xc + vx * 10, yc + vy * 10, r1, d1);
        radec_to_std(r0, d0, r0, d0, x0, y0);
        radec_to_std(r1, d1, r0, d0, x1, y1);
        rate = std::hypot(x1 - x0, y1 - y0) / 10;
        pa = std::fmod(std::atan2(x1 - x0, y1 - y0) * 180 / PI + 360, 360);
    };

    // predicted position in frame n and forced measurement; fixed-source test (vacant position) and track emission
    struct FM {
        int n;
        Detection d;
        double snr;
    };
    auto measure = [&](const Cand &c, std::vector<FM> &got) {
        got.clear();
        for (int n = 0; n < N; ++n) {
            double px = c.x + c.vx * (S[n].t - tref), py = c.y + c.vy * (S[n].t - tref), u, v;
            S[n].A.apply(px, py, u, v);
            FM m{n, Detection(), 0};
            if (forced_measure(R.frames[S[n].fi], S[n].M, u, v, dp.gain, m.d, m.snr, fixed_fb) && m.snr >= min_snr)
                got.push_back(m);
        }
    };
    // each measurement's position, seen in frames where the object has moved ≥ 2.5 FWHM away, must show no source:
    // combined SNR (weights a_m) < 3; no verifiable frame (unresolved motion) → reject
    auto vacated = [&](const Cand &c, const std::vector<FM> &got) {
        for (auto &g : got) {
            const Frame &fn = R.frames[S[g.n].fi];
            double ra, dec;
            fn.plate.pix_to_radec(g.d.x, g.d.y, ra, dec);
            double num = 0, den = 0;
            for (int m = 0; m < N; ++m) {
                if (m == g.n)
                    continue;
                const Frame &fm = R.frames[S[m].fi];
                double x = g.d.x - fn.offx + fm.offx, y = g.d.y - fn.offy + fm.offy;
                if (!plate_radec_to_pix(fm.plate, ra, dec, x, y))
                    continue;
                double px = c.x + c.vx * (S[m].t - tref), py = c.y + c.vy * (S[m].t - tref), u, v;
                S[m].A.apply(px, py, u, v);
                if (std::hypot(u - x, v - y) < 2.5 * std::max(S[m].M.fwhm, S[g.n].M.fwhm))
                    continue;
                // masked there (known star/static source within the mask radius): the measurement lies on a
                // fixed source — reject
                if (x >= 0 && y >= 0 && x <= fm.im.w - 1 && y <= fm.im.h - 1 && S[m].M.at(x, y) == 0.0)
                    return false;
                num += S[m].a * S[m].M.at(x, y);
                den += S[m].a * S[m].a;
            }
            if (den <= 0 || num / std::sqrt(den) >= 3.0)
                return false;
        }
        return true;
    };
    std::vector<Cand> accepted;
    auto duplicate = [&](const Cand &c) {
        for (auto &a : accepted) {
            bool same = true;
            for (int n = 0; n < N && same; ++n) {
                double t = S[n].t - tref;
                same = std::hypot(c.x + c.vx * t - a.x - a.vx * t, c.y + c.vy * t - a.y - a.vy * t) < 2 * fwhm_px;
            }
            if (same)
                return true;
        }
        return false;
    };
    auto emit = [&](const Cand &c, const std::vector<FM> &got, const char *how) {
        Track T;
        T.idx.assign(R.frames.size(), -1);
        for (auto &g : got) {
            if (!(g.d.flux > 0)) // non-positive flux: no defined magnitude — measurement discarded
                continue;
            Frame &f = R.frames[S[g.n].fi];
            f.det.push_back(g.d);
            if (f.typ.size() + 1 == f.det.size())
                f.typ.push_back(0);
            int k = int(f.det.size()) - 1;
            T.idx[S[g.n].fi] = k;
            Measurement me;
            me.frame = S[g.n].fi;
            me.x = g.d.x;
            me.y = g.d.y;
            f.plate.pix_to_radec(g.d.x, g.d.y, me.ra, me.dec);
            me.flux = g.d.flux;
            me.mag = f.plate.zp - 2.5 * std::log10(g.d.flux);
            me.snr = g.d.snr;
            me.snr_tot = g.d.snr_tot;
            me.fwhm_as = g.d.fwhm * 0.5 * (f.plate.scale_x() + f.plate.scale_y());
            me.jd = f.im.jd_mid;
            me.note = "";
            me.remark = g.d.converged ? how : std::string(how) + ", forced";
            T.meas.push_back(me);
        }
        to_sky_v(c.vx, c.vy, T.rate, T.pa);
        out.push_back(T);
        accepted.push_back(c);
    };

    // diagnostics: stack_probe = "ra dec jd dra_h ddec_h" (position and motion ″/h of a known object) → u per frame
    if (!s.str("stack_probe").empty()) {
        double pra, pde, pjd, pvr, pvd;
        if (std::sscanf(s.str("stack_probe").c_str(), "%lf %lf %lf %lf %lf", &pra, &pde, &pjd, &pvr, &pvd) == 5) {
            std::string l = "stack probe:";
            for (int n = 0; n < N; ++n) {
                const Frame &f = R.frames[S[n].fi];
                double dt = (f.im.jd_mid - pjd) * 24, ra, dec;
                std_to_radec(pvr * dt, pvd * dt, pra, pde, ra, dec);
                double x = f.im.w / 2.0, y = f.im.h / 2.0;
                plate_radec_to_pix(f.plate, ra, dec, x, y);
                double best = -99;
                for (int dy = -3; dy <= 3; ++dy)
                    for (int dx = -3; dx <= 3; ++dx)
                        best = std::max(best, S[n].M.at(x + dx, y + dy));
                std::snprintf(b, sizeof b, " [t=%.1f (%.0f,%.0f) u=%.1f max7=%.1f a=%.2g]", S[n].t, x, y, S[n].M.at(x, y), best, S[n].a);
                l += b;
            }
            R.log.push_back(l);
        }
    }
    // ---- 1. 2-detection tracks + forced confirmation in the remaining frames
    if (do_pairs) {
        // free detections (not fixed, not in tracks) in ref coordinates at the frame epoch
        std::vector<std::vector<int>> fi(N);
        std::vector<std::vector<Pt>> fp(N);
        for (int n = 0; n < N; ++n) {
            const Frame &f = R.frames[S[n].fi];
            // inverse of the affine ref → frame
            const Affine &A = S[n].A;
            double det = A.a[1] * A.b[2] - A.a[2] * A.b[1];
            for (size_t i = 0; i < f.det.size(); ++i) {
                if (i < f.typ.size() && f.typ[i] != 0)
                    continue;
                bool intrack = false;
                for (auto &q : trk[S[n].fi])
                    intrack |= std::hypot(q.x - f.det[i].x, q.y - f.det[i].y) < 1.0;
                if (intrack || S[n].M.at(f.det[i].x, f.det[i].y) == 0.0)
                    continue;
                double ex = f.det[i].x - A.a[0], ey = f.det[i].y - A.b[0];
                fi[n].push_back(int(i));
                fp[n].push_back({(A.b[2] * ex - A.a[2] * ey) / det, (-A.b[1] * ex + A.a[1] * ey) / det});
            }
        }
        if (get("stack_debug", 0) != 0) {
            std::string l = "stack dbg free:";
            for (int n = 0; n < N; ++n) {
                l += " " + std::to_string(fp[n].size());
                for (int i : fi[n]) {
                    char q[64];
                    std::snprintf(q, sizeof q, " (%.0f,%.0f)", R.frames[S[n].fi].det[i].x, R.frames[S[n].fi].det[i].y);
                    l += q;
                }
            }
            R.log.push_back(l);
        }
        int npairs = 0;
        for (int a = 0; a < N; ++a)
            for (int c = a + 1; c < N; ++c)
                npairs += (S[c].t - S[a].t) * vmax / sc > 2 * fwhm_px;
        double sig_pos = std::max(0.15 / sc, 0.1 * fwhm_px); // typical position σ of a faint detection (px)
        for (int a = 0; a < N; ++a)
            for (int c = a + 1; c < N; ++c) {
                double dt = S[c].t - S[a].t;
                if (dt * vmax / sc <= 2 * fwhm_px)
                    continue;
                double area = PI * (std::pow(vmax * dt / sc, 2) - std::pow(vmin * dt / sc, 2));
                double rho = double(fp[c].size()) / (double(fr.im.w) * fr.im.h);
                double Nb = std::max(1.0, std::min(double(fp[c].size()), rho * area));
                double Na = std::max(1.0, double(fp[a].size()));
                for (size_t i = 0; i < fp[a].size(); ++i)
                    for (size_t j = 0; j < fp[c].size(); ++j) {
                        double vx = (fp[c][j].x - fp[a][i].x) / dt, vy = (fp[c][j].y - fp[a][i].y) / dt;
                        double v = std::hypot(vx, vy) * sc;
                        if (v < vmin || v > vmax)
                            continue;
                        const Detection &da = R.frames[S[a].fi].det[fi[a][i]], &dc = R.frames[S[c].fi].det[fi[c][j]];
                        double ma = R.frames[S[a].fi].plate.zp - 2.5 * std::log10(std::max(da.flux, 1e-3));
                        double mc = R.frames[S[c].fi].plate.zp - 2.5 * std::log10(std::max(dc.flux, 1e-3));
                        if (std::fabs(ma - mc) > 1.5)
                            continue;
                        Cand cd{fp[a][i].x + vx * (tref - S[a].t), fp[a][i].y + vy * (tref - S[a].t), vx, vy, 0};
                        // combined forced SNR in other frames with an accurate prediction (σ_pred ≤ 0.35 FWHM). Expected
                        // brightness in each frame: u_m ≈ a_m·Â, Â fitted on the two detections; a forced value far above
                        // (> 2.5·expected + 3) indicates a fixed source at the predicted position, not the object.
                        const double ua = S[a].M.at(da.x, da.y), uc = S[c].M.at(dc.x, dc.y);
                        const double Ahat = (S[a].a * ua + S[c].a * uc) / (S[a].a * S[a].a + S[c].a * S[c].a);
                        bool contaminated = false;
                        std::vector<int> used_m;
                        double num = 0, den = 0;
                        for (int m = 0; m < N; ++m) {
                            // only frames from another epoch (object ≥ 1.5 FWHM away from both detections): a nearly
                            // simultaneous frame confirms the detection, not the velocity
                            double vpx = std::hypot(vx, vy);
                            if (m == a || m == c || std::fabs(S[m].t - S[a].t) * vpx < 1.5 * fwhm_px ||
                                std::fabs(S[m].t - S[c].t) * vpx < 1.5 * fwhm_px)
                                continue;
                            double sp = sig_pos * std::hypot(S[m].t - S[c].t, S[m].t - S[a].t) / dt;
                            if (sp > 0.35 * fwhm_px)
                                continue;
                            double px = cd.x + vx * (S[m].t - tref), py = cd.y + vy * (S[m].t - tref), u, w;
                            S[m].A.apply(px, py, u, w);
                            double uu = S[m].M.at(u, w);
                            if (uu == 0.0)
                                continue;
                            contaminated |= uu > 2.5 * S[m].a * Ahat + 3;
                            used_m.push_back(m);
                            num += S[m].a * std::min(uu, clip);
                            den += S[m].a * S[m].a;
                        }
                        if (den <= 0 || contaminated)
                            continue;
                        double sf = num / std::sqrt(den);
                        // chance alignment: E = N_a·N_b(velocity annulus)·P(forced SNR ≥ sf | noise)·number of frame pairs
                        double E = Na * Nb * gauss_sf(sf) * npairs;
                        if (get("stack_debug", 0) != 0) {
                            std::snprintf(b, sizeof b, "stack dbg pair %d-%d: forced SNR %.1f, N_a %.0f N_b %.1f E %.2g", S[a].fi + 1,
                                          S[c].fi + 1, sf, Na, Nb, E);
                            R.log.push_back(b);
                        }
                        if (E > get("false_alarm", 0.05))
                            continue;
                        cd.s = stat(cd.x, cd.y, cd.vx, cd.vy, treal);
                        if (duplicate(cd))
                            continue;
                        std::vector<FM> got;
                        measure(cd, got);
                        // both detections re-measured; a confirmation frame with SNR ≥ 5 at the prediction (position then
                        // known to ≲ FWHM/4) must accept the PSF fit — otherwise the position is contaminated
                        bool ha = false, hc = false, fitbad = false;
                        for (auto &g : got) {
                            ha |= g.n == a;
                            hc |= g.n == c;
                        }
                        for (int m : used_m) {
                            double px = cd.x + vx * (S[m].t - tref), py = cd.y + vy * (S[m].t - tref), u, w;
                            S[m].A.apply(px, py, u, w);
                            bool in = false;
                            for (auto &g : got)
                                in |= g.n == m;
                            fitbad |= S[m].M.at(u, w) >= std::max(min_snr, 5.0) && !in;
                        }
                        if (fitbad)
                            continue;
                        if (!ha || !hc || !vacated(cd, got))
                            continue;
                        std::snprintf(b, sizeof b,
                                      "2-detection track (frames %d and %d) confirmed: forced SNR %.1f, E %.2g, %zu measurements",
                                      S[a].fi + 1, S[c].fi + 1, sf, E, got.size());
                        R.log.push_back(b);
                        emit(cd, got, "forced phot");
                    }
            }
    }
    if (!do_stack)
        return out;

    // ---- 2. stacking along a velocity grid
    // grid step Δv = √2·FWHM/arc: the nearest node misses the shift by ≤ FWHM/2 at the arc ends
    // (Shao et al. 2014 §3); pixel grid with binning b (sampling of the PSF-smoothed map) if the cost exceeds
    // stack_budget operations
    const int Wr = fr.im.w, Hr = fr.im.h;
    const double vmin_px = vmin / sc, vmax_px = vmax / sc;
    if (!(vmax_px > vmin_px))
        return out;
    double dv = std::sqrt(2.0) * fwhm_px / span;
    const double budget = get("stack_budget", 4e9);
    int bin = 1;
    auto nvel = [&](double d) { return PI * (vmax_px * vmax_px - vmin_px * vmin_px) / (d * d); };
    while (nvel(dv) * double(Wr) * Hr / (bin * bin) * N > budget && bin < 4)
        ++bin;
    while (nvel(dv) * double(Wr) * Hr / (bin * bin) * N > budget)
        dv *= 1.25;
    const int w = (Wr + bin - 1) / bin, h = (Hr + bin - 1) / bin;
    // velocities → integer shifts (on the binned grid) per frame, without repetition
    auto shifts_for = [&](const std::vector<double> &tt, std::vector<std::vector<int>> &tup, std::vector<std::pair<double, double>> &vel) {
        std::map<std::vector<int>, int> seen;
        int nv = int(std::ceil(vmax_px / dv));
        for (int j = -nv; j <= nv; ++j)
            for (int i = -nv; i <= nv; ++i) {
                double vx = i * dv, vy = j * dv, vv = std::hypot(vx, vy);
                if (vv < vmin_px || vv > vmax_px)
                    continue;
                std::vector<int> k(size_t(2 * N));
                for (int n = 0; n < N; ++n) {
                    k[size_t(2 * n)] = int(std::lround(vx * (tt[n] - tref) / bin));
                    k[size_t(2 * n + 1)] = int(std::lround(vy * (tt[n] - tref) / bin));
                }
                if (seen.emplace(k, int(tup.size())).second) {
                    tup.push_back(k);
                    vel.push_back({vx, vy});
                }
            }
    };
    // maps U_n = (a_n/‖a‖)·min(u_n, clip) resampled onto the ref (binned) grid, with zero border P
    int P = int(std::ceil(vmax_px * std::max(tref, span - tref) / bin)) + 2;
    const int Wp = w + 2 * P, Hp = h + 2 * P;
    std::vector<std::vector<float>> U(N, std::vector<float>(size_t(Wp) * Hp, 0.f));
    for (int n = 0; n < N; ++n)
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                double u, v;
                S[n].A.apply(x * bin, y * bin, u, v);
                U[n][size_t(y + P) * Wp + x + P] = float(S[n].a / anorm * std::min(S[n].M.at_s(u, v), clip));
            }
    struct Peak {
        int x, y, t;
        float s;
    };
    auto search = [&](const std::vector<double> &tt, int keep, size_t &ntup) {
        std::vector<std::vector<int>> tup;
        std::vector<std::pair<double, double>> vel;
        shifts_for(tt, tup, vel);
        ntup = tup.size();
        std::vector<float> best(size_t(w) * h, -1e9f), acc(static_cast<size_t>(w));
        std::vector<int> bt(size_t(w) * h, -1);
        std::vector<const float *> rows(N);
        for (size_t q = 0; q < tup.size(); ++q) {
            const auto &k = tup[q];
            for (int y = 0; y < h; ++y) {
                for (int n = 0; n < N; ++n)
                    rows[n] = &U[n][size_t(y + P + k[size_t(2 * n + 1)]) * Wp + P + k[size_t(2 * n)]];
                std::fill(acc.begin(), acc.end(), 0.f);
                for (int n = 0; n < N; ++n) {
                    const float *r = rows[n];
                    for (int x = 0; x < w; ++x)
                        acc[x] += r[x];
                }
                float *bb = &best[size_t(y) * w];
                int *bi = &bt[size_t(y) * w];
                for (int x = 0; x < w; ++x)
                    if (acc[x] > bb[x]) {
                        bb[x] = acc[x];
                        bi[x] = int(q);
                    }
            }
        }
        // local maxima (radius ~ FWHM) in decreasing order
        std::vector<Peak> pk;
        int rr = std::max(1, int(std::lround(fwhm_px / bin)));
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                float v = best[size_t(y) * w + x];
                if (v < 3.5f)
                    continue;
                bool mx = true;
                for (int dy = -rr; dy <= rr && mx; ++dy)
                    for (int dx = -rr; dx <= rr && mx; ++dx) {
                        int xx = x + dx, yy = y + dy;
                        if ((dx || dy) && xx >= 0 && yy >= 0 && xx < w && yy < h && best[size_t(yy) * w + xx] > v)
                            mx = false;
                    }
                if (mx)
                    pk.push_back({x, y, bt[size_t(y) * w + x], v});
            }
        std::sort(pk.begin(), pk.end(), [](const Peak &p, const Peak &q) { return p.s > q.s; });
        if (int(pk.size()) > keep)
            pk.resize(size_t(keep));
        std::vector<Cand> cs;
        for (auto &p : pk)
            cs.push_back(refine(Cand{double(p.x * bin), double(p.y * bin), vel[size_t(p.t)].first, vel[size_t(p.t)].second, 0}, tt, dv));
        std::sort(cs.begin(), cs.end(), [](const Cand &p, const Cand &q) { return p.s > q.s; });
        return cs;
    };
    // control: permuted times (the permutation that shifts the epochs the most)
    std::vector<int> perm(N), bestp;
    std::iota(perm.begin(), perm.end(), 0);
    double bd = -1;
    int np = 0;
    do {
        double d = 0;
        for (int n = 0; n < N; ++n)
            d += std::fabs(treal[perm[n]] - treal[n]);
        if (d > bd + 1e-9) {
            bd = d;
            bestp = perm;
        }
    } while (std::next_permutation(perm.begin(), perm.end()) && ++np < 5040);
    std::vector<double> tperm(N);
    for (int n = 0; n < N; ++n)
        tperm[n] = treal[bestp[n]];
    size_t ntup = 0, ntup_c = 0;
    std::vector<Cand> ctrl = search(tperm, 10, ntup_c);
    std::vector<Cand> cands = search(treal, 60, ntup);
    // effective number of hypotheses: independent positions (useful area / FWHM²) × distinct velocities (Heinze et al. 2015)
    double usable = 0;
    for (int y = 0; y < Hr; y += 2)
        for (int x = 0; x < Wr; x += 2)
            usable += S[0].M.at_s(x, y) != 0.0 ? 4 : 0;
    double neff = std::max(1.0, usable / (fwhm_px * fwhm_px)) * double(std::max<size_t>(ntup, 1));
    double k_gauss = gauss_isf(alpha / neff), k_ctrl = ctrl.empty() ? 0 : ctrl[0].s;
    double thr = std::max({k_gauss, k_ctrl, get("stack_min_snr", 0.0)});
    if (get("stack_debug", 0) != 0)
        for (int which = 0; which < 2; ++which)
            for (size_t q = 0; q < std::min<size_t>(5, which ? cands.size() : ctrl.size()); ++q) {
                const Cand &c = which ? cands[q] : ctrl[q];
                const std::vector<double> &tt = which ? treal : tperm;
                std::string l = which ? "stack dbg real" : "stack dbg ctrl";
                for (int n = 0; n < N; ++n) {
                    double px = c.x + c.vx * (tt[n] - tref), py = c.y + c.vy * (tt[n] - tref), u, v;
                    S[n].A.apply(px, py, u, v);
                    std::snprintf(b, sizeof b, " [%.0f,%.0f u=%.1f s=%.1f]", u, v, S[n].M.at(u, v), S[n].M.at_s(u, v));
                    l += b;
                }
                std::snprintf(b, sizeof b, " S=%.1f", c.s);
                R.log.push_back(l + b);
            }
    // epochs: groups of frames in which the object (at the candidate velocity) moves < 1.5 FWHM. With < 3 epochs the
    // velocity is not verified (any pair of sources across 2 epochs defines a "track"); and the SNR without the
    // strongest epoch must stay ≥ stack_loo (3): a single epoch (a near-simultaneous pair, a transient) does not make a track.
    const double loo_min = get("stack_loo", 3.0);
    auto epochs_ok = [&](const Cand &c) {
        double vv = std::max(std::hypot(c.vx, c.vy), 1e-9), te = 1.5 * fwhm_px / vv;
        std::vector<int> ep(N, 0);
        int ne = 0;
        for (int n = 1; n < N; ++n)
            ep[n] = S[n].t - S[n - 1].t > te ? ++ne : ne;
        ++ne;
        if (ne < 3)
            return false;
        std::vector<double> num(ne, 0), den(ne, 0);
        for (int n = 0; n < N; ++n) {
            double px = c.x + c.vx * (S[n].t - tref), py = c.y + c.vy * (S[n].t - tref), u, v;
            S[n].A.apply(px, py, u, v);
            num[ep[n]] += S[n].a * std::min(S[n].M.at_s(u, v), clip);
            den[ep[n]] += S[n].a * S[n].a;
        }
        int kb = 0;
        for (int k = 1; k < ne; ++k)
            if (num[k] / std::sqrt(den[k]) > num[kb] / std::sqrt(den[kb]))
                kb = k;
        double sn = 0, sd = 0;
        for (int k = 0; k < ne; ++k)
            if (k != kb) {
                sn += num[k];
                sd += den[k];
            }
        return sd > 0 && sn / std::sqrt(sd) >= loo_min;
    };
    int nacc = 0;
    for (auto &c : cands) {
        if (c.s < thr)
            break;
        const bool sdbg = get("stack_debug", 0) != 0;
        if (duplicate(c) || !epochs_ok(c)) {
            if (sdbg) {
                std::snprintf(b, sizeof b, "stack dbg reject S=%.1f: %s", c.s, duplicate(c) ? "duplicate" : "epochs/LOO");
                R.log.push_back(b);
            }
            continue;
        }
        std::vector<FM> got;
        measure(c, got);
        if (int(got.size()) >= min_frames && vacated(c, got) && snr_extra > 0 && snr_extra < min_snr) {
            // extra measurements after acceptance: in the frames left out, forced flux with position FIXED at the
            // prediction if the matched-filter SNR there is ≥ stack_frame_snr_extra
            std::vector<char> in(N, 0);
            for (auto &g : got)
                in[g.n] = 1;
            for (int n = 0; n < N; ++n) {
                if (in[n])
                    continue;
                double px = c.x + c.vx * (S[n].t - tref), py = c.y + c.vy * (S[n].t - tref), u, v;
                S[n].A.apply(px, py, u, v);
                FM m{n, Detection(), 0};
                if (forced_measure(R.frames[S[n].fi], S[n].M, u, v, dp.gain, m.d, m.snr, 2) && m.snr >= snr_extra)
                    got.push_back(m);
            }
            std::sort(got.begin(), got.end(), [](const FM &a, const FM &b) { return a.n < b.n; });
        }
        // at most one fixed-position (fallback) measurement among those forming the track: the rest must converge
        int nfix = 0;
        for (auto &g : got)
            nfix += !g.d.converged && g.snr >= min_snr;
        if (nfix > 1) {
            std::vector<FM> keep;
            for (auto &g : got)
                if (g.d.converged || g.snr < min_snr)
                    keep.push_back(g);
            got.swap(keep);
        }
        if (int(got.size()) < min_frames || !vacated(c, got)) {
            if (sdbg) {
                std::snprintf(b, sizeof b, "stack dbg reject S=%.1f: %s (%zu measurements)", c.s,
                              int(got.size()) < min_frames ? "too few forced measurements" : "position not free", got.size());
                std::string l = b;
                for (int n = 0; n < N; ++n) {
                    double px = c.x + c.vx * (S[n].t - tref), py = c.y + c.vy * (S[n].t - tref), u, v, sn = 0;
                    S[n].A.apply(px, py, u, v);
                    Detection d;
                    bool ok = forced_measure(R.frames[S[n].fi], S[n].M, u, v, dp.gain, d, sn);
                    std::snprintf(b, sizeof b, " [q%d ok %d snr %.2f Δ %.2f px]", n + 1, int(ok), sn, std::hypot(d.x - u, d.y - v));
                    l += b;
                }
                R.log.push_back(l);
            }
            continue;
        }
        double rate, pa;
        to_sky_v(c.vx, c.vy, rate, pa);
        std::snprintf(b, sizeof b, "stack: candidate SNR %.1f, %.2f\"/min PA %.0f°, %zu forced measurements (SNR ≥ %.0f)", c.s, rate, pa,
                      got.size(), min_snr);
        R.log.push_back(b);
        emit(c, got, "synthetic tracking");
        ++nacc;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::snprintf(b, sizeof b,
                  "stack: %d frames, %zu rates (%.2f–%.2f\"/min, step %.3f\"/min, binning %d), N_eff %.2g, threshold %.2f "
                  "(gauss %.2f, controle %.2f), %d trilha(s), %.1f s",
                  N, ntup, vmin, vmax, dv * sc, bin, neff, thr, k_gauss, k_ctrl, nacc, secs);
    R.log.push_back(b);
    return out;
}

} // namespace internal
} // namespace astro
