// Trailed-PSF fitting for moving objects and joint fitting with neighbouring stars.
//
// References (independent implementation):
//   - Vereš et al. (2012), PASP 124, 1197, "Improved asteroid astrometry and photometry with trail fitting" — trail
//     model: Gaussian PSF convolved with a line segment of length L = v·t_exp along the motion direction; the
//     segment centre corresponds to mid-exposure. Here the Gaussian is elliptical, aligned with the detector axes
//     (non-square pixels), and the integral along the segment has a closed form using the error function.
//   - Stetson (1987), PASP 99, 191 — DAOPHOT: simultaneous fitting of neighbouring sources (groups) with the frame PSF,
//     target photometry from the fitted flux; PSF → aperture correction from the median over stars.
//   - Levenberg (1944), Marquardt (1963) — non-linear least squares; covariance matrix (JᵀJ)⁻¹·s².
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace astro {
namespace internal {

double psf_profile(const PsfShape &P, double u, double v) {
    const double sx2 = P.sx * P.sx, sy2 = P.sy * P.sy;
    const double ax = 1 / sx2, ay = 1 / sy2, norm = 1 / (2 * PI * std::sqrt(sx2 * sy2));
    if (P.len < 1e-3)
        return norm * std::exp(-0.5 * (ax * u * u + ay * v * v));
    // ∫_{−L/2}^{L/2} N(p − s·d; Σ) ds / L, with A = Σ⁻¹: a = dᵀAd, b = dᵀAp
    const double a = ax * P.dx * P.dx + ay * P.dy * P.dy, b = ax * P.dx * u + ay * P.dy * v;
    const double q = ax * u * u + ay * v * v - b * b / a, k = std::sqrt(a / 2), h = P.len / 2, s0 = b / a;
    return norm * std::exp(-0.5 * q) * std::sqrt(PI / (2 * a)) * (std::erf(k * (h - s0)) - std::erf(k * (-h - s0))) / P.len;
}

double psf_pixel(const PsfShape &P, double u, double v) {
    if (P.fast) {
        PsfShape Q = P;
        Q.sx = std::sqrt(P.sx * P.sx + 1.0 / 12);
        Q.sy = std::sqrt(P.sy * P.sy + 1.0 / 12);
        Q.fast = false;
        return psf_profile(Q, u, v);
    }
    // pixel integral by 3×3 midpoint rule (the point profile has σ ≈ 1 px: the box matters)
    double s = 0;
    for (int j = -1; j <= 1; ++j)
        for (int i = -1; i <= 1; ++i)
            s += psf_profile(P, u + i / 3.0, v + j / 3.0);
    return s / 9;
}

namespace {

// Solves M·x = r (n×n) by Gauss–Jordan with partial pivoting; inv receives M⁻¹ if non-null.
bool solve(std::vector<double> M, std::vector<double> r, int n, std::vector<double> &x, std::vector<double> *inv) {
    std::vector<double> I(size_t(n) * n, 0);
    for (int i = 0; i < n; ++i)
        I[size_t(i) * n + i] = 1;
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int i = c + 1; i < n; ++i)
            if (std::fabs(M[size_t(i) * n + c]) > std::fabs(M[size_t(piv) * n + c]))
                piv = i;
        if (std::fabs(M[size_t(piv) * n + c]) < 1e-300)
            return false;
        for (int k = 0; k < n; ++k) {
            std::swap(M[size_t(c) * n + k], M[size_t(piv) * n + k]);
            std::swap(I[size_t(c) * n + k], I[size_t(piv) * n + k]);
        }
        std::swap(r[c], r[piv]);
        double d = M[size_t(c) * n + c];
        for (int k = 0; k < n; ++k) {
            M[size_t(c) * n + k] /= d;
            I[size_t(c) * n + k] /= d;
        }
        r[c] /= d;
        for (int i = 0; i < n; ++i)
            if (i != c) {
                double f = M[size_t(i) * n + c];
                if (f == 0)
                    continue;
                for (int k = 0; k < n; ++k) {
                    M[size_t(i) * n + k] -= f * M[size_t(c) * n + k];
                    I[size_t(i) * n + k] -= f * I[size_t(c) * n + k];
                }
                r[i] -= f * r[c];
            }
    }
    x = r;
    if (inv)
        *inv = I;
    return true;
}

} // namespace

bool fit_scene(const Image &im, double cx, double cy, int R, const PsfShape &target, const PsfShape &star, std::vector<SceneSource> &src,
               double &bkg, double *rms_out, int max_iter, double capsule) {
    // square window around (cx, cy)
    std::vector<double> X, Y, Z;
    int x0 = int(std::lround(cx)), y0 = int(std::lround(cy));
    for (int y = y0 - R; y <= y0 + R; ++y)
        for (int x = x0 - R; x <= x0 + R; ++x)
            if (x >= 0 && y >= 0 && x < im.w && y < im.h) {
                if (capsule > 0) {
                    double u = x - cx, v = y - cy, h = target.len / 2 + capsule, t = std::clamp(u * target.dx + v * target.dy, -h, h);
                    if (std::hypot(u - t * target.dx, v - t * target.dy) > capsule)
                        continue;
                }
                X.push_back(x);
                Y.push_back(y);
                Z.push_back(im.pix[size_t(y) * im.w + x]);
            }
    const int npix = int(Z.size());
    // parameters: background, then per source: flux [, x, y if free] [, σx, σy if the shape is free]
    std::vector<double> p{bkg};
    struct Map {
        int f, x, y, sx, sy, len = -1, ang = -1;
    };
    std::vector<Map> mp;
    for (auto &s : src) {
        Map m{int(p.size()), -1, -1, -1, -1};
        p.push_back(s.flux);
        if (s.free_pos) {
            m.x = int(p.size());
            p.push_back(s.x);
            m.y = int(p.size());
            p.push_back(s.y);
        }
        if (s.free_shape) {
            const PsfShape &P = s.trail ? target : star;
            m.sx = int(p.size());
            p.push_back(P.sx);
            m.sy = int(p.size());
            p.push_back(P.sy);
        }
        if (s.trail && s.free_len) {
            m.len = int(p.size());
            p.push_back(target.len);
            m.ang = int(p.size());
            p.push_back(std::atan2(target.dy, target.dx));
        }
        mp.push_back(m);
    }
    const int np = int(p.size());
    if (npix < np + 4)
        return false;
    // extra rows: length prior (residual (L0 − L)·pix_noise/len_sigma, same scale as the pixels)
    std::vector<int> prior_of;
    for (size_t k = 0; k < src.size(); ++k)
        if (mp[k].len >= 0)
            prior_of.push_back(int(k));
    const int nrow = npix + int(prior_of.size());
    for (int k : prior_of)
        Z.push_back(target.len * src[k].pix_noise / std::max(src[k].len_sigma, 1e-3));
    auto model = [&](const std::vector<double> &q, int i) {
        if (i >= npix) {
            int k = prior_of[size_t(i - npix)];
            return q[mp[k].len] * src[k].pix_noise / std::max(src[k].len_sigma, 1e-3);
        }
        double v = q[0];
        for (size_t k = 0; k < src.size(); ++k) {
            PsfShape P = src[k].trail ? target : star;
            if (mp[k].sx >= 0) {
                P.sx = q[mp[k].sx];
                P.sy = q[mp[k].sy];
            }
            if (mp[k].len >= 0) {
                P.len = q[mp[k].len];
                P.dx = std::cos(q[mp[k].ang]);
                P.dy = std::sin(q[mp[k].ang]);
            }
            double sx = mp[k].x >= 0 ? q[mp[k].x] : src[k].x, sy = mp[k].y >= 0 ? q[mp[k].y] : src[k].y;
            v += q[mp[k].f] * psf_pixel(P, X[i] - sx, Y[i] - sy);
        }
        return v;
    };
    auto chi2 = [&](const std::vector<double> &q) {
        double c = 0;
        for (int i = 0; i < nrow; ++i) {
            double r = Z[i] - model(q, i);
            c += r * r;
        }
        return c;
    };
    auto valid = [&](const std::vector<double> &q) {
        for (size_t k = 0; k < src.size(); ++k) {
            if (mp[k].x >= 0 && (std::fabs(q[mp[k].x] - src[k].x) > R || std::fabs(q[mp[k].y] - src[k].y) > R))
                return false;
            if (mp[k].sx >= 0 && (q[mp[k].sx] < 0.2 || q[mp[k].sy] < 0.2 || q[mp[k].sx] > R || q[mp[k].sy] > R))
                return false;
            if (mp[k].len >= 0 && (q[mp[k].len] < 0 || q[mp[k].len] > 2 * R))
                return false;
        }
        return true;
    };
    double lambda = 1e-3, c0 = chi2(p);
    std::vector<double> JtJ, Jtr, J(size_t(nrow) * np);
    bool conv = false;
    for (int it = 0; it < max_iter && !conv; ++it) {
        // Jacobian by central differences (exact for the linear parameters)
        for (int a = 0; a < np; ++a) {
            double h = (a == 0 || std::find_if(mp.begin(), mp.end(), [&](const Map &m) { return m.f == a; }) != mp.end()) ? 1.0 : 1e-4;
            std::vector<double> qp = p, qm = p;
            qp[a] += h;
            qm[a] -= h;
            for (int i = 0; i < nrow; ++i)
                J[size_t(i) * np + a] = (model(qp, i) - model(qm, i)) / (2 * h);
        }
        JtJ.assign(size_t(np) * np, 0);
        Jtr.assign(np, 0);
        for (int i = 0; i < nrow; ++i) {
            double r = Z[i] - model(p, i);
            for (int a = 0; a < np; ++a) {
                Jtr[a] += J[size_t(i) * np + a] * r;
                for (int b = 0; b < np; ++b)
                    JtJ[size_t(a) * np + b] += J[size_t(i) * np + a] * J[size_t(i) * np + b];
            }
        }
        bool improved = false;
        for (int tries = 0; tries < 12 && !improved; ++tries) {
            std::vector<double> M = JtJ, d;
            for (int a = 0; a < np; ++a)
                M[size_t(a) * np + a] *= 1 + lambda;
            if (!solve(M, Jtr, np, d, nullptr))
                return false;
            std::vector<double> q = p;
            double step = 0;
            for (int a = 0; a < np; ++a)
                q[a] += d[a];
            for (int a = 0; a < np; ++a) // relative step for all parameters
                step = std::max(step, std::fabs(d[a]) / (std::fabs(p[a]) + 1e-3));
            double c1 = valid(q) ? chi2(q) : INFINITY;
            if (c1 <= c0) {
                conv = step < 1e-7 || (c0 - c1) < 1e-12 * c0;
                p = q;
                c0 = c1;
                lambda = std::max(lambda * 0.3, 1e-9);
                improved = true;
            } else
                lambda *= 10;
        }
        if (!improved)
            conv = true;
    }
    // covariance: s²·(JᵀJ)⁻¹ with s² = χ²/(n − p)
    std::vector<double> inv, dummy;
    bool have_cov = solve(JtJ, Jtr, np, dummy, &inv);
    double s2 = c0 / std::max(1, npix - np);
    bkg = p[0];
    for (size_t k = 0; k < src.size(); ++k) {
        src[k].flux = p[mp[k].f];
        src[k].eflux = have_cov ? std::sqrt(std::max(0.0, s2 * inv[size_t(mp[k].f) * np + mp[k].f])) : NAN;
        if (mp[k].x >= 0) {
            src[k].x = p[mp[k].x];
            src[k].y = p[mp[k].y];
            src[k].ex = have_cov ? std::sqrt(std::max(0.0, s2 * inv[size_t(mp[k].x) * np + mp[k].x])) : NAN;
            src[k].ey = have_cov ? std::sqrt(std::max(0.0, s2 * inv[size_t(mp[k].y) * np + mp[k].y])) : NAN;
            src[k].exy = have_cov ? s2 * inv[size_t(mp[k].x) * np + mp[k].y] : NAN;
            if (have_cov && src[k].gain > 0) {
                // sandwich: Cov(p) = A⁻¹ (Σᵢ Vᵢ Jᵢ Jᵢᵀ) A⁻¹; prior row with variance pix_noise² (its own scale).
                // Only the x and y rows of A⁻¹ are needed: u_a = Σ_b A⁻¹[a][b] J_ib.
                const int ix = mp[k].x, iy = mp[k].y;
                const double sky2 = src[k].pix_noise * src[k].pix_noise;
                double vxx = 0, vyy = 0, vxy = 0;
                for (int i = 0; i < nrow; ++i) {
                    double ux = 0, uy = 0;
                    for (int b = 0; b < np; ++b) {
                        ux += inv[size_t(ix) * np + b] * J[size_t(i) * np + b];
                        uy += inv[size_t(iy) * np + b] * J[size_t(i) * np + b];
                    }
                    const double V = i < npix ? sky2 + std::max(model(p, i) - p[0], 0.0) / src[k].gain : sky2;
                    vxx += V * ux * ux;
                    vyy += V * uy * uy;
                    vxy += V * ux * uy;
                }
                src[k].ex = std::sqrt(vxx);
                src[k].ey = std::sqrt(vyy);
                src[k].exy = vxy;
            }
        }
        if (mp[k].sx >= 0) {
            src[k].fit_sx = p[mp[k].sx];
            src[k].fit_sy = p[mp[k].sy];
        }
        if (mp[k].len >= 0) {
            src[k].fit_len = p[mp[k].len];
            src[k].fit_dx = std::cos(p[mp[k].ang]);
            src[k].fit_dy = std::sin(p[mp[k].ang]);
        }
    }
    if (rms_out)
        *rms_out = std::sqrt(c0 / npix);
    return conv;
}

// Frame PSF and aperture correction: fits free σx, σy on unsaturated, isolated references.
bool frame_psf(const Frame &f, double saturation, int aperture, PsfShape &psf, double &apcorr) {
    std::vector<double> vx, vy, vr;
    for (size_t k = 0; k < f.plate.ref_idx.size() && vx.size() < 60; ++k) {
        if (!f.plate.ref_ok[k])
            continue;
        const Detection &d = f.det[f.plate.ref_idx[k]];
        if (d.peak >= saturation || d.snr_tot < 20 || d.flux <= 0)
            continue;
        // isolated: no other detection within 3·FWHM + aperture
        bool iso = true;
        for (auto &e : f.det)
            if (&e != &d && std::hypot(e.x - d.x, e.y - d.y) < 3 * d.fwhm + aperture) {
                iso = false;
                break;
            }
        if (!iso)
            continue;
        double s0 = std::max(0.5, d.fwhm / 2.3548);
        PsfShape P{s0, s0, 0, 1, 0};
        std::vector<SceneSource> src(1);
        src[0] = SceneSource{d.x, d.y, d.flux, true, true, true};
        double b = d.bkg;
        int R = std::max(3, int(std::ceil(3 * s0)) + 1);
        if (!fit_scene(f.im, d.x, d.y, R, P, P, src, b, nullptr) || !(src[0].flux > 0))
            continue;
        vx.push_back(src[0].fit_sx);
        vy.push_back(src[0].fit_sy);
        vr.push_back(d.flux / src[0].flux);
    }
    if (vx.size() < 5)
        return false;
    psf = PsfShape{median(vx), median(vy), 0, 1, 0};
    // aperture correction with the fixed PSF (as it will be used on the target)
    std::vector<double> vr2;
    for (size_t k = 0; k < f.plate.ref_idx.size() && vr2.size() < 60; ++k) {
        if (!f.plate.ref_ok[k])
            continue;
        const Detection &d = f.det[f.plate.ref_idx[k]];
        if (d.peak >= saturation || d.snr_tot < 20 || d.flux <= 0)
            continue;
        std::vector<SceneSource> src(1);
        src[0] = SceneSource{d.x, d.y, d.flux, true, true, false};
        double b = d.bkg;
        int R = std::max(3, int(std::ceil(3 * std::max(psf.sx, psf.sy))) + 1);
        if (fit_scene(f.im, d.x, d.y, R, psf, psf, src, b, nullptr) && src[0].flux > 0)
            vr2.push_back(d.flux / src[0].flux);
    }
    apcorr = vr2.size() >= 5 ? clipped_median(vr2, 3.0, 10) : median(vr);
    return std::isfinite(apcorr) && apcorr > 0;
}

// Plate standard coordinates (″) at pixel (x, y).
static void pix_to_std(const Plate &p, double x, double y, double &xi, double &eta) {
    auto t = p.terms(x, y);
    xi = eta = 0;
    for (size_t k = 0; k < t.size() && k < p.cx.size(); ++k) {
        xi += t[k] * p.cx[k];
        eta += t[k] * p.cy[k];
    }
}

// Plate inverse by Newton (numerical Jacobian); false if it does not converge.
static bool std_to_pix(const Plate &p, double xi, double eta, double &x, double &y) {
    for (int it = 0; it < 8; ++it) {
        double a, b, ax, bx, ay, by;
        pix_to_std(p, x, y, a, b);
        pix_to_std(p, x + 0.5, y, ax, bx);
        pix_to_std(p, x, y + 0.5, ay, by);
        double j11 = (ax - a) / 0.5, j21 = (bx - b) / 0.5, j12 = (ay - a) / 0.5, j22 = (by - b) / 0.5, det = j11 * j22 - j12 * j21;
        if (std::fabs(det) < 1e-12)
            return false;
        double ex = xi - a, ey = eta - b, dx = (j22 * ex - j12 * ey) / det, dy = (-j21 * ex + j11 * ey) / det;
        x += dx;
        y += dy;
        if (std::fabs(dx) + std::fabs(dy) < 1e-6)
            return true;
    }
    return false;
}

void refine_movers(Reduction &R, const Catalog &cat, const Settings &s) {
    auto get = [&](const char *k, double def) { return s.str(k).empty() ? def : s.num(k); };
    // trail_fit: −1 automatic (default: trail when L = v·t_exp ≥ 1.5·FWHM), 0 off, 1 always trail, 2 star PSF
    const int mode = int(get("trail_fit", -1));
    const bool deblend = get("trail_deblend", 1) != 0, phot = get("psf_phot", 1) != 0;
    if ((mode == 0 && !phot) || R.movers.empty())
        return;
    const DetectParams dp = detect_params(s);
    std::vector<int> have(R.frames.size(), 0); // 0 not computed, 1 ok, −1 failed
    std::vector<PsfShape> psf(R.frames.size());
    std::vector<double> apc(R.frames.size(), 1);
    char b[256];
    for (auto &T : R.movers) {
        int nok = 0;
        for (auto &m : T.meas) {
            if (m.frame < 0 || m.frame >= int(R.frames.size()))
                continue;
            Frame &f = R.frames[m.frame];
            if (!f.solved || f.stacked)
                continue;
            if (have[m.frame] == 0)
                have[m.frame] = frame_psf(f, dp.saturation, dp.aperture, psf[m.frame], apc[m.frame]) ? 1 : -1;
            if (have[m.frame] < 0)
                continue;
            // trail direction and length in pixels: sky velocity → pixels via the plate Jacobian
            PsfShape tr = psf[m.frame];
            if (mode != 2) {
                double a, c, ax, cx, ay, cy;
                pix_to_std(f.plate, m.x, m.y, a, c);
                pix_to_std(f.plate, m.x + 1, m.y, ax, cx);
                pix_to_std(f.plate, m.x, m.y + 1, ay, cy);
                double j11 = ax - a, j21 = cx - c, j12 = ay - a, j22 = cy - c, det = j11 * j22 - j12 * j21;
                double vxi = T.rate * std::sin(T.pa * PI / 180), veta = T.rate * std::cos(T.pa * PI / 180); // ″/min
                double px = (j22 * vxi - j12 * veta) / det, py = (-j21 * vxi + j11 * veta) / det;           // px/min
                double vp = std::hypot(px, py);
                tr.len = vp * f.im.exptime / 60;
                if (vp > 0) {
                    tr.dx = px / vp;
                    tr.dy = py / vp;
                }
            }
            const double smax = std::max(tr.sx, tr.sy), fw = 2.3548 * smax;
            // automatic mode: trail-model positions only when the trail is resolved (L ≥ 1.5·FWHM); below that
            // the point centroid is just as good and the fit serves photometry only
            const bool use_trail = mode == 1 || (mode < 0 && tr.len >= 1.5 * fw);
            const int W = int(std::ceil(3 * smax + (use_trail ? 0.6 : 0.5) * tr.len)) + 1;
            std::vector<SceneSource> src(1);
            src[0] = SceneSource{m.x, m.y, std::max(m.flux, 1.0) / apc[m.frame], true, true, false};
            if (use_trail && tr.len >= 1.5 * fw) {
                // free length and direction, length with prior L = v·t_exp (σ = 0.5 px ⊕ 15 %: linking velocity,
                // shutter and nominal exposure time)
                src[0].free_len = true;
                src[0].len_sigma = std::hypot(0.5, 0.15 * tr.len);
                const int id = T.idx[m.frame];
                src[0].pix_noise = id >= 0 && id < int(f.det.size()) && f.det[id].sigma > 0 ? f.det[id].sigma : 1.0;
            }
            if (use_trail && get("trail_cov_poisson", 1) != 0) {
                const int id = T.idx[m.frame];
                if (id >= 0 && id < int(f.det.size()) && f.det[id].sigma > 0) {
                    src[0].pix_noise = f.det[id].sigma;
                    src[0].gain = R.gain > 0 ? R.gain : get("gain", 1.0);
                }
            }
            // neighbours: catalogue stars (position fixed by the plate) and stationary detections with no catalogue match
            int nnb = 0;
            if (deblend) {
                const double rnb = std::max(2 * fw + tr.len / 2, double(W));
                std::vector<std::pair<double, double>> nb;
                double xi0, eta0;
                pix_to_std(f.plate, m.x, m.y, xi0, eta0);
                for (size_t i = 0; i < cat.size(); ++i) {
                    double xi, eta;
                    radec_to_std(cat.ra[i], cat.dec[i], f.plate.ra0, f.plate.dec0, xi, eta);
                    if (std::hypot(xi - xi0, eta - eta0) > 2 * rnb * std::max(f.plate.scale_x(), f.plate.scale_y()))
                        continue;
                    double x = m.x, y = m.y;
                    if (std_to_pix(f.plate, xi, eta, x, y) && std::hypot(x - m.x, y - m.y) < rnb)
                        nb.push_back({x, y});
                }
                for (size_t i = 0; i < f.det.size(); ++i) {
                    if (i >= f.typ.size() || f.typ[i] == 0 || f.typ[i] == 3 || std::hypot(f.det[i].x - m.x, f.det[i].y - m.y) >= rnb)
                        continue;
                    bool dup = false;
                    for (auto &q : nb)
                        dup |= std::hypot(q.first - f.det[i].x, q.second - f.det[i].y) < 0.5 * fw;
                    if (!dup)
                        nb.push_back({f.det[i].x, f.det[i].y});
                }
                for (auto &q : nb) {
                    // star almost on top of the target: degenerate positions, not separable
                    if (std::hypot(q.first - m.x, q.second - m.y) < 0.5 * fw)
                        continue;
                    src.push_back(SceneSource{q.first, q.second, 0.1 * src[0].flux, false, false, false});
                    ++nnb;
                }
            }
            double bk = f.det.empty() ? 0 : f.det[T.idx[m.frame] >= 0 ? T.idx[m.frame] : 0].bkg;
            const SceneSource ini = src[0];
            if (!fit_scene(f.im, m.x, m.y, W, tr, psf[m.frame], src, bk, nullptr) || !(src[0].flux > 0) ||
                std::hypot(src[0].x - ini.x, src[0].y - ini.y) > std::max({1.0, 0.5 * fw, use_trail ? 0.25 * tr.len : 0.0}))
                continue;
            if (mode > 0 || use_trail) {
                m.x = src[0].x;
                m.y = src[0].y;
                f.plate.pix_to_radec(m.x, m.y, m.ra, m.dec);
                ++nok;
                // centroid covariance from the fit itself, s²·(JᵀJ)⁻¹ (anisotropic on a trail: the error along it
                // is ~L/√12 larger than across); propagated to the ADES σ via J·C·Jᵀ (report.cpp)
                if (use_trail && std::isfinite(src[0].ex) && std::isfinite(src[0].ey) && std::isfinite(src[0].exy) && src[0].ex > 0 &&
                    src[0].ey > 0 && get("trail_cov", 1) != 0) {
                    m.cxx = src[0].ex * src[0].ex;
                    m.cyy = src[0].ey * src[0].ey;
                    m.cxy = src[0].exy;
                }
                m.trail_len_px = src[0].free_len ? src[0].fit_len : tr.len;
                m.trail_dx = src[0].free_len ? src[0].fit_dx : tr.dx;
                m.trail_dy = src[0].free_len ? src[0].fit_dy : tr.dy;
            }
            if (phot) {
                m.flux = src[0].flux * apc[m.frame];
                m.mag = f.plate.zp - 2.5 * std::log10(m.flux);
            }
            std::snprintf(b, sizeof b, "trail image %d: L %.2f px%s, σ %.2f/%.2f px, %d neighbour(s), Δ %+.2f %+.2f px, mag %.2f",
                          m.frame + 1, src[0].free_len ? src[0].fit_len : tr.len, use_trail ? " (trail position)" : "", tr.sx, tr.sy, nnb,
                          src[0].x - ini.x, src[0].y - ini.y, f.plate.zp - 2.5 * std::log10(src[0].flux * apc[m.frame]));
            R.log.push_back(b);
        }
        // velocity and PA recomputed from the remeasured positions (linear fit of ξ(t), η(t))
        if (nok > 0 && T.meas.size() >= 2) {
            double ra0 = T.meas[0].ra, dec0 = T.meas[0].dec, t0 = T.meas[0].jd;
            double n = 0, st = 0, stt = 0, sx = 0, sy = 0, stx = 0, sty = 0;
            for (auto &m : T.meas) {
                double xi, eta, t = (m.jd - t0) * 1440;
                radec_to_std(m.ra, m.dec, ra0, dec0, xi, eta);
                n += 1;
                st += t;
                stt += t * t;
                sx += xi;
                sy += eta;
                stx += t * xi;
                sty += t * eta;
            }
            double den = n * stt - st * st;
            if (den > 0) {
                double vx = (n * stx - st * sx) / den, vy = (n * sty - st * sy) / den;
                T.rate = std::hypot(vx, vy);
                T.pa = std::fmod(std::atan2(vx, vy) * 180 / PI + 360, 360);
            }
        }
    }
}

} // namespace internal
} // namespace astro
