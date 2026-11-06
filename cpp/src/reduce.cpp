// Plate reduction (least-squares plate constants) and differential photometry against the catalog.
//
// References (independent implementation):
//   - Green (1985), "Spherical Astronomy", ch. 13 — plate constants: standard coordinates (ξ, η) as a polynomial
//     in (x, y), least-squares fitted to the reference stars.
//   - König (1962), "Astrometry with astrographs", in Astronomical Techniques — higher-order plate models.
//   - Iterative sigma rejection (3σ until convergence) and median/MAD: Hoaglin, Mosteller & Tukey (1983).
//   - Golub & Van Loan, "Matrix Computations", §5.3 — QR least squares (Householder).
//   - Photometric zero point: sigma-clipped median of (m_cat − m_inst), m_inst = −2.5·log10(flux)
//     (e.g. Howell 2006, "Handbook of CCD Astronomy", ch. 5).
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>

namespace astro {
namespace internal {

double median(std::vector<double> &v) {
    if (v.empty())
        return NAN;
    size_t h = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + h, v.end());
    double m = v[h];
    if (v.size() % 2 == 0) {
        double lo = *std::max_element(v.begin(), v.begin() + h);
        m = 0.5 * (m + lo);
    }
    return m;
}

double mad_sigma(std::vector<double> v) {
    double m = median(v);
    for (auto &x : v)
        x = std::fabs(x - m);
    return 1.4826 * median(v);
}

double clipped_median(std::vector<double> v, double k, int iters, double *sigma, int *n) {
    double m = NAN, s = 0;
    for (int it = 0; it < iters && !v.empty(); ++it) {
        std::vector<double> t = v;
        m = median(t);
        s = mad_sigma(v);
        if (!(s > 0)) {
            // degenerate distribution (nearly all values equal): plain standard deviation
            double a = 0;
            for (double x : v)
                a += (x - m) * (x - m);
            s = std::sqrt(a / v.size());
            break;
        }
        std::vector<double> keep;
        for (double x : v)
            if (std::fabs(x - m) <= k * s)
                keep.push_back(x);
        if (keep.size() == v.size())
            break;
        v.swap(keep);
    }
    if (!v.empty()) {
        std::vector<double> t = v;
        m = median(t);
        double a = 0;
        for (double x : v)
            a += (x - m) * (x - m);
        s = v.size() > 1 ? std::sqrt(a / (v.size() - 1)) : s;
    }
    if (sigma)
        *sigma = s;
    if (n)
        *n = int(v.size());
    return m;
}

bool lstsq(const std::vector<double> &A0, int m, int n, const std::vector<double> &b0, std::vector<double> &c) {
    if (m < n || n <= 0)
        return false;
    // columns normalized for stability; Householder QR on A (column-major working copy)
    std::vector<double> A(A0), b(b0), sc(n, 1.0);
    for (int j = 0; j < n; ++j) {
        double s = 0;
        for (int i = 0; i < m; ++i)
            s += A[size_t(i) * n + j] * A[size_t(i) * n + j];
        s = std::sqrt(s);
        if (s <= 0)
            return false;
        sc[j] = s;
        for (int i = 0; i < m; ++i)
            A[size_t(i) * n + j] /= s;
    }
    for (int k = 0; k < n; ++k) {
        double nrm = 0;
        for (int i = k; i < m; ++i)
            nrm += A[size_t(i) * n + k] * A[size_t(i) * n + k];
        nrm = std::sqrt(nrm);
        if (nrm < 1e-12)
            return false;
        double alpha = A[size_t(k) * n + k] > 0 ? -nrm : nrm;
        std::vector<double> v(m - k);
        for (int i = k; i < m; ++i)
            v[i - k] = A[size_t(i) * n + k];
        v[0] -= alpha;
        double vn = 0;
        for (double x : v)
            vn += x * x;
        if (vn < 1e-300)
            continue;
        for (int j = k; j < n; ++j) {
            double d = 0;
            for (int i = k; i < m; ++i)
                d += v[i - k] * A[size_t(i) * n + j];
            d = 2 * d / vn;
            for (int i = k; i < m; ++i)
                A[size_t(i) * n + j] -= d * v[i - k];
        }
        double d = 0;
        for (int i = k; i < m; ++i)
            d += v[i - k] * b[i];
        d = 2 * d / vn;
        for (int i = k; i < m; ++i)
            b[i] -= d * v[i - k];
    }
    c.assign(n, 0);
    for (int k = n - 1; k >= 0; --k) {
        double s = b[k];
        for (int j = k + 1; j < n; ++j)
            s -= A[size_t(k) * n + j] * c[j];
        if (std::fabs(A[size_t(k) * n + k]) < 1e-12)
            return false;
        c[k] = s / A[size_t(k) * n + k];
    }
    for (int j = 0; j < n; ++j)
        c[j] /= sc[j];
    return true;
}

double centroid_sigma(const Detection &d, double scale, double gain) {
    double g = gain > 0 ? gain : 1.0;
    double sp = std::max(d.fwhm, 0.5) / 2.3548200450309493; // σ_PSF (px)
    double Ne = std::max(d.flux, 1e-9) * g, B = d.sigma * d.sigma * g * g;
    double v = sp * sp / Ne * (1 + 8 * PI * sp * sp * B / Ne);
    return std::sqrt(v) * scale;
}

bool invert(std::vector<double> &A, int n) {
    std::vector<double> I(size_t(n) * n, 0.0);
    for (int i = 0; i < n; ++i)
        I[size_t(i) * n + i] = 1;
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r)
            if (std::fabs(A[size_t(r) * n + c]) > std::fabs(A[size_t(piv) * n + c]))
                piv = r;
        if (!(std::fabs(A[size_t(piv) * n + c]) > 1e-300))
            return false;
        for (int k = 0; k < n; ++k) {
            std::swap(A[size_t(c) * n + k], A[size_t(piv) * n + k]);
            std::swap(I[size_t(c) * n + k], I[size_t(piv) * n + k]);
        }
        double d = A[size_t(c) * n + c];
        for (int k = 0; k < n; ++k) {
            A[size_t(c) * n + k] /= d;
            I[size_t(c) * n + k] /= d;
        }
        for (int r = 0; r < n; ++r)
            if (r != c) {
                double f = A[size_t(r) * n + c];
                if (f == 0)
                    continue;
                for (int k = 0; k < n; ++k) {
                    A[size_t(r) * n + k] -= f * A[size_t(c) * n + k];
                    I[size_t(r) * n + k] -= f * I[size_t(c) * n + k];
                }
            }
    }
    A.swap(I);
    return true;
}

// Plate-constant fit with iterative nearest-neighbor matching (shrinking radius), weighted least squares
// (weights 1/σ²) and 3σ rejection; fixed order or chosen by k-fold cross-validation (Stone 1974, JRSS B 36,
// 111; Hastie, Tibshirani & Friedman 2009, "The Elements of Statistical Learning", §7.10), which also gives the
// out-of-sample rms reported in the ADES errors.
bool fit_plate(Frame &f, const std::vector<Pt> &cat_std, const Catalog &cat, const SolveCtx &c, const Affine &T0) {
    Plate &p = f.plate;
    p.ra0 = c.ra0;
    p.dec0 = c.dec0;
    p.order = 1;
    p.cx = {T0.a[0], T0.a[1], T0.a[2]};
    p.cy = {T0.b[0], T0.b[1], T0.b[2]};
    const double sc = std::sqrt(std::fabs(T0.a[1] * T0.b[2] - T0.a[2] * T0.b[1]));
    const double rads[] = {std::max(4 * c.match_rad, 4 * sc), std::max(2 * c.match_rad, 2 * sc), c.match_rad, c.match_rad};
    const bool has_err = cat.err.size() == cat.size();
    std::vector<int> di, ci;
    std::vector<char> ok;
    std::vector<double> wt; // 1/σ² weights per pair
    // color term: extra column (BP − RP − c☉) in the system (0 for references without color), only with ≥ 10 references
    // with color and color spread ≥ 0.2 mag (Pier et al. 2003, AJ 125, 1559, §3: differential chromatic refraction and
    // lateral color); the target (solar color) is evaluated with a zero term
    std::vector<double> ck;
    bool use_col = false;
    auto ncols = [&](int order) { return plate_nterms(order) + (use_col ? 1 : 0); };
    auto design = [&](int order, size_t k) {
        auto t = p.terms(f.det[di[k]].x, f.det[di[k]].y);
        if (use_col)
            t.push_back(ck[k]);
        return t;
    };
    // weighted fit of order `order` on pairs flagged in `use`; returns residuals of all pairs
    auto solve_order = [&](int order, const std::vector<char> &use, std::vector<double> &qx, std::vector<double> &qy) {
        p.order = order;
        int nt = ncols(order), m = 0;
        std::vector<double> A, bx, by;
        for (size_t k = 0; k < di.size(); ++k)
            if (use[k]) {
                auto t = design(order, k);
                double sw = std::sqrt(wt[k]);
                for (int q = 0; q < nt; ++q)
                    A.push_back(t[q] * sw);
                bx.push_back(cat_std[ci[k]].x * sw);
                by.push_back(cat_std[ci[k]].y * sw);
                ++m;
            }
        return m >= nt + 1 && lstsq(A, m, nt, bx, qx) && lstsq(A, m, nt, by, qy);
    };
    auto resid = [&](int order, const std::vector<double> &qx, const std::vector<double> &qy, size_t k, double &rx, double &ry) {
        p.order = order;
        auto t = design(order, k);
        double a = 0, b = 0;
        for (int q = 0; q < ncols(order); ++q) {
            a += t[q] * qx[q];
            b += t[q] * qy[q];
        }
        rx = a - cat_std[ci[k]].x;
        ry = b - cat_std[ci[k]].y;
    };
    // k-fold cross-validation (k = 5) on good pairs: per-axis out-of-sample rms (″)
    auto cross_val = [&](int order, double &cvx, double &cvy) {
        const int K = 5;
        double sx = 0, sy = 0, sw = 0;
        int n = 0, j = 0;
        std::vector<int> fold(di.size(), -1);
        for (size_t k = 0; k < di.size(); ++k)
            if (ok[k])
                fold[k] = j++ % K;
        for (int g = 0; g < K; ++g) {
            std::vector<char> use(di.size());
            for (size_t k = 0; k < di.size(); ++k)
                use[k] = ok[k] && fold[k] != g;
            std::vector<double> qx, qy;
            if (!solve_order(order, use, qx, qy))
                return false;
            for (size_t k = 0; k < di.size(); ++k)
                if (ok[k] && fold[k] == g) {
                    double rx, ry;
                    resid(order, qx, qy, k, rx, ry);
                    sx += wt[k] * rx * rx;
                    sy += wt[k] * ry * ry;
                    sw += wt[k];
                    ++n;
                }
        }
        // weighted rms (equals plain rms without weights): typical error of a reference left out of the fit
        cvx = std::sqrt(sx / std::max(sw, 1e-300));
        cvy = std::sqrt(sy / std::max(sw, 1e-300));
        return n > 0;
    };
    // fit with iterative 3σ rejection (residuals normalized by each pair's σ when weighted)
    std::vector<double> rx, ry;
    auto set_coef = [&](int order, const std::vector<double> &qx, const std::vector<double> &qy) {
        int nt = plate_nterms(order);
        p.cx.assign(qx.begin(), qx.begin() + nt);
        p.cy.assign(qy.begin(), qy.begin() + nt);
        p.col_ra = use_col ? qx[nt] : 0;
        p.col_dec = use_col ? qy[nt] : 0;
    };
    auto fit_clip = [&](int order) {
        int nt = ncols(order);
        ok.assign(di.size(), 1);
        rx.assign(di.size(), 0);
        ry.assign(di.size(), 0);
        for (int it = 0; it < 20; ++it) {
            std::vector<double> qx, qy;
            if (!solve_order(order, ok, qx, qy))
                return false;
            set_coef(order, qx, qy);
            double sx = 0, sy = 0, sw = 0;
            int m = 0;
            for (size_t k = 0; k < di.size(); ++k) {
                resid(order, qx, qy, k, rx[k], ry[k]);
                if (ok[k]) {
                    sx += wt[k] * rx[k] * rx[k];
                    sy += wt[k] * ry[k] * ry[k];
                    sw += wt[k];
                    ++m;
                }
            }
            // weighted σ with degrees of freedom (m − nt): σ_k = σ·sqrt(w̄/w_k)
            double f_dof = double(m) / std::max(1, m - nt), s1 = std::sqrt(sx / sw * f_dof), s2 = std::sqrt(sy / sw * f_dof);
            double wbar = sw / std::max(1, m);
            bool changed = false;
            int kept = 0;
            for (size_t k = 0; k < di.size(); ++k) {
                double g = std::sqrt(wbar / wt[k]);
                char nk = std::fabs(rx[k]) <= 3 * s1 * g + 1e-9 && std::fabs(ry[k]) <= 3 * s2 * g + 1e-9;
                if (nk != ok[k])
                    changed = true;
                ok[k] = nk;
                kept += nk;
            }
            if (kept < nt + 2)
                return false;
            if (!changed)
                break;
        }
        return true;
    };
    for (int pass = 0; pass < 4; ++pass) {
        double r = rads[pass];
        // matching: for each catalog star, the nearest unsaturated detection (and vice versa, one-to-one)
        std::vector<double> pxi(f.det.size()), peta(f.det.size());
        for (size_t i = 0; i < f.det.size(); ++i) {
            auto t = p.terms(f.det[i].x, f.det[i].y);
            double a = 0, b = 0;
            for (size_t k = 0; k < t.size(); ++k) {
                a += t[k] * p.cx[k];
                b += t[k] * p.cy[k];
            }
            pxi[i] = a;
            peta[i] = b;
        }
        std::vector<int> best_cat(f.det.size(), -1);
        std::vector<double> best_d(f.det.size(), r * r);
        for (size_t j = 0; j < cat_std.size(); ++j) {
            int bi = -1;
            double bd = r * r;
            for (size_t i = 0; i < f.det.size(); ++i) {
                double dx = pxi[i] - cat_std[j].x, dy = peta[i] - cat_std[j].y;
                double d = dx * dx + dy * dy;
                if (d < bd) {
                    bd = d;
                    bi = int(i);
                }
            }
            if (bi >= 0 && bd < best_d[bi]) {
                best_d[bi] = bd;
                best_cat[bi] = int(j);
            }
        }
        di.clear();
        ci.clear();
        wt.clear();
        for (size_t i = 0; i < f.det.size(); ++i)
            if (best_cat[i] >= 0 && f.det[i].peak < c.saturation) {
                di.push_back(int(i));
                ci.push_back(best_cat[i]);
                // σ² = σ_centroid² + σ_cat² + floor² (″²); σ_centroid = FWHM / (2.355·SNR) (e.g. Neuschaefer & Windhorst
                // 1995, ApJS 96, 371; Mighell 2005, MNRAS 361, 861 for the Cramér–Rao bound)
                const Detection &d = f.det[i];
                double scen = d.fwhm * sc / (2.3548 * std::max(d.snr_tot, 1.0));
                double scat = has_err && std::isfinite(cat.err[best_cat[i]]) ? cat.err[best_cat[i]] : c.cat_err;
                wt.push_back(c.weighted ? 1.0 / (scen * scen + scat * scat + c.floor * c.floor) : 1.0);
            }
        ck.assign(di.size(), 0.0);
        use_col = false;
        if (c.color && cat.has_color() && pass >= 2) {
            double sc1 = 0, sc2 = 0;
            int n = 0;
            for (size_t k = 0; k < di.size(); ++k) {
                double cc = cat.color[ci[k]] - SOLAR_BP_RP;
                if (std::isfinite(cc) && std::fabs(cc) < 3) {
                    ck[k] = cc;
                    sc1 += cc;
                    sc2 += cc * cc;
                    ++n;
                }
            }
            use_col = n >= 10 && sc2 / n - (sc1 / n) * (sc1 / n) >= 0.04;
            p.ncol = use_col ? n : 0;
        }
        // polynomial order: only with enough references (≥ 3 per term); 0 = automatic via cross-validation
        int order = pass >= 2 ? (c.order > 0 ? c.order : 3) : 1;
        while (order > 1 && int(di.size()) < 3 * plate_nterms(order))
            --order;
        int nt = ncols(order);
        if (int(di.size()) < std::max(plate_nterms(1) + 2, 4))
            return false;
        if (pass >= 2 && c.order == 0) {
            // pick the order with the lowest out-of-sample rms (requires ≥ 3 % gain to increase the order)
            if (!fit_clip(1))
                return false;
            int best = 1;
            double bcv = INFINITY;
            for (int o = 1; o <= order; ++o) {
                double cx2, cy2;
                if (cross_val(o, cx2, cy2) && std::hypot(cx2, cy2) < bcv * (o == 1 ? 1.0 : 0.97)) {
                    bcv = std::hypot(cx2, cy2);
                    best = o;
                }
            }
            order = best;
            nt = ncols(order);
        }
        if (int(di.size()) < std::max(nt + 2, 4) || !fit_clip(order))
            return false;
        p.order = order;
        if (pass == 3) {
            // weighted rms with degrees of freedom: sqrt(Σw·r² / Σw · m/(m − nt))
            double sx = 0, sy = 0, sw = 0;
            int m = 0;
            for (size_t k = 0; k < di.size(); ++k)
                if (ok[k]) {
                    sx += wt[k] * rx[k] * rx[k];
                    sy += wt[k] * ry[k] * ry[k];
                    sw += wt[k];
                    ++m;
                }
            double fd = double(m) / std::max(1, m - nt);
            p.rms_ra = std::sqrt(sx / sw * fd);
            p.rms_dec = std::sqrt(sy / sw * fd);
            p.dra_mean = p.rms_ra;
            p.dde_mean = p.rms_dec;
            p.nref = m;
            if (!cross_val(order, p.rms_cv_ra, p.rms_cv_dec))
                p.rms_cv_ra = p.rms_cv_dec = 0;
            p.order = order;
            // error model: per-axis σ_sys such that Σ r²/(σ_c² + σ_cat² + σ_sys²) = m − nt (reduced χ² = 1), with σ_c
            // from the Cramér–Rao bound; covariance of the weighted-fit coefficients (weights w) via the sandwich form
            // (AᵀWA)⁻¹·AᵀW·V·WA·(AᵀWA)⁻¹ with V = σ_c² + σ_cat² + σ_sys² (Huber 1967; White 1980)
            std::vector<double> vk(di.size(), 0);
            for (size_t k = 0; k < di.size(); ++k) {
                double sce = centroid_sigma(f.det[di[k]], sc, c.gain);
                double scat = has_err && std::isfinite(cat.err[ci[k]]) ? cat.err[ci[k]] : c.cat_err;
                vk[k] = sce * sce + scat * scat;
            }
            auto sys_of = [&](const std::vector<double> &r) {
                auto chi = [&](double s2) {
                    double q = 0;
                    for (size_t k = 0; k < di.size(); ++k)
                        if (ok[k])
                            q += r[k] * r[k] / (vk[k] + s2);
                    return q - (m - nt);
                };
                if (chi(0) <= 0)
                    return 0.0;
                double lo = 0, hi = 1.0;
                while (chi(hi) > 0 && hi < 1e4)
                    hi *= 4;
                for (int it = 0; it < 60; ++it) {
                    double mid = 0.5 * (lo + hi);
                    (chi(mid) > 0 ? lo : hi) = mid;
                }
                return std::sqrt(0.5 * (lo + hi));
            };
            p.sys_ra = sys_of(rx);
            p.sys_dec = sys_of(ry);
            // the target is not a reference: the systematic error affecting it is the out-of-sample prediction error. With few
            // references or a high order the in-sample χ² underestimates σ_sys (the fit absorbs part of the error; degrees of
            // freedom m − nt only correct this on average). Alternative estimate from k-fold cross-validation:
            // σ_sys,cv² = rms_cv² − ⟨σ_c² + σ_cat²⟩_w (out-of-sample rms minus the references' own noise);
            // the larger of the two is adopted.
            if (p.rms_cv_ra > 0 && p.rms_cv_dec > 0) {
                double sv = 0, sw = 0;
                for (size_t k = 0; k < di.size(); ++k)
                    if (ok[k]) {
                        sv += wt[k] * vk[k];
                        sw += wt[k];
                    }
                double mv = sw > 0 ? sv / sw : 0;
                p.sys_ra = std::max(p.sys_ra, std::sqrt(std::max(0.0, p.rms_cv_ra * p.rms_cv_ra - mv)));
                p.sys_dec = std::max(p.sys_dec, std::sqrt(std::max(0.0, p.rms_cv_dec * p.rms_cv_dec - mv)));
            }
            p.nx0 = 0.5 * (f.im.w - 1);
            p.ny0 = 0.5 * (f.im.h - 1);
            p.ns = std::max(1.0, 0.5 * std::max(f.im.w, f.im.h));
            std::vector<double> Bn(size_t(nt) * nt, 0), Mx(size_t(nt) * nt, 0), My(size_t(nt) * nt, 0);
            Plate q;
            q.order = order;
            for (size_t k = 0; k < di.size(); ++k)
                if (ok[k]) {
                    auto t = q.terms((f.det[di[k]].x - p.nx0) / p.ns, (f.det[di[k]].y - p.ny0) / p.ns);
                    if (use_col)
                        t.push_back(ck[k]);
                    double w = wt[k];
                    for (int i = 0; i < nt; ++i)
                        for (int j = 0; j < nt; ++j) {
                            double tt = t[i] * t[j];
                            Bn[size_t(i) * nt + j] += w * tt;
                            Mx[size_t(i) * nt + j] += w * w * (vk[k] + p.sys_ra * p.sys_ra) * tt;
                            My[size_t(i) * nt + j] += w * w * (vk[k] + p.sys_dec * p.sys_dec) * tt;
                        }
                }
            p.cov_ra.clear();
            p.cov_dec.clear();
            if (invert(Bn, nt)) {
                auto sand = [&](const std::vector<double> &Mm) {
                    std::vector<double> T(size_t(nt) * nt, 0), O(size_t(nt) * nt, 0);
                    for (int i = 0; i < nt; ++i)
                        for (int j = 0; j < nt; ++j)
                            for (int l = 0; l < nt; ++l)
                                T[size_t(i) * nt + j] += Bn[size_t(i) * nt + l] * Mm[size_t(l) * nt + j];
                    for (int i = 0; i < nt; ++i)
                        for (int j = 0; j < nt; ++j)
                            for (int l = 0; l < nt; ++l)
                                O[size_t(i) * nt + j] += T[size_t(i) * nt + l] * Bn[size_t(l) * nt + j];
                    return O;
                };
                p.cov_ra = sand(Mx);
                p.cov_dec = sand(My);
                if (use_col) {
                    // the target has solar color (zero term): only the polynomial-term block enters the prediction
                    const int np = plate_nterms(order);
                    std::vector<double> cr(size_t(np) * np), cd(size_t(np) * np);
                    for (int i = 0; i < np; ++i)
                        for (int j = 0; j < np; ++j) {
                            cr[size_t(i) * np + j] = p.cov_ra[size_t(i) * nt + j];
                            cd[size_t(i) * np + j] = p.cov_dec[size_t(i) * nt + j];
                        }
                    p.cov_ra = cr;
                    p.cov_dec = cd;
                }
            }
        }
    }
    p.nmatch = int(di.size());
    p.ref_idx = di;
    p.cat_idx = ci;
    p.ref_ok = ok;
    return p.nref >= 4;
}

void fit_zero_point(Frame &f, const Catalog &cat, double saturation) {
    // m_cat + 2.5·log10(F) = zp + k·(c − c☉) per reference, with c = BP − RP (Gaia). The color term k converts the
    // catalog band (RP) to the frame filter (photometric transformation equation; Harris, FitzGerald & Reed 1981,
    // PASP 93, 507): without it, the median zero point of the references (field stars redder than the Sun)
    // makes a solar-color asteroid look ~0.5 mag brighter in g and ~0.3 mag fainter in z than in r, and the
    // magnitude series of a track across filters (SDSS r, i, z, g) varies more than the real brightness. Least-squares
    // fit with iterative 3σ rejection (MAD); only with ≥ 10 references with color and color spread
    // ≥ 0.2 mag (otherwise k = 0 and zp = median, as before). zp is referred to solar color (SOLAR_BP_RP).
    Plate &p = f.plate;
    std::vector<double> dm, cc;
    for (size_t k = 0; k < p.ref_idx.size(); ++k) {
        const Detection &d = f.det[p.ref_idx[k]];
        const int ci = p.cat_idx[k];
        if (!p.ref_ok[k] || d.flux <= 0 || d.peak >= saturation || !std::isfinite(cat.mag[ci]))
            continue;
        dm.push_back(cat.mag[ci] + 2.5 * std::log10(d.flux));
        cc.push_back(cat.has_color() ? cat.color[ci] - SOLAR_BP_RP : NAN);
    }
    if (dm.empty())
        return;
    double s = 0;
    int n = 0;
    p.zp = clipped_median(dm, 3.0, 10, &s, &n);
    p.zp_rms = s;
    p.nphot = n;
    p.dmag_mean = s;
    p.zp_color = 0;
    std::vector<char> use(dm.size(), 0);
    int nc = 0;
    for (size_t i = 0; i < dm.size(); ++i)
        if (std::isfinite(cc[i]) && std::fabs(cc[i]) < 3 && std::fabs(dm[i] - p.zp) < 3 * std::max(s, 0.02))
            use[i] = 1, ++nc;
    if (nc < 10)
        return;
    double a = p.zp, b = 0, rs = s;
    for (int it = 0; it < 5; ++it) {
        double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
        for (size_t i = 0; i < dm.size(); ++i)
            if (use[i]) {
                S += 1;
                Sx += cc[i];
                Sy += dm[i];
                Sxx += cc[i] * cc[i];
                Sxy += cc[i] * dm[i];
            }
        double var = Sxx / S - (Sx / S) * (Sx / S);
        if (S < 10 || var < 0.04)
            return;
        b = (S * Sxy - Sx * Sy) / (S * Sxx - Sx * Sx);
        a = (Sy - b * Sx) / S;
        std::vector<double> ad;
        for (size_t i = 0; i < dm.size(); ++i)
            if (std::isfinite(cc[i]))
                ad.push_back(std::fabs(dm[i] - a - b * cc[i]));
        rs = 1.4826 * median(ad);
        for (size_t i = 0; i < dm.size(); ++i)
            use[i] = std::isfinite(cc[i]) && std::fabs(cc[i]) < 3 && std::fabs(dm[i] - a - b * cc[i]) < 3 * std::max(rs, 0.02);
    }
    p.zp = a;
    p.zp_color = b;
    p.zp_rms = rs;
    p.dmag_mean = rs;
}

} // namespace internal

using namespace internal;

// Deferred trail search (after the plate): per-frame static-source mask — plate references, projected catalog
// stars and detections with a counterpart (within max(3″, seeing)) in another solved frame — and detect_trails
// outside it. Point fragments of each trail get typ = 3 (excluded from linking) and trails are appended to the end
// of Frame::det (reference indices preserved). Each frame's maps are released afterwards.
static void late_trails(Reduction &R, const Catalog &cat, const DetectParams &dp, std::vector<DetectScratch> &scratch, bool dbg) {
    const size_t NF = R.frames.size();
    // sky positions (standard coordinates around a common point, ″) of all detections in solved frames
    double ra0 = NAN, dec0 = NAN;
    for (auto &f : R.frames)
        if (f.solved && !std::isfinite(ra0)) {
            ra0 = f.plate.ra0;
            dec0 = f.plate.dec0;
        }
    if (!std::isfinite(ra0))
        return;
    std::vector<std::vector<Pt>> X(NF);
    double cell = 3;
    for (size_t k = 0; k < NF; ++k) {
        const Frame &f = R.frames[k];
        if (!f.solved)
            continue;
        cell = std::max(cell, f.seeing);
        for (auto &d : f.det) {
            double ra, dec;
            Pt p;
            f.plate.pix_to_radec(d.x, d.y, ra, dec);
            radec_to_std(ra, dec, ra0, dec0, p.x, p.y);
            X[k].push_back(p);
        }
    }
    auto key = [&](double x, double y) { return std::make_pair(long(std::floor(x / cell)), long(std::floor(y / cell))); };
    std::vector<std::map<std::pair<long, long>, std::vector<int>>> G(NF);
    for (size_t k = 0; k < NF; ++k)
        for (size_t i = 0; i < X[k].size(); ++i)
            G[k][key(X[k][i].x, X[k][i].y)].push_back(int(i));
    std::vector<Pt> C(cat.size());
    std::map<std::pair<long, long>, std::vector<int>> GC;
    for (size_t j = 0; j < cat.size(); ++j) {
        radec_to_std(cat.ra[j], cat.dec[j], ra0, dec0, C[j].x, C[j].y);
        GC[key(C[j].x, C[j].y)].push_back(int(j));
    }
    char b[256];
    for (size_t k = 0; k < NF; ++k) {
        Frame &f = R.frames[k];
        if (!f.solved || !scratch[k].ok)
            continue;
        const double sc = 0.5 * (f.plate.scale_x() + f.plate.scale_y()), r = std::max(3.0, f.seeing);
        std::vector<double> mask;
        for (size_t i = 0; i < f.det.size(); ++i) {
            const Pt &p = X[k][i];
            bool st = f.typ[i] == 2;
            auto kk = key(p.x, p.y);
            for (long a = kk.first - 1; a <= kk.first + 1 && !st; ++a)
                for (long c = kk.second - 1; c <= kk.second + 1 && !st; ++c) {
                    auto it = GC.find({a, c});
                    if (it != GC.end())
                        for (int j : it->second)
                            st = st || std::hypot(C[j].x - p.x, C[j].y - p.y) <= r;
                    for (size_t m = 0; m < NF && !st; ++m) {
                        if (m == k || !R.frames[m].solved)
                            continue;
                        auto jt = G[m].find({a, c});
                        if (jt != G[m].end())
                            for (int q : jt->second)
                                st = st || std::hypot(X[m][q].x - p.x, X[m][q].y - p.y) <= r;
                    }
                }
            if (st) {
                mask.push_back(f.det[i].x);
                mask.push_back(f.det[i].y);
                mask.push_back(r / sc);
            }
        }
        DetectStats ds;
        std::vector<char> frag;
        std::vector<Detection> T = detect_trails(f.im, dp, scratch[k], f.det, mask, frag, &ds);
        for (size_t i = 0; i < frag.size() && i < f.typ.size(); ++i)
            if (frag[i] && f.typ[i] != 2)
                f.typ[i] = 3;
        for (auto &t : T) {
            f.det.push_back(t);
            f.typ.push_back(0);
        }
        if (dbg) {
            std::snprintf(b, sizeof b,
                          "trails image %zu: deferred search %.3f s (%zu masked); kernels %.3f s (%d), mask/components %.3f s "
                          "(%d comp., %d cand.), fits %.3f s (%d, %d in prescreen), replacement %.3f s; %d trail(s)",
                          k + 1, ds.t_total, mask.size() / 3, ds.t_kernels, ds.nkern, ds.t_mask, ds.ncomp, ds.ncand, ds.t_fits, ds.nfit,
                          ds.nskip, ds.t_replace, ds.ntrail);
            R.log.push_back(b);
        }
        scratch[k] = DetectScratch{};
    }
}

Reduction reduce(std::vector<Image> images, const Catalog &cat_in, double ra0, double dec0, const Settings &s_in, bool moving) {
    Reduction R;
    // gain: Settings Gain when set; otherwise GAIN/EGAIN from the FITS header (median over images); otherwise 1 e⁻/ADU
    Settings s = s_in;
    if (s.str("gain").empty() || !(s.num("gain") > 0)) {
        std::vector<double> g;
        for (auto &im : images)
            if (im.gain > 0)
                g.push_back(im.gain);
        if (!g.empty()) {
            char gb[32];
            std::snprintf(gb, sizeof gb, "%.6g", median(g));
            s.set("gain", gb);
        }
    }
    auto get = [&](const char *k, double def) { return s.str(k).empty() ? def : s.num(k); };
    // reference catalog: excluding stars brighter than catalog_bright_mag (saturated); indices in Plate::cat_idx
    Catalog cat;
    for (size_t i = 0; i < cat_in.size(); ++i)
        if (cat_in.mag[i] >= get("catalog_bright_mag", -99))
            cat.push_from(cat_in, i);
    DetectParams dp = detect_params(s);
    R.gain = dp.gain;
    double focal = get("focal_length_mm", 1000);
    SolveCtx c;
    c.sx = get("pixel_width_mm", 0.009) / focal * 206264.806;
    c.sy = get("pixel_height_mm", 0.009) / focal * 206264.806;
    c.pa = get("position_angle", 0);
    c.pointing = get("pointing_tolerance", 5);
    c.match_rad = std::max(get("match_radius", 2.0), 1.0 * std::max(c.sx, c.sy));
    c.scale_tol = std::clamp(get("focal_length_tolerance", 1.0) / 100 * 5, 0.02, 0.2);
    // plate_order 0 = automatic (cross-validation); plate_weights 1 = 1/σ² weights; σ_cat (no catalog error) and floor in ″
    c.order = std::clamp(int(get("plate_order", 1)), 0, 3);
    c.weighted = get("plate_weights", 1) != 0;
    c.cat_err = get("plate_cat_err", 0.01);
    c.floor = get("plate_floor", 0.10);
    c.color = get("plate_color", 1) != 0;
    c.saturation = dp.saturation;
    c.gain = dp.gain;
    auto project = [&](double a0, double d0) {
        std::vector<Pt> q(cat.size());
        for (size_t i = 0; i < cat.size(); ++i)
            radec_to_std(cat.ra[i], cat.dec[i], a0, d0, q[i].x, q[i].y);
        return q;
    };
    std::vector<Pt> cat_nom = project(ra0, dec0);
    char b[256];
    const bool tdbg = get("trail_debug", 0) != 0;
    // trail_late 1 (default): trail search runs after the plate, only outside static sources (detect_trails)
    const bool late = dp.trails && get("trail_late", 1) != 0;
    std::vector<DetectScratch> scratch(images.size());
    for (auto &im : images) {
        Frame f;
        f.im = std::move(im);
        DetectStats ds;
        f.det = detect_sources(f.im, dp, tdbg ? &ds : nullptr, late ? &scratch[R.frames.size()] : nullptr);
        if (tdbg && !late) {
            std::snprintf(b, sizeof b,
                          "trails image %zu: detection %.3f s; kernels %.3f s (%d), mask/components %.3f s (%d comp., %d cand.), "
                          "fits %.3f s (%d), replacement %.3f s; %d trail(s)",
                          R.frames.size() + 1, ds.t_total, ds.t_kernels, ds.nkern, ds.t_mask, ds.ncomp, ds.ncand, ds.t_fits, ds.nfit,
                          ds.t_replace, ds.ntrail);
            R.log.push_back(b);
        }
        R.frames.push_back(std::move(f));
    }
    int first_solved = -1;
    // two passes: a frame whose field cannot be identified against the catalog (e.g. bright stars saturated
    // in one filter but not another) is aligned to an already solved frame; the 2nd pass makes this order-independent
    // (previously, a failure on the 1st frame had nothing to align to)
    for (int pass = 0; pass < 2; ++pass)
        for (size_t k = 0; k < R.frames.size(); ++k) {
            Frame &f = R.frames[k];
            if (f.solved || (pass == 1 && first_solved < 0))
                continue;
            const bool last = pass == 1 || first_solved >= 0;
            Affine T;
            c.ra0 = ra0;
            c.dec0 = dec0;
            bool ok = initial_solution(f, cat_nom, cat, c, T);
            if (!ok && first_solved >= 0) {
                // alignment with an already solved frame (triangles between detection lists, in pixels)
                const Frame &g = R.frames[first_solved];
                std::vector<Pt> P, Q;
                for (auto &d : f.det)
                    if (P.size() < 60)
                        P.push_back({d.x, d.y});
                for (auto &d : g.det)
                    if (Q.size() < 60)
                        Q.push_back({d.x, d.y});
                Affine A;
                if (match_triangles(P, Q, 3.0, 0.05, A) >= 6) {
                    // pixel(f) → pixel(g) → standard coordinates around (ra0, dec0) via g's linear plate
                    double xi0, eta0;
                    radec_to_std(g.plate.ra0, g.plate.dec0, ra0, dec0, xi0, eta0);
                    const Plate &q = g.plate;
                    T.a[0] = q.cx[0] + q.cx[1] * A.a[0] + q.cx[2] * A.b[0] + xi0;
                    T.a[1] = q.cx[1] * A.a[1] + q.cx[2] * A.b[1];
                    T.a[2] = q.cx[1] * A.a[2] + q.cx[2] * A.b[2];
                    T.b[0] = q.cy[0] + q.cy[1] * A.a[0] + q.cy[2] * A.b[0] + eta0;
                    T.b[1] = q.cy[1] * A.a[1] + q.cy[2] * A.b[1];
                    T.b[2] = q.cy[1] * A.a[2] + q.cy[2] * A.b[2];
                    ok = true;
                }
            }
            if (!ok) {
                if (last) {
                    std::snprintf(b, sizeof b, "image %zu: field not identified (%zu detections)", k + 1, f.det.size());
                    R.log.push_back(b);
                }
                continue;
            }
            if (!fit_plate(f, cat_nom, cat, c, T)) {
                std::snprintf(b, sizeof b, "image %zu: plate fit failed", k + 1);
                R.log.push_back(b);
                continue;
            }
            // reproject the catalog with the tangent point at the image center and redo the fit
            double rc, dc;
            f.plate.pix_to_radec((f.im.w - 1) / 2.0, (f.im.h - 1) / 2.0, rc, dc);
            std::vector<Pt> cat_c = project(rc, dc);
            Affine T2;
            {
                // linear part of the current plate expressed around the new tangent point
                std::vector<Pt> P, Q;
                for (size_t i = 0; i < f.plate.ref_idx.size(); ++i)
                    if (f.plate.ref_ok[i]) {
                        P.push_back({f.det[f.plate.ref_idx[i]].x, f.det[f.plate.ref_idx[i]].y});
                        Q.push_back(cat_c[f.plate.cat_idx[i]]);
                    }
                fit_affine(P, Q, T2);
            }
            c.ra0 = rc;
            c.dec0 = dc;
            if (!fit_plate(f, cat_c, cat, c, T2)) {
                std::snprintf(b, sizeof b, "image %zu: plate fit failed", k + 1);
                R.log.push_back(b);
                continue;
            }
            fit_zero_point(f, cat, dp.saturation);
            if (get("plate_debug", 0) != 0) {
                // diagnostics: residual of each reference (plate − catalog, ″) with position, magnitude, color, SNR and FWHM,
                // to study biases by magnitude, color, chip position and pixel phase ("ref|..." lines in the log)
                const double scd = 0.5 * (f.plate.scale_x() + f.plate.scale_y());
                for (size_t i = 0; i < f.plate.ref_idx.size(); ++i) {
                    const Detection &d = f.det[f.plate.ref_idx[i]];
                    const int ci = f.plate.cat_idx[i];
                    auto t = f.plate.terms(d.x, d.y);
                    double xi = 0, eta = 0;
                    for (size_t q = 0; q < t.size(); ++q) {
                        xi += t[q] * f.plate.cx[q];
                        eta += t[q] * f.plate.cy[q];
                    }
                    const double cc = cat.has_color() && std::isfinite(cat.color[ci]) && std::fabs(cat.color[ci] - SOLAR_BP_RP) < 3
                                          ? cat.color[ci] - SOLAR_BP_RP
                                          : 0;
                    xi += f.plate.col_ra * cc;
                    eta += f.plate.col_dec * cc;
                    std::snprintf(b, sizeof b, "ref|%zu|%.3f|%.3f|%.3f|%.3f|%.1f|%.2f|%+.4f|%+.4f|%d|%.4f|%.4f", k + 1, d.x, d.y,
                                  cat.mag[ci], cat.has_color() ? cat.color[ci] : NAN, d.snr_tot, d.fwhm, xi - cat_c[ci].x,
                                  eta - cat_c[ci].y, int(f.plate.ref_ok[i]), centroid_sigma(d, scd, dp.gain),
                                  cat.err.size() == cat.size() ? cat.err[ci] : c.cat_err);
                    R.log.push_back(b);
                }
            }
            f.solved = true;
            if (first_solved < 0) {
                first_solved = int(k);
                R.npattern = f.plate.nmatch;
            }
            std::vector<double> fw;
            for (size_t i = 0; i < f.plate.ref_idx.size(); ++i)
                if (f.plate.ref_ok[i])
                    fw.push_back(f.det[f.plate.ref_idx[i]].fwhm);
            f.seeing = fw.empty() ? 0 : median(fw) * 0.5 * (f.plate.scale_x() + f.plate.scale_y());
            std::snprintf(b, sizeof b, "image %zu: %d/%d references, rms %.2f\"/%.2f\", ZP %.3f ± %.3f (%d), color %+.4f/%+.4f \"/mag (%d)",
                          k + 1, f.plate.nref, f.plate.nmatch, f.plate.rms_ra, f.plate.rms_dec, f.plate.zp, f.plate.zp_rms, f.plate.nphot,
                          f.plate.col_ra, f.plate.col_dec, f.plate.ncol);
            R.log.push_back(b);
        }
    // offset of each frame relative to the first solved one (pixels)
    if (first_solved >= 0) {
        const Frame &g = R.frames[first_solved];
        double rc, dc;
        g.plate.pix_to_radec((g.im.w - 1) / 2.0, (g.im.h - 1) / 2.0, rc, dc);
        for (auto &f : R.frames)
            if (f.solved) {
                const Plate &p = f.plate;
                double xi, eta;
                radec_to_std(rc, dc, p.ra0, p.dec0, xi, eta);
                double a = p.cx[1], bb = p.cx[2], cc = p.cy[1], d = p.cy[2], det = a * d - bb * cc;
                double px = (d * (xi - p.cx[0]) - bb * (eta - p.cy[0])) / det, py = (-cc * (xi - p.cx[0]) + a * (eta - p.cy[0])) / det;
                f.offx = px - (f.im.w - 1) / 2.0;
                f.offy = py - (f.im.h - 1) / 2.0;
            }
    }
    for (auto &f : R.frames) {
        f.typ.assign(f.det.size(), 0);
        for (size_t i = 0; i < f.plate.ref_idx.size(); ++i)
            if (f.plate.ref_ok[i])
                f.typ[f.plate.ref_idx[i]] = 2;
    }
    if (late)
        late_trails(R, cat, dp, scratch, tdbg);
    if (moving) {
        R.movers = link_movers(R, cat_in, s);
        for (auto &T : stack_search(R, cat_in, s))
            R.movers.push_back(std::move(T));
        refine_movers(R, cat, s);
        drop_incoherent_movers(R);
    }
    if (get("plate_debug", 0) != 0)
        // diagnostics: all detections in solved frames with plate positions ("det|..." lines)
        for (size_t k = 0; k < R.frames.size(); ++k) {
            const Frame &f = R.frames[k];
            if (!f.solved)
                continue;
            for (size_t i = 0; i < f.det.size(); ++i) {
                const Detection &d = f.det[i];
                double ra, dec;
                f.plate.pix_to_radec(d.x, d.y, ra, dec);
                std::snprintf(b, sizeof b, "det|%zu|%.3f|%.3f|%.6f|%+.6f|%.1f|%.2f|%.2f|%d", k + 1, d.x, d.y, ra, dec, d.snr_tot, d.fwhm,
                              f.plate.zp - 2.5 * std::log10(std::max(d.flux, 1e-9)), i < f.typ.size() ? f.typ[i] : -1);
                R.log.push_back(b);
            }
        }
    return R;
}

} // namespace astro
