// Gnomonic projection, plate constants and field identification (detection ↔ catalogue matching).
//
// References (independent implementation):
//   - Green (1985), "Spherical Astronomy", ch. 13 — standard coordinates (gnomonic projection) and plate constants.
//   - Calabretta & Greisen (2002), A&A 395, 1077 — TAN projection.
//   - Groth (1986), AJ 91, 1244 — point-list matching by triangles.
//   - Valdes, Campusano, Velasquez & Stetson (1995), PASP 107, 1119 — FOCAS: space of invariant triangle ratios
//     (b/a, c/a) and voting on vertex correspondences.
//   - Fischler & Bolles (1981), Comm. ACM 24, 381 — RANSAC (hypothesis consensus).
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <tuple>

namespace astro {

static constexpr double D2R = PI / 180, AS2R = PI / 180 / 3600;

void radec_to_std(double ra, double dec, double ra0, double dec0, double &xi, double &eta) {
    double a = ra * D2R, d = dec * D2R, a0 = ra0 * D2R, d0 = dec0 * D2R;
    double cosc = std::sin(d0) * std::sin(d) + std::cos(d0) * std::cos(d) * std::cos(a - a0);
    xi = std::cos(d) * std::sin(a - a0) / cosc / AS2R;
    eta = (std::cos(d0) * std::sin(d) - std::sin(d0) * std::cos(d) * std::cos(a - a0)) / cosc / AS2R;
}

void std_to_radec(double xi, double eta, double ra0, double dec0, double &ra, double &dec) {
    double x = xi * AS2R, y = eta * AS2R, a0 = ra0 * D2R, d0 = dec0 * D2R;
    double den = std::cos(d0) - y * std::sin(d0);
    double a = a0 + std::atan2(x, den);
    dec = std::atan2(std::sin(d0) + y * std::cos(d0), std::hypot(x, den)) / D2R;
    ra = std::fmod(a / D2R + 360.0, 360.0);
}

// Polynomial terms: 1, x, y | x², xy, y² | x³, x²y, xy², y³ (order of the cx/cy coefficients).
std::vector<double> Plate::terms(double x, double y) const {
    std::vector<double> t{1, x, y};
    if (order >= 2) {
        t.push_back(x * x);
        t.push_back(x * y);
        t.push_back(y * y);
    }
    if (order >= 3) {
        t.push_back(x * x * x);
        t.push_back(x * x * y);
        t.push_back(x * y * y);
        t.push_back(y * y * y);
    }
    return t;
}

void Plate::pred_sigma(double x, double y, double &sra, double &sdec) const {
    sra = sdec = 0;
    Plate q;
    q.order = order;
    auto t = q.terms((x - nx0) / ns, (y - ny0) / ns);
    const size_t n = t.size();
    if (cov_ra.size() != n * n || cov_dec.size() != n * n)
        return;
    double a = 0, b = 0;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            a += t[i] * cov_ra[i * n + j] * t[j];
            b += t[i] * cov_dec[i * n + j] * t[j];
        }
    sra = std::sqrt(std::max(a, 0.0));
    sdec = std::sqrt(std::max(b, 0.0));
}

void Plate::pix_to_radec(double x, double y, double &ra, double &dec) const {
    auto t = terms(x, y);
    double xi = 0, eta = 0;
    for (size_t i = 0; i < t.size() && i < cx.size() && i < cy.size(); ++i) {
        xi += t[i] * cx[i];
        eta += t[i] * cy[i];
    }
    std_to_radec(xi, eta, ra0, dec0, ra, dec);
}

void plate_jacobian(const Plate &p, double x, double y, double J[4], double ra_t, double dec_t) {
    double r0, d0, r1, d1, r2, d2, u0, v0, u1, v1, u2, v2;
    p.pix_to_radec(x, y, r0, d0);
    p.pix_to_radec(x + 1, y, r1, d1);
    p.pix_to_radec(x, y + 1, r2, d2);
    if (!std::isfinite(ra_t) || !std::isfinite(dec_t)) {
        ra_t = r0;
        dec_t = d0;
    }
    radec_to_std(r0, d0, ra_t, dec_t, u0, v0);
    radec_to_std(r1, d1, ra_t, dec_t, u1, v1);
    radec_to_std(r2, d2, ra_t, dec_t, u2, v2);
    J[0] = u1 - u0;
    J[1] = u2 - u0;
    J[2] = v1 - v0;
    J[3] = v2 - v0;
}

double Plate::scale_x() const {
    return cx.size() >= 3 ? std::hypot(cx[1], cy[1]) : 0;
}
double Plate::scale_y() const {
    return cy.size() >= 3 ? std::hypot(cx[2], cy[2]) : 0;
}

namespace internal {

int plate_nterms(int order) {
    return order >= 3 ? 10 : order == 2 ? 6 : 3;
}

bool fit_affine(const std::vector<Pt> &P, const std::vector<Pt> &Q, Affine &T) {
    int m = int(std::min(P.size(), Q.size()));
    if (m < 3)
        return false;
    std::vector<double> A(size_t(m) * 3), bu(m), bv(m), cu, cv;
    for (int i = 0; i < m; ++i) {
        A[size_t(i) * 3] = 1;
        A[size_t(i) * 3 + 1] = P[i].x;
        A[size_t(i) * 3 + 2] = P[i].y;
        bu[i] = Q[i].x;
        bv[i] = Q[i].y;
    }
    if (!lstsq(A, m, 3, bu, cu) || !lstsq(A, m, 3, bv, cv))
        return false;
    for (int k = 0; k < 3; ++k) {
        T.a[k] = cu[k];
        T.b[k] = cv[k];
    }
    return true;
}

namespace {

struct Tri {
    int v[3];      // vertices: A (longest ∩ middle side), B (longest ∩ shortest), C (middle ∩ shortest)
    double r1, r2; // b/a, c/a
    double a;      // lado maior
    int sgn;       // orientation (sign of the cross product)
};

static std::vector<Tri> triangles(const std::vector<Pt> &P, int n, double min_side) {
    std::vector<Tri> out;
    n = std::min<int>(n, int(P.size()));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j + 1; k < n; ++k) {
                int id[3] = {i, j, k};
                // side s[t] is opposite vertex id[t]
                double s[3] = {std::hypot(P[j].x - P[k].x, P[j].y - P[k].y), std::hypot(P[i].x - P[k].x, P[i].y - P[k].y),
                               std::hypot(P[i].x - P[j].x, P[i].y - P[j].y)};
                int o[3] = {0, 1, 2};
                std::sort(o, o + 3, [&](int a, int b) { return s[a] > s[b]; });
                double a = s[o[0]], b = s[o[1]], c = s[o[2]];
                if (c < min_side || b / a > 0.97 || c / a < 0.12)
                    continue;
                Tri t;
                // the vertex opposite the shortest side is shared by the longest and middle ones, etc.
                t.v[0] = id[o[2]];
                t.v[1] = id[o[1]];
                t.v[2] = id[o[0]];
                t.r1 = b / a;
                t.r2 = c / a;
                t.a = a;
                const Pt &A = P[t.v[0]], &B = P[t.v[1]], &C = P[t.v[2]];
                double cr = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x);
                t.sgn = cr > 0 ? 1 : -1;
                out.push_back(t);
            }
    return out;
}

// Similarity P→Q from two pairs; mirrored if mirror. false if degenerate.
static bool sim2(const Pt &p1, const Pt &p2, const Pt &q1, const Pt &q2, bool mirror, Affine &T, double &scale) {
    using C = std::complex<double>;
    C z1(p1.x, mirror ? -p1.y : p1.y), z2(p2.x, mirror ? -p2.y : p2.y), w1(q1.x, q1.y), w2(q2.x, q2.y);
    if (std::abs(z2 - z1) < 1e-9)
        return false;
    C m = (w2 - w1) / (z2 - z1), t = w1 - m * z1;
    scale = std::abs(m);
    double mr = m.real(), mi = m.imag(), sg = mirror ? -1 : 1;
    T.a[0] = t.real();
    T.a[1] = mr;
    T.a[2] = -mi * sg;
    T.b[0] = t.imag();
    T.b[1] = mi;
    T.b[2] = mr * sg;
    return true;
}

// Agreeing pairs: for each transformed P, the nearest Q within tol (one Q per P, the nearest).
static int consensus(const std::vector<Pt> &P, const std::vector<Pt> &Q, const Affine &T, double tol,
                     std::vector<std::pair<int, int>> *pairs) {
    std::vector<int> best(Q.size(), -1);
    std::vector<double> bd(Q.size(), tol * tol);
    for (size_t i = 0; i < P.size(); ++i) {
        double u, v;
        T.apply(P[i].x, P[i].y, u, v);
        int bj = -1;
        double b = tol * tol;
        for (size_t j = 0; j < Q.size(); ++j) {
            double d = (Q[j].x - u) * (Q[j].x - u) + (Q[j].y - v) * (Q[j].y - v);
            if (d < b) {
                b = d;
                bj = int(j);
            }
        }
        if (bj >= 0 && b < bd[bj]) {
            bd[bj] = b;
            best[bj] = int(i);
        }
    }
    int n = 0;
    if (pairs)
        pairs->clear();
    for (size_t j = 0; j < Q.size(); ++j)
        if (best[j] >= 0) {
            ++n;
            if (pairs)
                pairs->push_back({best[j], int(j)});
        }
    return n;
}

} // namespace

int match_triangles(const std::vector<Pt> &P, const std::vector<Pt> &Q, double tol, double scale_tol, Affine &T) {
    const int nP = std::min<int>(int(P.size()), 25), nQ = std::min<int>(int(Q.size()), 90);
    if (nP < 4 || nQ < 4)
        return 0;
    auto tp = triangles(P, nP, 3 * tol), tq = triangles(Q, nQ, 3 * tol);
    // index of tq by cell of the (r1, r2) space
    const double cell = 0.01, eps = 0.006;
    std::map<std::pair<int, int>, std::vector<int>> grid;
    for (size_t i = 0; i < tq.size(); ++i)
        grid[{int(tq[i].r1 / cell), int(tq[i].r2 / cell)}].push_back(int(i));
    // votes per parity (0: same orientation, 1: mirrored)
    std::vector<int> votes[2] = {std::vector<int>(size_t(nP) * nQ, 0), std::vector<int>(size_t(nP) * nQ, 0)};
    for (auto &a : tp) {
        int c1 = int(a.r1 / cell), c2 = int(a.r2 / cell);
        for (int d1 = -1; d1 <= 1; ++d1)
            for (int d2 = -1; d2 <= 1; ++d2) {
                auto it = grid.find({c1 + d1, c2 + d2});
                if (it == grid.end())
                    continue;
                for (int k : it->second) {
                    const Tri &b = tq[k];
                    if (std::fabs(a.r1 - b.r1) > eps || std::fabs(a.r2 - b.r2) > eps)
                        continue;
                    double sc = b.a / a.a;
                    if (std::fabs(sc - 1) > scale_tol)
                        continue;
                    int par = a.sgn == b.sgn ? 0 : 1;
                    for (int v = 0; v < 3; ++v)
                        ++votes[par][size_t(a.v[v]) * nQ + b.v[v]];
                }
            }
    }
    int best_n = 0;
    Affine bestT;
    for (int par = 0; par < 2; ++par) {
        // most-voted correspondences (best Q for each P)
        std::vector<std::tuple<int, int, int>> cand; // votos, i, j
        for (int i = 0; i < nP; ++i) {
            int bj = -1, bv = 0;
            for (int j = 0; j < nQ; ++j)
                if (votes[par][size_t(i) * nQ + j] > bv) {
                    bv = votes[par][size_t(i) * nQ + j];
                    bj = j;
                }
            if (bv >= 2)
                cand.push_back({bv, i, bj});
        }
        std::sort(cand.begin(), cand.end(), [](auto &x, auto &y) { return std::get<0>(x) > std::get<0>(y); });
        if (cand.size() > 20)
            cand.resize(20);
        // two-pair hypotheses → consensus over all points
        for (size_t u = 0; u < cand.size(); ++u)
            for (size_t w = u + 1; w < cand.size(); ++w) {
                int i1 = std::get<1>(cand[u]), j1 = std::get<2>(cand[u]), i2 = std::get<1>(cand[w]), j2 = std::get<2>(cand[w]);
                if (j1 == j2)
                    continue;
                Affine H;
                double sc;
                if (!sim2(P[i1], P[i2], Q[j1], Q[j2], par == 1, H, sc) || std::fabs(sc - 1) > scale_tol)
                    continue;
                int n = consensus(P, Q, H, tol, nullptr);
                if (n > best_n) {
                    best_n = n;
                    bestT = H;
                }
            }
    }
    if (best_n < 4)
        return 0;
    // refine with a least-squares affine fit over the consensus
    std::vector<std::pair<int, int>> pr;
    consensus(P, Q, bestT, tol, &pr);
    std::vector<Pt> a, b;
    for (auto &q : pr) {
        a.push_back(P[q.first]);
        b.push_back(Q[q.second]);
    }
    Affine R = bestT;
    if (a.size() >= 6 && fit_affine(a, b, R) && consensus(P, Q, R, tol, nullptr) >= best_n)
        bestT = R;
    T = bestT;
    return consensus(P, Q, T, tol, nullptr);
}

bool initial_solution(const Frame &f, const std::vector<Pt> &cat_std, const Catalog &cat, const SolveCtx &c, Affine &T) {
    (void)cat;
    const Image &im = f.im;
    double cxp = (im.w - 1) / 2.0, cyp = (im.h - 1) / 2.0;
    // detections in ″ (nominal scale), centred; already sorted by flux
    std::vector<Pt> P;
    for (auto &d : f.det) {
        if (d.peak >= c.saturation)
            continue;
        P.push_back({(d.x - cxp) * c.sx, (d.y - cyp) * c.sy});
        if (P.size() >= 60)
            break;
    }
    // catalogue stars in a box covering the field plus the pointing uncertainty (cat_std already sorted by mag)
    double hw = 0.5 * std::hypot(im.w * c.sx, im.h * c.sy) + c.pointing * 60;
    double fov = im.w * c.sx * im.h * c.sy, box = 4 * hw * hw;
    std::vector<Pt> Q;
    for (auto &q : cat_std)
        if (std::fabs(q.x) < hw && std::fabs(q.y) < hw)
            Q.push_back(q);
    if (P.size() < 4 || Q.size() < 4)
        return false;
    // expected fraction of catalogue stars in the field → how many to use for triangles
    size_t nq = size_t(std::ceil(25 * box / fov * 0.8));
    std::vector<Pt> Qb(Q.begin(), Q.begin() + std::min(Q.size(), std::max<size_t>(nq, 25)));
    double tol = std::max(3 * c.match_rad, 3 * std::max(c.sx, c.sy));
    // matching: first with the brightest; consensus check with more points
    Affine S;
    int n = match_triangles(P, Qb, tol, c.scale_tol, S);
    if (n < 6)
        return false;
    // compose: pixel → centred ″ → standard
    T.a[0] = S.a[0] - S.a[1] * cxp * c.sx - S.a[2] * cyp * c.sy;
    T.a[1] = S.a[1] * c.sx;
    T.a[2] = S.a[2] * c.sy;
    T.b[0] = S.b[0] - S.b[1] * cxp * c.sx - S.b[2] * cyp * c.sy;
    T.b[1] = S.b[1] * c.sx;
    T.b[2] = S.b[2] * c.sy;
    return true;
}

} // namespace internal
} // namespace astro
