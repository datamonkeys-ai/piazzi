// Moving-object search: separation of fixed stars and linking of detections by uniform rectilinear motion.
//
// References (independent implementation):
//   - Kubica et al. (2007), Icarus 189, 151 — linking detections across nights/frames by linear-motion hypotheses
//     and verification in the intermediate frames ("tracklets").
//   - Denneau et al. (2013), PASP 125, 357 — Pan-STARRS MOPS: intra-night tracklets from detection pairs with
//     bounded velocity and confirmation by linear prediction; removal of stationary sources (static sky) before
//     linking.
//   - Linear least-squares fit of ξ(t), η(t) (constant position and velocity).
//   - Forced aperture photometry (Stetson 1987, PASP 99, 191; CCD equation, Howell 2006 §4.4) at the position of each
//     detection in the frames where the object has already moved away: a fixed source is still there, a mover is not.
//   - Position angle of motion measured from north through east.
//   - Trails: length and direction of each elongated detection consistent with v·t_exp and the track direction
//     (Vereš et al. 2012, PASP 124, 1197; Dawson et al. 2016, AJ 152, 112).
#include "piazzi/reduce.hpp"
#include "reduce_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>

namespace astro {
namespace internal {

namespace {

// Simple spatial index (grid of cells of side `cell`) for neighbor search in standard coordinates (″).
struct Grid {
    double cell = 1;
    std::unordered_map<long long, std::vector<int>> m;
    // shift as unsigned: (i << 32) with negative i is UB in C++17
    static long long key(long long i, long long j) {
        return (long long)((static_cast<unsigned long long>(i) << 32) ^ (static_cast<unsigned long long>(j) & 0xffffffffULL));
    }
    void build(const std::vector<Pt> &P, double c) {
        cell = c;
        m.clear();
        for (size_t k = 0; k < P.size(); ++k)
            m[key(std::llround(std::floor(P[k].x / cell)), std::llround(std::floor(P[k].y / cell)))].push_back(int(k));
    }
    // indices within r of (x, y) (r ≤ cell)
    template <class Fn> void near(double x, double y, Fn fn) const {
        long long i0 = std::llround(std::floor(x / cell)), j0 = std::llround(std::floor(y / cell));
        for (long long i = i0 - 1; i <= i0 + 1; ++i)
            for (long long j = j0 - 1; j <= j0 + 1; ++j) {
                auto it = m.find(key(i, j));
                if (it != m.end())
                    for (int k : it->second)
                        fn(k);
            }
    }
};

// (α, δ) → pixel by numerical (Newton) inversion of the plate; (x, y) holds the initial guess.
bool radec_to_pix(const Plate &p, double ra, double dec, double &x, double &y) {
    for (int it = 0; it < 8; ++it) {
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
            return std::isfinite(x) && std::isfinite(y);
    }
    return std::isfinite(x) && std::isfinite(y);
}

// Flux fraction of a circular Gaussian (σ in px) centered `d` px from the center of an aperture of radius `rap`,
// relative to the same aperture centered on the source (numerical integration on a 0.1·rap grid).
double offset_aperture_frac(double d, double rap, double sg) {
    auto frac = [&](double off) {
        double sum = 0;
        const int n = 10;
        for (int j = -n; j <= n; ++j)
            for (int i = -n; i <= n; ++i) {
                double u = i * rap / n, v = j * rap / n;
                if (u * u + v * v <= rap * rap)
                    sum += std::exp(-((u - off) * (u - off) + v * v) / (2 * sg * sg));
            }
        return sum;
    };
    double f0 = frac(0);
    return f0 > 0 ? frac(d) / f0 : 0;
}

} // namespace

bool flux_coherent(const std::vector<double> &mag, const std::vector<double> &smag, bool color_terms) {
    // A Solar System object has, in magnitudes referred to the solar color (per-frame color term), the same flux
    // in all filters up to its light curve and intrinsic color: σ_int = 0.2 mag (rotational variation in ≲ 10 min and
    // spread of asteroid colors around solar, ~0.1 mag in a*, i − z; Ivezić et al. 2001, AJ 122, 2749),
    // 0.7 mag without color terms (magnitudes then carry the color difference between object and references).
    // Pearson χ² test in flux (not magnitude: at low SNR the magnitude error is asymmetric):
    // χ² = Σ (F_m − F̄)²/(σ_m² + (σ_int·F̄)²), F̄ weighted mean; rejected if P(χ²_{n−1}) < 1 %. The level is
    // stricter than the usual 10⁻³ because the measured σ excludes selection bias (Eddington 1913, MNRAS 73, 359): a
    // noise track at SNR 2–7 only passes where the threshold is high if the noise peak there is high, so the least
    // sensitive frame (z in SDSS) carries the brightest detection — the typical "red" false-positive track.
    const double sint = (color_terms ? 0.2 : 0.7) / 1.0857;
    std::vector<double> fl, sf;
    for (size_t i = 0; i < mag.size() && i < smag.size(); ++i)
        if (std::isfinite(mag[i]) && mag[i] < 90 && std::isfinite(smag[i])) {
            double f = std::pow(10.0, -0.4 * (mag[i] - 20));
            fl.push_back(f);
            sf.push_back(f * smag[i] / 1.0857);
        }
    const int k = int(fl.size()) - 1;
    if (k < 1)
        return true;
    double fb = 0;
    for (int it = 0; it < 3; ++it) {
        double sw = 0, swf = 0;
        for (size_t i = 0; i < fl.size(); ++i) {
            double w = 1 / (sf[i] * sf[i] + std::pow(sint * fb, 2));
            sw += w;
            swf += w * fl[i];
        }
        fb = swf / sw;
    }
    double chi = 0;
    for (size_t i = 0; i < fl.size(); ++i)
        chi += std::pow(fl[i] - fb, 2) / (sf[i] * sf[i] + std::pow(sint * fb, 2));
    // 99 % quantile of χ²_k by the Wilson & Hilferty approximation (1931, PNAS 17, 684), z = 2.326
    double q = 2.0 / (9 * k), lim = k * std::pow(1 - q + 2.326 * std::sqrt(q), 3);
    return chi <= lim;
}

void drop_incoherent_movers(Reduction &R) {
    bool all_color = true;
    for (auto &f : R.frames)
        if (f.solved && !f.stacked)
            all_color = all_color && f.plate.zp_color != 0;
    std::vector<Track> keep;
    for (auto &T : R.movers) {
        std::vector<double> m, sm;
        for (auto &me : T.meas) {
            m.push_back(me.mag);
            sm.push_back(1.0857 / std::max(me.snr_tot > 0 ? me.snr_tot : me.snr, 1.0));
        }
        if (flux_coherent(m, sm, all_color))
            keep.push_back(std::move(T));
    }
    R.movers.swap(keep);
}

bool forced_flux(const Image &im, double x, double y, double rap, double gain, double &flux, double &sigma) {
    const double rin = rap + 3, rout = rap + 8;
    const int R = int(std::ceil(rout)), px = int(std::lround(x)), py = int(std::lround(y));
    if (px - R < 0 || py - R < 0 || px + R >= im.w || py + R >= im.h)
        return false;
    std::vector<double> sky;
    for (int v = py - R; v <= py + R; ++v)
        for (int u = px - R; u <= px + R; ++u) {
            double r2 = (u - x) * (u - x) + (v - y) * (v - y);
            if (r2 >= rin * rin && r2 <= rout * rout)
                sky.push_back(im.pix[size_t(v) * im.w + u]);
        }
    if (sky.size() < 10)
        return false;
    double sb = 0;
    int nsky = 0;
    double b = clipped_median(sky, 3.0, 10, &sb, &nsky);
    double f = 0, npix = 0;
    for (int v = py - int(rap) - 1; v <= py + int(rap) + 1; ++v)
        for (int u = px - int(rap) - 1; u <= px + int(rap) + 1; ++u) {
            int in = 0;
            for (int sy = 0; sy < 5; ++sy)
                for (int sx = 0; sx < 5; ++sx) {
                    double dx = u - 0.4 + 0.2 * sx - x, dy = v - 0.4 + 0.2 * sy - y;
                    in += dx * dx + dy * dy <= rap * rap;
                }
            if (!in)
                continue;
            f += in / 25.0 * (im.pix[size_t(v) * im.w + u] - b);
            npix += in / 25.0;
        }
    flux = f;
    sigma = std::sqrt(std::max(f, 0.0) / (gain > 0 ? gain : 1) + npix * (1 + npix / std::max(nsky, 1)) * sb * sb);
    return sigma > 0;
}

std::vector<Track> link_movers(Reduction &R, const Catalog &cat, const Settings &s) {
    std::vector<Track> out;
    std::vector<int> F; // resolved (non-stacked) frames, in time order
    for (size_t k = 0; k < R.frames.size(); ++k)
        if (R.frames[k].solved && !R.frames[k].stacked)
            F.push_back(int(k));
    std::stable_sort(F.begin(), F.end(), [&](int a, int b) { return R.frames[a].im.jd_mid < R.frames[b].im.jd_mid; });
    const int N = int(F.size());
    if (N < 3)
        return out;
    const DetectParams dp = detect_params(s);
    // overlapping-source separation (Stetson 1987, PASP 99, 191): a free detection with a neighbor ≥ 1.5× brighter
    // within 2.5·FWHM has its centroid contaminated by the neighbor's wing; the position is redone by a joint
    // two-PSF fit in the frame (common background) and adopted if it converges near the original
    for (int n = 0; n < N; ++n) {
        Frame &f = R.frames[F[n]];
        PsfShape psf;
        double apc;
        if (!frame_psf(f, dp.saturation, dp.aperture, psf, apc))
            continue;
        const double fw = 2.3548 * 0.5 * (psf.sx + psf.sy);
        for (size_t i = 0; i < f.det.size(); ++i) {
            if (i < f.typ.size() && (f.typ[i] == 2 || f.typ[i] == 3))
                continue;
            Detection &d = f.det[i];
            int nb = -1;
            double bd = 2.5 * fw;
            for (size_t j = 0; j < f.det.size(); ++j) {
                double dd = std::hypot(f.det[j].x - d.x, f.det[j].y - d.y);
                if (j != i && dd < bd && f.det[j].flux > 1.5 * d.flux) {
                    bd = dd;
                    nb = int(j);
                }
            }
            if (nb < 0 || d.flux <= 0)
                continue;
            std::vector<SceneSource> src(2);
            src[0] = SceneSource{d.x, d.y, d.flux, false, true, false};
            src[1] = SceneSource{f.det[nb].x, f.det[nb].y, f.det[nb].flux, false, true, false};
            double bk = d.bkg;
            int Rw = int(std::ceil(3 * 0.5 * (psf.sx + psf.sy) + bd));
            if (!fit_scene(f.im, d.x, d.y, Rw, psf, psf, src, bk, nullptr) || !(src[0].flux > 0) ||
                std::hypot(src[0].x - d.x, src[0].y - d.y) > 1.5 || std::hypot(src[1].x - f.det[nb].x, src[1].y - f.det[nb].y) > 1.0)
                continue;
            d.x = src[0].x;
            d.y = src[0].y;
        }
    }
    // common standard coordinates (tangent point of the first frame), in ″; magnitude and magnitude σ per detection
    double ra0 = R.frames[F[0]].plate.ra0, dec0 = R.frames[F[0]].plate.dec0;
    std::vector<std::vector<Pt>> X(N);
    std::vector<std::vector<double>> MAG(N), SMAG(N), SIG(N);
    double prms = 0, seeing = 0;
    for (int n = 0; n < N; ++n) {
        Frame &f = R.frames[F[n]];
        // astrometric uncertainty of each detection (″, per axis): centroid σ_psf/SNR (King 1983, PASP 95, 163;
        // Mendez et al. 2014, PASP 126, 798 — Cramér-Rao bound for a Gaussian PSF) with σ_psf = FWHM/2.3548 of the
        // detection itself (≥ frame seeing: extended sources have worse centroids), in quadrature with the plate rms
        // per axis (model, catalog and differential chromatic refraction errors between filters)
        const double sc = 0.5 * (f.plate.scale_x() + f.plate.scale_y());
        const double sp = std::hypot(f.plate.rms_ra, f.plate.rms_dec) / std::sqrt(2.0);
        for (auto &d : f.det) {
            double fw = std::max(d.fwhm * sc, f.seeing);
            SIG[n].push_back(std::hypot(fw / 2.3548 / std::max(d.snr_tot, 1.0), sp));
        }
        for (auto &d : f.det) {
            double ra, dec;
            Pt p;
            f.plate.pix_to_radec(d.x, d.y, ra, dec);
            radec_to_std(ra, dec, ra0, dec0, p.x, p.y);
            X[n].push_back(p);
            MAG[n].push_back(d.flux > 0 ? f.plate.zp - 2.5 * std::log10(d.flux) : 99);
            SMAG[n].push_back(1.0857 / std::max(d.snr_tot, 1.0));
        }
        prms = std::max(prms, std::hypot(f.plate.rms_ra, f.plate.rms_dec));
        seeing = std::max(seeing, f.seeing);
    }
    const double tmin = R.frames[F[0]].im.jd_mid, span_min = (R.frames[F[N - 1]].im.jd_mid - tmin) * 1440;
    auto tmin_of = [&](int n) { return (R.frames[F[n]].im.jd_mid - tmin) * 1440; };
    // compatible brightness between two measurements: |Δm| < 1 mag (inter-filter color, variability) + 3σ of the
    // difference
    auto same_mag = [](double m1, double s1, double m2, double s2) { return std::fabs(m1 - m2) < 1.0 + 3 * std::hypot(s1, s2); };

    auto get = [&](const char *k, double def) { return s.str(k).empty() ? def : s.num(k); };
    // ---- fixed sources
    // (a) plate references; (b) match with a catalog star of compatible brightness (static sky,
    // Denneau et al. 2013 §4); (c) match with a compatible-brightness detection in frames far apart in time
    // (|Δt| ≥ arc/4: in near-simultaneous frames a slow object also looks fixed) in most of those frames.
    // "Match" is statistical: separation < k·√(σ_i² + σ_j²) with k = motion_sigma (default 5; for a fixed source the
    // separation is Rayleigh-distributed and P(> 5σ) = e^−12.5 ≈ 4·10⁻⁶). The same k defines resolved motion below, so
    // "fixed" and "moving" are two faces of a single test — no floors in ″, which on short arcs with sharp images
    // (SDSS: 5 min, 1.4″) erased objects moving 1–4″.
    const double ksig = std::max(get("motion_sigma", 5.0), 1.0);
    double sigmax = 0;
    for (int n = 0; n < N; ++n)
        for (double v : SIG[n])
            sigmax = std::max(sigmax, v);
    const double rfix = std::max(ksig * std::sqrt(2.0) * sigmax, 1e-3); // largest possible radius (grid cell)
    std::vector<Pt> C(cat.size());
    for (size_t j = 0; j < cat.size(); ++j)
        radec_to_std(cat.ra[j], cat.dec[j], ra0, dec0, C[j].x, C[j].y);
    Grid gc;
    gc.build(C, rfix);
    std::vector<Grid> gx(N);
    for (int n = 0; n < N; ++n)
        gx[n].build(X[n], rfix);
    const double tau = span_min / 4;
    std::vector<std::vector<char>> star(N);
    for (int n = 0; n < N; ++n) {
        Frame &f = R.frames[F[n]];
        star[n].assign(X[n].size(), 0);
        for (size_t i = 0; i < X[n].size(); ++i) {
            const Pt &p = X[n][i];
            // typ 3: point fragment of a trail (deferred search) — excluded from linking, like a fixed source
            bool st = i < f.typ.size() && (f.typ[i] == 2 || f.typ[i] == 3);
            if (!st)
                gc.near(p.x, p.y, [&](int j) {
                    if (!st && std::hypot(C[j].x - p.x, C[j].y - p.y) < ksig * SIG[n][i] && std::isfinite(cat.mag[j]) &&
                        same_mag(MAG[n][i], SMAG[n][i], f.plate.star_mag(cat.mag[j], cat.has_color() ? cat.color[j] : NAN), 0.0))
                        st = true;
                });
            if (!st) {
                int elig = 0, cnt = 0;
                for (int m = 0; m < N; ++m) {
                    if (m == n || std::fabs(tmin_of(m) - tmin_of(n)) < tau)
                        continue;
                    ++elig;
                    bool hit = false;
                    gx[m].near(p.x, p.y, [&](int k) {
                        if (!hit && std::hypot(X[m][k].x - p.x, X[m][k].y - p.y) < ksig * std::hypot(SIG[n][i], SIG[m][k]) &&
                            same_mag(MAG[n][i], SMAG[n][i], MAG[m][k], SMAG[m][k]))
                            hit = true;
                    });
                    cnt += hit;
                }
                st = elig > 0 && 2 * cnt >= elig;
            }
            star[n][i] = st;
            if (st && i < f.typ.size() && f.typ[i] == 0)
                f.typ[i] = 1;
        }
    }

    // ---- fixed-source test by forced photometry: at the position of each track detection, in the other frames,
    // measure the aperture flux; if there is significant flux (> 3σ) compatible with the detection's (0.4–2.5× after
    // the zero-point difference), the detection is a fixed source. On short arcs the object is still near in those
    // frames (less than one or two FWHM): the fraction g of its light falling in the aperture (shifted Gaussian, from
    // the fitted track) is subtracted using the flux measured at the predicted position, and the frame only counts when
    // the two hypotheses separate — g ≤ 1/2 and a fixed source of 0.4× the flux would be detected at 3σ. Without any
    // such frame, the decision falls to the significant-motion test (|v| ≥ k·σ_v) and the fixed-source flagging above.
    // link_stat_blend (default 1): with ≥ 3 detections, a detection ≥ 0.5 mag brighter than the track median is
    // treated as a blend with a fixed source (object passing over a star) and is not evidence of a stationary
    // source — the star is there in every frame, the object is not; the other detections decide
    const bool blend_skip = get("link_stat_blend", 1) != 0;
    std::unordered_map<long long, std::array<double, 5>> fcache;
    auto stationary = [&](const std::vector<int> &pick, double ax, double bx, double ay, double by) {
        std::vector<double> mg;
        for (int m = 0; m < N; ++m)
            if (pick[m] >= 0)
                mg.push_back(MAG[m][pick[m]]);
        const double mmed = mg.size() >= 3 ? median(mg) : -99;
        for (int m = 0; m < N; ++m) {
            if (pick[m] < 0)
                continue;
            const Frame &fm = R.frames[F[m]];
            const Detection &d = fm.det[pick[m]];
            if (!(d.flux > 0) || (blend_skip && MAG[m][pick[m]] < mmed - 0.5))
                continue;
            double ra, dec;
            fm.plate.pix_to_radec(d.x, d.y, ra, dec);
            const Pt &p = X[m][pick[m]];
            for (int k = 0; k < N; ++k) {
                if (k == m)
                    continue;
                const Frame &fk = R.frames[F[k]];
                double t = tmin_of(k), ox = ax + bx * t, oy = ay + by * t;
                double x = d.x - fm.offx + fk.offx, y = d.y - fm.offy + fk.offy;
                if (!radec_to_pix(fk.plate, ra, dec, x, y))
                    continue;
                // small aperture (≈ 0.67·FWHM, near-optimal for a point source) to reduce contamination by neighbors;
                // expected flux of a fixed source: the one measured in m with the same aperture, corrected by the
                // zero-point difference and by the Gaussian encircled energy 1 − exp(−r²/2σ²) in each frame
                double sxm = 0.5 * (fm.plate.scale_x() + fm.plate.scale_y()), sxk = 0.5 * (fk.plate.scale_x() + fk.plate.scale_y());
                double fwm = std::max(fm.seeing / sxm, 1.0), fwk = std::max(fk.seeing / sxk, 1.0);
                double rap = std::max(1.5, 0.67 * std::max(fwm, fwk));
                double g = offset_aperture_frac(std::hypot(ox - p.x, oy - p.y) / sxk, rap, fwk / 2.3548);
                if (g > 0.5)
                    continue;
                double fl, sg, f0, s0, ft = 0, st = 0;
                // the two forced photometries do not depend on the track: cache keyed by (frame m, detection, frame k)
                const long long fkey = ((long long)(m * N + k) << 32) | (long long)(unsigned)pick[m];
                auto fit = fcache.find(fkey);
                if (fit == fcache.end()) {
                    std::array<double, 5> v{0, 0, 0, 0, 0};
                    v[4] = forced_flux(fk.im, x, y, rap, dp.gain, v[0], v[1]) && forced_flux(fm.im, d.x, d.y, rap, dp.gain, v[2], v[3]);
                    fit = fcache.emplace(fkey, v).first;
                }
                fl = fit->second[0];
                sg = fit->second[1];
                f0 = fit->second[2];
                s0 = fit->second[3];
                (void)s0;
                if (fit->second[4] == 0 || !(f0 > 0))
                    continue;
                if (g > 1e-3) {
                    double ora, odec, u = x, v = y;
                    std_to_radec(ox, oy, ra0, dec0, ora, odec);
                    if (!radec_to_pix(fk.plate, ora, odec, u, v) || !forced_flux(fk.im, u, v, rap, dp.gain, ft, st))
                        continue;
                    fl -= g * std::max(ft, 0.0);
                    sg = std::hypot(sg, g * st);
                }
                auto ee = [&](double fw) {
                    double sg2 = fw / 2.3548;
                    return 1 - std::exp(-rap * rap / (2 * sg2 * sg2));
                };
                double fexp = f0 * std::pow(10.0, 0.4 * (fk.plate.zp - fm.plate.zp)) * ee(fwk) / ee(fwm);
                if (0.4 * fexp < 3 * sg)
                    continue;
                if (fl > 3 * sg && fl > 0.4 * fexp && fl < 2.5 * fexp)
                    return true;
            }
        }
        return false;
    };
    // velocity limits (″/min): MinRate/MaxRate if set
    // (the physical minimum is not a fixed velocity: each anchor pair requires displacement ≥ k·σ of the pair, and the
    // fitted track requires |v| ≥ k·σ_v)
    double vmin = get("min_rate", 0), vmax = get("max_rate", 100.0);
    const double tol = std::max(1.5, 3 * prms);
    std::vector<std::vector<char>> used(N);
    for (int n = 0; n < N; ++n)
        used[n].assign(X[n].size(), 0);
    // anchor pairs (a, b) with at least one frame between them, from the largest interval to the smallest
    std::vector<std::pair<int, int>> anchors;
    for (int a = 0; a < N; ++a)
        for (int b = a + 2; b < N; ++b)
            anchors.push_back({a, b});
    std::stable_sort(anchors.begin(), anchors.end(), [](auto &p, auto &q) { return (p.second - p.first) > (q.second - q.first); });
    // tracks with ≥ 3 detections; with ≥ 4 frames one miss between the anchors is allowed (source on a star, noise)
    const int min_frames = 3, max_missing = N >= 4 ? 1 : 0;
    std::vector<Grid> gfree(N);
    for (int n = 0; n < N; ++n)
        gfree[n].build(X[n], tol);
    // density of free (non-fixed) detections per frame (″⁻²) and number of free ones
    std::vector<double> rho(N), nfree(N);
    for (int n = 0; n < N; ++n) {
        const Frame &f = R.frames[F[n]];
        nfree[n] = double(std::count(star[n].begin(), star[n].end(), 0));
        rho[n] = nfree[n] / std::max(1.0, f.im.w * f.plate.scale_x() * f.im.h * f.plate.scale_y());
    }
    auto f_area = [&](int n) {
        const Frame &f = R.frames[F[n]];
        return f.im.w * f.plate.scale_x() * f.im.h * f.plate.scale_y();
    };
    // diagnostic: link_probe = "ra dec jd dra_h ddec_h" (known object, ″/h) → per frame the detection closest to the
    // predicted position (within 10″) and, in the linking test chain, the first test that rejects it ("link probe"
    // lines)
    std::vector<int> probe(N, -1);
    bool probing = false;
    char pb[400];
    if (!s.str("link_probe").empty()) {
        double pra, pde, pjd, pvr, pvd;
        if (std::sscanf(s.str("link_probe").c_str(), "%lf %lf %lf %lf %lf", &pra, &pde, &pjd, &pvr, &pvd) == 5) {
            probing = true;
            for (int n = 0; n < N; ++n) {
                const Frame &f = R.frames[F[n]];
                double dt = (f.im.jd_mid - pjd) * 24, ra, dec;
                Pt q;
                std_to_radec(pvr * dt, pvd * dt, pra, pde, ra, dec);
                radec_to_std(ra, dec, ra0, dec0, q.x, q.y);
                double bd = 10;
                for (size_t i = 0; i < X[n].size(); ++i) {
                    double d = std::hypot(X[n][i].x - q.x, X[n][i].y - q.y);
                    if (d < bd) {
                        bd = d;
                        probe[n] = int(i);
                    }
                }
                if (probe[n] >= 0) {
                    const Detection &d = f.det[probe[n]];
                    std::snprintf(
                        pb, sizeof pb,
                        "link probe: frame %d: det %d a %.2f\" (fixed %d, typ %d), mag %.2f, SNR %.1f, FWHM %.2f px, trail %.1f px; "
                        "free %.0f, ρ %.2e/″²",
                        n + 1, probe[n], bd, int(star[n][probe[n]]), probe[n] < int(f.typ.size()) ? f.typ[probe[n]] : -1, MAG[n][probe[n]],
                        d.snr_tot, d.fwhm, d.trail_len, nfree[n], rho[n]);
                } else
                    std::snprintf(pb, sizeof pb, "link probe: frame %d: no detection within 10\"; free %.0f", n + 1, nfree[n]);
                R.log.push_back(pb);
            }
        }
    }
    auto why = [&](bool on, int a, int b, const char *msg, double v1 = NAN, double v2 = NAN) {
        if (!on)
            return;
        std::snprintf(pb, sizeof pb, "link probe: anchors %d–%d: rejected by %s (%.3g, %.3g)", a + 1, b + 1, msg, v1, v2);
        R.log.push_back(pb);
    };
    // expected number of chance tracks as good as the candidate (Kubica et al. 2007 §3; Denneau et al. 2013 §5):
    // N_a (free anchors) × N_b (partners in the velocity annulus) × Π_m ρ_m·π·r_m² in the verification frames, with
    // r_m = fit residual (minimum 3σ of the detection: the region where the prediction would accept it). The track is
    // rejected if this number exceeds false_alarm (default 0.05).
    const double alpha = get("false_alarm", 0.05);
    // Conditioning (link_mag_e, default 0 — on the NEA set it found 2019 AN5 but brought in a chance track): N_a, N_b and
    // ρ_m count only free detections with brightness compatible (same_mag) with the anchor — the brightness test is
    // already applied to each hypothesis, so a chance track needs compatible partners; and (link_trail_e, default 1) with
    // an elongated anchor, the velocity annulus shrinks to the sector the trail allows (|ΔL| ≤ max(2 px, L/2), |Δθ| ≤
    // 20°, both directions), as trail_coherent enforces.
    const bool ext_rule = get("link_extended", 1) != 0;
    const bool mag_e = get("link_mag_e", 0) != 0, trail_e = get("link_trail_e", 1) != 0;
    // cache by (frame, anchor): the count depends only on the anchor magnitude (O(N_free) per key)
    std::unordered_map<long long, double> ncache;
    auto ncompat = [&](int m, int a, int i) {
        if (!mag_e)
            return nfree[m];
        const long long key = ((long long)(m * N + a) << 32) | (long long)(unsigned)i;
        auto it = ncache.find(key);
        if (it != ncache.end())
            return it->second;
        const double mg = MAG[a][i], sg = SMAG[a][i];
        double c = 0;
        for (size_t k = 0; k < X[m].size(); ++k)
            c += !star[m][k] && same_mag(mg, sg, MAG[m][k], SMAG[m][k]);
        ncache[key] = c;
        return c;
    };
    auto chance = [&](const std::vector<int> &pick, int a, int b, double ax, double bx, double ay, double by) {
        double dt = tmin_of(b) - tmin_of(a), dmin = std::max(vmin * dt, ksig * std::hypot(SIG[a][pick[a]], SIG[b][pick[b]]));
        double v1 = vmax, frac = 1;
        const Frame &fa = R.frames[F[a]];
        const Detection &da = fa.det[pick[a]];
        if (trail_e && da.trail_len > 0 && fa.im.exptime > 0) {
            const double sc = 0.5 * (fa.plate.scale_x() + fa.plate.scale_y()), tl = std::max(2.0, 0.5 * da.trail_len);
            dmin = std::max(dmin, std::max(da.trail_len - tl, 0.0) * sc * 60 / fa.im.exptime * dt);
            v1 = std::min(v1, (da.trail_len + tl) * sc * 60 / fa.im.exptime);
            frac = 2 * 40.0 / 360;
        }
        const double area = std::max(1.0, f_area(b));
        double nb = std::min(nfree[b], ncompat(b, a, pick[a]) / area * frac * PI * std::max(std::pow(v1 * dt, 2) - dmin * dmin, 0.0));
        double e = std::max(ncompat(a, a, pick[a]), 1.0) * std::max(nb, 1.0);
        for (int m = 0; m < N; ++m)
            if (pick[m] >= 0 && m != a && m != b) {
                double t = tmin_of(m);
                double r = std::max(3 * SIG[m][pick[m]], std::hypot(X[m][pick[m]].x - (ax + bx * t), X[m][pick[m]].y - (ay + by * t)));
                e *= std::min(1.0, ncompat(m, a, pick[a]) / std::max(1.0, f_area(m)) * PI * r * r);
            }
        return e;
    };
    // trail coherence (Vereš et al. 2012; Dawson et al. 2016): for an elongated detection, the measured length and
    // direction must agree with L = v·t_exp and with the track velocity direction in the frame (pixels via the plate).
    // Tolerances: |ΔL| ≤ max(2 px, L/2) and |Δθ| ≤ 20° (measured length is noisy at low SNR). A chance track
    // of three noise detections rarely has trails aligned with its own motion.
    auto trail_coherent = [&](const std::vector<int> &pick, double bx, double by) {
        for (int m = 0; m < N; ++m) {
            if (pick[m] < 0)
                continue;
            const Frame &f = R.frames[F[m]];
            const Detection &d = f.det[pick[m]];
            if (!(d.trail_len > 0) || !(f.im.exptime > 0))
                continue;
            // Jacobian (ξ, η)/(x, y) by finite differences at the detection position
            double J[4];
            plate_jacobian(f.plate, d.x, d.y, J, ra0, dec0);
            double j11 = J[0], j12 = J[1], j21 = J[2], j22 = J[3], det = j11 * j22 - j12 * j21;
            if (!(std::fabs(det) > 0))
                continue;
            double px = (j22 * bx - j12 * by) / det, py = (-j21 * bx + j11 * by) / det; // px/min
            double lp = std::hypot(px, py) * f.im.exptime / 60;
            if (std::fabs(d.trail_len - lp) > std::max(2.0, 0.5 * lp))
                return false;
            double da = std::fmod(std::fabs(std::atan2(py, px) - d.trail_ang), PI);
            if (std::min(da, PI - da) > 20 * PI / 180)
                return false;
        }
        return true;
    };
    long n_chance = 0, n_stat = 0;
    for (auto [a, b] : anchors) {
        double ta = tmin_of(a), tb = tmin_of(b);
        if (tb - ta <= 0)
            continue;
        for (size_t i = 0; i < X[a].size(); ++i) {
            const bool pa = probing && int(i) == probe[a];
            if (star[a][i] || used[a][i]) {
                why(pa && probe[b] >= 0, a, b, star[a][i] ? "fixed anchor" : "anchor already used");
                continue;
            }
            struct Best {
                std::vector<int> pick;
                int nf = 0;
                double rms = 1e30;
                double vx = 0, vy = 0;
            } best;
            for (size_t j = 0; j < X[b].size(); ++j) {
                const bool pp = pa && int(j) == probe[b];
                if (star[b][j] || used[b][j]) {
                    why(pp, a, b, star[b][j] ? "fixed 2nd anchor" : "2nd anchor already used");
                    continue;
                }
                double vx = (X[b][j].x - X[a][i].x) / (tb - ta), vy = (X[b][j].y - X[a][i].y) / (tb - ta);
                double v = std::hypot(vx, vy);
                // motion resolved between the anchors: displacement ≥ k·σ of the pair
                // extended source (FWHM > 1.3·seeing: galaxy, blend): the centroid shifts with filter and depth —
                // additionally requires displacement ≥ 0.5·its FWHM (link_extended, default 1)
                double dext = 0;
                if (ext_rule) {
                    const Frame &fa = R.frames[F[a]], &fb = R.frames[F[b]];
                    const double wa = fa.det[i].fwhm * 0.5 * (fa.plate.scale_x() + fa.plate.scale_y());
                    const double wb = fb.det[j].fwhm * 0.5 * (fb.plate.scale_x() + fb.plate.scale_y());
                    if (wa > 1.3 * fa.seeing || wb > 1.3 * fb.seeing)
                        dext = 0.5 * std::max(wa, wb);
                }
                if (v < vmin || v > vmax || v * (tb - ta) < std::max(dext, ksig * std::hypot(SIG[a][i], SIG[b][j]))) {
                    why(pp, a, b, "motion unresolved/out of range", v, v * (tb - ta));
                    continue;
                }
                // verification in all other frames by linear prediction
                std::vector<int> pick(N, -1);
                pick[a] = int(i);
                pick[b] = int(j);
                int missing = 0;
                for (int m = 0; m < N; ++m) {
                    if (m == a || m == b)
                        continue;
                    double t = tmin_of(m), px = X[a][i].x + vx * (t - ta), py = X[a][i].y + vy * (t - ta);
                    int bk = -1;
                    double bd = tol;
                    gfree[m].near(px, py, [&](int k) {
                        if (star[m][k] || used[m][k])
                            return;
                        double d = std::hypot(X[m][k].x - px, X[m][k].y - py);
                        if (d < bd) {
                            bd = d;
                            bk = k;
                        }
                    });
                    pick[m] = bk;
                    if (bk < 0 && m > a && m < b)
                        ++missing;
                }
                int nf = 0;
                for (int m = 0; m < N; ++m)
                    nf += pick[m] >= 0;
                if (missing > max_missing || nf < min_frames) {
                    if (pp) {
                        std::string l = "link probe: verification:";
                        for (int m = 0; m < N; ++m) {
                            double t = tmin_of(m), px = X[a][i].x + vx * (t - ta), py = X[a][i].y + vy * (t - ta);
                            std::snprintf(pb, sizeof pb, " [q%d pick %d probe %d Δ %.2f\"]", m + 1, pick[m], probe[m],
                                          probe[m] >= 0 ? std::hypot(X[m][probe[m]].x - px, X[m][probe[m]].y - py) : NAN);
                            l += pb;
                        }
                        R.log.push_back(l);
                    }
                    why(pp, a, b, "misses/frames within tolerance", missing, tol);
                    continue;
                }
                // linear fit ξ(t), η(t)
                double st = 0, sx = 0, sy = 0, stt = 0, stx = 0, sty = 0;
                for (int m = 0; m < N; ++m)
                    if (pick[m] >= 0) {
                        double t = tmin_of(m);
                        const Pt &p = X[m][pick[m]];
                        st += t;
                        stt += t * t;
                        sx += p.x;
                        sy += p.y;
                        stx += t * p.x;
                        sty += t * p.y;
                    }
                double den = nf * stt - st * st;
                if (den <= 0)
                    continue;
                double bx = (nf * stx - st * sx) / den, ax = (sx - bx * st) / nf;
                double by = (nf * sty - st * sy) / den, ay = (sy - by * st) / nf;
                double r2 = 0;
                for (int m = 0; m < N; ++m)
                    if (pick[m] >= 0) {
                        double t = tmin_of(m);
                        const Pt &p = X[m][pick[m]];
                        r2 += std::pow(p.x - (ax + bx * t), 2) + std::pow(p.y - (ay + by * t), 2);
                    }
                double rms = std::sqrt(r2 / std::max(1, nf - 2));
                const double rms_max = std::max(1.0, 1.5 * prms);
                // significant motion: |v| ≥ k·σ_v, σ_v = σ/√Σ(t − t̄)² with σ = max(mean detection σ, fit
                // rms) — likelihood-ratio test motion vs. fixed source (χ² with 2 d.o.f. for the velocity)
                {
                    double s2 = 0;
                    for (int m = 0; m < N; ++m)
                        if (pick[m] >= 0)
                            s2 += SIG[m][pick[m]] * SIG[m][pick[m]];
                    double sv = std::max(std::sqrt(s2 / nf), rms) / std::sqrt(den / nf);
                    if (std::hypot(bx, by) < ksig * sv) {
                        why(pp, a, b, "motion not significant", std::hypot(bx, by), ksig * sv);
                        continue;
                    }
                }
                // consistent brightness across frames: max and min magnitude differ by less than 1 mag (rotational
                // light curve over minutes–hours and asteroid color spread around solar, ≲ 0.3 mag; Ivezić et al. 2001,
                // AJ 122, 2749 — magnitudes are already referred to the solar color by each frame's color term) plus 3σ
                // of the difference (photometric noise; same form as same_mag)
                double mlo = 1e9, mhi = -1e9, slo = 0, shi = 0;
                for (int m = 0; m < N; ++m)
                    if (pick[m] >= 0) {
                        if (MAG[m][pick[m]] < mlo)
                            mlo = MAG[m][pick[m]], slo = SMAG[m][pick[m]];
                        if (MAG[m][pick[m]] > mhi)
                            mhi = MAG[m][pick[m]], shi = SMAG[m][pick[m]];
                    }
                if (!same_mag(mlo, slo, mhi, shi) || rms > rms_max) {
                    why(pp, a, b, !same_mag(mlo, slo, mhi, shi) ? "inconsistent brightness (mlo, mhi)" : "rms",
                        !same_mag(mlo, slo, mhi, shi) ? mlo : rms, !same_mag(mlo, slo, mhi, shi) ? mhi : rms_max);
                    continue;
                }
                if (!trail_coherent(pick, bx, by)) {
                    why(pp, a, b, "inconsistent trail");
                    continue;
                }
                // competing hypothesis: more detections wins; on a tie, lower rms
                if (nf < best.nf || (nf == best.nf && rms >= best.rms))
                    continue;
                // Bonferroni correction: the same track can arise from any of the tested anchor pairs
                double e = chance(pick, a, b, ax, bx, ay, by) * double(anchors.size());
                ++n_chance;
                if (e <= alpha)
                    ++n_stat;
                if (e > alpha || stationary(pick, ax, bx, ay, by)) {
                    why(pp, a, b, e > alpha ? "chance alignment E (E, α)" : "stationary source", e, alpha);
                    continue;
                }
                if (pp) {
                    std::snprintf(pb, sizeof pb, "link probe: anchors %d–%d: accepted (nf %d, rms %.2f\", E %.2g)", a + 1, b + 1, nf, rms,
                                  e);
                    R.log.push_back(pb);
                }
                best.pick = pick;
                best.nf = nf;
                best.rms = rms;
                best.vx = bx;
                best.vy = by;
            }
            if (best.pick.empty())
                continue;
            Track T;
            T.idx.assign(R.frames.size(), -1);
            for (int m = 0; m < N; ++m) {
                int k = best.pick[m];
                if (k < 0)
                    continue;
                used[m][k] = 1;
                Frame &f = R.frames[F[m]];
                const Detection &d = f.det[k];
                T.idx[F[m]] = k;
                Measurement me;
                me.frame = F[m];
                me.x = d.x;
                me.y = d.y;
                f.plate.pix_to_radec(d.x, d.y, me.ra, me.dec);
                me.flux = d.flux;
                me.mag = d.flux > 0 ? f.plate.zp - 2.5 * std::log10(d.flux) : NAN;
                me.snr = d.snr;
                me.snr_tot = d.snr_tot;
                me.fwhm_as = d.fwhm * 0.5 * (f.plate.scale_x() + f.plate.scale_y());
                me.jd = f.im.jd_mid;
                me.note = "";
                if (d.trail_len > 0) {
                    me.trail_len_px = d.trail_len;
                    me.trail_dx = std::cos(d.trail_ang);
                    me.trail_dy = std::sin(d.trail_ang);
                }
                T.meas.push_back(me);
            }
            T.rate = std::hypot(best.vx, best.vy);
            T.pa = std::fmod(std::atan2(best.vx, best.vy) * 180 / PI + 360, 360);
            out.push_back(T);
        }
    }
    if (get("trail_debug", 0) != 0) {
        char lb[160];
        std::snprintf(lb, sizeof lb, "linking trails: %ld hypotheses in the E test, %ld in the stationary test", n_chance, n_stat);
        R.log.push_back(lb);
    }
    return out;
}

} // namespace internal
} // namespace astro
