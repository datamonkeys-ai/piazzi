// Internal reduction utilities (robust statistics, least squares, pattern matching).
// Independent implementation from the literature:
//   - Robust statistics (median, MAD, sigma clipping): Hoaglin, Mosteller & Tukey (1983), "Understanding Robust and
//     Exploratory Data Analysis"; MAD factor 1.4826 for Gaussian σ.
//   - QR least squares (Householder): Golub & Van Loan, "Matrix Computations", §5.3.
//   - Triangle matching: Groth (1986), AJ 91, 1244; Valdes et al. (1995), PASP 107, 1119.
//   - RANSAC: Fischler & Bolles (1981), Comm. ACM 24, 381.
#pragma once
#include "piazzi/reduce.hpp"
#include <vector>

namespace astro {
namespace internal {

// Median of v (reorders v). NaN if empty.
double median(std::vector<double> &v);
// Sigma-clipped mean/median: iteratively drops |v - median| > k·σ (σ from MAD) until convergence.
// Returns the final median; sigma (std. dev. of kept values) and n (kept count) are optional.
double clipped_median(std::vector<double> v, double k, int iters, double *sigma = nullptr, int *n = nullptr);
// σ robusto pelo MAD (1,4826·MAD).
double mad_sigma(std::vector<double> v);

// Solves min ||A·c - b|| (A has m rows and n columns, row-major) by Householder QR with normalized columns.
// false if the system is singular or m < n.
bool lstsq(const std::vector<double> &A, int m, int n, const std::vector<double> &b, std::vector<double> &c);

// Per-axis centroid σ (″) of a Gaussian-profile source: Cramér–Rao bound
//   σ² = (σ_PSF²/N)·(1 + 8π·σ_PSF²·B/N)   (N = source electrons, B = per-pixel sky variance in e⁻²)
// (King 1983, PASP 95, 163; Mighell 2005, MNRAS 361, 861, eq. 46). scale = ″/px; gain = e⁻/ADU.
double centroid_sigma(const Detection &d, double scale, double gain);
// Inverse of an n×n symmetric positive matrix (Gauss–Jordan with pivoting); false if singular.
bool invert(std::vector<double> &A, int n);

// Number of terms of the plate polynomial of order `order` (1: 3, 2: 6, 3: 10).
int plate_nterms(int order);

// 2D affine transform: u = a0 + a1·x + a2·y, v = b0 + b1·x + b2·y.
struct Affine {
    double a[3] = {0, 1, 0}, b[3] = {0, 0, 1};
    void apply(double x, double y, double &u, double &v) const {
        u = a[0] + a[1] * x + a[2] * y;
        v = b[0] + b[1] * x + b[2] * y;
    }
};

struct Pt {
    double x, y;
};

// Triangle matching between points P (sorted by brightness, already in units comparable to Q) and Q.
// Finds a similarity transform (rotation + scale ≈ 1, with or without mirroring) P→Q by voting in the space
// of invariant ratios and RANSAC consensus; tol = acceptance radius of a pair (Q units).
// Returns the number of agreeing pairs (0 on failure) and the transform in T.
int match_triangles(const std::vector<Pt> &P, const std::vector<Pt> &Q, double tol, double scale_tol, Affine &T);

// Fits an affine transform by least squares to the pairs (P[i] → Q[i]).
bool fit_affine(const std::vector<Pt> &P, const std::vector<Pt> &Q, Affine &T);

// Solves the plate of a frame: catalog matching and plate-constant fit (plate_solve.cpp + reduce.cpp).
struct SolveCtx {
    double ra0, dec0;      // ponto tangente (centro nominal)
    double sx, sy;         // nominal scale ″/px in x and y
    double pa;             // nominal position angle (degrees)
    double pointing;       // pointing uncertainty (′)
    double match_rad;      // final match radius (″)
    double scale_tol;      // relative scale tolerance
    int order;             // plate polynomial order (0 = automatic via k-fold cross-validation)
    bool weighted = true;  // weights 1/σ², σ² = σ_centroid² + σ_cat² + floor²
    double cat_err = 0.01; // σ_cat (″) when the catalog has no errors
    double floor = 0.10;   // piso (″)
    bool color = true;     // per-axis color term in the plate fit (plate_color; Pier et al. 2003)
    double saturation;     // ADU
    double gain = 1.0;     // e-/ADU (centroid error model)
};
// Matches frame detections to the catalog (cat already projected to xi/eta ″) and returns the initial transform
// pixel→(xi, eta). false if there is no solution.
bool initial_solution(const Frame &f, const std::vector<Pt> &cat_std, const Catalog &cat, const SolveCtx &c, Affine &T);
// Plate fit (order c.order) with iterative nearest-neighbor matching and 3σ rejection; fills f.plate.
bool fit_plate(Frame &f, const std::vector<Pt> &cat_std, const Catalog &cat, const SolveCtx &c, const Affine &T0);
// Photometry: zero point from the sigma-clipped median of the references.
void fit_zero_point(Frame &f, const Catalog &cat, double saturation);
// Moving objects (linking.cpp). cat: full catalog (static sky) used to discard known stars.
std::vector<Track> link_movers(Reduction &R, const Catalog &cat, const Settings &s);
// Coherent-flux χ² test (magnitudes referred to solar color, with their σ) across a track's measurements;
// color_terms = all frames have a photometric color term. Implemented in linking.cpp.
bool flux_coherent(const std::vector<double> &mag, const std::vector<double> &smag, bool color_terms);
// Removes from R.movers the tracks whose final measurements (after refinement) fail flux_coherent.
void drop_incoherent_movers(Reduction &R);
// Forced circular-aperture photometry (radius rap px) at (x, y) with background annulus rap+3…rap+8; σ from CCD equation.
bool forced_flux(const Image &im, double x, double y, double rap, double gain, double &flux, double &sigma);

// ---- PSF/trail fitting (psf_fit.cpp)
// PSF shape: axis-aligned elliptical Gaussian (σx, σy in px) convolved with a segment of length len (px)
// along unit direction (dx, dy). len = 0: star PSF.
struct PsfShape {
    double sx = 1, sy = 1, len = 0, dx = 1, dy = 0;
    // pixel integral approximated by a single sample with σ² + 1/12 (1 px box ≈ Gaussian of variance 1/12): ~9× faster,
    // accurate to < 1 % for σ ≥ 1.2 px; used only in the trail search
    bool fast = false;
};
// Unit-flux point profile at offset (u, v) px from the center (segment center = mid-exposure).
double psf_profile(const PsfShape &P, double u, double v);
// The same profile integrated over the pixel centered at (u, v) (3×3 midpoint rule).
double psf_pixel(const PsfShape &P, double u, double v);
// A scene source: target (trail = true, shape `target`) or neighbor star (shape `star`). Input: guess; output: fit.
struct SceneSource {
    double x = 0, y = 0, flux = 0;
    bool trail = false, free_pos = true, free_shape = false;
    double ex = NAN, ey = NAN, exy = NAN, eflux = NAN, fit_sx = 0, fit_sy = 0; // fit σ (and x–y covariance); shape (free_shape)
    // trail (trail = true) with free length and direction: length has a Gaussian prior around target.len
    // (σ len_sigma px; prior residual scaled by per-pixel noise pix_noise), direction has no prior.
    bool free_len = false;
    double len_sigma = 1, pix_noise = 1, fit_len = 0, fit_dx = 1, fit_dy = 0;
    // gain > 0 (e⁻/ADU): ex/ey/exy covariance via the "sandwich" A⁻¹·(Jᵀ V J)·A⁻¹ (A = JᵀJ) of the unweighted estimator,
    // with V the per-pixel model variance (σ_sky² = pix_noise² + signal/gain, Poisson) instead of s²·A⁻¹ with
    // uniform s² = χ²/(n − p) — the latter treats sky and trail pixels as equally noisy
    double gain = 0;
};
// Joint LM fit of a constant background + sources in a (2R+1)² window centered at (cx, cy). bkg: guess/output.
// max_iter: Levenberg–Marquardt iterations (the trail search uses fewer). capsule > 0: only pixels ≤ capsule px from
// the target segment (center (cx, cy), direction of `target`, half-length target.len/2 + capsule, slack for the
// free length) — a tight window for long trails.
bool fit_scene(const Image &im, double cx, double cy, int R, const PsfShape &target, const PsfShape &star, std::vector<SceneSource> &src,
               double &bkg, double *rms = nullptr, int max_iter = 80, double capsule = 0);
// Frame PSF (median σx, σy of isolated references) and aperture correction (aperture flux / PSF flux).
bool frame_psf(const Frame &f, double saturation, int aperture, PsfShape &psf, double &apcorr);
// Re-measures linked objects with the trail model (Vereš et al. 2012) and a joint fit of neighbor stars
// (Stetson 1987). Settings: trail_fit = fit positions (−1 automatic [default]: trail if v·t_exp ≥ 1.5·FWHM;
// 0 off, 1 trail, 2 star PSF);
// psf_phot = magnitude from fitted flux (1 [default]) or aperture (0); trail_deblend = fit neighbors (1/0).
void refine_movers(Reduction &R, const Catalog &cat, const Settings &s);

// ---- stack search and 2-detection tracks (stack_search.cpp)
// Mesh background (32 px, sigma-clipped median, 3×3 median filter, bilinear) and per-pixel sky σ.
void background_mesh(const Image &im, std::vector<float> &B, std::vector<float> &S);
// Matched-filter SNR map of a frame: u = ψ/√φ, ψ = K ⊗ (D·m/σ²), φ = K² ⊗ (m/σ²) (Whidden et al. 2019), with
// static sources masked (m = 0); sqrt_phi = typical √φ (sensitivity); frame psf/apc.
struct MFMap {
    int w = 0, h = 0;
    std::vector<float> u, us, bkg, sig; // us: u with the additional stacking mask (empty = u)
    PsfShape psf;
    double apc = 1, fwhm = 2.5, noise = 1, sqrt_phi = 0;
    bool ok = false;
    double at(double x, double y) const;   // u bilinear (0 fora / mascarado)
    double at_s(double x, double y) const; // us bilinear
};
bool build_mf(const Frame &f, const Catalog &cat, const std::vector<Pt> &extra_mask, double saturation, int aperture, MFMap &M);
// Forced measurement at (x, y): fit of the frame PSF with free position (offset ≤ max(1 px, FWHM/2)); snr = u at
// the predicted position. Fills a Detection (flux on the aperture scale via the aperture correction).
// fixed_fallback: 1 = if the free-position fit fails, fit flux only with position fixed at (x, y); 2 = always
// fixed position (d.converged = false in these cases).
bool forced_measure(const Frame &f, const MFMap &M, double x, double y, double gain, Detection &d, double &snr, int fixed_fallback = 0);
// New tracks via forced confirmation of detection pairs (pair_tracks) and stacking along velocities
// (stack_search). Measurements are appended to Frame::det. Settings: stack_min_rate/stack_max_rate (″/min, default
// MinRate/1.0), stack_far (FA per sequence, 0.01), stack_frame_snr (3), stack_clip (6), stack_budget (4e9).
std::vector<Track> stack_search(Reduction &R, const Catalog &cat, const Settings &s);

} // namespace internal
} // namespace astro
