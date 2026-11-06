// Piazzi: C++ core for asteroid astrometry and photometry.
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace astro {

// 16-bit pixel image, row by row (pix[y*w + x]); outside the frame at() returns 0.
struct Image {
    int w = 0, h = 0;
    std::vector<uint16_t> pix;
    double jd_mid = 0;  // mid-exposure (UTC, JD)
    double exptime = 0; // seconds
    double gain = 0;    // e-/ADU from the FITS header (GAIN or EGAIN); 0 = unknown
    std::string path;
    int at(int x, int y) const { return (x < 0 || y < 0 || x >= w || y >= h) ? 0 : pix[size_t(y) * w + x]; }
};

// Reads SBIG ST-x (compressed or not) and integer FITS 8/16/32 bit; flip as flip_horizontal/flip_vertical.
Image load_image(const std::string &path, bool flip_h, bool flip_v, bool time_is_start = true);

// A detection (detect.cpp).
struct Detection {
    double x = 0, y = 0; // centroid (Gaussian fit or moments), 0-based pixels
    double flux = 0;     // aperture photometry (ADU, annulus background subtracted)
    double fwhm = 0;     // FWHM in pixels
    double snr = 0;      // peak SNR (amplitude / sky σ)
    double snr_tot = 0;  // aperture photometry SNR (CCD equation)
    double fit = 0;      // rms of the Gaussian fit residual / amplitude
    double peak = 0, bkg = 0, sigma = 0, det_snr = 0;
    int ix = 0, iy = 0;
    bool converged = false;
    // elongated source (fast-object trail): segment length (px, 0 = point source) and direction (rad, image x axis
    // → y axis); the position is the segment center (mid-exposure)
    double trail_len = 0, trail_ang = 0;
};

struct DetectParams {
    int aperture = 4;
    double limit = 4.0;
    double saturation = 60000;
    double min_fwhm = 0.7;
    double fit_rms = 0.20;
    bool bkground = true;
    int binning = 1;
    int margins[4] = {0, 0, 0, 0}; // left, right, top, bottom (px ignored)
    double gain = 1.0;             // e-/ADU (CCD equation); 1 if unknown
    bool matched = true;           // 2nd pass with PSF-matched filter (Gaussian of the 1st-pass median FWHM)
    double beta = 0;               // psf_beta: > 1 fits the centroid with a Moffat profile of exponent β; 0 = Gaussian
    bool trails = true;            // trail_detect: trail search with an elongated matched filter (segment ⊗ PSF)
    double trail_max = 30;         // trail_max_len: largest trail kernel (px; the fitted length may exceed it)
    bool trail_prescreen = true;   // trail_prescreen: aperture SNR in the kernel capsule before fitting
    int trail_iter = 60;           // trail_fit_iter: LM iterations in the search fits (60 up to v7)
    double trail_gate = 1.0;       // trail_gate: only pixels with S ≥ this enter the kernel search
};

// Trail search diagnostics (trail_debug): time of each phase (s) and counts.
struct DetectStats {
    double t_total = 0, t_kernels = 0, t_mask = 0, t_fits = 0, t_replace = 0;
    int nkern = 0, ncomp = 0, ncand = 0, nfit = 0, ntrail = 0, nskip = 0; // nskip: discarded by prescreening
};

// Maps kept by detect_sources for the deferred trail search (detect_trails): per-pixel sky background and σ, image
// filtered by the PSF (F), PSF σ (ks) and effective filter noise gain (ng). ok = maps present (matched filter ran).
struct DetectScratch {
    std::vector<float> bkg, rms, F;
    double ks = 0, ng = 1;
    bool ok = false;
};

// Source detection (mesh background, k·σ threshold, connected components, centroid, PSF and aperture).
// Output sorted by decreasing flux. stats (optional) receives the trail search diagnostics.
// scratch (optional): no trail search here; keeps the maps for detect_trails (after plate solving and classification).
std::vector<Detection> detect_sources(const Image &im, const DetectParams &p, DetectStats *stats = nullptr,
                                      DetectScratch *scratch = nullptr);
// Deferred trail search: same search as detect_sources, but ignoring components centered within ≤ r of a position in
// `mask` (x, y, r triplets in px: static sources — references, catalog stars, detections with a counterpart in another
// frame). Does not reorder D: returns the trails (to append at the end) and marks in frag the point fragments of each.
std::vector<Detection> detect_trails(const Image &im, const DetectParams &p, const DetectScratch &scratch, const std::vector<Detection> &D,
                                     const std::vector<double> &mask, std::vector<char> &frag, DetectStats *stats = nullptr);

// INI file settings (key=value, sections ignored).
struct Settings {
    std::map<std::string, std::string> kv;
    static Settings defaults();
    void load(const std::string &path);
    double num(const std::string &k) const;
    std::string str(const std::string &k) const;
    void set(const std::string &k, const std::string &v) { kv[k] = v; }
};

} // namespace astro
