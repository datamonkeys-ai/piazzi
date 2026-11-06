// Astrometric and photometric reduction (plate_solve.cpp, reduce.cpp, linking.cpp).
#pragma once
#include "piazzi/core.hpp"
#include <cmath>
#include <string>
#include <vector>

namespace astro {

constexpr double PI = 3.14159265358979323846;

// Reference catalog (at image epoch), sorted by magnitude.
struct Catalog {
    std::vector<double> ra, dec, mag;
    // position σ at image epoch (″, per axis), when the catalog provides position and proper-motion errors;
    // empty (or of a different size) = unknown.
    std::vector<double> err;
    // Gaia BP − RP color (mag) when the catalog provides BPmag and RPmag; empty (or NaN) = unknown
    std::vector<double> color;
    size_t size() const { return ra.size(); }
    bool has_color() const { return color.size() == ra.size() && !ra.empty(); }
    // copy of star i of `o` to the end of this catalog (keeps err/color when `o` has them)
    void push_from(const Catalog &o, size_t i) {
        ra.push_back(o.ra[i]);
        dec.push_back(o.dec[i]);
        mag.push_back(o.mag[i]);
        if (o.err.size() == o.size())
            err.push_back(o.err[i]);
        if (o.has_color())
            color.push_back(o.color[i]);
    }
};
// Solar BP − RP color in the Gaia DR3 system (Casagrande & VandenBerg 2018, MNRAS 479, L102: 0.82). Asteroids have
// near-solar colors (nearly flat reflectance in the visible; Ivezić et al. 2001, AJ 122, 2749), so the magnitude of
// a Solar System object is referred to this color.
constexpr double SOLAR_BP_RP = 0.82;
// Reads a VizieR asu-tsv response (Gaia DR3: RAJ2000 DEJ2000 pmRA pmDE ... RPmag), applies proper motion
// from J2000 to epoch_jyear and discards RPmag > maglim. If error columns exist (e_RA_ICRS/e_DE_ICRS/e_pmRA/e_pmDE or
// ra_error/dec_error/pmra_error/pmdec_error, in mas and mas/yr at the Gaia DR3 epoch 2016.0), propagates position σ to
// the epoch (Catalog::err); if RUWE exists, discards stars with RUWE > ruwe_max (1.4; 0 = no cut; Lindegren et
// al. 2021, A&A 649, A2).
Catalog load_vizier_tsv(const std::string &path, double epoch_jyear, double maglim, double ruwe_max = 1.4);
// Builds the URL of a VizieR box query (Gaia DR3, I/355/gaiadr3) with positions, proper motions, magnitudes,
// parallax, position and proper-motion errors (e_RA_ICRS…) and RUWE: with them each reference's catalog σ
// is the real one (≈ 1 mas for bright stars 10–20 years from the Gaia epoch, 10–15 mas at G ≈ 20) instead of the
// constant plate_cat_err, and the excess scatter of bright references goes into the plate σ_sys, where it belongs.
std::string gaia_box_url(double ra, double dec, double w_arcmin, double maglim, const std::string &host = "vizier.cds.unistra.fr");

// Box query URL for catalog `catalog` (4 UCAC4, 5 CMC-15, 6 PPMXL, 7 Gaia DR3, 8 ATLAS-RefCat2); "" if unknown.
std::string catalog_box_url(int catalog, double ra, double dec, double w_arcmin, double maglim,
                            const std::string &host = "vizier.cds.unistra.fr");

void radec_to_std(double ra, double dec, double ra0, double dec0, double &xi, double &eta); // arcsec
void std_to_radec(double xi, double eta, double ra0, double dec0, double &ra, double &dec); // degrees

struct Plate {
    double ra0 = 0, dec0 = 0;
    int order = 1;
    std::vector<double> cx, cy;
    double rms_ra = 0, rms_dec = 0, dra_mean = 0, dde_mean = 0;
    // out-of-sample rms (k-fold cross-validation on references, ″); 0 if not computed
    double rms_cv_ra = 0, rms_cv_dec = 0;
    int nref = 0, nmatch = 0, nphot = 0;
    std::vector<int> ref_idx, cat_idx;
    std::vector<char> ref_ok;
    // photometry: m_cat = zp + zp_color·(c − SOLAR_BP_RP) − 2.5·log10(flux), with c = star BP − RP; zp is the
    // zero point for a solar-colored object and zp_color the frame color term (frame filter × catalog band)
    double zp = 0, zp_rms = 0, dmag_mean = 0, zp_color = 0;
    // predicted instrumental magnitude in this frame (on the zp scale) of a catalog star of magnitude m and color c
    double star_mag(double m, double c) const { return std::isfinite(c) ? m - zp_color * (c - SOLAR_BP_RP) : m; }
    // error model (″): per-axis systematic component (unmodeled turbulence/distortion, estimated from the
    // references after removing centroid noise) and coefficient covariance (sandwich) in the normalized basis
    // u = (x − nx0)/ns, v = (y − ny0)/ns, for the variance of the plate prediction at a point
    double sys_ra = 0, sys_dec = 0, nx0 = 0, ny0 = 0, ns = 1;
    std::vector<double> cov_ra, cov_dec;
    // per-axis astrometric color term (″/mag): reference ξ, η get + col·(BP − RP − SOLAR_BP_RP) (differential
    // chromatic refraction and optical lateral color; Pier et al. 2003, AJ 125, 1559, §3); the plate polynomial is
    // referred to the solar color, that of the target. 0 when the catalog has no color or the term was not fitted.
    double col_ra = 0, col_dec = 0;
    int ncol = 0; // references with color used in the term (0 = no term)
    // plate prediction σ (″) at (x, y) per axis; 0 if the covariance was not computed
    void pred_sigma(double x, double y, double &sra, double &sdec) const;
    std::vector<double> terms(double x, double y) const;
    void pix_to_radec(double x, double y, double &ra, double &dec) const;
    double scale_x() const; // "/px along x
    double scale_y() const;
};

struct Measurement {
    int frame;
    double x, y, ra, dec, mag, flux, snr, snr_tot, fwhm_as, jd;
    std::string note;
    std::string remark; // ADES remarks (e.g. forced measurement of a track found by stacking)
    // centroid covariance (px²; NAN = no model, isotropic σ from the detection) and trail (length in px and
    // unit direction in pixels; 0 = point source)
    double cxx = NAN, cxy = NAN, cyy = NAN, trail_len_px = 0, trail_dx = 1, trail_dy = 0;
}; // note: MPC note (1 character) / ADES notes
struct Track {
    std::vector<int> idx;
    std::vector<Measurement> meas;
    double rate = 0, pa = 0;
    std::string desig;
};

struct Frame {
    Image im;
    std::vector<Detection> det;
    Plate plate;
    bool solved = false;
    double seeing = 0, offx = 0, offy = 0;
    std::vector<int> typ; // per detection: 0 unmatched, 1 linked star, 2 reference (empty in stacked frames)
    bool stacked = false; // Track & Stack frame
};

struct Reduction {
    std::vector<Frame> frames;
    std::vector<Track> movers;
    int npattern = 0;
    double gain = 1.0; // effective e-/ADU (Settings Gain; else header GAIN/EGAIN; else 1)
    std::vector<std::string> log;
};

Reduction reduce(std::vector<Image> images, const Catalog &cat, double ra0, double dec0, const Settings &s, bool moving);

// ADES PSV and MPC 80-column reports (public MPC specifications), CRLF.
std::string ades_report(const Reduction &r, const std::vector<Track> &tracks, const Settings &s);
std::string phot_report(const Reduction &r, const std::vector<Track> &tracks, const Settings &s);
std::string mpc_report(const Reduction &r, const std::vector<Track> &tracks, const Settings &s);
// per-axis σ (″) of a measurement (the ADES rmsRA/rmsDec); cra/cde (optional) receive only the centroid term.
void measurement_sigma(const Reduction &R, const Track &T, const Measurement &m, const Settings &s, double &era, double &ede,
                       double *cra = nullptr, double *cde = nullptr);

// trkSub/temporary designation of an unidentified track (idx = position among the report tracks)
std::string temp_id(const Track &T, size_t idx);

// Utilities
std::string jd_to_isot(double jd);                  // "YYYY-MM-DDTHH:MM:SS" (truncated, like astropy isot[:19])
std::string fmt_ra(double ra, bool extra = false);  // "HH MM SS.ss"
std::string fmt_dec(double de, bool extra = false); // "+DD MM SS.s"
std::string packed_from(const std::string &desig);

// Manual measurement (like detection, without the threshold/FWHM filters): finds the brightest pixel within the
// aperture radius around (x, y) and does background + centroid + PSF + aperture. false if there is no signal.
bool measure_at(const Image &im, double x, double y, const DetectParams &p, Detection &out);
DetectParams detect_params(const Settings &s);
// (α, δ) → pixel by numerical (Newton) inversion of the plate; (x, y) is the initial guess.
bool plate_radec_to_pix(const Plate &p, double ra, double dec, double &x, double &y);
// Jacobian ∂(ξ, η)/∂(x, y) of the plate at (x, y), ″/px, J = {∂ξ/∂x, ∂ξ/∂y, ∂η/∂x, ∂η/∂y}, by 1 px differences; ξ, η in
// the tangent plane at (ra_t, dec_t) — NAN: at the point itself (ξ = east = RA·cos δ, η = north).
void plate_jacobian(const Plate &p, double x, double y, double J[4], double ra_t = NAN, double dec_t = NAN);

} // namespace astro
