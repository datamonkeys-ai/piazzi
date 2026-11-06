// C API of the Piazzi core, used by the iOS app (Swift) and by any other language.
//
// Concurrency: an AstroSession is NOT thread-safe. All calls on the same session must be serialized
// by the caller (one thread at a time, or an external lock). Distinct sessions are independent.
//
// Null pointer: every function accepts s == NULL and returns the failure value (0, -1, "" or "[]"), with no effect.
//
// Returned strings (const char *): owned by the session and valid until the next call that rewrites the SAME buffer
// (or until astro_free). Copy them if you need to keep them. Shared buffers:
//   - astro_get                                     -> own buffer (getv)
//   - astro_catalog_url, astro_catalog_url_for      -> same buffer (buf)
//   - astro_result_json                             -> json buffer, also rewritten by astro_reduce,
//     astro_set_mover_desig, astro_add_manual, astro_set_note, astro_set_excluded and astro_stack
//   - astro_report                                  -> own buffer (report)
//   - astro_image_header                            -> own buffer (hdr)
//   - astro_identify_offline, astro_image_catalog   -> same buffer (buf3)
//   - astro_known_offline                           -> own buffer (buf2)
//   (with s == NULL, astro_identify_offline/astro_known_offline use a static global buffer, also not thread-safe)
//
// Return codes: unless stated otherwise, int = 1 success / 0 failure (e.g. missing frame or unsolved plate);
// astro_reduce, astro_stack and astro_add_manual return -1 on error; counts (astro_add_image, astro_image_count,
// astro_frame_count, astro_load_orbits) return the number of items.
//
// Output buffer sizes: rgba (astro_render, astro_render_levels) must hold w*h*4 bytes, with w/h from
// astro_frame_size; out of astro_measure must hold 12 doubles; counts of astro_histogram, `bins` doubles, and minmax2,
// 2 doubles (minimum and maximum).
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
typedef struct AstroSession AstroSession;
// Piazzi version ("0.2.0"; PIAZZI_VERSION in piazzi/version.h); static string
const char *astro_version(void);
// schema of the settings exposed to the apps (static JSON): {"groups":[{id,en,pt,collapsed}],"settings":[{key, group,
// type (number/int/bool/enum/text), en, pt, hint_en, hint_pt, default, min, max, options:[{value,en,pt}]}]}
const char *astro_settings_schema(void);
AstroSession *astro_new(void);
void astro_free(AstroSession *s);
void astro_load_cfg(AstroSession *s, const char *path);
void astro_set(AstroSession *s, const char *key, const char *value);
const char *astro_get(AstroSession *s, const char *key);
int astro_add_image(AstroSession *s, const char *path); // returns the number of images
void astro_clear_images(AstroSession *s);
const char *astro_catalog_url(AstroSession *s, const char *ra, const char *dec);
// Reduces the images: returns the number of moving objects, or -1 on error. On error, astro_result_json holds
// {"error": English message, "code": stable code, "detail": optional}; apps translate by code:
//   no_images, image_unreadable (detail = path), catalog_empty, internal (detail = exception text)
int astro_reduce(AstroSession *s, const char *catalog_path, const char *ra, const char *dec, const char *object);
const char *astro_result_json(AstroSession *s);
const char *astro_report(AstroSession *s, int which); // 0 = ADES, 1 = PhotReport, 2 = MPC
int astro_frame_count(AstroSession *s);
int astro_frame_size(AstroSession *s, int frame, int *w, int *h);
void astro_render(AstroSession *s, int frame, unsigned char *rgba); // w*h*4 bytes, automatic contrast
// RA/Dec (degrees) -> pixel in image `frame` via the fitted plate; returns 0 if there is no plate.
int astro_sky_to_pix(AstroSession *s, int frame, double ra, double dec, double *x, double *y);
// Field center (degrees) of image `frame`, via the plate.
int astro_frame_center(AstroSession *s, int frame, double *ra, double *dec);
int astro_pix_to_sky(AstroSession *s, int frame, double x, double y, double *ra, double *dec);
void astro_set_mover_desig(AstroSession *s, int index, const char *desig);

// ===================== full parity (see session_internal.hpp) =====================
// --- calibration: dark and flat applied to the images on the next astro_reduce (empty path = remove)
int astro_set_dark(AstroSession *s, const char *path); // 1 = ok
int astro_set_flat(AstroSession *s, const char *path);
// --- time / header (index = astro_add_image order)
int astro_image_count(AstroSession *s);
const char *astro_image_header(AstroSession *s, int index);     // header as text (FITS: cards; SBIG: key = value)
double astro_image_jd(AstroSession *s, int index);              // instant used (including time_offset_s etc.)
void astro_set_image_jd(AstroSession *s, int index, double jd); // overrides the instant (0 = revert to header)
// --- manual measurement (centroid + PSF like detection) on reduced frame `frame`
//     out[12] = ra, dec, mag, snr_peak, fwhm_arcsec, x, y, fit_rms, flux, snr_total, dra_arcsec, ddec_arcsec (plate)
int astro_measure(AstroSession *s, int frame, double x, double y, double *out);
//     creates/extends a hand-measured object (track = -1: new); returns the object index (included in reports)
int astro_add_manual(AstroSession *s, int track, int frame, double x, double y);
void astro_set_note(AstroSession *s, int track, int frame, const char *note); // 1-character MPC note (e.g. "I")
void astro_set_excluded(AstroSession *s, int track, int frame, int excluded); // excludes the measurement from reports
// --- Track & Stack: sum/mean/median shifting by rate (″/min) and PA (°); mode 0 sum, 1 mean, 2 median
//     appends a "stacked" frame to the result and returns its index (or -1); detections/measurements work on it
int astro_stack(AstroSession *s, double rate_arcsec_min, double pa_deg, int mode);
// --- catalogs: 4 UCAC4, 5 CMC-15, 6 PPMXL, 7 Gaia DR3, 8 ATLAS-RefCat2 (same numbering as the .cfg Catalog key)
const char *astro_catalog_url_for(AstroSession *s, int catalog, const char *ra, const char *dec);
// --- display: levels as percentiles (e.g. 50 and 99.5), mode 0 linear 1 sqrt 2 log, invert 0/1
void astro_render_levels(AstroSession *s, int frame, double lo_pct, double hi_pct, int mode, int invert, unsigned char *rgba);
int astro_histogram(AstroSession *s, int frame, int bins, double *counts, double *minmax2);
// --- save 16-bit FITS with WCS (TAN, from the fitted plate)
int astro_save_fits(AstroSession *s, int frame, const char *path);
// --- offline orbits: MPCORB.DAT (or a subset) and COMET.DAT; returns how many orbits were read
int astro_load_orbits(AstroSession *s, const char *mpcorb_path, const char *comet_path);
//     objects near (ra, dec) at instant jd, topocentric from [ObservingSite] (Longitude/Latitude/Height); JSON
const char *astro_identify_offline(AstroSession *s, double jd, double ra, double dec, double radius_deg);
//     known objects in the field of the reduced frame, with x/y; JSON
const char *astro_known_offline(AstroSession *s, int frame, double vmax);
// --- frame detection list: JSON with x, y, ra, dec, mag, snr, fwhm, fit and type: R used reference, r rejected
//     reference, S star, T point fragment of a trail, L attached trail (trail_len > 0), M mover, ? unmatched
const char *astro_image_catalog(AstroSession *s, int frame);

#ifdef __cplusplus
}
#endif
