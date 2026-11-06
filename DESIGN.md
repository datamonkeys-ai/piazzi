# Piazzi — design and algorithms

Piazzi's core (`cpp/src`) is an independent implementation built from published literature. This document lists
what each stage does and where the method comes from. Source files cite the same references in their headers.

## 1. Image input (`image_io.cpp`)
FITS (BITPIX 8/16/32/−32/−64, BZERO/BSCALE) and SBIG ST-x (raw or compressed). Dimensions and data sizes are
validated before allocation. Times are converted to mid-exposure (start + exposure/2 when the header gives the start).
Dark and flat frames are applied when configured. Gain: the `gain` setting when given, otherwise the FITS `GAIN`
(or `EGAIN`) keyword (median over the images; ZTF: 6.2 e⁻/ADU), otherwise 1 e⁻/ADU. It enters the Poisson term of
the CCD equation (aperture SNR) and the source-electron count of the centroid error (ADES σ).

## 2. Source detection (`detect.cpp`)
- Background mesh (32–64 px cells): sigma-clipped median per cell, 3×3 median filter on the mesh, bilinear
  interpolation; noise from the MAD. — *Bertin & Arnouts 1996, A&AS 117, 393 (SExtractor)*
- Detection: k·σ threshold on a 3×3 binomial-smoothed image, 8-connected components with a minimum area;
  blended sources split at 5×5 local maxima. The threshold is in units of the filtered noise, σ·√(Σw²)/Σw.
- Optional matched filter (`matched_filter=1`, default 0): a second pass convolves with a separable Gaussian of the
  median FWHM of the well-measured first-pass sources (white noise ⇒ maximal point-source SNR). On synthetic data
  (FWHM 5 px) it finds a source at optimal SNR ≈ 5.5 that the 3×3 kernel misses, and on the ZTF benchmark it
  raised the source count by ~50 % in bad seeing, but it found no additional target and added one false track, so
  it is off by default. — *Irwin 1985, MNRAS 214, 575; Bertin & Arnouts 1996 §6*
- Trail detection (fast movers, `trail_detect`, default 1; `trail_max_len`, longest kernel, default 30 px; runs after
  the matched-filter pass): a point-source kernel loses
  SNR on a trail as ≈ √(FWHM/L), and segmentation splits a 2–5 FWHM trail into 2–3 point detections whose centroids
  sit anywhere along it. After the matched-filter pass, S = F/σ_F (point-filter SNR map) is averaged along segments
  of length L = 2·FWHM·2^j (≤ `trail_max_len`) at angles spaced so the segment end moves ≤ σ_PSF — the matched filter to
  a trail (segment ⊗ PSF), with the noise of each kernel the larger of the analytic value (Gaussian autocorrelation
  of S) and an empirical MAD over source-free segments (ZTF cutouts have correlated noise: the analytic value alone
  put 8 % of the pixels above threshold). Samples are clipped at 1.5× the threshold so a bright star crossed by a
  segment does not make a "trail". A pixel enters the trail mask if the clipped response is ≥ k + 1 (several hundred
  kernels per pixel), the unclipped response beats the point filter everywhere on the segment (≥ 1.1·max S − 3σ),
  and the samples are uniformly high (median ≥ ½ mean); components with fewer than 6 pixels or ¼ of the kernel length
  are noise peaks, and components centred on a star-shaped point detection or within max(3, 2·log₁₀ SNR*)·FWHM + L/2
  of a point source with SNR* > 3× the component's response are skipped (real, non-Gaussian PSF wings give arc-like
  responses around bright stars on ZTF). The 40 most elongated components (unclipped elongated/point response ratio)
  are fitted — deep SDSS frames have hundreds of elongated galaxies — with the trail model (below), inside a capsule of
  radius 3σ + 3 px around the kernel segment, extended by the same amount at both ends (pixel integral approximated by
  one sample with σ² + 1/12 when σ ≥ 1.2 px), with
  position, flux, length (weak prior) and angle free; it is kept if L ≥ 1.5·FWHM, if the profile along it has no
  dip (two-star pairs), and if it beats the competing model of the point detections lying on it (two PSFs when
  there is only one) by an information criterion χ²/σ² + 4·n_par (penalty between Akaike's 2 and Schwarz's ln n;
  chosen on the synthetic tests), and its capsule-aperture SNR is ≥ k + 0.5 (fainter elongated sources keep their point
  detections). The trail replaces its point fragments; position = segment centre = mid-exposure,
  flux from a capsule aperture (radius = `aperture_radius` around the segment), `Detection::trail_len/trail_ang`
  carry the shape. — *Dawson et al. 2016, AJ 152, 112; Nir, Zackay & Ofek 2018, AJ 156, 229; Vereš et al. 2012;
  Akaike 1974, IEEE TAC 19, 716; Schwarz 1978, Ann. Stat. 6, 461*
- Deferred trail search (v8, `trail_late`, default 1): `detect_sources` keeps its background/σ/filtered maps
  (`DetectScratch`) and the search runs after the plate solution (`late_trails` in `reduce.cpp`), skipping components
  centred within max(3″, seeing) of a static source — plate references, projected catalog stars and detections with a
  counterpart in another solved frame (galaxies and blends are static; a trail is in one frame only). Point fragments
  of each trail are marked `typ = 3` (excluded from linking and from the neighbour list of the PSF refit) and the trail
  detections are appended at the end of `Frame::det`, so reference indices do not move. A cheap capsule-aperture SNR
  pre-screen (`trail_prescreen`, the same k + 0.5 cut applied after the fit) runs before each Levenberg–Marquardt fit,
  and only pixels with S ≥ `trail_gate` (1.0) enter the kernel loop. `trail_debug=1` logs the time of each phase.
- Centroid by intensity-weighted moments, refined by a circular 2D Gaussian fit with Levenberg–Marquardt.
  — *Stetson 1987, PASP 99, 191; Auer & van Altena 1978, AJ 83, 531*
- Aperture photometry with area weighting; sky from a sigma-clipped annulus; SNR from the CCD equation.
  — *Howell 2006, Handbook of CCD Astronomy*
- Trail fitting for moving objects (`psf_fit.cpp`, after linking): the model is the frame PSF — an axis-aligned
  elliptical Gaussian whose σx, σy are the medians fitted on isolated reference stars — convolved with a segment of
  length v·t_exp along the motion direction (rate/PA mapped to pixels through the plate Jacobian), integrated over the
  pixel (3×3 midpoint rule). The segment integral has a closed form with erf. Levenberg–Marquardt fits
  (x0, y0, flux, background); for resolved trails also length and angle, the length with a Gaussian prior
  L = v·t_exp (σ = 0.5 px ⊕ 15 %, a pseudo-residual scaled by the pixel noise); the segment centre is mid-exposure.
  — *Vereš et al. 2012, PASP 124, 1197*
- Neighbour separation: catalog stars (positions fixed by the plate) and stationary detections within
  max(2·FWHM + L/2, fit window) are fitted simultaneously with the target (fluxes free, frame PSF). The target
  magnitude comes from its fitted flux times the frame aperture correction (median aperture/PSF flux ratio of
  the reference stars), so it stays on the zero point of the aperture photometry. — *Stetson 1987, PASP 99, 191*
- Settings: `psf_phot` (default 1: magnitude from the joint fit; 0: aperture), `trail_deblend` (default 1),
  `trail_fit` (positions from the fit: −1 = automatic, default — trail positions when the predicted trail
  L = v·t_exp ≥ 1.5·FWHM, detection centroid otherwise; 0 = never; 1 = always trail; 2 = star PSF). Below 1.5 FWHM the
  trail fit did not improve astrometry (trails ≈ 1 FWHM, SNR 7–15), above it a Gaussian centroid of an
  elongated source is noisier and, when segmentation fragments the trail, biased along it.

## 3. Plate solution (`plate_solve.cpp`, `reduce.cpp`)
- Gnomonic projection to standard coordinates. — *Calabretta & Greisen 2002, A&A 395, 1077*
- Triangle matching in invariant ratio space between the brightest detections and catalog stars, with the
  nominal scale from focal length and pixel size, both parities, RANSAC consensus and affine refinement.
  — *Groth 1986, AJ 91, 1244; Valdes et al. 1995, PASP 107, 1119; Fischler & Bolles 1981*
- Plate constants: least-squares polynomial of order 1–3 (Householder QR), iterative nearest-neighbour matching
  with a shrinking radius and 3σ rejection until convergence. — *Green, Spherical Astronomy; Golub & Van Loan*
- Weights 1/σ² per reference (`plate_weights`, default 1), σ² = σ_centroid² + σ_cat² + floor², σ_centroid =
  FWHM/(2.355·SNR), σ_cat from the catalog position and proper-motion errors propagated from the Gaia DR3 epoch
  (2016.0) to the image epoch when the TSV has `e_RA_ICRS/e_DE_ICRS/e_pmRA/e_pmDE` (or `ra_error`…), otherwise
  `plate_cat_err` (0.01″); floor `plate_floor` (0.10″). Rejection uses residuals scaled by each σ. References with
  RUWE > `plate_ruwe` (1.4; 0 = keep all) are dropped when the column exists. The Gaia box query printed by `piazzi` (`gaia_box_url`) asks for
  these error columns and RUWE, so σ_cat is the real one per star: ≈ 1 mas for G < 17 ten to twenty years from
  2016.0, 10–15 mas at G ≈ 20. With the constant 0.01″ the bright SDSS references (SNR > 150, own noise 0.002″)
  scattered 0.010–0.018″ about the plate and that excess was booked as "catalogue error", leaving σ_sys at
  0.003–0.008″; with the real σ_cat it goes to σ_sys, where it belongs (benchmark e20 below).
  — *Lindegren et al. 2021, A&A 649, A2 (RUWE)*
- Astrometric colour term (`plate_color`, default 1): the design matrix gets one extra column per axis,
  (BP − RP) − (BP − RP)☉, so ξ_cat = P(x, y) + k_ξ·(c − c☉) and likewise η; the target (solar colour) is evaluated with the
  polynomial alone. Needs ≥ 10 references with colour and a colour spread ≥ 0.2 mag, else k = 0. The term absorbs
  differential chromatic refraction (direction: parallactic; size ∝ tan z × Δn between the object's and each
  reference's effective wavelength) and the lateral colour of the optics without a model of either — the two
  coefficients per frame are fitted on the references like the plate constants, and enter the degrees of freedom of
  σ_sys and the covariance. On the 32 SDSS sequences (15 159 references) the residual-versus-colour slopes before the
  term were, in the scan direction (RA, image y) / cross-scan (Dec, x): r −0.0054/+0.0037, i −0.0030/+0.0012,
  z −0.0017/−0.0007, g −0.0202/−0.0122 ″ per mag of BP − RP (±0.0005), i.e. a solar-coloured object measured against
  references of mean colour ≈ 1.4 was displaced by +0.005″ (r) to +0.015″ (g); with the term every slope is
  < 0.0013 ″/mag. The expected atmospheric DCR (Edlén dispersion, Δλ_eff ≈ 10 nm in g and 6 nm in r between the Sun
  and a BP − RP = 2 star) is 0.06·tan z ″ in g and 0.016·tan z in r, i.e. 0.01–0.04″ at the SDSS tan z = 0.4–0.9,
  the same order as measured; but most SDSS fields were within 30° of the meridian (parallactic angle |q| < 30°, DCR
  along Dec) while the larger slope is along the scan, so the term is partly instrumental — SDSS's own astrometric
  pipeline fits colour terms in both axes for the same reason. The unit test (`test_linking` §8) injects
  +0.03/−0.02 ″/mag and recovers the coefficients to 0.001 ″/mag and the solar-coloured target to < 0.002″.
  — *Pier et al. 2003, AJ 125, 1559, §3; Filippenko 1982, PASP 94, 715 (DCR); Edlén 1966, Metrologia 2, 71*
- `plate_debug=1` prints every reference (position, magnitude, colour, SNR, FWHM, residual per axis, σ_c, σ_cat) and every
  detection with its plate position in the log, for the systematics studies below.
- Order selection: `plate_order = 0` (the default when no configuration fixes it) picks order 1–3 by 5-fold cross
  validation of the clipped reference set (a higher order must lower the out-of-sample rms by ≥ 3 %). The
  out-of-sample (cross-validated) rms is stored in `Plate::rms_cv_*` for every order and is the plate term of the
  ADES uncertainties. Reported plate rms values are weighted rms. — *Stone 1974, JRSS B 36, 111; Hastie et al.
  2009, Elements of Statistical Learning §7.10*
- Error model of the plate at a point (used for the ADES uncertainties): a systematic term σ_sys per axis
  (atmospheric image motion, unmodelled distortion) chosen so that Σ r²/(σ_c² + σ_cat² + σ_sys²) = m − n_terms over
  the reference residuals (reduced χ² = 1 once their own centroid noise σ_c and catalogue error σ_cat are removed;
  n_terms counts the colour term), and the covariance of the
  plate coefficients by the sandwich form (AᵀWA)⁻¹AᵀW·V·WA(AᵀWA)⁻¹ in a normalised (u, v) basis, giving the
  prediction variance tᵀC t at the target. The target is not a reference, so what reaches it is the out-of-sample
  error: σ_sys is the larger of the in-sample χ² estimate and the cross-validated one, σ_sys,cv² = rms_cv² −
  ⟨σ_c² + σ_cat²⟩_w (with ~70–90 SDSS references and order 3 the in-sample χ² often returned σ_sys = 0 on one
  axis, leaving a 0.004″ ADES σ for a target whose repeat scatter is 0.015″). — *Huber 1967; White 1980,
  Econometrica 48, 817; Stone 1974*
- Frames whose field is not identified against the catalog are aligned (triangles between detection lists) to an
  already solved frame; a second pass makes this independent of frame order. Before, a failure in the first frame
  had nothing to align to: in 4 SDSS sequences the r frame (first in time) failed because its brightest stars are
  saturated and skipped, so the 60 brightest unsaturated detections overlapped poorly with the brightest catalog
  stars (4 triangle matches, < 6 needed) while i/z/g matched.
- Photometric zero point with a colour term: m_cat + 2.5·log F = zp + k·(c − c☉) per reference, c = Gaia BP − RP
  (`Catalog::color`, read when the TSV has `BPmag` and `RPmag`), least squares with 3σ (MAD) clipping; zp is the
  zero point for a solar-coloured object, (BP − RP)☉ = 0.82, and `Plate::star_mag` gives a catalogue star's expected
  magnitude on the frame (used by the fixed-star test and the stack mask). Needs ≥ 10 references with colour and a
  colour spread ≥ 0.2 mag; otherwise k = 0 and zp = sigma-clipped median, as before. Without it, the zero point of
  the (redder than solar) field stars made a solar-coloured asteroid ≈ 0.5 mag brighter in g and ≈ 0.3 mag fainter
  in z than in r: on SDSS r, i, z, g sequences the target series spanned 1.0–1.5 mag for every detected target and
  1.6 mag for (1030) Vitja (14.29 / 14.82 / 15.17 / 13.55), above the linker's amplitude limit; with the colour term
  it is 14.9 / 14.8 / 14.9 / 14.8. Asteroids are close to solar colour (*Ivezić et al. 2001, AJ 122, 2749*), the
  transformation equation is the standard one (*Harris, FitzGerald & Reed 1981, PASP 93, 507*), solar colour from
  *Casagrande & VandenBerg 2018, MNRAS 479, L102*.

## 4. Moving objects (`linking.cpp`)
- Per-detection astrometric uncertainty (per axis): σ = (FWHM/2.355)/SNR ⊕ σ_plate, with FWHM the detection's own
  (≥ frame seeing, so extended sources get larger σ) and σ_plate the per-axis plate rms. — *King 1983, PASP 95, 163;
  Mendez et al. 2014, PASP 126, 798 (Cramér–Rao bound for a Gaussian PSF)*
- Static sky: a detection is stationary if it is a plate reference, if it lies within k·σ of a catalog star of
  compatible brightness (|Δm| < 1 + 3σ_Δm), or if it has a brightness-compatible counterpart within k·√(σ_i² + σ_j²)
  in at least half of the frames that are ≥ ¼ of the arc away in time. k = `motion_sigma` (default 5): for a fixed
  source the separation is Rayleigh-distributed and P(> 5σ) = e^−12.5. There is no fixed floor in arcseconds (the
  former r_fix = max(2″, seeing/2) and v_min = 2·r_fix/arc made every object moving < 4″ over a 5-min SDSS arc
  invisible). — *Denneau et al. 2013, PASP 125, 357, §4*
- Deblending before linking: a free detection with a ≥ 1.5× brighter neighbour closer than 2.5·FWHM is re-centred by
  a joint two-PSF fit (frame PSF, common background). — *Stetson 1987, PASP 99, 191*
- Linking: anchor pairs (largest time baseline first), linear prediction in the other frames (tolerance
  max(1.5″, 3·plate rms)), least-squares ξ(t), η(t), residual and magnitude-spread (≤ 1.5 mag) checks; ≥ 3
  detections, one missing frame allowed when N ≥ 4. Competing hypotheses for an anchor: more detections wins, then
  lower rms. — *Kubica et al. 2007, Icarus 189, 151*
- Trail coherence: on every elongated detection of a candidate track, the measured length and direction must match
  L = v·t_exp and the track direction mapped to pixels through the plate Jacobian (|ΔL| ≤ max(2 px, L/2),
  |Δθ| ≤ 20°). — *Vereš et al. 2012; Dawson et al. 2016*
- Resolved motion (replaces the minimum rate): each anchor pair must be displaced by ≥ k·√(σ_a² + σ_b²), and the fitted
  track must have |v| ≥ k·σ_v, σ_v = σ/√Σ(t − t̄)² with σ = max(mean detection σ, fit rms) — a likelihood-ratio test
  of linear motion against a fixed source. `min_rate` still applies if set.
- Chance-alignment control: expected number of random tracks as good as the candidate,
  E = N_a · N_b(annulus between the resolved displacement and v_max·Δt) · Π_m ρ_m π r_m² · N_anchor-pairs
  (ρ_m = density of free detections, r_m = fit residual, ≥ 3σ of the detection; the last factor is a Bonferroni
  correction because the same track can arise from any anchor pair); rejected if E > `false_alarm` (default 0.05).
- Conditioned chance alignment (v8): N_a, N_b and ρ_m count only free detections whose brightness is compatible
  (`same_mag`) with the anchor (`link_mag_e`) — the brightness test is applied to every hypothesis, so a chance track
  needs compatible partners; when the anchor is a trail, the velocity annulus shrinks to the sector its length and
  direction admit (|ΔL| ≤ max(2 px, L/2), |Δθ| ≤ 20°, both senses; `link_trail_e`), the constraint `trail_coherent`
  imposes anyway. *Kubica et al. 2007; Denneau et al. 2013*
- `link_probe = "ra dec jd dRA/h dDec/h"` logs, for the detection nearest the prediction in every frame, its state and
  the first linking test that rejects it (`link probe`); `piazzi --inject-probe I` sets it for an injected object.
- Vacated-position test (forced photometry): at the position of every track detection, in the other frames, the flux
  in a small aperture (0.67·FWHM, energy-corrected and scaled by the zero points) must not show a source of comparable
  brightness (> 3σ and 0.4–2.5× the expected flux). In short arcs the object is still close: the fraction g of its
  light falling in the aperture (offset Gaussian, numerically integrated, offset from the fitted track) is subtracted
  using the flux measured at the predicted position, and a frame only counts when the hypotheses separate (g ≤ ½ and a
  fixed source at 0.4× would be a 3σ detection). With no such frame the decision rests on the resolved-motion test.
  With ≥ 3 detections, one ≥ 0.5 mag brighter than the track median is treated as a blend with a fixed star (object
  crossing it) and gives no evidence (`link_stat_blend`): the star is there in all frames, the object is not.

Brightness tests on a track (magnitudes on the solar-colour scale above): (a) during linking, largest minus smallest
magnitude < 1 mag + 3σ of the difference (the form of the fixed-star test; was a flat 1.5 mag, which rejected
SNR-5 detections of real objects); (b) after the final PSF refit (`drop_incoherent_movers`), a Pearson χ² of the
fluxes about their weighted mean with an intrinsic σ_int = 0.2 mag (rotation over ≲ 10 min and asteroid colour
scatter about solar; 0.7 mag if some frame has no colour term), rejected at P(χ²_{n−1}) < 1 % (quantile by Wilson &
Hilferty 1931, PNAS 17, 684). The 1 % level is stricter than the usual 10⁻³ because the measured σ ignores the
selection bias of noise tracks (*Eddington 1913, MNRAS 73, 359*): an SNR 2–7 chance alignment passes only where the
threshold is high if the noise peak there is high, so its least-sensitive frame (z in SDSS) carries its brightest
point (typical SDSS false track 21.8 / 20.9 / 20.5 / 22.1 in r / i / z / g). The test runs on the refined fluxes
only: at detection level a contaminated aperture (ZTF 2004 TV295, r 1.2 mag too bright) would reject a real object.

## 4b. Faint movers: shift-and-stack and 2-detection tracks (`stack_search.cpp`)
Runs after the linking of §4, on the solved frames (≥ 3), and adds tracks for objects too faint to be detected in
enough single frames.
- Per-frame matched-filter map (*Whidden et al. 2019, AJ 157, 119, "KBMOD"*): background mesh as in §2, sky σ per
  pixel, ψ = K ⊗ (D·m/σ²), φ = K² ⊗ (m/σ²), u = ψ/√φ (SNR of a point source at each pixel, unit variance under
  noise), K = the frame PSF (§2 fit on reference stars). The mask m removes the static sky: plate references and
  stationary detections (§4), catalog stars projected with the plate (including those brighter than `catalog_bright_mag`) and
  saturated pixels, each to the radius where its Gaussian wing falls to 0.3σ; detections of tracks already linked.
  For the stack only, also: every detection with a positional counterpart (any brightness, < max(3″, seeing)) in
  another frame; pixels whose zero-velocity stack is > 4σ while u > 1 in every frame (static sources the detector
  rejected: galaxies, blends); and connected u > 3 regions larger than twice the area of a point source of the same
  peak (halos, background residuals).
- Stack statistic for a hypothesis (x₀, v): S = Σ aₙ·min(uₙ(x₀ + v(tₙ − t_ref)), c) / √Σ aₙ², with aₙ = gₙ√φₙ the
  frame sensitivity (gₙ = 10^{0.4(ZPₙ − ZP_ref)}: same magnitude in every frame) and c = `stack_clip` (6) so that one
  frame (cosmic ray, transient) cannot carry a hypothesis. Maps are resampled on the first frame's pixel grid through
  an affine fitted to the plate mapping.
- Velocity grid (*Shao et al. 2014, ApJ 782, 1; Zhai et al. 2014, ApJ 792, 60*): square grid in (vx, vy) with step
  √2·FWHM/arc (the nearest node misplaces the end frames by ≤ FWHM/2), from max(`min_rate`, 1.5·FWHM/arc) to
  `stack_max_rate` (default 1.0″/min, main belt near opposition; faster objects trail and are bright enough for §4).
  Velocities giving the same integer shifts are evaluated once; the pixel grid is binned (sub-sampled) when
  pixels × velocities × frames exceeds `stack_budget` (4·10⁹). A running maximum over velocities keeps the best S per
  pixel; local maxima are refined in (x₀, y₀, vx, vy) by pattern search on the full-resolution maps.
- False-alarm control (*Heinze et al. 2015, AJ 150, 125*): threshold = max(k_G, k_C). k_G is the Gaussian tail for
  `stack_far` (0.01) expected false alarms per sequence given N_eff = (unmasked area/FWHM²) × (distinct velocities);
  k_C is the maximum refined S of the same search with the frame times permuted (the permutation that moves the times
  most): the hypotheses no longer correspond to straight-line motion, so the control measures the real null tail
  (noise, subtraction residuals, real objects in only some frames) on the same data.
- Acceptance: S ≥ threshold; ≥ 3 epochs (groups of frames between which the object moves ≥ 1.5 FWHM — with two
  epochs any pair of sources defines a velocity, so nothing is verified); S without the strongest epoch ≥
  `stack_loo` (3); ≥ `stack_min_frames` (3) frames individually measured with SNR ≥ `stack_frame_snr` (3); vacated
  position: at each measured position, the frames where the object is ≥ 2.5 FWHM away must show no source
  (combined u < 3) and must not be masked there.
- Per-frame measurement (forced photometry/astrometry, *Stetson 1987*): at the predicted position, a fit of the
  frame PSF with free position (shift ≤ max(1 px, FWHM/2)) and flux; the SNR is u at the predicted (forced) position;
  the flux goes to the aperture scale with the aperture correction. Only frames with SNR ≥ `stack_frame_snr` are
  written. These ADES lines carry `remarks` = "synthetic tracking" (the ADES `notes` field only takes MPC note codes).
- 2-detection tracks (`pair_tracks`, default 0): two free detections (not stationary, not in a track) at different
  epochs with a velocity in the stack range; forced u in the frames of other epochs (only those where the predicted
  position error ≤ 0.35 FWHM) combined as above. Accepted if the chance-alignment number E = N_a·N_b(velocity
  annulus)·P(S_forced | noise)·(number of frame pairs) ≤ `false_alarm`, the forced values are not far above the
  brightness predicted from the two detections (≤ 2.5× + 3σ: otherwise a static source sits there), a forced frame
  with SNR ≥ 5 accepts the PSF fit, and the vacated-position test passes. Remarks "forced phot". Off by default: on
  the development benchmark it found no target and added one unidentified track (see Accuracy).
- Settings: `stack_search` (default 1), `pair_tracks` (0), `stack_min_rate`, `stack_max_rate`, `stack_far`, `stack_clip`,
  `stack_frame_snr`, `stack_min_frames`, `stack_loo`, `stack_budget`; diagnostics `stack_debug=1` (top candidates of the
  search and of the control in the log) and `stack_probe="ra dec jd dRA/h dDec/h"` (u along a known ephemeris).

## 5. Identification (`api_orbits.cpp`, online services)
Offline: two-body propagation of MPCORB/COMET elements with an Earth ephemeris fitted to ERFA `epv00`.
Online: IMCCE SkyBoT cone search.

## 6. Reports (`report.cpp`)
ADES PSV and 80-column MPC formats as published by the Minor Planet Center. Uncertainties, per axis (the final error
model, every term from the data of the frame itself, no empirical factor):

    rmsRA, rmsDec = √(σ_c² + σ_sys² + σ_pred²)
    σ_c²    = (σ_PSF²/N)·(1 + 8π·σ_PSF²·B/N)          Cramér–Rao bound of the Gaussian centroid, σ_PSF from the
                                                        detection's own fit (a trailed object has a larger one), N source
                                                        electrons, B sky variance per pixel (e⁻²)
    σ_sys   = max(σ_χ², σ_cv) per axis                 σ_χ²: Σ r_i²/(σ_c,i² + σ_cat,i² + σ_sys²) = m − n_terms over the
                                                        references; σ_cv² = rms_cv² − ⟨σ_c² + σ_cat²⟩_w (5-fold cross
                                                        validation); σ_cat,i the Gaia error of each reference propagated
                                                        to the epoch
    σ_pred² = tᵀ C t                                    prediction variance of the plate polynomial at the target
                                                        (sandwich covariance of the weighted fit; the colour term is
                                                        zero for the solar-coloured target)

Trails (v8, `trail_cov`, default 1): σ_c is replaced by the anisotropic centroid covariance projected on the sky,
C_sky = J·C_px·Jᵀ (J = local plate Jacobian, `plate_jacobian`), with C_px from the Levenberg–Marquardt trail fit
(s²·(JᵀJ)⁻¹ of `refine_movers`) or, without that fit, the Cramér–Rao model of segment ⊗ PSF: along
(σ_PSF² + L²/12)/N·(1 + 2AB/N), across σ_PSF²/N·(1 + 2AB/N), A = 4πσ_PSF²(1 + L/(2√π σ_PSF)). The v7 isotropic
σ inflated by L²/12 also in the background term overstated trail errors ~4× (injection: |O−C|/σ = 0.17).
v10 (`trail_cov_poisson`, default 1): the trail fit is unweighted, so its covariance is the sandwich
A⁻¹(Σ V_i J_i J_iᵀ)A⁻¹ (A = JᵀJ) with V_i = σ_sky² + signal_i/gain per pixel, not s²·A⁻¹ with a uniform s² = χ²/(n−p)
that weighs sky and trail pixels alike (positions unchanged; only trails with the trail-model position use it). SDSS trail
injection (2 seeds, 539 trailed measurements): |O−C|/σ_c per axis RA/Dec 0.69/0.82 → 0.64/0.77 (all scenario
measurements 0.76/0.83 → 0.69/0.73; combined 0.795 → 0.706). The low |O−C|/σ against the *total* ADES σ in injections
(≈ 0.4 per axis for trails) is not an inflated σ: the injected truth is drawn through the plate of the first
reduction, so σ_sys ⊕ σ_pred (≈ 0.02″ on SDSS, dominant for bright trails) is absent from the injection O−C by
construction; on the sky (NEA dev vs Horizons) the same σ gives 0.6–0.75 per axis. `run.py inject` reports the
per-axis ratio against σ_c (`z_c_ra/dec_median`, trails `z_c_trail_*`). Remaining: along-track 0.75 vs cross-track 0.46
(the cross-track σ of long trails is still conservative).

(*King 1983, PASP 95, 163; Mighell 2005, MNRAS 361, 861; Huber 1967; White 1980; Stone 1974*). What is **not** in it,
deliberately: the ephemeris error of whatever the measurement is compared with, and the per-sequence offsets that such
an error produces (see Benchmark: the Horizons formal σ is 0.001–0.005″ for numbered asteroids at 1998–2008 epochs,
yet their 4-filter mean offsets scatter by ≈ 0.02″ per axis). Manual measurements without a
detection fall back to √(rms_cv² + (FWHM/(2.355·SNR))²). rmsMag = √(zp_rms² + (1.0857/SNR)²), logSNR = log10(SNR),
seeing = median FWHM of reference stars × scale. Precision: rmsRA/rmsDec with 3 decimals (2 decimals rounded the
0.01–0.03″ SDSS σ by up to ±35 %); `obsTime` with milliseconds (whole seconds shift the epoch by 0.5 s on average,
v/120 ″ along the motion); RA/Dec with 6 decimals in degrees when σ < 0.1″ or `ades_six_digit_positions` (5 decimals are 0.036″
steps, rounding rms 0.010″ per axis — comparable to the SDSS σ).
Unidentified tracks get no provID/permID; every track gets a unique `trkSub` (also the 80-column temporary
designation, columns 6–12): 'P' + 4 base-36 digits of the UTC minute of its first measurement + 2 of the track
index (7 characters).
Header `# telescope`: ADES requires `design`, `aperture` and `detector`; without `telescope_design`/`telescope_aperture_m` they are
derived from the free-text `telescope` setting (design keyword Schmidt/Refractor/Cassegrain(Ritchey)/Reflector(Newton),
else "Unknown"; aperture from "1.2-m"/"30 cm", else focal length / f-ratio "f/N"). Sexagesimal formatting clamps the
rounding carry (59.97″ → 00.0, never "-0.0", which broke the 80-column Dec field). Format check (offline):
`run.py check-ades --fields <set>` (or `--ades/--mpc <file>`) applies ADES 2022 PSV rules (version line, required
header blocks, required columns, trkSub ≤ 8 alphanumerics and unique per track, provID/permID pattern, ISO obsTime
with Z, RA/Dec decimals vs rmsRA/rmsDec — rounding step/√12 ≤ σ/3 and step ≥ σ/1000 —, rms > 0, mag/band/photCat,
astCat list, stn = mpcCode, mode CCD) and 80-column rules (exact length, date 16–32, RA 33–44, Dec 45–56 with
in-range fields, blank 57–65, mag 66–70, band 71, station 78–80 = COD); ctest `report_format` runs it on the bundled example.

## Accuracy
Test set: the bundled example (`cpp/gui/assets/example/`, see its README): three 30 s ZTF exposures (r, i, g) of
(5157) Hindemith, 2021-02-06, 356 × 356 px crops (1.01″/px), with the Gaia DR3 excerpt `gaia_dr3.tsv` and
`example.cfg`. Reference positions from JPL Horizons for station I41 at mid-exposure
(`cpp/tests/data/horizons_hindemith.txt`, query in the header); offline identification uses osculating elements from
Horizons in MPCORB/COMET layout (`cpp/tests/data/orbits_*.dat`). Current defaults (`plate_order=1` from the cfg):

| Image | O−C RA·cosδ | O−C Dec | Reference stars | Plate rms RA / Dec | Out-of-sample rms | Mag |
|---|---|---|---|---|---|---|
| 1 | −0.018″ | +0.004″ | 87 | 0.021″ / 0.022″ | 0.021″ / 0.023″ | 17.03 |
| 2 | +0.001″ | +0.007″ | 89 | 0.032″ / 0.022″ | 0.032″ / 0.023″ | 17.01 |
| 3 | −0.014″ | −0.015″ | 77 | 0.053″ / 0.032″ | 0.054″ / 0.033″ | 16.93 |

Max residual 0.021″, rms 0.016″; motion 0.353″/min at PA 276.5° (Horizons 0.353″/min, 276.6°); `plate_order=0` picks
order 1 on all three plates (same residuals). Checked by `cpp/tests/accuracy.sh` (≤ 0.06″ per ADES line) and
`test_accuracy` (each ≤ 0.06″, rms ≤ 0.045″, plate rms < 0.15″ with ≥ 50 references, rate within 0.02″/min and PA
within 3° of Horizons, magnitude spread ≤ 0.25) in CI — limits about 2–3× the measured values. The example is a bright
(SNR 40–60), well-sampled, slow target, so it guards against regressions rather than measuring the error model; the
benchmarks below (ZTF, SDSS, blind set) do that. `test_centroid` checks the pixel-phase bias and model floor of the
centroid on synthetic stars (§ Systematics).

On synthetic data (`test_psf_fit`, L = 6 px, SNR ≈ 15) the trail model halves the centre scatter (0.03 px vs 0.06 px
for a star PSF) and recovers the centre to < 0.001 px at high SNR; a star 2.7 px away biases a lone fit by 2 px and 3×
in flux, while the joint fit recovers both to < 0.02 px and < 1 %.

ADES uncertainty calibration on the ZTF development benchmark (§ Benchmark, 20 sequences, 51 target measurements):
median |O−C|/√(rmsRA² + rmsDec²) = 0.78 (0.83 expected for correctly scaled Gaussian errors; 0.52 with the previous
out-of-sample-rms model, i.e. σ about 1.6× too large), 0.83 for SNR < 10 and 0.75 for SNR 10–30.

Faint-mover search (§4b) on the ZTF development benchmark (20 sequences; FP = tracks without a known object per
sequence; piazzi time for the 20 sequences):

| label | change | target det. | \|O−C\| med / rms / p95 | ratio | FP | time |
|---|---|---|---|---|---|---|
| v2-dev | before | 15/20 | 0.058 / 0.140 / 0.307″ | 0.78 | 0.10 | 16 s |
| v3-gain | FITS GAIN (6.2 e⁻/ADU) | 15/20 | 0.059 / 0.140 / 0.307″ | 0.79 | 0.10 | 16 s |
| v3-s1 | + stack (catalog/stationary mask) + pairs | 16/20 | 0.058 / 0.137 / 0.302″ | 0.78 | 0.20 | 34 s |
| v3-s2 | + counterpart/static/extended masks, ≥ 3 epochs, LOO | 16/20 | same | 0.78 | 0.50 | 31 s |
| v3-s3 | + ≥ 3 measured frames, masked vacated position rejects | 16/20 | 0.058 / 0.137 / 0.302″ | 0.78 | 0.15 | 30 s |
| v3-s4 | pairs: confirmation only from other epochs | 15/20 | 0.059 / 0.140 / 0.307″ | 0.79 | 0.10 | 30 s |
| v3-s5 | stack: no mask of bright single detections | 16/20 | 0.064 / 0.229 / 0.502″ | 0.87 | 0.10 | 30 s |
| v3-s6 | ≥ 2 measured frames (rejected) | 16/20 | same | 0.87 | 0.45 | 30 s |
| v3-s9 | pairs without 3rd measured frame + contamination tests | 16/20 | same | 0.87 | 0.15 | 30 s |
| **v3-dev** | pair_tracks = 0 (default), stack_search = 1 | **16/20** | 0.064 / 0.229 / 0.502″ | 0.87 | 0.10 | 31 s |

The gained target is 2000 SB142 (V 19.7, three forced measurements at SNR 3.5 / 4.9 / 9.5, |O−C| 1.1 / 0.5 / 0.6″
with ADES σ 0.63 / 0.35 / 0.33″, i.e. consistent with their errors). It alone moves the O−C rms and p95 (the other 51
measurements are unchanged); the stack also recovered known non-target objects (2008 TP133, 2005 EO341) and added no
unidentified track. The 2-detection path found no target and one unidentified track (v3-s9), so it is off by default.

Injection–recovery (`piazzi --inject 8`, 3 seeds × 20 sequences, PSF = frame seeing, magnitude on the frame zero
point (Gaia RP), 0.15–0.8″/min, recovery = a track with ≥ 2 measurements within 2″):

| mag | stack off | stack on |
|---|---|---|
| 18.0–18.5 | 85 % | 85 % |
| 18.5–19.0 | 88 % | 88 % |
| 19.0–19.5 | 78 % | 83 % |
| 19.5–20.0 | 42 % | 55 % |
| 20.0–20.5 | 10 % | 18 % |
| 20.5–21.0 | 1 % | 4 % |
| 21.0–21.5 | 0 % | 0 % |

Losses at bright magnitudes are objects on stars, at the edges, or in sequences with only two epochs.

### Fast movers (v7, synthetic — `test_trail`)
Trails rendered as segment ⊗ Gaussian PSF (σ = 1.27 px), 40 random angles/positions per cell, 40 bright stars per
160² frame; "SNR" is the optimal (matched-filter) SNR. Without the trail search (detection as in HEAD 063a20e) no trail was recognised:
10–15 px trails were split into 2–3 point detections (single detection in 4–36 of 40) and the nearest fragment was off
by an rms of 0.45–2.2 px along the trail.

| L (px) | SNR | detected | single detection | along / across bias (px) | along / across rms (px) | Cramér–Rao |
|---|---|---|---|---|---|---|
| 5 | 10 / 20 / 50 | 39 / 40 / 40 | 38 / 40 / 40 | ≤ 0.02 / ≤ 0.01 | 0.33–0.05 / 0.17–0.04 | 0.27–0.05 |
| 10 | 10 / 20 / 50 | 40 / 40 / 40 | 34 / 40 / 40 | +0.09, 0.00, +0.03 / ≤ 0.02 | 0.92, 0.19, 0.08 / 0.18, 0.09, 0.04 | 0.45, 0.22, 0.09 |
| 15 | 10 / 20 / 50 | 40 / 40 / 40 | 29 / 40 / 40 | −0.09, +0.03, +0.03 / ≤ 0.03 | 1.54, 0.32, 0.11 / 0.30, 0.10, 0.03 | 0.64, 0.32, 0.13 |

Full reductions (1″/px, 30 s): 12″/min over 4 frames (6-px trails, SNR 40) worst |O−C| 0.11″; 20″/min over 3 frames
(10 px, SNR 40) 0.15″; 19.8″/min at SNR 15 linked in all 4 frames, worst 0.79″ (≈ 2× the along-track Cramér–Rao
bound); trails perpendicular to their own motion are not linked; star-only fields: 1 spurious "trail" (a bright
unresolved star pair) in 10 dense frames, detection count changed by 1 of 1921. The trail centre is mid-exposure
(`Image::jd_mid`). Real-data NEA evaluation is pending a development set of fast movers (the blind ZTF set is held
out). On the slow-mover sets nothing changed: ZTF dev (v7-dev) and SDSS (sdss-v7) have the same detections and O−C
as v6, SDSS false tracks 1 → 0; run time ×2.8 (ZTF) and ×2.9 (SDSS) from the trail search.

## Injection tests
`piazzi --inject N` draws N synthetic movers into the real development images (ZTF dev, SDSS, NEA dev — never the
blind `fields_test.json`): the frame PSF (σx, σy of isolated references, `frame_psf`) convolved with the exposure
segment (`psf_pixel`), total flux = zero-point aperture flux / aperture correction, Poisson noise by the gain; truth is
exact (mid-exposure position through the first reduction's plate). Options: `--inject-rate`, `--inject-mag`,
`--inject-snr` (per-frame optimal SNR instead of magnitude), `--inject-color` (odd frames fainter), `--inject-cross`
(over a star in one frame), `--inject-edge` (leaves the cutout), `--inject-probe I`. `run.py inject --scenario
slow|trail|faint|cross|color|edge|all` aggregates completeness by rate and magnitude/SNR, O−C against the truth and
|O−C|/σ (results/inject-*.json). 6 objects per sequence, seed 1; ZTF dev 20 sequences, SDSS first 8:

| scenario | ZTF v7 → v8 | SDSS v7 → v8 |
|---|---|---|
| slow (0.15–0.8″/min, V 18–21.5) | 65.0 → 64.2 % | 58.3 → 58.3 % |
| trail (1–30″/min, V 16–19.5) | 82.5 → 82.5 % | 87.5 → 87.5 % |
| faint (SNR 2–4 per frame) | 2.5 → 2.5 % | 0 → 0 % |
| cross (over a star) | 17.5 → 20.8 % | 12.5 → 41.7 % |
| color (±0.5 mag between frames) | 80.0 → 79.2 % | 95.8 → 95.8 % |
| edge (leaves the cutout, ≥ 3 frames in) | 71.4 → 71.4 % | 36.1 → 58.3 % |

With `link_mag_e=1` (off by default) cross reached 65 % and edge 78 % on SDSS, and the NEA set gained 2019 AN5, but
one chance track appeared on the NEA set. SDSS trail |O−C|/σ (per axis, 0.674 calibrated): 0.17 (v7) → 0.44 (v8, J·C·Jᵀ) — still conservative. ZTF injected
movers at ≤ 10″/min are < 1.5 FWHM long in 30 s and stay point sources. Fixed-position extra forced measurements
(`stack_frame_snr_extra` = 2) did not change completeness at V 20–20.5 (3/20) and are off by default.

## Benchmark
Reproducible accuracy benchmark on public data (`tools/benchmark/`, Python 3 standard library only; details in
`tools/benchmark/README.md`). Images: ZTF science-image cutouts (12′) from IRSA; reference catalog: Gaia DR3 box printed
by `piazzi` (VizieR, falling back to the same box through the NOIRLab Astro Data Lab TAP when VizieR is down); known
objects: IMCCE SkyBoT (site I41); truth: JPL Horizons astrometric RA/Dec (`CENTER='I41@399'`, `QUANTITIES='1,36'`,
`EXTRA_PREC`) at the ADES `obsTime`s; quantity 36 gives the formal 3σ uncertainty of the ephemeris per axis
(`RA_3sigma`, `DEC_3sigma`, plane-of-sky arcseconds), stored per measurement as σ_eph = 3σ/3. When VizieR is down the
Gaia box comes, with the same columns, from the NOIRLab Astro Data Lab TAP and, failing that, from the ESA Gaia archive
TAP. Settings: `focal_length_mm=3063 pixel_width_mm=pixel_height_mm=0.015 mpc_code=I41 time_reference=0 catalog_faint_mag=19`. Timing: IRSA
cutouts carry no `DATE-OBS`; `OBSJD` (= UTC-OBS) is the exposure start and `EXPTIME`=30 s, so with time_reference=0 the
reported time is OBSJD + 15 s (mid-exposure).

```sh
python3 -I tools/benchmark/run.py select --n 30     # regenerate tools/benchmark/fields.json (versioned; fixed seed)
python3 -I tools/benchmark/run.py fetch             # optional: prefetch cutouts (IRSA often answers 503)
python3 -I tools/benchmark/run.py run --label base --piazzi cpp/build/piazzi
python3 -I tools/benchmark/run.py compare tools/benchmark/results/base.json tools/benchmark/results/new.json
```

Metrics: detection rate of the SkyBoT target (overall and by V), |O−C| median / rms / p95 and per axis, median ADES
σ (rmsRA⊕rmsDec), median |O−C|/σ (0.83 for well-calibrated 2-D Gaussian errors: median of a Rayleigh variable
1.177σ over √2σ), tracks without a known object (false positives) per sequence, run time. To separate Piazzi's error
from the error of the truth, three more views: (b) the same metrics on the *clean* measurements only — targets whose
Horizons σ_eph < 0.05″ on both axes (numbered or long-arc objects) — with the per-axis median |O−C|/σ (0.674
calibrated); (c) the *corrected* ratio |O−C|/√(σ² + σ_eph²); and the *within-sequence* ratio: residual minus the
sequence mean, scaled by √(n/(n−1)), over σ — the ephemeris error is constant over a 5-minute arc, so this isolates the
measurement error at the cost of also removing any Piazzi error common to the frames (the plates are independent, so
that is ≈ σ_pred/√n ≈ 0.002″). `run.py rescore <label>` recomputes all of it from the saved ADES files and the Horizons
cache. A track is associated
with the target (or a `known` object of the field file) when the median distance to its linearly predicted position
is < 10″. A track that matches none of them is then checked against **every** SkyBoT object of the field (cone of
0.30° for SDSS frames — two half-diagonals, the target can be anywhere in the frame — and 0.15° for 12′ ZTF cutouts,
at the sequence mid-time, site = sequence MPC code), with radius max(10″, min(3·ErrPos, 20″)) from SkyBoT's own
position uncertainty; only tracks matching nothing count as false positives (`fp_per_seq`; the old count is kept as
`fp_per_seq_target_only`). `run.py rescore <label>` re-applies this to an old result from its saved ADES files
(`cache/out/<label>`), and `run.py selftest` checks the association offline.

Base results (commit `ebc257f`, 21 sequences, 2019–2023, targets V 17.6–20.0, 0.1–0.7″/min main-belt):

| metric | value |
|---|---|
| sequences run / data failures | 20 / 1 (cutout at quadrant edge) |
| plates not solved | 0 |
| target detected | 13 / 20 (65 %): V 17–18 2/3, 18–19 4/6, 19–20 5/7, 20–21 2/4 |
| target measurements | 45 |
| \|O−C\| median / rms / p95 | 0.058″ / 0.201″ / 0.388″ |
| O−C mean RA·cosδ / Dec | −0.025″ / −0.001″ |
| O−C rms RA·cosδ / Dec | 0.177″ / 0.096″ |
| median ADES σ | 0.127″ |
| median \|O−C\|/σ | 0.52 (σ overestimated ≈ 2×) |
| unknown tracks per sequence | 0.80 (14 of 16 in one crowded galactic-plane field) |
| piazzi time (20 sequences) | 12.5 s |

Findings: astrometry of detected objects is excellent (median 0.06″) and the ADES σ is conservative; misses are
detection/linking failures (no moving track found for 7 targets, all plates solved), and a dense Milky-Way field with
two exposures 80 s apart produced 14 spurious links.

### Independent test set: SDSS (`tools/benchmark/fields_sdss.json`) — held out

`fields.json` (ZTF) is the development set; `fields_sdss.json` is a held-out test set from a different instrument,
reserved for the final evaluation (nothing in `piazzi` is tuned on it). Data: SDSS DR17 calibrated frames (2.5-m,
Apache Point, MPC 645, 0.396″/px, full 2048×1489 frames). The drift-scan camera images each run/camcol/field in the
order r, i, u, z, g, 71.7 s apart (53.9 s effective exposure); a sequence uses **r, i, z, g** (4 epochs over 287 s; u is
dropped — shallow and ≈1.5–2 mag fainter for asteroids, too shallow to help). Drift-scan timing:
mid-exposure of a point at row `y` is `TAI − EXPTIME/2 + y·EXPTIME/2048` (header TAI = first row; checked against the
per-object `TAI_<f>`/`rowc_<f>` in SkyServer photoObj to 0.08 s), converted to UTC with TAI−UTC = 32/33/34 s; each
image is stamped with the time at the target's row in that filter. Frames are converted (stdlib) from nanomaggies to
electrons (per-column calibration, median ALLSKY sky, approximate DR7 gain) as 16-bit FITS with `DATE-OBS` = UTC
mid-exposure and `GAIN=1`. Settings: `focal_length_mm=12501 pixel_width_mm=pixel_height_mm=0.024 mpc_code=645 time_reference=1 catalog_faint_mag=20`.
Targets: SkyBoT (`-loc 645`) objects V 15–21 in random SkyServer `Field` rows with |β| ≤ 12°; truth: Horizons
`CENTER='645@399'`. Final evaluation:

```sh
python3 -I tools/benchmark/run.py run --fields tools/benchmark/fields_sdss.json --label sdss-final --piazzi cpp/build/piazzi
```

Results (32 sequences):

| run | change | plates failed | detected | \|O−C\| med / rms / p95 | median \|O−C\|/σ | unknown tracks/seq |
|---|---|---|---|---|---|---|
| sdss-final | fixed floors (r_fix ≥ 2″, displacement ≥ 4″) | 4 | 1/32 | 0.015 / 0.018 / 0.025″ (1 object) | 0.68 | 0.25 |
| e1 | σ-based stationary/motion tests, two-pass plate, g-corrected vacated test | 0 | 26/32 | 0.034 / 0.079 / 0.183″ | 1.60 | 0.56 |
| **sdss-v4** | + Bonferroni factor (anchor pairs) in the chance test | **0** | **26/32** | 0.034 / 0.079 / 0.183″ | 1.60 | 0.38 |

v5 experiments (FP = tracks matching no SkyBoT object of the whole field; sdss-v4 rescored with this rule: 0.156, old
rule 0.375; v4-dev: 0.05 either way; the ratio for calibrated errors is 0.83):

| run | change | SDSS det. | SDSS \|O−C\| med / rms / p95 | σ med | ratio | FP | ZTF dev det. / ratio / FP |
|---|---|---|---|---|---|---|---|
| sdss-v4 | (rescored) | 26/32 | 0.034 / 0.079 / 0.183″ | 0.028 | 1.59 | 0.156 | 16/20 / 0.87 / 0.05 |
| e9 | photometric colour term (zp at solar colour) | 29/32 | 0.029 / 0.075 / 0.158″ | 0.028 | 1.55 | 0.094 | — |
| e10 | + rmsRA/rmsDec with 3 decimals | 29/32 | same | 0.021 | 1.63 | 0.094 | — |
| e11 | + obsTime ms, RA/Dec 6 decimals when σ < 0.1″ | 29/32 | 0.028 / 0.076 / 0.167″ | 0.021 | 1.33 | 0.094 | — |
| e12/e13 | + astrometric colour (DCR) term, always / if > 3σ (rejected) | 29/32 | 0.030 / 0.075 / 0.158″ | 0.021 | 1.39–1.40 | 0.094 | — |
| e14 | e11 + σ_sys from cross validation | 29/32 | 0.028 / 0.076 / 0.167″ | 0.024 | 1.15 | 0.094 | 15/20 / 0.77 / 0.05 |
| e15 | + amplitude test 1 mag + 3σ (was 1.5 mag flat) | 29/32 | same | 0.024 | 1.15 | 0.219 | 16/20 / 0.77 / 0.05 |
| e17 | + flux χ² (P < 10⁻³, σ_int 0.3) | 29/32 | same | 0.024 | 1.15 | 0.219 | — |
| e18 | χ² P < 1 %, σ_int 0.2, in linking and after refit | 29/32 | same | 0.024 | 1.15 | 0.031 | 15/20 / 0.77 / 0.05 |
| **sdss-v5 / v5-dev** | χ² only after refit | **29/32** | 0.028 / 0.076 / 0.167″ | 0.024 | 1.15 | **0.031** | **16/20** / 0.77 / 0.05 |

v5: detection by V 15–16 2/2, 16–17 6/6 (Vitja, Onizaki recovered by the colour term; 2000 WG as well); misses: the
two slow objects (0.17, 0.24″/min) and 2000 CZ68 (V 20.9). The remaining SDSS false positive is an SNR ≈ 3–6 track.
ZTF dev: same targets, |O−C| median 0.064 → 0.060″ (output precision), ratio 0.87 → 0.77 (σ slightly conservative
there, 0.83 expected). By SNR (SDSS v5, median ratio): log SNR < 1.5 1.24, 1.5–2 1.46, 2–2.4 0.99, > 2.4 1.06.

Detection by V (sdss-v4): 15–16 1/2, 16–17 5/6, 17–18 6/6, 18–19 6/6, 19–20 5/6, 20–21 3/6. ZTF dev set unchanged
(16/20, same O−C), unknown tracks 0.10 → 0.05/seq. Some "unknown" SDSS tracks are real asteroids missed by the
benchmark's association (e.g. sdss_004145_2_0087 track 1 = (162241) 1999 TZ217, SNR ≈ 100, 4.5″ from its geocentric
SkyBoT position); the rest are SNR 2–7 chance alignments. Misses: two slow objects (0.17, 0.24″/min ⇒ 0.8–1.2″ over
the arc, below 5σ for their SNR), two at V ≈ 20.3–20.9, two bright (V 15.6, 16.0).

### Systematics (SDSS held-out set, 32 sequences × 4 frames, `plate_debug=1`, 15 159 accepted references)

Residual of the references about the plate, before the v6 changes, binned by every quantity that could bias the target:

| test | result |
|---|---|
| pixel phase (x − round x, 5 bins per axis) | mean residual < 0.001″ in every bin (SDSS FWHM 2.5–4 px is well sampled). Synthetic injection (`test_centroid`: pixel-integrated Gaussian and Moffat β = 2.5 stars at known sub-pixel phases, Poisson noise, SNR ≈ 500): bias ≤ 0.0015 px at FWHM 2.5–3.5 px, ≤ 0.006 px for a Moffat star at FWHM 1.8 px; rms 0.003–0.006 px. The Gaussian fit evaluated at pixel centres is adequate; a pixel-integrated model was not needed. — *Anderson & King 2000, PASP 112, 1360* |
| position on the chip (6 bins in x and in y) | ≤ 0.003″ (significant at 4σ with n ≈ 2000 per bin, i.e. real distortion below the order-3 polynomial, but 10× below the measurement noise) |
| magnitude (10–21) | ≤ 0.003″ at any brightness: no magnitude equation (no flux-dependent CTE trailing along the scan) |
| colour (BP − RP 0.3–2.5) | +0.013″ (bluest) to −0.006″ (reddest) along the scan, half that cross-scan; removed by the plate colour term (§3) |
| SNR of the reference | bright references (SNR > 150, own noise 0.002″) scatter 0.010–0.013″ (x) / 0.012–0.018″ (y): the floor of the data (differential atmospheric image motion over 54 s on a 13′ field is ≈ 0.01–0.02″, *Lindegren 1980, A&A 89, 41*), now booked in σ_sys through the real σ_cat |
| timing (drift scan) | along-track mean offset of the bright targets vs rate: slope 0.007″ per ″/min ⇒ δt = 0.4 ± 1 s; the row-time formula is right. Within-sequence scatter along/across the motion 0.021/0.020″: the trail (0.7–2.3 px) does not add noise beyond its larger fitted FWHM |
| annual aberration | Gaia positions and Horizons astrometric positions are both unaberrated; the plate absorbs the common term; the first-order differential term over the 8.4′ semi-diagonal is 0.05″ but linear (absorbed by the affine terms), the second order is 1.2·10⁻⁴″; diurnal 8·10⁻⁴″ — nothing to correct |
| parallax of the references | Gaia RP < 20 references are > 100 pc as a rule (< 10 mas annual parallax); with 80 references and 3σ rejection the plate effect is ≲ 1 mas; not applied (`Plx` is read but unused) |
| proper motion | propagated from the catalogue epoch to the image epoch (Julian year of the first frame); Gaia errors propagated with (t − 2016.0) |
| DCR | quantified above (§3): 0.01–0.04″ expected in g, ≤ 0.01″ in r/i/z at the SDSS airmasses; handled by the colour term rather than by a refraction model (no observatory coordinates are needed) |
| PSF model | Moffat centroid (`psf_beta=2.5`, e22) instead of Gaussian: O−C unchanged (median 0.030 → 0.030″, rms 0.075 → 0.075″), one detection lost and one false track more; rejected. The synthetic test shows why: a Gaussian fit of a Moffat star costs 0.004–0.006 px, below everything else |

### v6: error model and systematics (both sets; SDSS held out, ZTF dev)

All v5 results were rescored with the Horizons uncertainties (no re-run). Columns: raw ratio |O−C|/σ (0.83 calibrated);
clean = σ_eph < 0.05″ (SDSS 92 measurements / 23 sequences; ZTF all 54 / 51); per-axis clean ratio RA / Dec (0.674);
within = within-sequence ratio and per axis.

| run | change | SDSS det. | \|O−C\| med / rms / p95 (clean) | σ med | ratio raw / clean / RA / Dec | within / RA / Dec | FP | ZTF det. / \|O−C\| med / rms / p95 / ratio / within / FP |
|---|---|---|---|---|---|---|---|---|
| sdss-v5 / v5-dev | rescored | 29/32 | 0.022 / 0.040 / 0.100 | 0.024 | 1.15 / 1.01 / 1.27 / 0.81 | 0.86 / 0.88 / 0.71 | 0.031 | 16/20 / 0.060 / 0.229 / 0.502 / 0.77 / 0.81 / 0.05 |
| e20 | real Gaia σ_cat + RUWE in the query | 29/32 | 0.024 / 0.040 / 0.101 | 0.025 | 1.09 / 0.95 / 1.01 / 0.71 | 0.80 / 0.85 / 0.67 | 0.062 | 15/20 / 0.056 / 0.139 / 0.308 / 0.75 / 0.76 / 0.05 |
| e21 | + plate colour term | 29/32 | 0.023 / 0.039 / 0.091 | 0.023 | 1.18 / 1.03 / 1.03 / 0.75 | 0.82 / 0.72 / 0.69 | 0.031 | 15/20 / 0.059 / 0.140 / 0.307 / 0.75 / 0.82 / 0.05 |
| e22 | + Moffat centroid β = 2.5 (rejected) | 28/32 | 0.022 / 0.040 / 0.093 | 0.024 | 1.09 / 0.94 / 0.97 / 0.70 | 0.80 / 0.72 / 0.71 | 0.062 | — |
| **sdss-v6 / v6-dev** | = e21 | **29/32** | 0.023 / 0.039 / 0.091 | 0.023 | 1.18 / 1.03 / 1.03 / 0.75 | **0.82** / 0.72 / 0.69 | **0.031** | **15/20** / 0.059 / 0.140 / 0.307 / 0.75 / **0.82** / 0.05 |

Reading: (1) the ephemeris explains the raw SDSS excess — the Horizons formal σ_eph is 0.001–0.005″ for the numbered
targets, 0.04–0.08″ for the short-arc ones, yet the sequence-mean offsets of numbered objects scatter by ≈ 0.02″ per
axis (10× the formal value: orbits at 1998–2008 epochs rest on pre-Gaia astrometry; Horizons says its values "can be
optimistic" with catalogue systematics), so even the clean ratio stays ≈ 1.0 while the within-sequence ratio, which
removes the constant, is 0.82 on both instruments (0.83 expected) — the model is consistent between SDSS (σ ≈ 0.023″,
FWHM 1.4″, 54 s) and ZTF (σ ≈ 0.07″, 2″, 30 s) with no adjustment factor. (2) The real σ_cat (e20) fixed the
SNR-dependence: with the constant 0.01″, targets with σ < 0.012″ had a within ratio of 1.74 and those with σ 0.03–0.05″
0.69; now 0.90 (σ 0.012–0.02″) and 0.76 (0.03–0.05″), and the per-axis RA/Dec asymmetry (1.27/0.81 → 1.03/0.75)
is gone — the "RA excess" of v5 was the under-estimated σ_sys of bright targets, not a bias. (3) The colour term (e21)
moves the clean mean O−C from +0.008/+0.007″ to +0.006/+0.001″ and lowers p95 0.101 → 0.091″; it makes σ slightly
smaller (part of σ_sys was colour), hence the raw ratio rises. (4) On ZTF, e20 loses 2000 SB142 (V 19.7, the
synthetic-tracking detection at stack SNR 8.3 against a 6.7 threshold, per-frame SNR 3.5/4.9/9.5): the per-star σ_cat
changes the plates by a few mas and the candidate falls under the threshold; `plate_ruwe=0` restores the 46 stars with
RUWE > 1.4 and the reference counts but not the detection. The ZTF rms/p95 drop (0.229 → 0.140″) is only that lost
noisy target; the other 51 measurements are unchanged to 0.003″. Kept, because the σ_cat change is correct physics and
the detection was at the edge of the false-alarm control.

Remaining misses, checked with `plate_debug` (all are limits of the data or of the false-alarm control, not fixed by
tuning): SDSS 2000 VX33 (0.17″/min, 0.9″ over 287 s; detected at SNR 50/37/11/42 but the 72-s displacement, 0.2″, is
≈ 5·√2·σ for the SNR-50 frames and far below it for z, so two of the four detections are classed stationary and 2
free detections cannot make a track); 2006 EK66 (V 20.2, 0.24″/min; SNR 32/22/5/31, the z detection stationary at
its σ, three free detections remain, and a 3-point track with one verification frame among 300–440 free detections per
frame has an expected chance count E ≈ 4 ≫ `false_alarm` = 0.05 at the default `max_rate` of 100″/min); 2000 CZ68
(V 20.9, SNR 14/13.5/–/11.6, same E ≈ 9). ZTF: 2000 YY39 and 2008 UC72 (two epochs: exposure pairs 80 s apart;
UC72 also SNR 4.4/3.0 in two frames only), 1996 OA (cutout at the quadrant edge, data error), 2001 SX8 (at x = 7–19 px
of the cutout in frames 2–3, outside it in frame 1: the SkyBoT position used to centre the cutout was off by ≈ 6′),
2006 WX108 (V 20.0 at the 19.7–20.1 limiting magnitude: SNR 2.4 and 5.6 in two frames, nothing in the third).

### v8 (development sets and the blind set)
| set | found | median / p95 |O−C| | false/seq | piazzi time |
|---|---|---|---|---|
| ZTF dev v7 (catalog_faint_mag 19) | 15/20 | 0.059 / 0.307″ | 0.05 | 90 s |
| ZTF dev v8 (catalog_faint_mag 19) | 15/20 | 0.059 / 0.307″ | 0.05 | 49 s |
| SDSS v7 → v8 | 29/32 → 29/32 | 0.030 / 0.163″ (same) | 0 → 0.031 | 457 → 362 s |
| NEA dev (13, 1–12.5″/min) eced3c4 → v8 | 11/13 → 11/13 | 0.101 / 0.331″ (same) | 0.077 → 0.077 | 77 → 42 s |
| **ZTF blind** be8aa23 → eced3c4 → v8 | 16 → 15 → 16/21 | 0.128/0.898 → 0.106/0.317 → 0.108/0.323″ (well-known orbits: 0.087″) | 0 → 0 → 0 | 33 → 99 → 57 s |

Querying Gaia to catalog_faint_mag 20.5 (the ZTF image depth) found no new targets and added one 3-point false track
(ztf_20201216), so the benchmark keeps catalog_faint_mag 19. The blind set (`fields_test.json`, 23 sequences including 10 fast
NEAs) was run only three times, each after development was frozen: be8aa23 (before trail work), eced3c4 (trails) and
v8. The trail work cut the blind p95 from 0.90″ to 0.32″ (fast NEAs no longer biased by fragmented trails) at
unchanged false detections (0). The SDSS v8 "false" track is 4 coherent trails at 5.3″/min, V ≈ 21, fit rms 0.03″, absent
from SkyBoT — probably an unknown asteroid found by the deferred trail search, counted as false. On the NEA set the
deferred search first let a galaxy whose centroid shifts between g and r (0.6″) link as a 0.01″/min "mover"; extended
detections (FWHM > 1.3·seeing) now must move ≥ 0.5·FWHM between anchors (`link_extended`). 2000 SB142 (lost since v6)
is recovered with `stack_fixed_pos=1` (fixed-position fallback when the free forced fit fails in one frame) but that
option makes 12 chance tracks out of 150 random transients per frame in `test_linking`, so it is off.
`select`/`select-targets` now reject cutouts that would cross the quadrant edge.

### Comparison with ZTF's own MPC astrometry (`run.py vs-mpc`)
`vs-mpc` pairs each Piazzi target measurement with the observation that the survey itself (I41 = ZTF, 645 = SDSS) sent
to the Minor Planet Center (public `get-obs` API, ADES_DF, cached in `cache/mpcobs/`), nearest in time within 60 s.
Both are reduced against Horizons at their own times, Δ = (P−H)(t_P) − (MPC−H)(t_MPC); a positions-only Δ (both
against Horizons at Piazzi's mid-exposure time) is also kept.

| set | pairs (seq.) | |P−MPC| med/rms/p95 | positions only | mean bias RA·cosδ/Dec | |O−C| med Piazzi / MPC | Piazzi closer |
|---|---|---|---|---|---|---|
| ZTF dev (v9) | 7 (4) | 0.069/0.096/0.151″ | 0.054/0.058/0.084″ | +0.054/−0.010″ | 0.020 / 0.064″ | 7/7 |
| NEA dev (v9) | 9 (3) | 0.484/0.509/0.697″ | 0.067/0.068/0.095″ | −0.464/+0.039″ | 0.095 / 0.512″ | 9/9 |
| SDSS (v8, held out, read only) | 44 (22) | 0.054/0.098/0.194″ | 0.121/0.138/0.257″ | −0.003/−0.015″ | 0.033 / 0.083″ | 36/44 |

The sample is small: ZTF reported only 7 of the 31 Piazzi-measured ZTF/NEA sequences to the MPC (the other nights have
no I41 observation of the target within ±1 h). ZTF's 2019–2020 MPC `obstime`s are 15.0 s early (exposure start)
while the positions are mid-exposure: their |O−C| falls from 0.138″ to 0.071″ (main belt) and from 0.512″ to 0.102″
(NEAs) when evaluated at mid-exposure; the 2021 submission (Hindemith) agrees with Piazzi's time to < 0.3 s. The
positions-only comparison is therefore the fair one for ZTF: Piazzi and ZTF agree to 0.05–0.07″ median, with a
+0.05″ RA offset (ZTF's records name catalog Gaia1; Piazzi uses Gaia DR3). As reported, Piazzi is closer to Horizons
in every pair; with ZTF's time corrected to mid-exposure Piazzi is still closer in 5/7 (main belt, median 0.020″ vs
0.071″) and 6/9 (NEAs, 0.095″ vs 0.102″) pairs. SDSS reported to the MPC as well (UCAC1, 0.001 s/0.01″ precision, no rms); its times are
right, Piazzi and SDSS agree to 0.054″ median and Piazzi is closer to Horizons in 36 of 44 pairs (0.033″ vs 0.083″).

## Known limitations
- Single-frame depth: on ZTF 30 s images the shift-and-stack search (§4b) raises completeness at V 19.5–20 from 42 to
  55 % (injection test), but it emits only frames with SNR ≥ 3 and needs ≥ 3 such frames, so objects at SNR 2–3 per
  frame are detected in the stack but usually not reported; their low-SNR measurements also have large (correctly
  stated) errors.
- Sequences with only two epochs (two exposures minutes apart plus one more, e.g. 2000 YY39, 2008 UC72) cannot
  confirm a velocity: any pair of sources across two epochs is a valid line. Neither the stack nor the 2-detection
  path accepts them.
- Objects that leave the field (2001 SX8) or are not detected and invisible in a frame (2006 WX108) are not linked.
- The 2-detection path (`pair_tracks`) is off by default; the stack search costs ~0.5–1.5 s per 12′ ZTF cutout, and
  `stack_max_rate` (1″/min) excludes fast movers from it.
- SDSS: the raw median |O−C|/σ against JPL stays above 0.83 (1.18 in v6) because of the ephemeris: Horizons' formal σ
  is 0.001–0.005″ for the numbered targets at 1998–2008 epochs, but their four-filter mean offsets scatter by ≈ 0.02″
  per axis (up to 0.08″ for (2) Jillpernell-class objects, 0.18–0.26″ for the V ≈ 20.8 short-arc ones), consistent
  across four independent plates with σ_pred ≈ 0.003″ — so the offsets are in the truth, which the ADES σ must not
  contain. Within-sequence the ratio is 0.82 (0.83 expected) and per axis 0.72/0.69 (0.674), see Benchmark.
- Remaining SDSS systematics after the v6 changes (§ Benchmark, "Systematics"): nothing above the 0.003″ level is
  left in the references; the target's O−C mean on the clean set is +0.006″/+0.001″.
- The motion threshold is statistical (k = 5σ): objects that move less than ~5·√2·σ over the arc (≈ 1″ at SNR 10
  with 1.4″ seeing) are indistinguishable from fixed sources and are not reported.
- The shift-and-stack search keeps its 1.5·FWHM minimum displacement: its static-sky mask removes anything present
  at the same pixel in all frames, so sub-FWHM movers are left to the detection-level linker.
- Plate matching ignores the expected position angle (any rotation is accepted).
- Aperture photometry of stars is not deblended (only the moving target is fitted jointly with its neighbours).
- Trail fitting does not improve astrometry for trails ≲ 1 FWHM, so the automatic `trail_fit` uses it only from
  1.5 FWHM. The stack and pair searches (`stack_max_rate` 1″/min) use a point PSF; at their rates v·t_exp ≪ FWHM, and
  raising `stack_max_rate` to NEA rates would multiply the velocity grid by (v_max/1″/min)² (search cost and look-elsewhere
  threshold) for objects that are bright enough to be detected per frame — fast movers are left to the trail
  detection + linker.
- v8 time target (≤ 1.3× the pre-trail time) is not met on SDSS (≈ 2.4×): with galaxies masked the fits shrink, but
  the kernel loop (0.6–0.8 s per 2048×1489 frame) remains; `trail_fit_iter` 25 (instead of 60) failed `test_trail`.
- Not done in v8: curvature tolerance for fast movers, brightness slack for blended trail detections; SNR 2–4 objects
  (injection ~0 %); 2019 AN5 and 2018 VW3 (NEA dev, 11–12″/min) are still missed with the default settings.
- Trail detection limits: trails crossing a star are usually not recognised (the point model of the star wins) and
  fall back to point fragments; at optimal SNR ≈ 10 about 10 % of 10-px trails stay fragmented (the two-PSF model
  ties, or the capsule SNR is below k + 0.5); trails much longer than `trail_max_len` (30 px) are fitted inside a capsule
  that may cut their ends; only the 40 most elongated components per frame are fitted, so a faint trail in a field
  crowded with galaxies can be skipped; the search costs ~0.3–1.5 s per frame (ZTF dev run time ×3.7, SDSS ×2, mostly
  the Levenberg–Marquardt fits of elongated galaxies and blends, which also become "trail" detections — static, so
  they do not form tracks). Unresolved bright star pairs/trios can still pass as a short trail (1 in 10 dense synthetic
  star fields). The along-track centroid error grows as √(σ² + L²/12)/SNR (Cramér–Rao): a 10-px trail at SNR 15 has
  ≈ 0.3–0.4 px per-frame scatter, which no centroiding removes.
- The PSF model is an axis-aligned elliptical Gaussian (no rotation, no wings); the aperture correction absorbs the
  flux scale but not shape errors.
