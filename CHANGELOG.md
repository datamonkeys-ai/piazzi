# Changelog
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow [SemVer](https://semver.org/).

## [0.2.1] — 2026-11-06
### Known issues
- The Windows build is cross-compiled and packaged but not yet tested on real Windows hardware.

### Changed
- **Settings format (breaking):** every setting key now has its own snake_case name, for example
  `focal_length_mm`, `pixel_width_mm`, `flip_horizontal`, `catalog_faint_mag`, `plate_order` and `trail_detect`. Older
  key names are no longer read. Update existing `.cfg` files and saved profiles (the bundled example already uses the
  new names).
- Code comments, log messages, CLI help, CI and benchmark reports are in English.
- Added `CITATION.cff`, a README screenshot and badges. `tools/benchmark/results/` keeps only the final runs (see its
  README).
- Tests no longer need any external download: image tests use the bundled ZTF example ((5157) Hindemith,
  `cpp/gui/assets/example/`) with a new JPL Horizons reference (`cpp/tests/data/horizons_hindemith.txt`, station I41,
  query in the header). Limits are ~2–3× the measured accuracy (O−C max 0.021″, rms 0.016″ → ≤ 0.06″ / 0.045″; plate
  rms ≤ 0.053″ → < 0.15″ with ≥ 50 references; rate/PA vs Horizons within 0.02″/min and 3°; magnitude spread
  0.09 → ≤ 0.25). `test_orbits` uses versioned osculating elements from Horizons in MPCORB/COMET layout
  (`cpp/tests/data/orbits_*.dat`, origin in `cpp/tests/data/README.md`); `test_image_io` adds a synthetic
  uncompressed SBIG file (pixels, exposure, mid-exposure time).
- Tests no longer download third-party sample data: the external download script, its CI step and the old reference
  files were removed. DESIGN.md "Accuracy" now describes the bundled ZTF example.
- README/THIRD_PARTY_NOTICES describe `.cfg` support neutrally (INI-style configuration files).
- `run.py` writes paths relative to the repository root (or just the file name when outside it) in results and
  messages (`run`, `inject`, `vs-mpc`, `check-ades`); existing results under `tools/benchmark/results/` were rewritten
  without local absolute paths (`<local>/…`).
- Core error messages in the C API are now English and actionable, with a stable `code` (`no_images`,
  `image_unreadable`, `catalog_empty`, `internal`) and optional `detail` in the error JSON; the apps translate by code.
- `astro_reduce` honours `plate_ruwe` like the CLI (default 1.4, unchanged behaviour).
### Fixed
- Without `--object`/object name the API no longer labels every mover "Piazzi": the designation is empty, ADES carries
  only the `trkSub` and the MPC report uses the per-track temporary designation; apps show the `trkSub` or "objeto N".
- Result JSON never contains NaN/Inf (written as `null`); added per measurement `trksub`, `sig_ra`, `sig_dec` (″, the
  ADES σ), `trail_len` (″), `remarks`, and per mover `trksub` (existing fields unchanged).
- Undefined behaviour in the linking grid hash (left shift of a negative value).
- Shift-and-stack search: forced measurements with non-positive flux are discarded (no `log10` of ≤ 0).
- `astro_image_catalog`: trail fragments get type `T` and attached trails `L` (`R` stays "reference").
- Dockerfile header no longer claims that the image runs the tests.
- ADES `# telescope` no longer written with an empty `design` and `aperture 0.0` when `telescope_design`/`telescope_aperture_m` are
  unset: both are derived from the `telescope` description (all ZTF/SDSS/NEA benchmark reports were affected).
- Trail centroid covariance (ADES σ of trailed movers): Poisson-weighted sandwich covariance of the unweighted trail fit
  instead of a uniform χ²/(n−p) scale (`trail_cov_poisson`, default 1); SDSS trail injection |O−C|/σ_c 0.795 → 0.706.
  `run.py inject` adds per-axis ratios against σ_c (the injection O−C has no plate error by construction).
- 80-column report: a seconds field rounding to 60 no longer prints "-0.0" (e.g. `-00 58 -0.0`), which broke the Dec
  columns 45–56 (2 SDSS measurements); same clamp in RA.
### Added
- `astro_settings_schema()`: JSON list of the exposed settings (key, group camera/detection/plate/motion/reports/
  advanced, English and Portuguese label and hint, type number/int/bool/enum/text with options, min/max, core default),
  including the new `gain`, `trail_fit` (−1 auto, 0 off, 1 trail, 2 star PSF), `trail_detect`, `stack_search`,
  `plate_weights`, `plate_color`, `motion_sigma` and `plate_ruwe`. Desktop, iOS and web settings screens are generated from
  it: collapsible groups ("Advanced" closed), checkboxes, menus, range-checked numeric fields, hints, restore defaults.
- English user interface by default with optional Portuguese in the three apps (desktop/web remember the choice; iOS
  follows the system language with an in-app override). Strings centralised per app.
- "Open example": bundled offline example (three 6′ ZTF crops of (5157) Hindemith, 2021-02-06, Gaia DR3 excerpt,
  `example.cfg`, ~1.7 MB) in `cpp/gui/assets/example/`, `web/example/` and the iOS bundle; MPC submission is blocked
  for example data and its observatory code is `XXX`. Credits (ZTF/IRSA, ESA/Gaia/DPAC) in the example README and
  THIRD_PARTY_NOTICES.
- App icon (`packaging/icon/piazzi.svg`, PNG 16–1024, `.icns`, `.ico`, Windows `.rc`, iOS AppIcon, macOS bundle icon,
  `packaging/icon.png` for the AppImage, SDL window icon).
- `run.py check-ades`: offline ADES 2022 (PSV) and MPC 80-column format checker for the benchmark outputs or single
  files; ctest `report_format` applies it to the reports of the bundled example.
- `astro_version()`, `piazzi --version`, single version source (`cpp/include/piazzi/version.h`) used by CMake and
  shown in the desktop, iOS and web home screens.

## [0.2.0] — 2026-11-05
### Added
- Trailed fast movers (NEAs): trail detection and fitting (segment ⊗ PSF), deferred trail search after the plate
  solution, anisotropic trail covariance in the ADES σ. Blind ZTF p95 |O−C| 0.90″ → 0.32″.
- Calibrated ADES uncertainties (median |O−C|/σ 0.78 on ZTF dev vs 0.83 expected; previously σ ≈ 1.6× too large).
- Shift-and-stack search for faint slow movers (injection: V 19.5–20 completeness 42 → 55 %).
- `piazzi --inject`: synthetic movers drawn into real images for injection–recovery tests.
- Benchmarks (`tools/benchmark`): ZTF development set, held-out SDSS set and a blind ZTF set run only after
  development was frozen (16/21 found, no false tracks).
### Known limitations
- SNR 2–4 objects are rarely reported; two-epoch sequences cannot confirm a velocity; some 11–12″/min NEAs are still
  missed; trail search costs ≈ 2.4× run time on SDSS frames. See DESIGN.md.

## [0.1.0] — 2026-11-02
First public release.
### Added
- Dependency-free C++17 core with a stable C API: FITS/SBIG reading, dark/flat, source detection, plate reduction
  against Gaia DR3 (and UCAC4, CMC-15, PPMXL, ATLAS-RefCat2), photometry, moving-object linking, Track & Stack,
  manual measurement, offline identification (MPCORB/COMET) and ADES/MPC/photometry reports.
- Desktop app (macOS, Windows, Linux), iOS/iPadOS/Mac Catalyst app, web app (WebAssembly) and `piazzi` CLI.
- Independent algorithms from published literature (see DESIGN.md); accuracy tests against JPL Horizons.
- SkyBoT identification, "What is this?" lookup (SIMBAD, VSX, Gaia, SkyBoT) and MPC submission with test mode.
