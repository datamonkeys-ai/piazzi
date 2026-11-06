# Third-party notices

Piazzi is MIT-licensed. It bundles or builds against the following third-party components:

| Component | Where | License |
|---|---|---|
| Source Serif 4 (Adobe) | `cpp/gui/assets/fonts/serif.ttf` | SIL Open Font License 1.1 — `cpp/gui/assets/fonts/OFL-serif.txt` |
| JetBrains Mono | `cpp/gui/assets/fonts/mono.ttf` | SIL Open Font License 1.1 — `cpp/gui/assets/fonts/OFL-mono.txt` |
| Dear ImGui (fetched at build time) | desktop app | MIT |
| SDL2 (system or fetched at build time) | desktop app | zlib |
| Earth ephemeris table | `cpp/src/earth_cheb.inc` | Chebyshev fit to ERFA `epv00` (BSD-3-Clause, derived from IAU SOFA); the table is a numerical fit, not ERFA code |
| JPL Horizons ephemeris of (5157) Hindemith | `cpp/tests/data/horizons_hindemith.txt` (test data) | NASA/JPL Solar System Dynamics, Horizons API (public data; query recorded in the file header) |
| Osculating elements of 6 bodies (MPCORB/COMET layout) | `cpp/tests/data/orbits_*.dat` (test data) | NASA/JPL Solar System Dynamics, Horizons API (public data; see `cpp/tests/data/README.md`) |
| ZTF image crops of (5157) Hindemith, 2021-02-06 | `cpp/gui/assets/example/*.fits` (copied to `web/example/`, `ios/Piazzi/Example/`) | Zwicky Transient Facility public data via NASA/IPAC IRSA; credit text in `cpp/gui/assets/example/README.md` ("Based on observations obtained with the Samuel Oschin Telescope 48-inch and the 60-inch Telescope at the Palomar Observatory as part of the Zwicky Transient Facility project…") |
| Gaia DR3 star catalog excerpt for the example | `cpp/gui/assets/example/gaia_dr3.tsv` (and copies) | ESA/Gaia/DPAC (CC BY-SA 3.0 IGO), served by the NSF NOIRLab Astro Data Lab TAP |
| Piazzi icon | `packaging/icon/` | original artwork, MIT (same as Piazzi) |

Online services queried at run time (not redistributed): ESA Gaia DR3 via CDS VizieR, IMCCE SkyBoT, CDS SIMBAD and
Sesame, AAVSO VSX and the IAU Minor Planet Center. Please acknowledge them in publications as their terms request.
