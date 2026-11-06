#!/bin/sh
# Report format: generates ADES (PSV) and MPC 80-column reports from the bundled example ((5157) Hindemith, ZTF) and
# applies the rules of tools/benchmark/run.py check-ades (ADES 2022 + 80-column format). Telescope holds only the free
# description: the # telescope design/aperture must be inferred from it.
# Usage: sh cpp/tests/report_format.sh [build-dir]   (from the repository root)
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd); BUILD=${1:-$ROOT/cpp/build}; OUT=$(mktemp -d)
E="$ROOT/cpp/gui/assets/example"
"$BUILD/piazzi" --cfg "$E/example.cfg" --set observer=Test --set measurer=Test --set submitter=Test \
  --set mpc_code=568 --set "telescope=0.30-m f/8 Reflector + CCD" \
  --catalog "$E/gaia_dr3.tsv" --center "07 11 43.534" "+22 24 17.78" --object "5157" --out "$OUT" \
  "$E/ztf_hindemith_0.fits" "$E/ztf_hindemith_1.fits" "$E/ztf_hindemith_2.fits" > "$OUT/log" 2>&1
grep -q '^! design Reflector' "$OUT/ADESReport.psv" && grep -q '^! aperture 0.30' "$OUT/ADESReport.psv" ||
  { echo "FAIL: # telescope not inferred from the description"; exit 1; }
python3 -I "$ROOT/tools/benchmark/run.py" check-ades --ades "$OUT/ADESReport.psv" --mpc "$OUT/MPCReport.txt"
