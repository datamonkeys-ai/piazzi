#!/bin/sh
# End-to-end accuracy on the bundled example (cpp/gui/assets/example: 3 ZTF frames of (5157) Hindemith):
#  1. runs the CLI on the 3 FITS with the example Gaia DR3 catalog and checks that the ADES/MPC reports contain 3 lines
#     of the object with positions within ≤ 0.06" of JPL Horizons (cpp/tests/data/horizons_hindemith.txt);
#  2. runs test_accuracy (plates, rate, magnitudes, identification).
# Usage: sh cpp/tests/accuracy.sh [build-dir]   (from the repository root)
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd); BUILD=${1:-$ROOT/cpp/build}; OUT=$(mktemp -d)
E="$ROOT/cpp/gui/assets/example"
"$BUILD/piazzi" --cfg "$E/example.cfg" --set observer=Test --set measurer=Test --set submitter=Test \
  --catalog "$E/gaia_dr3.tsv" --center "07 11 43.534" "+22 24 17.78" --object "5157" --out "$OUT" \
  "$E/ztf_hindemith_0.fits" "$E/ztf_hindemith_1.fits" "$E/ztf_hindemith_2.fits"
n=$(grep -c '^5157 ' "$OUT/ADESReport.psv" || true); m=$(grep -c '^05157 ' "$OUT/MPCReport.txt" || true)
[ "$n" = 3 ] && [ "$m" = 3 ] || { echo "FAIL: expected 3 object lines (ADES $n, MPC $m)"; exit 1; }
# ADES O−C (ra/dec in degrees, 6 decimals) against Horizons, in image order
grep -v '^#' "$ROOT/cpp/tests/data/horizons_hindemith.txt" | awk '{for(i=3;i<=NF;i++) if ($i ~ /^[0-9.]+$/) {print $i, $(i+1); break}}' > "$OUT/hz"
grep '^5157 ' "$OUT/ADESReport.psv" | awk -F'|' '{gsub(/ /,"",$7); gsub(/ /,"",$8); print $7, $8}' > "$OUT/obs"
paste -d' ' "$OUT/obs" "$OUT/hz" | awk 'BEGIN{bad=0} {c=cos($2*3.14159265358979/180); dra=($1-$3)*3600*c; dde=($2-$4)*3600; d=sqrt(dra*dra+dde*dde);
  printf "ADES image %d: dRA* %+.3f\" dDec %+.3f\" |%.3f\"|\n", NR, dra, dde, d; if (d>0.06) bad=1} END{exit bad}'
"$BUILD/test_accuracy" "$ROOT"
