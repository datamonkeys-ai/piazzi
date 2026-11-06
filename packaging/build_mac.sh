#!/usr/bin/env bash
# Universal Release build (arm64+x86_64) of Piazzi.app with static SDL2; signs it and produces dist/Piazzi-mac.dmg.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$ROOT/cpp}"
BUILD="${BUILD:-$ROOT/build-mac}"
DIST="$ROOT/dist"
IDENTITY="${IDENTITY:--}"   # "Developer ID Application: Name (TEAMID)"; without IDENTITY signs ad hoc (local use only)
ARCHS="${ARCHS:-arm64;x86_64}"
mkdir -p "$DIST"
cmake -S "$SRC" -B "$BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release -DASTRO_FETCH_SDL2=ON \
  -DCMAKE_OSX_ARCHITECTURES="$ARCHS" -DCMAKE_OSX_DEPLOYMENT_TARGET=11.0
cmake --build "$BUILD" --target piazzi_app
APP="$BUILD/Piazzi.app"
lipo -info "$APP/Contents/MacOS/Piazzi"
if otool -L "$APP/Contents/MacOS/Piazzi" | grep -qE '/opt/homebrew|/usr/local'; then
  echo "ERROR: non-system dependency found"; otool -L "$APP/Contents/MacOS/Piazzi"; exit 1
fi
if [ "$IDENTITY" = "-" ]; then
  codesign --force --deep -s - "$APP"   # ad hoc: runs on this Mac; set IDENTITY to distribute
elif security find-identity -v -p codesigning | grep -qF "$IDENTITY"; then
  if ! codesign --force --deep --options runtime --timestamp -s "$IDENTITY" "$APP"; then
    echo "WARNING: timestamp failed; signing without --timestamp"
    codesign --force --deep --options runtime -s "$IDENTITY" "$APP"
  fi
  codesign --verify --deep --strict --verbose=2 "$APP"
else
  echo "WARNING: identity '$IDENTITY' not found; app NOT signed"
fi
STAGE="$(mktemp -d)"; cp -R "$APP" "$STAGE/"; ln -s /Applications "$STAGE/Applications"
rm -f "$DIST/Piazzi-mac.dmg"
hdiutil create -volname Piazzi -srcfolder "$STAGE" -ov -format UDZO "$DIST/Piazzi-mac.dmg"
rm -rf "$STAGE"
codesign -s "$IDENTITY" "$DIST/Piazzi-mac.dmg" 2>/dev/null || true
# Optional notarization with an App Store Connect API key:
#   NOTARY_KEY=path/AuthKey_XXXX.p8 NOTARY_KEY_ID=XXXX NOTARY_ISSUER=uuid IDENTITY="Developer ID Application: ..." sh packaging/build_mac.sh
if [ -n "$NOTARY_KEY" ] && [ -n "$NOTARY_KEY_ID" ] && [ -n "$NOTARY_ISSUER" ]; then
  xcrun notarytool submit "$DIST/Piazzi-mac.dmg" --key "$NOTARY_KEY" --key-id "$NOTARY_KEY_ID" --issuer "$NOTARY_ISSUER" --wait
  xcrun stapler staple "$DIST/Piazzi-mac.dmg"
  spctl -a -t open --context context:primary-signature -v "$DIST/Piazzi-mac.dmg"
fi
ls -lh "$DIST/Piazzi-mac.dmg"
