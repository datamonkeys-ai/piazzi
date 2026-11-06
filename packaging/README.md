# Desktop app packaging (Piazzi)

Everything runs locally on this Mac; nothing is published or uploaded. Output goes to `dist/`.

| Script | Output | How |
|---|---|---|
| `packaging/build_mac.sh` | `dist/Piazzi-mac.dmg` | Universal Release (arm64+x86_64), static SDL2, signed with Developer ID + hardened runtime + timestamp |
| `packaging/build_linux.sh [amd64\|arm64]` | `dist/Piazzi-linux-x86_64.AppImage` / `-aarch64.AppImage` | Docker `ubuntu:22.04` (glibc ≥ 2.35), appimagetool; `.tar.gz` fallback |
| `packaging/build_windows.sh` | `dist/Piazzi-windows-x64.zip` | Docker `debian:bookworm` + mingw-w64 (posix); single static exe, GUI subsystem |

Requirements: Xcode CLT, cmake, ninja, Docker Desktop. The first run downloads the SDL2 2.30.9 source (FetchContent).

## CMake
- `cpp/CMakeLists.txt` compiles all `cpp/gui/*.cpp` (glob with `CONFIGURE_DEPENDS`).
- `-DASTRO_FETCH_SDL2=ON` (used by the scripts) builds static SDL2 from source → self-contained app.
  Without it, the system SDL2 is used (`find_package`), as before.
- `cpp/gui/assets/` is copied to `Piazzi.app/Contents/Resources/assets/` (macOS) or to `assets/` next to
  the executable (Windows/Linux). In the AppImage it lives in `usr/bin/assets`, and `AppRun` does `cd` to `usr/bin`
  before running (so relative `assets/...` paths work; file paths passed on the command line must be
  absolute).
- Windows (MinGW): `-static -static-libgcc -static-libstdc++ -mwindows`; `_popen/_pclose` are available in msvcrt.

## Variables
- mac: `IDENTITY=...` (e.g. `Developer ID Application: Name (TEAMID)`; without it the app is signed ad hoc, for local use only), `ARCHS="arm64"` for a single architecture.
  If `--timestamp` fails (no network), it signs again without a timestamp and warns.
- Windows: `WINE_TEST=1` runs the exe under wine+xvfb in a container just to check that it starts.
- All: `SRC=/path/to/repo` to build another tree (e.g. `git archive HEAD | tar -x -C /tmp/x`).

## Notes for users
- **macOS**: release DMGs are signed with a Developer ID and notarized by Apple, so they open without warnings.
  `build_mac.sh` notarizes and staples when `NOTARY_KEY`, `NOTARY_KEY_ID` and `NOTARY_ISSUER` (App Store Connect API
  key) are set; without them the DMG is only signed and Gatekeeper asks for right-click → Open on first launch.
- **Windows**: unsigned exe → SmartScreen shows "Windows protected your PC": More info → Run anyway.
  The app uses the system `curl.exe` (Windows 10 1803+).
- **Linux**: `chmod +x Piazzi-linux-*.AppImage && ./Piazzi-linux-*.AppImage`. Requires libGL/X11 or Wayland
  (loaded dynamically by SDL) and FUSE 2 (`libfuse2`); without FUSE: `./Piazzi-*.AppImage --appimage-extract-and-run`.

## Technical note (AppImage under emulation)
appimagetool is itself an AppImage; the magic bytes `AI\x02` at offset 8 make binfmt (Docker's Rosetta/qemu on Mac)
refuse to execute it ("Exec format error"). The script zeroes those 3 bytes in a copy and uses `--appimage-extract`.
The generated AppImage is normal (its bytes are left intact).
