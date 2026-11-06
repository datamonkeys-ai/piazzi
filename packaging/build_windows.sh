#!/usr/bin/env bash
# Cross-build Windows x64 (mingw-w64, posix threads) via Docker → dist/Piazzi-windows-x64.zip
# SDL2 built statically from source; static libstdc++/libgcc/winpthread; GUI subsystem (-mwindows).
# WINE_TEST=1 → tries to start the exe under wine (headless, xvfb) just to check that it loads.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$ROOT}"
mkdir -p "$ROOT/dist"
IMG=piazzi-mingw-build
docker build --platform linux/amd64 -t "$IMG" - <<'DOCKER'
FROM debian:bookworm
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends g++-mingw-w64-x86-64-posix gcc-mingw-w64-x86-64-posix \
    mingw-w64-tools cmake ninja-build git ca-certificates file zip && rm -rf /var/lib/apt/lists/*
RUN printf 'set(CMAKE_SYSTEM_NAME Windows)\nset(CMAKE_SYSTEM_PROCESSOR x86_64)\nset(CMAKE_C_COMPILER x86_64-w64-mingw32-gcc-posix)\nset(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++-posix)\nset(CMAKE_RC_COMPILER x86_64-w64-mingw32-windres)\nset(CMAKE_FIND_ROOT_PATH /usr/x86_64-w64-mingw32)\nset(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)\nset(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)\nset(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)\n' > /mingw.cmake
DOCKER
docker run --rm --platform linux/amd64 -v "$SRC:/src:ro" -v "$ROOT/dist:/dist" "$IMG" bash -euxc '
  cmake -S /src/cpp -B /b -G Ninja -DCMAKE_TOOLCHAIN_FILE=/mingw.cmake -DCMAKE_BUILD_TYPE=Release -DASTRO_FETCH_SDL2=ON
  cmake --build /b --target piazzi_app
  x86_64-w64-mingw32-strip /b/Piazzi.exe
  file /b/Piazzi.exe
  echo "Imported DLLs:"; x86_64-w64-mingw32-objdump -p /b/Piazzi.exe | grep "DLL Name" || true
  rm -rf /tmp/Piazzi && mkdir /tmp/Piazzi && cp /b/Piazzi.exe /tmp/Piazzi/ && cp -r /b/assets /tmp/Piazzi/
  ls /b/*.dll 2>/dev/null && cp /b/*.dll /tmp/Piazzi/ || true
  rm -f /dist/Piazzi-windows-x64.zip && (cd /tmp && zip -qr /dist/Piazzi-windows-x64.zip Piazzi)
  chown '"$(id -u):$(id -g)"' /dist/Piazzi-windows-x64.zip
'
if [ "${WINE_TEST:-0}" = 1 ]; then
  docker run --rm --platform linux/amd64 -v "$ROOT/dist:/dist:ro" debian:bookworm bash -euc '
    dpkg --add-architecture i386 >/dev/null; apt-get update -qq && apt-get install -y -qq wine wine64 xvfb xauth unzip >/dev/null
    cd /tmp && unzip -q /dist/Piazzi-windows-x64.zip && cd Piazzi
    timeout 20 xvfb-run -a wine ./Piazzi.exe; rc=$?; echo "wine exit=$rc (124 = still running after 20 s → started)"' || true
fi
ls -lh "$ROOT/dist/Piazzi-windows-x64.zip"
