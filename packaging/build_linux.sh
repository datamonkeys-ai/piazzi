#!/usr/bin/env bash
# Linux build via Docker (Ubuntu 22.04 → glibc 2.35) and AppImage packaging (fallback: .tar.gz).
# Usage: packaging/build_linux.sh [amd64|arm64]   (default amd64)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$ROOT}"
ARCH="${1:-amd64}"
case "$ARCH" in amd64) AARCH=x86_64;; arm64) AARCH=aarch64;; *) echo "invalid arch"; exit 1;; esac
mkdir -p "$ROOT/dist"
IMG="piazzi-linux-build:$ARCH"
docker build --platform "linux/$ARCH" -t "$IMG" - <<'DOCKER'
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake ninja-build git ca-certificates curl file \
    libgl1-mesa-dev libx11-dev libxext-dev libxrandr-dev libxcursor-dev libxi-dev libxinerama-dev libxss-dev libxxf86vm-dev \
    libwayland-dev libxkbcommon-dev wayland-protocols libegl1-mesa-dev libdbus-1-dev libudev-dev libasound2-dev libpulse-dev \
    libfuse2 squashfs-tools && rm -rf /var/lib/apt/lists/*
DOCKER
docker run --rm --platform "linux/$ARCH" -e AARCH="$AARCH" -v "$SRC:/src:ro" -v "$ROOT/dist:/dist" "$IMG" bash -euxc '
  cmake -S /src/cpp -B /b -G Ninja -DCMAKE_BUILD_TYPE=Release -DASTRO_FETCH_SDL2=ON -DOpenGL_GL_PREFERENCE=LEGACY \
        -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc"
  cmake --build /b --target piazzi_app
  strip /b/piazzi_app; file /b/piazzi_app; ldd /b/piazzi_app || true
  # AppDir  (the dd above zeroes the AppImage "magic" in bytes 8-10: without it binfmt/Rosetta/qemu refuse to run it)
  A=/AppDir; mkdir -p $A/usr/bin $A/usr/share/applications $A/usr/share/icons/hicolor/256x256/apps
  cp /b/piazzi_app $A/usr/bin/piazzi; cp -r /b/assets $A/usr/bin/assets
  printf "[Desktop Entry]\nType=Application\nName=Piazzi\nExec=piazzi\nIcon=piazzi\nCategories=Science;Astronomy;\n" > $A/piazzi.desktop
  cp $A/piazzi.desktop $A/usr/share/applications/
  ICON=$(ls /src/packaging/icon.png 2>/dev/null || true)
  if [ -n "$ICON" ]; then cp "$ICON" $A/piazzi.png; else
    printf "\x89PNG\r\n\x1a\n\0\0\0\rIHDR\0\0\0\x01\0\0\0\x01\x08\x06\0\0\0\x1f\x15\xc4\x89\0\0\0\rIDATx\x9cc\xf8\x0f\0\0\x01\x01\0\x05\x18\xd8N\0\0\0\0IEND\xaeB\x60\x82" > $A/piazzi.png; fi
  cp $A/piazzi.png $A/usr/share/icons/hicolor/256x256/apps/
  printf "#!/bin/sh\nHERE=\"\$(dirname \"\$(readlink -f \"\$0\")\")\"\ncd \"\$HERE/usr/bin\" && exec ./piazzi \"\$@\"\n" > $A/AppRun; chmod +x $A/AppRun
  OUT=/dist/Piazzi-linux-$AARCH
  rm -f $OUT.AppImage $OUT.tar.gz
  if curl -fsSL -o /tmp/appimagetool https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-$AARCH.AppImage \
     && chmod +x /tmp/appimagetool \
     && dd if=/dev/zero of=/tmp/appimagetool bs=1 count=3 seek=8 conv=notrunc 2>/dev/null \
     && (cd /tmp && ./appimagetool --appimage-extract >/dev/null) \
     && ARCH=$AARCH /tmp/squashfs-root/AppRun $A $OUT.AppImage; then
    echo "AppImage OK"
  else
    echo "AppImage failed; producing tar.gz"
    T=/tmp/Piazzi; mkdir -p $T; cp -r $A/usr/bin/* $T/
    printf "#!/bin/sh\ncd \"\$(dirname \"\$0\")\" && exec ./piazzi \"\$@\"\n" > $T/run.sh; chmod +x $T/run.sh
    tar -C /tmp -czf $OUT.tar.gz Piazzi
  fi
  chown -R '"$(id -u):$(id -g)"' /dist
'
ls -lh "$ROOT/dist"/Piazzi-linux-$AARCH.*
