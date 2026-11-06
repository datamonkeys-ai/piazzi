#!/bin/sh
# Compiles the C++ core (same as iPhone/CLI) to WebAssembly.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
export EMSDK_PYTHON=${EMSDK_PYTHON:-/opt/homebrew/bin/python3.15}
FUNCS=$(grep -o 'astro_[a-z_]*(' "$ROOT/cpp/include/piazzi/astro_c.h" | tr -d '(' | sort -u | sed 's/^/_/' | paste -sd, -)
em++ -O2 -fexceptions -std=c++17 -I"$ROOT/cpp/include" \
  $(ls "$ROOT"/cpp/src/*.cpp | grep -v '/main_') \
  -o "$ROOT/web/astro.js" --no-entry -sMODULARIZE=1 -sEXPORT_NAME=AstroModule -sENVIRONMENT=web \
  -fexceptions -sALLOW_MEMORY_GROWTH=1 -sEXPORTED_FUNCTIONS="$FUNCS,_malloc,_free" \
  -sEXPORTED_RUNTIME_METHODS=FS,ccall,cwrap,UTF8ToString,HEAPU8,HEAPF64,getValue
echo "ok: web/astro.js + web/astro.wasm"
