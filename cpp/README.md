# Piazzi — C++ core (native on every platform)

Pure C++17, no dependencies: builds on macOS, Windows, Linux and WebAssembly (Emscripten).
Own algorithms, based on published literature — see [`../DESIGN.md`](../DESIGN.md).

```sh
cmake -S cpp -B cpp/build -G Ninja && cmake --build cpp/build
ctest --test-dir cpp/build     # all tests (use the example in cpp/gui/assets/example)
sh cpp/tests/accuracy.sh       # positions compared against JPL Horizons (versioned Gaia DR3 catalog)
# WebAssembly (browser):
emcmake cmake -S cpp -B cpp/build-wasm -G Ninja && cmake --build cpp/build-wasm
```

Usage: `piazzi --cfg config.cfg --catalog gaia.tsv --center "07 11 43.534" "+22 24 17.78" --object "5157" --out dir images...`
(without `--catalog`, prints the URL for a Gaia DR3 box query to VizieR; accepts asu-txt and asu-tsv).
