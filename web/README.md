# Piazzi web

The same C++ core (WebAssembly) as the iPhone app and the CLI, with the "pocket observatory" interface.
The reduction runs entirely in the browser; network access only for Gaia DR3, SkyBoT and "What is this?" (via a local proxy, because of CORS).

```sh
./web/build.sh          # builds cpp/ -> web/astro.js + astro.wasm (Emscripten)
python3 web/serve.py    # http://localhost:8777 and on the local network (same Wi-Fi)
```

## Example (offline)

"Open example" on the start screen (or `http://localhost:8000/?example=1`) loads `example/` — three ZTF crops of
asteroid (5157) Hindemith and a Gaia DR3 TSV (credits in `example/README.md`) — and runs the whole reduction without
any network access (SkyBoT is skipped). The session is marked as example data: the results show "Example data — not
your observation"; never submit these reports to the MPC.

## Language and settings

The interface is in English by default, with Portuguese in Settings → Language (saved in `localStorage`); all strings
live in `strings.js`. The parameter list is generated from the core's `astro_settings_schema()` (grouped, "Advanced"
collapsed, numbers validated against min/max, ↺ per field and "Restore defaults").
