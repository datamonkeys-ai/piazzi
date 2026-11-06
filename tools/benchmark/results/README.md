# Benchmark results

Each run is a pair: `<label>.json` (machine-readable, input to `run.py compare` / `summary`) and `<label>.md`
(human-readable report). All positions are scored against JPL Horizons ephemerides.

| Run | Data set | What it shows |
|---|---|---|
| `ztf-cego.*`, `ztf-cego-v7.*`, `ztf-cego-v8.*` | ZTF blind test set (`fields_test.json`) | The only runs on the blind set, each made once after freezing the code. Kept together so the history is visible, including the earlier, weaker results. |
| `v12-dev.*` | ZTF development set (`fields.json`) | Latest development run. |
| `sdss-v10.*` | SDSS DR17, station 645 (`fields_sdss.json`) | Independent instrument, 5-minute arcs. |
| `nea-v10.*` | Fast near-Earth asteroids on ZTF (`fields_nea_dev.json`) | Trailed and fast-moving objects. |
| `inject-v8-ztf.*`, `inject-v8-sdss.*` | Synthetic objects injected into real images | Completeness and O−C against an exactly known truth. |
| `v9-dev-vsmpc.*`, `sdss-v8-vsmpc.*`, `nea-v9-vsmpc.*` | Same images as above | Piazzi compared with the survey's own official astrometry at the MPC (I41 and 645). |

Intermediate development runs are not kept here. Reproduce any run with the commands in
[`../README.md`](../README.md).
