# Piazzi vs official astrometry at the MPC — nea-v9

Δ = (Piazzi − Horizons)(t_P) − (MPC − Horizons)(t_MPC); pair = observation from the same station with |Δt| < 60 s.

**all pairs** — 9 pairs in 3 sequences (median |Δt| 15.00 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.484 / 0.509 / 0.697 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.067 / 0.068 / 0.095 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | -0.464 / 0.039 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.489 / 0.141 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.095 / 0.512 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.126 / 0.562 |
| O−C mean RA·cosδ / Dec Piazzi (") | 0.061 / 0.019 |
| O−C mean RA·cosδ / Dec MPC (") | 0.525 / -0.020 |
| Piazzi closer to Horizons | 9 of 9 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 9, -15.00 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.512 / 0.102 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.092 / 0.160 (9 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | None/None |
| MPC catalog | Gaia1 |

**clean pairs (σ_eph < 0.05" per axis)** — 9 pairs in 3 sequences (median |Δt| 15.00 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.484 / 0.509 / 0.697 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.067 / 0.068 / 0.095 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | -0.464 / 0.039 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.489 / 0.141 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.095 / 0.512 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.126 / 0.562 |
| O−C mean RA·cosδ / Dec Piazzi (") | 0.061 / 0.019 |
| O−C mean RA·cosδ / Dec MPC (") | 0.525 / -0.020 |
| Piazzi closer to Horizons | 9 of 9 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 9, -15.00 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.512 / 0.102 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.092 / 0.160 (9 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | None/None |
| MPC catalog | Gaia1 |

| sequence | target | Piazzi measurements | pairs | Δ med (") | \|O−C\| med Piazzi / MPC (") | closer | note |
|---|---|---|---|---|---|---|---|
| ztf_20210905_f296_c6_q3_gr | Poseidon | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20200124_f551_c9_q3_gr | Geographos | 4 | 4 | 0.346 | 0.102 / 0.397 | Piazzi | |
| ztf_20201017_f692_c6_q2_r | Ryugu | 5 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20191102_f471_c10_q4_gr | Mithra | 3 | 2 | 0.484 | 0.036 / 0.502 | Piazzi | |
| ztf_20200718_f622_c1_q4_r | Atira | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20191029_f472_c16_q1_r | Aten | 4 | 3 | 0.694 | 0.110 / 0.748 | Piazzi | |
| ztf_20210110_f396_c13_q2_gir | Toutatis | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20180519_f790_c13_q2_r | 2018 JA | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20210816_f646_c4_q4_gir | 2021 CH2 | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20180604_f630_c6_q1_gr | 2018 EJ4 | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20200823_f740_c6_q1_gr | 2019 AN5 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20181031_f784_c6_q2_r | 2018 VW3 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20201029_f406_c2_q2_gr | 2020 UN6 | 4 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
