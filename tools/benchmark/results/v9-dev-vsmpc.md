# Piazzi vs official astrometry at the MPC — v9-dev

Δ = (Piazzi − Horizons)(t_P) − (MPC − Horizons)(t_MPC); pair = observation from the same station with |Δt| < 60 s.

**all pairs** — 7 pairs in 4 sequences (median |Δt| 14.90 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.069 / 0.096 / 0.151 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.054 / 0.058 / 0.084 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | 0.054 / -0.010 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.081 / 0.051 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.020 / 0.064 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.065 / 0.121 |
| O−C mean RA·cosδ / Dec Piazzi (") | -0.008 / 0.021 |
| O−C mean RA·cosδ / Dec MPC (") | -0.063 / 0.032 |
| Piazzi closer to Horizons | 7 of 7 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 4, -15.20 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.138 / 0.071 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.051 / 0.096 (7 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | None/None |
| MPC catalog | Gaia1 |

**clean pairs (σ_eph < 0.05" per axis)** — 7 pairs in 4 sequences (median |Δt| 14.90 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.069 / 0.096 / 0.151 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.054 / 0.058 / 0.084 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | 0.054 / -0.010 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.081 / 0.051 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.020 / 0.064 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.065 / 0.121 |
| O−C mean RA·cosδ / Dec Piazzi (") | -0.008 / 0.021 |
| O−C mean RA·cosδ / Dec MPC (") | -0.063 / 0.032 |
| Piazzi closer to Horizons | 7 of 7 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 4, -15.20 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.138 / 0.071 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.051 / 0.096 (7 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | None/None |
| MPC catalog | Gaia1 |

| sequence | target | Piazzi measurements | pairs | Δ med (") | \|O−C\| med Piazzi / MPC (") | closer | note |
|---|---|---|---|---|---|---|---|
| ztf_20210206_f565_c15_q2_gir | Hindemith | 3 | 3 | 0.054 | 0.017 / 0.060 | Piazzi | |
| ztf_20210912_f338_c4_q2_gir | 1999 JH79 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20230907_f340_c16_q2_gir | Kudo | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20200709_f333_c1_q1_gr | 2000 YY39 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20231101_f449_c10_q4_gir | 2000 QF79 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20191023_f451_c8_q1_r | 1999 GL16 | 5 | 1 | 0.148 | 0.019 / 0.168 | Piazzi | |
| ztf_20191101_f1553_c10_q4_gr | 1996 OA | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20211113_f517_c16_q4_gir | 2000 RA10 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20220907_f606_c3_q3_gir | 2001 SX8 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20201216_f565_c15_q2_gir | 1994 AL10 | 3 | 2 | 0.104 | 0.091 / 0.133 | Piazzi | |
| ztf_20191029_f504_c7_q2_gr | 2002 SP24 | 5 | 0 | — | — | — | no I41 observation at the MPC within 60 s (3 from the station on the same night ±1 h) |
| ztf_20191029_f504_c7_q2_r | 2002 SP24 | 5 | 1 | 0.082 | 0.051 / 0.108 | Piazzi | |
| ztf_20221119_f554_c1_q3_gir | 2004 TV295 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20191124_f558_c2_q3_gr | 2005 GV102 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20210331_f420_c13_q4_gr | 2008 UC72 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20220830_f556_c15_q1_gir | 2001 QO198 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20221119_f468_c6_q1_gir | 2000 SB142 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20200111_f556_c6_q1_gr | 2006 WX108 | 0 | 0 | — | — | — | target not measured by Piazzi |
| ztf_20221025_f500_c10_q1_gir | 2008 VR67 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20230417_f323_c9_q3_gir | 1999 TC155 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| ztf_20230615_f323_c9_q3_gir | 2000 AD40 | 3 | 0 | — | — | — | no I41 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
