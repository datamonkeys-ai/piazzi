# Piazzi vs official astrometry at the MPC — sdss-v8

Δ = (Piazzi − Horizons)(t_P) − (MPC − Horizons)(t_MPC); pair = observation from the same station with |Δt| < 60 s.

**all pairs** — 44 pairs in 22 sequences (median |Δt| 6.88 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.054 / 0.098 / 0.194 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.121 / 0.138 / 0.257 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | -0.003 / -0.015 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.085 / 0.049 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.033 / 0.083 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.072 / 0.116 |
| O−C mean RA·cosδ / Dec Piazzi (") | 0.020 / -0.010 |
| O−C mean RA·cosδ / Dec MPC (") | 0.023 / 0.005 |
| Piazzi closer to Horizons | 36 of 44 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 24, 13.21 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.087 / 0.148 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.024 / — (0 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | 0.001/0.01 |
| MPC catalog | UCAC1 |

**clean pairs (σ_eph < 0.05" per axis)** — 34 pairs in 17 sequences (median |Δt| 6.88 s)

| metric | value |
|---|---|
| \|Piazzi − MPC\| median / rms / p95 (") | 0.064 / 0.108 / 0.199 |
| \|Piazzi − MPC\| positions only (both at Piazzi's time) median / rms / p95 (") | 0.129 / 0.144 / 0.263 |
| mean bias Piazzi − MPC RA·cosδ / Dec (") | -0.011 / -0.019 |
| rms Piazzi − MPC RA·cosδ / Dec (") | 0.095 / 0.053 |
| \|O−C\| Horizons median Piazzi / MPC (") | 0.028 / 0.085 |
| \|O−C\| Horizons rms Piazzi / MPC (") | 0.039 / 0.115 |
| O−C mean RA·cosδ / Dec Piazzi (") | 0.002 / 0.001 |
| O−C mean RA·cosδ / Dec MPC (") | 0.012 / 0.020 |
| Piazzi closer to Horizons | 32 of 34 |
| pairs with \|Δt\| > 1 s: n, median Δt (s, MPC − Piazzi) | 19, 13.20 |
| &nbsp;&nbsp;\|O−C\| median MPC at MPC time / at Piazzi time (") | 0.090 / 0.163 |
| median ADES σ Piazzi / MPC (", quadrature) | 0.022 / — (0 with rms at the MPC) |
| MPC declared precision (precRA s / precDec ") | 0.001/0.01 |
| MPC catalog | UCAC1 |

| sequence | target | Piazzi measurements | pairs | Δ med (") | \|O−C\| med Piazzi / MPC (") | closer | note |
|---|---|---|---|---|---|---|---|
| sdss_002583_4_0207 | 2001 QP112 | 4 | 2 | 0.062 | 0.061 / 0.003 | MPC | |
| sdss_004927_5_0327 | 2000 VX33 | 0 | 0 | — | — | — | target not measured by Piazzi |
| sdss_004198_3_0287 | 2006 EK66 | 0 | 0 | — | — | — | target not measured by Piazzi |
| sdss_006794_6_0207 | 2001 FB171 | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_004198_4_0407 | 2000 CZ68 | 0 | 0 | — | — | — | target not measured by Piazzi |
| sdss_003360_1_0407 | 1990 SH5 | 4 | 2 | 0.056 | 0.032 / 0.070 | Piazzi | |
| sdss_004145_2_0087 | 2000 WG | 4 | 2 | 0.200 | 0.029 / 0.179 | Piazzi | |
| sdss_006004_3_0207 | 2006 BP162 | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_000756_1_0327 | 2003 HE38 | 4 | 2 | 0.038 | 0.030 / 0.061 | Piazzi | |
| sdss_000745_1_0047 | 2012 DM98 | 4 | 2 | 0.045 | 0.086 / 0.106 | Piazzi | |
| sdss_003644_3_0047 | 2005 TQ171 | 4 | 2 | 0.060 | 0.060 / 0.032 | MPC | |
| sdss_003644_6_0087 | 1989 SJ1 | 4 | 2 | 0.024 | 0.022 / 0.021 | MPC | |
| sdss_006004_1_0167 | Flavia | 4 | 2 | 0.051 | 0.026 / 0.071 | Piazzi | |
| sdss_001140_6_0127 | 2014 QQ366 | 4 | 2 | 0.036 | 0.265 / 0.236 | MPC | |
| sdss_002247_2_0207 | Lucretius Carus | 4 | 2 | 0.038 | 0.018 / 0.041 | Piazzi | |
| sdss_008095_2_0207 | Ringleb | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_000094_1_0167 | Marciaprill | 4 | 2 | 0.052 | 0.045 / 0.091 | Piazzi | |
| sdss_007913_1_0127 | 1998 XF93 | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_002125_1_0287 | 2000 AU57 | 4 | 2 | 0.031 | 0.028 / 0.046 | Piazzi | |
| sdss_005052_4_0087 | Damir | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_001755_1_0367 | 2000 QF84 | 4 | 2 | 0.123 | 0.027 / 0.101 | Piazzi | |
| sdss_004128_6_0127 | 1999 XL130 | 4 | 2 | 0.139 | 0.037 / 0.150 | Piazzi | |
| sdss_006121_4_0207 | Gotha | 4 | 2 | 0.171 | 0.014 / 0.170 | Piazzi | |
| sdss_000094_6_0287 | Jillpernell | 4 | 2 | 0.100 | 0.106 / 0.135 | Piazzi | |
| sdss_002986_4_0167 | Saale | 4 | 2 | 0.019 | 0.038 / 0.043 | Piazzi | |
| sdss_003635_4_0047 | Street | 4 | 2 | 0.032 | 0.013 / 0.037 | Piazzi | |
| sdss_003354_5_0047 | Moore | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_006354_6_0287 | Boynton | 4 | 2 | 0.020 | 0.009 / 0.026 | Piazzi | |
| sdss_005935_5_0247 | Bernardus | 4 | 2 | 0.109 | 0.024 / 0.104 | Piazzi | |
| sdss_000125_2_0087 | Dylanthomas | 4 | 2 | 0.214 | 0.043 / 0.246 | Piazzi | |
| sdss_004203_6_0447 | Vitja | 4 | 0 | — | — | — | no 645 observation at the MPC within 60 s (0 from the station on the same night ±1 h) |
| sdss_004263_5_0287 | Onizaki | 4 | 2 | 0.090 | 0.014 / 0.081 | Piazzi | |
