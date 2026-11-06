# Benchmark Piazzi — v12-dev

piazzi: `cpp/build/piazzi` (edb21af) · 2026-11-06 00:42 · fields.json: 21 sequences

| metric | value |
|---|---|
| sequences (run / data failures) | 21 (20 / 1) |
| plate not solved (some image < 6 refs) | 0 |
| target detected | 15 (75%) |
| target measurements | 51 |
| |O−C| median (") | 0.059 |
| |O−C| rms (") | 0.140 |
| |O−C| p95 (") | 0.307 |
| O−C mean RA·cosδ / Dec (") | 0.005 / -0.008 |
| O−C rms RA·cosδ / Dec (") | 0.100 / 0.097 |
| median ADES sigma (", RA/Dec in quadrature) | 0.071 |
| median |O−C|/sigma ratio | 0.752 |
| median ephemeris sigma (Horizons, 1σ, quadrature) (") | 0.001 (51 measurements) |
| median corrected ratio |O−C|/√(σ²+σ_eph²) | 0.744 |
| clean measurements (σ_eph < 0.05" per axis) | 51 in 15 sequences |
|   |O−C| median / rms / p95 (") | 0.059 / 0.140 / 0.307 |
|   O−C mean RA·cosδ / Dec (") | 0.005 / -0.008 |
|   O−C rms RA·cosδ / Dec (") | 0.100 / 0.097 |
|   median |O−C|/σ ratio (0.83 calibrated) / corrected | 0.752 / 0.744 |
|   median per-axis |O−C|/σ RA / Dec (0.674 calibrated) | 0.538 / 0.627 |
| internal consistency (residual − sequence mean) median |O−C|/σ / per axis RA / Dec | 0.822 / 0.594 / 0.920 |
| tracks with no known object (false positives) per sequence | 0.050 |
|   same, associating only target/known (old metric) | 0.050 |
| tracks of other SkyBoT objects in the field | 0 |
| total time (s) / piazzi only (s) | 50.0 / 49.5 |

Detection by V bin:

| V | n | detected | rate |
|---|---|---|---|
| 17–18 | 3 | 3 | 100% |
| 18–19 | 6 | 4 | 67% |
| 19–20 | 7 | 5 | 71% |
| 20–21 | 4 | 3 | 75% |

| sequence | target | V | "/min | refs | det | |O−C| med (") | σ_eph (") | FP | status |
|---|---|---|---|---|---|---|---|---|---|
| ztf_20210206_f565_c15_q2_gir | Hindemith | 17.6 | 0.34 | 249,255,233 | yes | 0.02 | 0.001 | 0 | ok |
| ztf_20210912_f338_c4_q2_gir | 1999 JH79 | 17.7 | 0.35 | 198,223,189 | yes | 0.06 | 0.002 | 0 | ok |
| ztf_20230907_f340_c16_q2_gir | Kudo | 17.9 | 0.44 | 147,118,144 | yes | 0.07 | 0.001 | 0 | ok |
| ztf_20200709_f333_c1_q1_gr | 2000 YY39 | 18.0 | 0.60 | 3775,3876,2329 | no | — | — | 0 | ok |
| ztf_20231101_f449_c10_q4_gir | 2000 QF79 | 18.4 | 0.48 | 71,79,54 | yes | 0.12 | 0.012 | 0 | ok |
| ztf_20191023_f451_c8_q1_r | 1999 GL16 | 18.6 | 0.51 | 57,55,56,58,54 | yes | 0.05 | 0.001 | 0 | ok |
| ztf_20191101_f1553_c10_q4_gr | 1996 OA | 18.6 | 0.44 | — | no | — | — | 0 | data_error (recorte pequeno (116x712, borda do quadrante)) |
| ztf_20211113_f517_c16_q4_gir | 2000 RA10 | 18.7 | 0.59 | 160,148,164 | yes | 0.04 | 0.002 | 0 | ok |
| ztf_20220907_f606_c3_q3_gir | 2001 SX8 | 18.7 | 0.66 | 49,47,44 | no | — | — | 0 | ok |
| ztf_20201216_f565_c15_q2_gir | 1994 AL10 | 18.9 | 0.54 | 332,360,378 | yes | 0.06 | 0.001 | 0 | ok |
| ztf_20191029_f504_c7_q2_gr | 2002 SP24 | 19.1 | 0.51 | 58,74,71,70,71 | yes | 0.04 | 0.001 | 0 | ok |
| ztf_20191029_f504_c7_q2_r | 2002 SP24 | 19.1 | 0.51 | 73,74,73,74,70 | yes | 0.02 | 0.001 | 0 | ok |
| ztf_20221119_f554_c1_q3_gir | 2004 TV295 | 19.5 | 0.43 | 75,71,58 | yes | 0.20 | 0.036 | 0 | ok |
| ztf_20191124_f558_c2_q3_gr | 2005 GV102 | 19.6 | 0.53 | 46,44,62 | yes | 0.09 | 0.032 | 0 | ok |
| ztf_20210331_f420_c13_q4_gr | 2008 UC72 | 19.6 | 0.61 | 54,54,67,68 | no | — | — | 1 | ok |
| ztf_20220830_f556_c15_q1_gir | 2001 QO198 | 19.6 | 0.52 | 66,58,66 | yes | 0.07 | 0.001 | 0 | ok |
| ztf_20221119_f468_c6_q1_gir | 2000 SB142 | 19.7 | 0.57 | 29,32,31 | no | — | — | 0 | ok |
| ztf_20200111_f556_c6_q1_gr | 2006 WX108 | 20.0 | 0.55 | 69,88,95 | no | — | — | 0 | ok |
| ztf_20221025_f500_c10_q1_gir | 2008 VR67 | 20.0 | 0.45 | 65,56,63 | yes | 0.19 | 0.044 | 0 | ok |
| ztf_20230417_f323_c9_q3_gir | 1999 TC155 | 20.0 | 0.66 | 153,156,136 | yes | 0.13 | 0.039 | 0 | ok |
| ztf_20230615_f323_c9_q3_gir | 2000 AD40 | 20.0 | 0.11 | 150,149,131 | yes | 0.20 | 0.003 | 0 | ok |
