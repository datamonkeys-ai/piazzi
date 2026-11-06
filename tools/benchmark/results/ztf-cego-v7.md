# Benchmark Piazzi — ztf-cego-v7

piazzi: `cpp/build/piazzi` (eced3c4) · 2026-11-05 17:38 · fields.json: 23 sequences

| metric | value |
|---|---|
| sequences (run / data failures) | 23 (21 / 2) |
| plate not solved (some image < 6 refs) | 0 |
| target detected | 15 (71%) |
| target measurements | 52 |
| |O−C| median (") | 0.106 |
| |O−C| rms (") | 0.182 |
| |O−C| p95 (") | 0.317 |
| O−C mean RA·cosδ / Dec (") | 0.024 / -0.026 |
| O−C rms RA·cosδ / Dec (") | 0.145 / 0.111 |
| median ADES sigma (", RA/Dec in quadrature) | 0.120 |
| median |O−C|/sigma ratio | 1.030 |
| median ephemeris sigma (Horizons, 1σ, quadrature) (") | 0.007 (52 measurements) |
| median corrected ratio |O−C|/√(σ²+σ_eph²) | 0.884 |
| clean measurements (σ_eph < 0.05" per axis) | 29 in 9 sequences |
|   |O−C| median / rms / p95 (") | 0.087 / 0.116 / 0.210 |
|   O−C mean RA·cosδ / Dec (") | 0.014 / -0.030 |
|   O−C rms RA·cosδ / Dec (") | 0.092 / 0.070 |
|   median |O−C|/σ ratio (0.83 calibrated) / corrected | 0.995 / 0.995 |
|   median per-axis |O−C|/σ RA / Dec (0.674 calibrated) | 0.736 / 0.600 |
| internal consistency (residual − sequence mean) median |O−C|/σ / per axis RA / Dec | 0.800 / 0.726 / 0.523 |
| tracks with no known object (false positives) per sequence | 0.000 |
|   same, associating only target/known (old metric) | 0.000 |
| tracks of other SkyBoT objects in the field | 0 |
| total time (s) / piazzi only (s) | 102.2 / 98.8 |

Detection by V bin:

| V | n | detected | rate |
|---|---|---|---|
| 16–17 | 2 | 1 | 50% |
| 17–18 | 1 | 0 | 0% |
| 18–19 | 7 | 6 | 86% |
| 19–20 | 6 | 4 | 67% |
| 20–21 | 4 | 3 | 75% |
| 21–22 | 1 | 1 | 100% |

| sequence | target | V | "/min | refs | det | |O−C| med (") | σ_eph (") | FP | status |
|---|---|---|---|---|---|---|---|---|---|
| ztf_20200405_f756_c15_q1_r | 2020 FW4 | 18.6 | 15.76 | 97,96,100 | yes | 0.13 | 0.135 | 0 | ok |
| ztf_20190325_f525_c12_q2_gir | 2019 FT | 21.0 | 13.43 | 141,162,170,177 | yes | 0.32 | 0.198 | 0 | ok |
| ztf_20230125_f465_c11_q3_gir | 2011 AG5 | 18.0 | 2.97 | 278,283,257 | no | — | — | 0 | ok |
| ztf_20230813_f544_c16_q2_gr | 2023 PC1 | 19.4 | 7.45 | 566,570,522 | no | — | — | 0 | ok |
| ztf_20190727_f383_c10_q2_r | 2019 NX6 | 18.0 | 9.83 | 2790,2795,2775,2751,2744 | yes | 0.10 | 0.191 | 0 | ok |
| ztf_20190428_f715_c16_q4_gr | 2019 FP2 | 18.9 | 6.19 | 54,55,62,62 | yes | 0.06 | 0.075 | 0 | ok |
| ztf_20200430_f821_c11_q4_gi | 2020 FA7 | 19.1 | 9.74 | 97,115,116 | yes | 0.28 | 0.166 | 0 | ok |
| ztf_20200123_f751_c6_q2_r | 2020 BO1 | 16.2 | 20.09 | 216,219,214 | no | — | — | 0 | ok |
| ztf_20201029_f404_c13_q2_gr | 2020 TY1 | 18.7 | 6.72 | 51,45,81,86 | yes | 0.16 | 0.092 | 0 | ok |
| ztf_20200725_f753_c8_q3_gir | 1990 MF | 17.1 | 11.41 | 68,60,87 | no | — | — | 0 | ok |
| ztf_20191030_f607_c2_q1_r | 1998 WF23 | 18.0 | 0.44 | 89,87,92,86,86 | yes | 0.10 | 0.000 | 0 | ok |
| ztf_20210922_f607_c1_q3_gir | 1999 JC47 | 18.0 | 0.16 | 114,175,159 | yes | 0.05 | 0.001 | 0 | ok |
| ztf_20191114_f447_c14_q3_gr | 1999 VP143 | 18.6 | 0.12 | — | no | — | — | 0 | data_error (recorte pequeno (145x712, borda do quadrante)) |
| ztf_20230615_f373_c3_q1_gir | 1999 JK58 | 18.8 | 0.21 | — | no | — | — | 0 | data_error (recorte pequeno (74x352, borda do quadrante)) |
| ztf_20200207_f447_c14_q3_r | 2000 HM29 | 19.1 | 1.10 | 72,70,74,72 | yes | 0.09 | 0.000 | 0 | ok |
| ztf_20210109_f499_c2_q3_gir | 3503 P-L | 19.6 | 0.83 | 52,44,58 | yes | 0.11 | 0.001 | 0 | ok |
| ztf_20220416_f519_c2_q4_gir | 2001 VV32 | 19.7 | 0.29 | 53,54,33 | no | — | — | 0 | ok |
| ztf_20220823_f556_c4_q3_gir | 2000 GP71 | 19.9 | 0.26 | 75,86,86 | yes | 0.07 | 0.001 | 0 | ok |
| ztf_20220409_f519_c2_q4_gir | 2005 ET101 | 20.0 | 0.14 | 27,32,32 | yes | 0.15 | 0.004 | 0 | ok |
| ztf_20220723_f496_c1_q4_gir | 1998 SC171 | 20.0 | 0.08 | 74,74,76 | yes | 0.08 | 0.009 | 0 | ok |
| ztf_20220823_f607_c1_q3_gir | 2001 SP75 | 20.0 | 0.47 | 140,173,173 | yes | 0.09 | 0.003 | 0 | ok |
| ztf_20230427_f373_c3_q1_gir | 1999 VR157 | 20.0 | 0.60 | 43,34,40 | no | — | — | 0 | ok |
| ztf_20211030_f784_c1_q2_gr | Schorria | 16.3 | 1.26 | 79,63,62 | yes | 0.06 | 0.001 | 0 | ok |
