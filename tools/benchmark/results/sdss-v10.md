# Benchmark Piazzi — sdss-v10

piazzi: `cpp/build/piazzi` (4146b8d) · 2026-11-05 21:00 · fields.json: 32 sequences

| metric | value |
|---|---|
| sequences (run / data failures) | 32 (32 / 0) |
| plate not solved (some image < 6 refs) | 0 |
| target detected | 29 (91%) |
| target measurements | 116 |
| |O−C| median (") | 0.030 |
| |O−C| rms (") | 0.075 |
| |O−C| p95 (") | 0.163 |
| O−C mean RA·cosδ / Dec (") | 0.025 / -0.011 |
| O−C rms RA·cosδ / Dec (") | 0.061 / 0.043 |
| median ADES sigma (", RA/Dec in quadrature) | 0.023 |
| median |O−C|/sigma ratio | 1.183 |
| median ephemeris sigma (Horizons, 1σ, quadrature) (") | 0.005 (116 measurements) |
| median corrected ratio |O−C|/√(σ²+σ_eph²) | 0.891 |
| clean measurements (σ_eph < 0.05" per axis) | 92 in 23 sequences |
|   |O−C| median / rms / p95 (") | 0.023 / 0.039 / 0.091 |
|   O−C mean RA·cosδ / Dec (") | 0.006 / 0.001 |
|   O−C rms RA·cosδ / Dec (") | 0.031 / 0.023 |
|   median |O−C|/σ ratio (0.83 calibrated) / corrected | 1.034 / 0.925 |
|   median per-axis |O−C|/σ RA / Dec (0.674 calibrated) | 1.031 / 0.751 |
| internal consistency (residual − sequence mean) median |O−C|/σ / per axis RA / Dec | 0.824 / 0.720 / 0.693 |
| tracks with no known object (false positives) per sequence | 0.031 |
|   same, associating only target/known (old metric) | 0.250 |
| tracks of other SkyBoT objects in the field | 7 |
| total time (s) / piazzi only (s) | 394.5 / 391.9 |

Detection by V bin:

| V | n | detected | rate |
|---|---|---|---|
| 15–16 | 2 | 2 | 100% |
| 16–17 | 6 | 6 | 100% |
| 17–18 | 6 | 6 | 100% |
| 18–19 | 6 | 6 | 100% |
| 19–20 | 6 | 5 | 83% |
| 20–21 | 6 | 4 | 67% |

| sequence | target | V | "/min | refs | det | |O−C| med (") | σ_eph (") | FP | status |
|---|---|---|---|---|---|---|---|---|---|
| sdss_002583_4_0207 | 2001 QP112 | 19.6 | 0.56 | 143,141,147,143 | yes | 0.06 | 0.076 | 0 | ok |
| sdss_004927_5_0327 | 2000 VX33 | 19.4 | 0.17 | 99,98,96,89 | no | — | — | 0 | ok |
| sdss_004198_3_0287 | 2006 EK66 | 20.2 | 0.24 | 90,92,99,82 | no | — | — | 0 | ok |
| sdss_006794_6_0207 | 2001 FB171 | 19.3 | 0.32 | 624,631,669,598 | yes | 0.02 | 0.005 | 0 | ok |
| sdss_004198_4_0407 | 2000 CZ68 | 20.9 | 0.42 | 81,81,81,82 | no | — | — | 0 | ok |
| sdss_003360_1_0407 | 1990 SH5 | 18.1 | 0.64 | 86,81,84,88 | yes | 0.03 | 0.018 | 0 | ok |
| sdss_004145_2_0087 | 2000 WG | 20.3 | 0.64 | 71,72,81,70 | yes | 0.03 | 0.046 | 0 | ok |
| sdss_006004_3_0207 | 2006 BP162 | 20.9 | 0.60 | 74,74,79,74 | yes | 0.17 | 0.088 | 0 | ok |
| sdss_000756_1_0327 | 2003 HE38 | 19.2 | 0.59 | 94,90,92,97 | yes | 0.04 | 0.065 | 0 | ok |
| sdss_000745_1_0047 | 2012 DM98 | 19.7 | 0.62 | 116,116,114,114 | yes | 0.09 | 0.096 | 0 | ok |
| sdss_003644_3_0047 | 2005 TQ171 | 20.0 | 0.61 | 235,234,245,214 | yes | 0.07 | 0.071 | 0 | ok |
| sdss_003644_6_0087 | 1989 SJ1 | 19.6 | 0.57 | 185,201,202,186 | yes | 0.05 | 0.005 | 0 | ok |
| sdss_006004_1_0167 | Flavia | 17.7 | 0.61 | 102,103,110,97 | yes | 0.03 | 0.001 | 0 | ok |
| sdss_001140_6_0127 | 2014 QQ366 | 20.8 | 0.42 | 80,84,89,79 | yes | 0.27 | 0.053 | 0 | ok |
| sdss_002247_2_0207 | Lucretius Carus | 16.4 | 0.56 | 91,92,94,90 | yes | 0.03 | 0.003 | 0 | ok |
| sdss_008095_2_0207 | Ringleb | 17.9 | 0.28 | 112,120,118,115 | yes | 0.02 | 0.001 | 0 | ok |
| sdss_000094_1_0167 | Marciaprill | 18.0 | 0.62 | 77,75,80,76 | yes | 0.05 | 0.025 | 1 | ok |
| sdss_007913_1_0127 | 1998 XF93 | 18.0 | 0.52 | 95,91,103,102 | yes | 0.02 | 0.001 | 0 | ok |
| sdss_002125_1_0287 | 2000 AU57 | 18.4 | 0.33 | 81,82,83,77 | yes | 0.02 | 0.020 | 0 | ok |
| sdss_005052_4_0087 | Damir | 18.3 | 0.90 | 69,68,73,68 | yes | 0.02 | 0.002 | 0 | ok |
| sdss_001755_1_0367 | 2000 QF84 | 18.3 | 0.59 | 84,79,85,85 | yes | 0.03 | 0.010 | 0 | ok |
| sdss_004128_6_0127 | 1999 XL130 | 17.3 | 0.56 | 70,73,74,71 | yes | 0.02 | 0.025 | 0 | ok |
| sdss_006121_4_0207 | Gotha | 15.9 | 0.38 | 88,81,86,80 | yes | 0.01 | 0.002 | 0 | ok |
| sdss_000094_6_0287 | Jillpernell | 17.4 | 0.42 | 67,68,71,67 | yes | 0.11 | 0.018 | 0 | ok |
| sdss_002986_4_0167 | Saale | 17.4 | 0.60 | 87,97,102,91 | yes | 0.02 | 0.003 | 0 | ok |
| sdss_003635_4_0047 | Street | 17.9 | 0.57 | 242,249,263,245 | yes | 0.01 | 0.003 | 0 | ok |
| sdss_003354_5_0047 | Moore | 16.2 | 0.64 | 78,79,83,78 | yes | 0.04 | 0.002 | 0 | ok |
| sdss_006354_6_0287 | Boynton | 16.6 | 0.52 | 76,78,85,81 | yes | 0.01 | 0.011 | 0 | ok |
| sdss_005935_5_0247 | Bernardus | 16.6 | 0.75 | 69,70,70,66 | yes | 0.02 | 0.002 | 0 | ok |
| sdss_000125_2_0087 | Dylanthomas | 16.9 | 0.65 | 66,64,67,65 | yes | 0.05 | 0.004 | 0 | ok |
| sdss_004203_6_0447 | Vitja | 15.6 | 0.43 | 79,81,70,82 | yes | 0.01 | 0.003 | 0 | ok |
| sdss_004263_5_0287 | Onizaki | 16.0 | 0.42 | 67,65,65,68 | yes | 0.02 | 0.002 | 0 | ok |
