# Benchmark Piazzi — nea-v10

piazzi: `cpp/build/piazzi` (4146b8d) · 2026-11-05 20:54 · fields.json: 13 sequences

| metric | value |
|---|---|
| sequences (run / data failures) | 13 (13 / 0) |
| plate not solved (some image < 6 refs) | 0 |
| target detected | 11 (85%) |
| target measurements | 43 |
| |O−C| median (") | 0.098 |
| |O−C| rms (") | 0.340 |
| |O−C| p95 (") | 0.331 |
| O−C mean RA·cosδ / Dec (") | -0.007 / 0.010 |
| O−C rms RA·cosδ / Dec (") | 0.140 / 0.310 |
| median ADES sigma (", RA/Dec in quadrature) | 0.076 |
| median |O−C|/sigma ratio | 1.072 |
| median ephemeris sigma (Horizons, 1σ, quadrature) (") | 0.020 (43 measurements) |
| median corrected ratio |O−C|/√(σ²+σ_eph²) | 0.899 |
| clean measurements (σ_eph < 0.05" per axis) | 27 in 7 sequences |
|   |O−C| median / rms / p95 (") | 0.072 / 0.173 / 0.226 |
|   O−C mean RA·cosδ / Dec (") | 0.011 / -0.029 |
|   O−C rms RA·cosδ / Dec (") | 0.138 / 0.105 |
|   median |O−C|/σ ratio (0.83 calibrated) / corrected | 0.843 / 0.829 |
|   median per-axis |O−C|/σ RA / Dec (0.674 calibrated) | 0.617 / 0.738 |
| internal consistency (residual − sequence mean) median |O−C|/σ / per axis RA / Dec | 0.983 / 0.736 / 1.079 |
| tracks with no known object (false positives) per sequence | 0.231 |
|   same, associating only target/known (old metric) | 0.231 |
| tracks of other SkyBoT objects in the field | 0 |
| total time (s) / piazzi only (s) | 44.6 / 43.8 |

Detection by V bin:

| V | n | detected | rate |
|---|---|---|---|
| 16–17 | 2 | 2 | 100% |
| 17–18 | 2 | 1 | 50% |
| 18–19 | 4 | 4 | 100% |
| 19–20 | 4 | 4 | 100% |
| 20–21 | 1 | 0 | 0% |

| sequence | target | V | "/min | refs | det | |O−C| med (") | σ_eph (") | FP | status |
|---|---|---|---|---|---|---|---|---|---|
| ztf_20210905_f296_c6_q3_gr | Poseidon | 18.0 | 1.06 | 32,33,48,48 | yes | 0.07 | 0.009 | 0 | ok |
| ztf_20200124_f551_c9_q3_gr | Geographos | 18.7 | 1.32 | 91,93,79,79 | yes | 0.10 | 0.000 | 0 | ok |
| ztf_20201017_f692_c6_q2_r | Ryugu | 17.4 | 1.49 | 219,217,218,219,213 | yes | 0.04 | 0.020 | 1 | ok |
| ztf_20191102_f471_c10_q4_gr | Mithra | 19.4 | 1.89 | 59,72,72 | yes | 0.06 | 0.009 | 0 | ok |
| ztf_20200718_f622_c1_q4_r | Atira | 18.5 | 2.12 | 42,47,47,50 | yes | 0.07 | 0.039 | 0 | ok |
| ztf_20191029_f472_c16_q1_r | Aten | 19.7 | 2.52 | 54,51,52,46 | yes | 0.15 | 0.004 | 0 | ok |
| ztf_20210110_f396_c13_q2_gir | Toutatis | 16.8 | 3.70 | 49,51,34 | yes | 0.07 | 0.002 | 0 | ok |
| ztf_20180519_f790_c13_q2_r | 2018 JA | 19.2 | 7.75 | 54,53,53,56 | yes | 0.09 | 0.118 | 0 | ok |
| ztf_20210816_f646_c4_q4_gir | 2021 CH2 | 18.1 | 8.28 | 570,497,513,572 | yes | 0.14 | 0.096 | 0 | ok |
| ztf_20180604_f630_c6_q1_gr | 2018 EJ4 | 16.6 | 9.98 | 171,171,190,188 | yes | 0.21 | 0.100 | 2 | ok |
| ztf_20200823_f740_c6_q1_gr | 2019 AN5 | 17.2 | 11.33 | 1190,1321,1309 | no | — | — | 0 | ok |
| ztf_20181031_f784_c6_q2_r | 2018 VW3 | 20.4 | 11.81 | 151,159,156 | no | — | — | 0 | ok |
| ztf_20201029_f406_c2_q2_gr | 2020 UN6 | 19.2 | 12.50 | 319,322,243,241 | yes | 0.30 | 0.281 | 0 | ok |
