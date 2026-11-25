# liao2005_sulfur

Oracle `liao2005-sulfur` (`../oracles.csv`), checked by `test-oracle-lli2005.R`.

- **Paper.** Liao CT, Lin TY, Iyer HK (2005). *Technometrics* 47(3):323-335.
  doi:10.1198/004017005000000102. Example 2 (section 5.4) and Table 2. The data
  come originally from Wang and Iyer (1994).
- **`data.csv`.** 10 rows: `bottle,replicate,sulfur` (sulfur content of SRM 2682,
  weight %; bottles 1227, 0410, 0241 and 0809 measured twice, 0398 and 1500
  once). Read `bottle` as character to keep the leading zeros.
- **`printed.csv`.** The printed values: the six bottle means, sse, and the
  two-sided (beta = 0.99, gamma = 0.95) tolerance interval for the true-value
  distribution N(mu, sigma_alpha^2), [4.473, 4.870]. The authors computed it by
  Monte Carlo with K = 10000 runs, so the printed limits carry their own
  Monte-Carlo error (about 0.005 per limit).
- **Transcription.** Typed from the printed Table 2 and section 5.4 on
  2026-11-25, for Perceval Sondag (AI-assisted transcription). The test
  checks the data against the printed bottle means (3 dp) and sse (4 dp).
- **SHA-256:**
  - `data.csv`: `8c38f0d72fca53545c2bf5d8cd0d3e9227f8f7f9ef760a1140d699653843f482`
  - `printed.csv`: `929283e9513cdcc28c234bda801f358d2f89bc2507098918a7af52eecfab5947`
