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
- **Transcription.** Transcribed from the printed Table 2 and section 5.4 on
  2026-11-25; verified against the printed bottle means (3 dp) and sse (4 dp),
  which the test recomputes from the data.
- **SHA-256:**
  - `data.csv`: `43a4a40b59cb344709f9c463fc6c0e6342eade64c280b2d0355d24f2071720ef`
  - `printed.csv`: `8f33ad690831c3c45ca6560a382198f68b03d3471d77472eaeed9bcf931fbb8e`

  (The provenance header lines were reworded on 2026-11-25; the data rows are
  unchanged.)
