# chow_shao1991

Oracle `chow-shao1991-t23` (`../oracles.csv`), checked by
`test-oracle-chow-shao.R`.

- **Source.** Chow S-C, Shao J (1991). Estimating drug shelf-life with random
  batches. *Biometrics* 47(3):1071-1079. doi:10.2307/2532659. Table 2, p. 1078:
  potency (% of label claim) of 24 batches at 0, 12, 24 and 36 months.
- **`data.csv`.** Long format (`batch`, `month`, `potency`), 96 rows, with a `#`
  header (read with `utils::read.csv(comment.char = "#")`). Transcribed in this
  repo by Perceval Sondag, 2026-11-22.
- **Transcription check.** Recomputing the paper's printed estimates of
  section 5 from the transcribed data reproduces them at the printed 3
  decimals: the WLS estimator with unstructured D (sample covariance of the 24
  batch vectors, divisor k - 1), a_w = 103.968 (SE .289), b_w = -.199 (SE .009),
  and the batch-21 OLS line 102.500 (SE 1.775), -.333 (SE .079).
- **Reference values.** None are stored: the exact t_23 limits are computed in
  closed form in the test (the mean and standard deviation of the per-batch
  OLS lines at each time point).
- **SHA-256** at transcription: `626c3922d91e68f047fe05471ea0550ede41fbc21b056f36a653a6e72936b087`.
