# chow_shao1991

Oracle `chow-shao1991-t23` (`../oracles.csv`), checked by
`test-oracle-chow-shao.R`.

- **Source.** Chow S-C, Shao J (1991). Estimating drug shelf-life with random
  batches. *Biometrics* 47(3):1071-1079. doi:10.2307/2532659. Table 2, p. 1078:
  potency (% of label claim) of 24 batches at 0, 12, 24 and 36 months.
- **`data.csv`.** Long format (`batch`, `month`, `potency`), 96 rows, with a `#`
  header (read with `utils::read.csv(comment.char = "#")`). Transcribed from
  Chow & Shao (1991), Table 2 (p. 1078) on 2026-11-22; verified against the
  printed section 5 estimates (below).
- **Transcription check.** Recomputing the paper's printed estimates of
  section 5 from the transcribed data reproduces them at the printed 3
  decimals: the WLS estimator with unstructured D (sample covariance of the 24
  batch vectors, divisor k - 1), a_w = 103.968 (SE .289), b_w = -.199 (SE .009),
  and the batch-21 OLS line 102.500 (SE 1.775), -.333 (SE .079).
- **Reference values.** None are stored: the exact t_23 limits are computed in
  closed form in the test (the mean and standard deviation of the per-batch
  OLS lines at each time point).
- **SHA-256:** `8f95199c092ec56008420c91d7d903f680a694b2b53e6d40b7c10c00dc2ab0de` (the provenance header line was reworded on
  2026-11-25; the data rows are unchanged).
