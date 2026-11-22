# oliva2025_pqri

Oracle `oliva2025-pqri` (`../oracles.csv`), checked by `test-oracle-oliva-pqri.R`.

- **Paper.** Oliva-Aviles C, Hauser P (2025). *Technometrics* 67(2):193-202.
  doi:10.1080/00401706.2024.2407324. Section 6. The paper prints shelf lives
  of 27.9 (unbalanced, 70/95), 23.3 (unbalanced, 90/95) and 30.1 (balanced,
  90/95) months.
- **Data origin.** The PQRI stability data of Schwenke J, Quinlan M, Stroup W,
  Forenzo P (2020). Estimating shelf life through tolerance intervals.
  *AAPS PharmSciTech* 21:290 (no DOI).
- **`data.csv`.** 63 rows: `lot,time,assay,unbalanced`. Lots 1-9 are PQRI
  batches b1, b2, b5, b6, b8, b9, b10, b12 and b20; times are 0, 3, 6, 9, 12,
  18 and 24 months. `unbalanced` is TRUE for the 54 rows kept in the paper's
  unbalanced dataset, which drops t >= 18 for lots 1, 4, 7 and t = 24 for
  lots 2, 5, 8. The values come from the data chunk of the authors'
  supplementary code (SupplementaryCode.Rmd, SHA-256
  `62e5da0e03d6d96a8dbdc8eb573f7b0ec8a40e0c933bd92f9fd79e11316815c9`). Only
  the data were taken; none of the authors' code is in the repo.
- **`reference.csv`.** GENERATED numbers, not printed ones. On 2026-11-22 the
  authors' functions (`degRates_sim`, `TI_sim`, extracted unmodified with
  `knitr::purl`) were run outside the repo. For both datasets, at
  t0 = 0, 3, ..., 30, they give:
  - the lower 70/95 and 90/95 tolerance limits of the batch-mean distribution,
    and the median;
  - at the code's defaults, 10000 runs and seed 2023 (`value_10000`);
  - over 5 seeds (2024-2028) of 100000 runs (`mean_100000`, `sd_100000`).

  It also gives the three shelf lives: the first point of the 0.01-month grid
  at which the limit is below 90. The `#` header records the R and package
  versions. At 10000 runs a shelf life has a Monte-Carlo SD of about
  0.2 month, so the test compares shelf lives with the 100000-run mean.
- **SHA-256:**
  - `data.csv`: `34f786df766dd98d574876d541af4758769d05decfe1f6c75c9a500b39a52cdf`
  - `reference.csv`: `cfecbf2aad7adecf2d38d61c47188b6e6a25f66fc8220cd111c7d928668255e6`
