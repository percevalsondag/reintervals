# francq2019_esm

Oracle `francq2019-s51-esm` (`../oracles.csv`), checked by
`test-oracle-francq-s51.R`.

- **Source.** Supporting Information, Table 1, of Francq BG, Lin D, Hoyer W
  (2019). *Statistics in Medicine* 38(30):5603-5622. doi:10.1002/sim.8386.
  Orthopedic surgery data: error on safe margin (ESM) of 23 surgeons
  (seniority 1 = senior, 2 = junior), each with both procedures (1 = free hand,
  2 = navigated) at 4 slices. 184 rows.
- **Extraction.** From the text of the supplement PDF by Perceval Sondag,
  2026-11-22. `data.csv` has a `#` header naming the source and the check,
  then `surgeon,seniority,procedure,slice,esm`. The PDF is not in the repo.
- **Check.** Per-procedure REML fits, `esm ~ 0 + factor(slice) + (1 | surgeon)`,
  give surgeon / residual variances 7.349 / 18.074 (free hand) and
  0.447 / 3.611 (navigated), as printed in section 5.1. The test also
  reproduces the corrected Table 2 (Statistics in Medicine 2020;39(14):2015).
- **SHA-256:** `data.csv`
  `41b826d623d58616c42f0dcc1c9bd5e8b0d9324ec42952effa128d0a623e5afd`
