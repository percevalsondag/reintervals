## Build data/hk_assay.rda: Hoffman & Kringle (2007), Table III (Example 1).
## 6 analytical runs x 3 replicates of an assay concentration.
##
## HK (2005), Table 6, prints run 5 / replicate 3 as 1.020; the correct value is
## 0.998 (HK 2007, Table III). Only 0.998 is consistent with the published summary
## statistics (mean 0.9798, MSB 0.003251, MSE 0.001167). Francq et al. (2019),
## Sections 5.2-5.3, used the misprinted 1.020.
##
## Run from the package root:  Rscript data-raw/hk_assay.R
conc <- list(
  c(0.969, 0.976, 0.938),   # run 1
  c(0.952, 0.993, 0.956),   # run 2
  c(0.989, 0.883, 0.981),   # run 3
  c(1.000, 0.969, 0.954),   # run 4
  c(0.959, 0.989, 0.998),   # run 5 (HK 2005 Table 6 misprints replicate 3 as 1.020)
  c(1.020, 1.090, 1.020)    # run 6
)
hk_assay <- data.frame(
  run       = factor(rep(seq_along(conc), lengths(conc))),
  replicate = unlist(lapply(lengths(conc), seq_len)),
  conc      = unlist(conc)
)
## guard: the published summary statistics
ybar_i <- tapply(hk_assay$conc, hk_assay$run, mean)
msb <- 3 * sum((ybar_i - mean(hk_assay$conc))^2) / 5
mse <- sum((hk_assay$conc - ybar_i[hk_assay$run])^2) / 12
stopifnot(round(mean(hk_assay$conc), 4) == 0.9798, round(msb, 6) == 0.003251,
          round(mse, 6) == 0.001167)
save(hk_assay, file = "data/hk_assay.rda", version = 2)
