# oracle: oliva2025-pqri (fixtures/oracles.csv; absorbs the former pending oliva2025-worked-example)
#
# Oliva-Aviles C, Hauser P (2025). Technometrics 67(2):193-202.
# doi:10.1080/00401706.2024.2407324. Section 6: the PQRI stability data
# (fixtures/oliva2025_pqri/data.csv; originally Schwenke et al. 2020, see the
# README), balanced (9 lots x 7 times) and unbalanced (54 rows). The target is the
# lower (beta, gamma) tolerance limit for the batch-mean distribution
# (ti_lmm(over = "group_mean")), beta = 0.70 and 0.90, gamma = 0.95.
#
# Reference: GENERATED numbers, not printed ones (fixtures/oliva2025_pqri/
# reference.csv). They were computed by running the authors' supplementary code
# outside the repo. None of that code is in the repo. The paper prints only the
# three shelf lives, 27.9 / 23.3 / 30.1 months.
#
# Same estimator, checked by reading both sources: per-lot OLS lines; pivots
# U ~ chi2(n - 2K), Q ~ chi2(K - 1), Z ~ N(0, 1); the tau2 equation on
# (-phi1 s2e, Inf) with the clamp to -phi1 s2e (phi1 = min_k leverage); a center
# weighted by G = diag(tau2 + s2e a_k), computed with the tau2 realization even
# when it is negative (both allow negative realizations); content variance
# max(0, tau2); quantile type 7. The only differences are Monte-Carlo: the
# authors reseed per draw (set.seed(runs * seed + j)), while the package draws
# vectors under one seed. So the comparison is statistical. It can show agreement
# to about 3 x the combined MC SD (about 0.03-0.08 % label claim here). It cannot
# detect a systematic difference smaller than that, nor show draw-by-draw identity.
#
# Gate: |package - reference mean| <= 3 sqrt(package mc_se^2 + reference SD^2),
# where the reference mean and SD are over 5 seeds of 100000 runs.

pqri_data <- function() {
  d <- utils::read.csv(test_path("fixtures", "oliva2025_pqri", "data.csv"), comment.char = "#")
  d$lot <- factor(d$lot)
  d
}
pqri_ref <- function() utils::read.csv(test_path("fixtures", "oliva2025_pqri", "reference.csv"),
                                       comment.char = "#")
pqri_fit <- function(d) suppressWarnings(suppressMessages(
  lme4::lmer(assay ~ time + (time | lot), data = d)))
pqri_ti <- function(fm, t0, M, P = c(0.70, 0.90), seed = 2023L) suppressWarnings(
  ti_lmm(fm, data.frame(time = t0), over = "group_mean", sides = "lower", P = P, level = 0.95,
         control = lmm_interval_control(M = M, seed = seed)))
pqri_sets <- function(d) list(balanced = d, unbalanced = d[d$unbalanced, ])

## first point of the 0.01-month grid at which the lower limit is below 90, by
## bisection inside [a, b] (the limit decreases in t0; common random numbers
## across t0 under one seed make it a smooth function of t0)
pqri_shelf <- function(fm, P, a, b, M = 100000L) {
  f <- function(tt) pqri_ti(fm, tt, M, P = P)$lower
  stopifnot(f(a) >= 90, f(b) < 90)
  lo <- round(a * 100); hi <- round(b * 100)
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (f(mid / 100) < 90) hi <- mid else lo <- mid
  }
  hi / 100
}

test_that("oliva2025-pqri: the data and the unbalanced deletion pattern", {
  d <- pqri_data()
  expect_identical(dim(d), c(63L, 4L))
  expect_true(all(table(d$lot) == 7L))
  expect_equal(sort(unique(d$time)), c(0, 3, 6, 9, 12, 18, 24))
  drop <- (d$lot %in% c(1, 4, 7) & d$time >= 18) | (d$lot %in% c(2, 5, 8) & d$time == 24)
  expect_identical(d$unbalanced, !drop)
  expect_identical(sum(d$unbalanced), 54L)
  r <- pqri_ref()
  expect_true(all(c("quantity", "dataset", "t0", "beta", "gamma", "value_10000",
                    "mean_100000", "sd_100000", "n_seeds") %in% names(r)))
})

test_that("oliva2025-pqri: package GPQ lower limits vs the authors' code (3 combined MC SD)", {
  skip_if_not_installed("lme4")
  skip_on_cran()                                   # M = 100000: about 40 s
  ref <- pqri_ref(); ref <- ref[ref$quantity == "lower_limit", ]
  t0 <- seq(0, 30, 3)
  tab <- NULL
  for (nm in c("balanced", "unbalanced")) {
    fm <- pqri_fit(pqri_sets(pqri_data())[[nm]])
    for (M in c(20000L, 100000L)) {
      r <- pqri_ti(fm, t0, M)
      dg <- attr(r, "diagnostics")
      expect_identical(unique(r$method), "gpq")
      e <- ref[match(paste(nm, r$time, 100 * r$P), paste(ref$dataset, ref$t0, ref$beta)), ]
      bound <- 3 * sqrt(dg$mc_se_lower^2 + e$sd_100000^2)
      expect_true(all(abs(r$lower - e$mean_100000) <= bound), info = paste(nm, M))
      tab <- rbind(tab, data.frame(dataset = nm, M = M, t0 = r$time, beta = 100 * r$P,
        ref_10000 = e$value_10000, ref_mean = e$mean_100000, ref_sd = e$sd_100000,
        package = r$lower, mc_se = dg$mc_se_lower, diff = r$lower - e$mean_100000,
        z = (r$lower - e$mean_100000) / sqrt(dg$mc_se_lower^2 + e$sd_100000^2)))
    }
  }
  cat("\noliva2025-pqri lower limits (package vs authors' code):\n")
  print(format(tab, digits = 6), row.names = FALSE)
})

test_that("oliva2025-pqri: shelf lives within 0.1 month of the printed values and of the authors' code", {
  skip_if_not_installed("lme4")
  skip_on_cran()                                   # about 30 GPQ calls at M = 100000
  ref <- pqri_ref(); ref <- ref[ref$quantity == "shelf_life", ]
  sets <- pqri_sets(pqri_data())
  cases <- data.frame(dataset = c("unbalanced", "unbalanced", "balanced"), beta = c(70, 90, 90),
                      printed = c(27.9, 23.3, 30.1), a = c(27, 21, 27), b = c(30, 24, 33))
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    sl <- pqri_shelf(pqri_fit(sets[[cs$dataset]]), cs$beta / 100, cs$a, cs$b)
    e <- ref[ref$dataset == cs$dataset & ref$beta == cs$beta, ]
    info <- paste(cs$dataset, cs$beta, "package", sl)
    expect_lte(abs(sl - cs$printed), 0.1 + 1e-9, label = info)
    expect_lte(abs(sl - e$mean_100000), 0.1 + 1e-9, label = info)
    cat(sprintf("\nshelf life %s (%d,95): package %.2f | printed %.1f | authors' code 10000 runs %.2f, 100000 runs mean %.2f (SD %.3f)",
                cs$dataset, cs$beta, sl, cs$printed, e$value_10000, e$mean_100000, e$sd_100000))
  }
  cat("\n")
})
