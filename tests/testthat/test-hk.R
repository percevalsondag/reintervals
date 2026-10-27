# Hoffman & Kringle (2005) one-way tolerance interval: the slope-free branch of
# the anova-mls core. Oracles at the printed precision:
#   O1-O4  Hoffman & Kringle (2005, 2007) on the correct data (hk_assay, HK 2007
#          Table III) and HK 2007 Example 2;
#   O5-O6  Francq, Lin & Hoyer (2019) Sections 5.2-5.3 (package ems-mls), which
#          used HK 2005 Table 6 with run 5 / replicate 3 misprinted as 1.020.

hk_printed <- function() {               # HK 2005 Table 6 as printed (Francq 5.2)
  d <- hk_assay
  d$conc[d$run == 5 & d$replicate == 3] <- 1.020
  d
}
hk_ex2 <- function() {                   # HK 2007 Table V (Example 2)
  conc <- list(c(0.1260, 0.1220, 0.1330), c(0.0969, 0.0963, 0.0931),
               c(0.0888, 0.0958, 0.0945), c(0.0991, 0.0921, 0.0982),
               c(0.1070, 0.0961, 0.1080), c(0.1150, 0.0989, 0.0961))
  data.frame(run = factor(rep(1:6, each = 3)), conc = unlist(conc))
}
fit_hk <- function(d) suppressMessages(lme4::lmer(conc ~ 1 + (1 | run), data = d))

## Independent HK implementation, written from the paper's formulas (HK 2007,
## Eqs. 1-2), sharing no code with the package.
hk_formula <- function(y, run, beta, gamma) {
  run <- factor(run); I <- nlevels(run); N <- length(y); Ji <- as.numeric(table(run))
  ybar_i <- tapply(y, run, mean); centre <- mean(ybar_i)
  JH <- I / sum(1 / Ji)
  MSBU <- JH * sum((ybar_i - centre)^2) / (I - 1)
  MSE <- sum((y - ybar_i[as.integer(run)])^2) / (N - I)
  sTOT2 <- MSBU / JH + (JH - 1) / JH * MSE
  Ne <- I * (MSBU + (JH - 1) * MSE) / MSBU
  H1 <- (I - 1) / qchisq(1 - gamma, I - 1) - 1
  H2 <- (N - I) / qchisq(1 - gamma, N - I) - 1
  hw <- qnorm((1 + beta) / 2) * sqrt(1 + 1 / Ne) *
    sqrt(sTOT2 + sqrt(H1^2 * (MSBU / JH)^2 + H2^2 * ((JH - 1) / JH)^2 * MSE^2))
  c(lower = centre - hw, upper = centre + hw)
}

test_that("hk_assay matches the published summary statistics (HK 2007 Table III)", {
  ybar_i <- tapply(hk_assay$conc, hk_assay$run, mean)
  expect_equal(round(mean(hk_assay$conc), 4), 0.9798)
  expect_equal(round(3 * sum((ybar_i - mean(hk_assay$conc))^2) / 5, 6), 0.003251)
  expect_equal(round(sum((hk_assay$conc - ybar_i[hk_assay$run])^2) / 12, 6), 0.001167)
})

test_that("O1: HK 2005, beta 0.95, gamma 0.90 -> (0.846, 1.114); default route is anova-mls", {
  skip_if_not_installed("lme4")
  x <- ti_lmm(fit_hk(hk_assay), P = 0.95, level = 0.90)
  expect_identical(x$method, "anova-mls")
  expect_equal(round(c(x$lower, x$upper), 3), c(0.846, 1.114))
  r <- ti_anova_raw(hk_assay$conc, time = NULL, lot = hk_assay$run, P = 0.95, level = 0.90)
  expect_equal(c(r$lower, r$upper), c(x$lower, x$upper), tolerance = 1e-14)
})

test_that("O2: HK 2007 Ex. 1, beta 0.667 -> (0.914, 1.046) with printed intermediates", {
  r <- ti_anova_raw(hk_assay$conc, time = NULL, lot = hk_assay$run, P = 0.667, level = 0.90)
  expect_equal(round(c(r$lower, r$upper), 3), c(0.914, 1.046))
  d <- attr(r, "diagnostics")
  expect_equal(round(d$Vhat_Y, 5), 0.00186)          # sTOT2
  expect_equal(round(d$C1, 4), 2.1050)               # H1
  expect_equal(round(d$C2, 4), 0.9036)               # H2
  # N_e: HK print 10.308, computed from their rounded MSB 0.003251 / MSE 0.001167
  # (6 * (0.003251 + 2 * 0.001167) / 0.003251 = 10.3076); unrounded it is 10.3096.
  expect_equal(d$n_E, 6 * (0.003251 + 2 * 0.001167) / 0.003251, tolerance = 2e-4)
  expect_equal(round(6 * (0.003251 + 2 * 0.001167) / 0.003251, 3), 10.308)
})

test_that("O3: HK 2007 Ex. 2, beta 0.667 -> (0.0799, 0.1265)", {
  d <- hk_ex2()
  r <- ti_anova_raw(d$conc, time = NULL, lot = d$run, P = 0.667, level = 0.90)
  expect_equal(round(c(r$lower, r$upper), 4), c(0.0799, 0.1265))
  expect_equal(round(attr(r, "diagnostics")$n_E, 3), 6.907)
})

test_that("O4: HK 2007 unbalanced (run 1 replicates 1-2 removed) -> (0.904, 1.048)", {
  d <- hk_assay[!(hk_assay$run == 1 & hk_assay$replicate %in% 1:2), ]
  r <- ti_anova_raw(d$conc, time = NULL, lot = d$run, P = 0.667, level = 0.90)
  expect_equal(round(c(r$lower, r$upper), 3), c(0.904, 1.048))
  expect_equal(round(r$estimate, 4), 0.9759)
  expect_equal(round(attr(r, "diagnostics")$n_E, 3), 9.270)
  expect_equal(round(attr(r, "diagnostics")$C2, 4), 1.0554)        # H2
})

test_that("O5: Francq 5.2 (printed data), ems-mls -> (0.845, 1.117)", {
  skip_if_not_installed("lme4")
  x <- ti_lmm(fit_hk(hk_printed()), P = 0.95, level = 0.90, method = "ems-mls")
  expect_identical(x$method, "ems-mls")
  expect_equal(round(c(x$lower, x$upper), 3), c(0.845, 1.117))
})

test_that("O6: Francq 5.3 (printed data minus replicate 1 of runs 1, 2, 6) -> (0.838, 1.127), df_PI 11.31", {
  skip_if_not_installed("lme4")
  # Closes the long-deferred Section 5.3 test: with the genuine (misprinted) raw
  # values the unbalanced Hessian df_PI reproduces the paper's 11.31 exactly.
  d <- hk_printed()
  fm <- fit_hk(d[!(d$run %in% c(1, 2, 6) & d$replicate == 1), ])
  x <- ti_lmm(fm, P = 0.95, level = 0.90, method = "ems-mls")
  expect_equal(round(c(x$lower, x$upper), 3), c(0.838, 1.127))
  expect_equal(round(lme4::fixef(fm)[[1]], 4), 0.9822)
  expect_equal(round(.design_of(fm)$n0, 2), 2.48)
  comp <- vc_extract(fm)
  expect_equal(round(comp$dfs[["pi"]], 2), 11.31)
  expect_false(attr(comp, "pi_df_fallback"))
})

test_that("HK and ems-mls coincide on balanced non-singular one-way data (tight bobyqa)", {
  skip_if_not_installed("lme4")
  tight <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12, maxfun = 1e5))
  for (d in list(hk_assay, hk_ex2())) {
    fm <- suppressMessages(lme4::lmer(conc ~ 1 + (1 | run), data = d, control = tight))
    expect_false(lme4::isSingular(fm))
    for (P in c(0.95, 0.667)) {
      h <- ti_lmm(fm, P = P, level = 0.90)
      e <- ti_lmm(fm, P = P, level = 0.90, method = "ems-mls")
      expect_equal(c(h$lower, h$upper), c(e$lower, e$upper), tolerance = 1e-6)
    }
  }
})

test_that("the slope-free core equals an independent HK implementation on 20 unbalanced datasets", {
  for (i in 1:20) {
    set.seed(1000 + i)
    I <- sample(3:8, 1)
    repeat { Ji <- sample(1:5, I, replace = TRUE); if (sum(Ji) > I && length(unique(Ji)) > 1) break }
    run <- factor(rep(seq_len(I), Ji))
    y <- rnorm(I, 0, runif(1, 0.2, 2))[run] + rnorm(length(run))
    beta <- runif(1, 0.6, 0.99); gamma <- runif(1, 0.6, 0.99)
    r <- ti_anova_raw(y, time = NULL, lot = run, P = beta, level = gamma)
    expect_equal(c(r$lower, r$upper), unname(hk_formula(y, run, beta, gamma)),
                 tolerance = 1e-12)
  }
})

test_that("HK route: method choices and the two-sided restriction", {
  skip_if_not_installed("lme4")
  fm <- fit_hk(hk_assay)
  expect_identical(ti_lmm(fm, method = "anova-mls")$method, "anova-mls")
  expect_identical(ti_lmm(fm, method = "ems-mls")$method, "ems-mls")
  expect_error(ti_lmm(fm, sides = "lower"), "two-sided only")
  expect_error(ti_lmm(fm, method = "gpq"), "not available")
  expect_error(ti_anova_raw(hk_assay$conc, time = NULL, lot = hk_assay$run, sides = "upper"),
               "two-sided only")
  expect_error(ti_anova_raw(hk_assay$conc, time = NULL, lot = hk_assay$run, t0 = 1),
               "only with a `time` covariate")
  # a one-way model WITH a covariate is not HK: here a clean-balanced fixed-slope
  # design, which stays on the EMS engine
  d <- hk_assay; d$x <- rep(c(-1, 0, 1), 6)
  expect_identical(ti_lmm(suppressMessages(lme4::lmer(conc ~ x + (1 | run), data = d)),
                          newdata = data.frame(x = 0))$method, "ems-mls")
})
