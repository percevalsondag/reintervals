# oracle: francq2019-s52-s53 (fixtures/oracles.csv)
#
# Francq BG, Lin D, Hoyer W (2019). Statistics in Medicine 38(30):5603-5622.
# doi:10.1002/sim.8386. Sections 5.2 (balanced) and 5.3 (unbalanced), end to end
# through the exported verbs, on the Hoffman D, Kringle R (2005), J Biopharm Stat
# 15:283-293, one-way data (6 runs x 3 replicates) as Francq et al. used them:
# hk_assay with run 5 / replicate 3 set to the HK 2005 Table 6 misprint 1.020
# (data-raw/hk_assay.R). Section 5.3 removes the first replicate of runs 1, 2
# and 6 (15 observations).
#
# The paper's exponents are typos: "x 10^7" is x 10^-7 throughout.
#
# Precision: each printed number is asserted at its printed digits (round() for
# fixed decimals, signif() for the variance-scale numbers). Where a printed
# number is not reproduced, the package value is asserted at its own value and
# the mismatch is diagnosed by a reconstruction of the paper's arithmetic; no
# tolerance is loosened.
#
# Section 5.3 mismatches (diagnosed below):
#   - Var(total) printed 8.316e-7; package 8.3173e-7. The print is the sum of
#     the printed, rounded entries: 8.605 + 6.151 - 2 x 3.22 = 8.316.
#   - df_PI printed 11.31448; package 11.31143. The print is 2 x 0.002169^2 /
#     8.316e-7, i.e. computed from the rounded total variance and the rounded
#     Var(total). At 2 dp both are 11.31.
#   - SE printed 0.01441; package (model-based GLS SE) 0.01429. The print is the
#     Kenward-Roger adjusted SE with W = the observed-information covariance of
#     the variance components (0.014414). pbkrtest, the package's
#     "kenward-roger" arm, uses the expected information (0.014387).
#   - CI printed [0.943, 1.022]; package default (Satterthwaite, df 4.18, model SE)
#     [0.943, 1.021]; package KR arm (pbkrtest, df 4.91) [0.945, 1.019]. The print
#     is mean +- t(0.975, 4.18) x 0.014414: the observed-information KR SE on the
#     Satterthwaite df.
#   - TI: the default route for an intercept-only one-way model is anova-mls
#     (Hoffman-Kringle unweighted-means MLS, ti_lmm / .is_oneway_intercept_only,
#     R/dispatch.R), giving [0.832, 1.136]. The printed interval is the Francq EMS
#     MLS (EMS_A = sigma_e^2 + 2.48 sigma_a^2, H_A on 5 df, H_e on 9 df) =
#     method = "ems-mls", which reproduces [0.838, 1.127]. On the balanced
#     Section 5.2 data the two engines coincide.

fr_data <- function(section = c("5.2", "5.3")) {
  section <- match.arg(section)
  d <- hk_assay
  d$conc[d$run == 5 & d$replicate == 3] <- 1.020
  if (section == "5.3") d <- d[!(d$run %in% c(1, 2, 6) & d$replicate == 1), ]
  d
}
fr_fit <- function(d) suppressMessages(lme4::lmer(conc ~ 1 + (1 | run), data = d, REML = TRUE))

## Kenward-Roger adjusted variance of the intercept with a given covariance W of
## the variance components (Kenward MG, Roger JH (1997), Biometrics 53(3):983-997,
## doi:10.2307/2533558; R_ij = 0 for a model linear in the components). Written
## from the formula; uses only the model matrices.
fr_kr_se <- function(fm, W) {
  y <- lme4::getME(fm, "y"); X <- lme4::getME(fm, "X"); N <- length(y)
  zt <- lme4::getME(fm, "Zt")
  A <- list(as.matrix(Matrix::crossprod(zt)), diag(N))
  vc <- as.data.frame(lme4::VarCorr(fm))$vcov
  V <- vc[1] * A[[1]] + vc[2] * A[[2]]
  Vi <- solve(V)
  Phi <- solve(t(X) %*% Vi %*% X)
  P <- lapply(A, function(a) -t(X) %*% Vi %*% a %*% Vi %*% X)
  S <- 0
  for (i in 1:2) for (j in 1:2)
    S <- S + W[i, j] * (t(X) %*% Vi %*% A[[i]] %*% Vi %*% A[[j]] %*% Vi %*% X -
                          P[[i]] %*% Phi %*% P[[j]])
  sqrt(drop(Phi + 2 * Phi %*% S %*% Phi))
}

## the replicate order: dropping replicate 1 (and not 2 or 3) of runs 1, 2, 6
## reproduces the printed mean and components
test_that("francq2019-s52-s53: Section 5.3 removes replicate 1 (not 2 or 3) of runs 1, 2 and 6", {
  skip_if_not_installed("lme4")
  d <- fr_data("5.2")
  hit <- vapply(1:3, function(r) {
    fm <- fr_fit(d[!(d$run %in% c(1, 2, 6) & d$replicate == r), ])
    vc <- as.data.frame(lme4::VarCorr(fm))$vcov
    round(lme4::fixef(fm)[[1]], 4) == 0.9822 && signif(vc[1], 3) == 0.000575 &&
      signif(vc[2], 4) == 0.001594
  }, logical(1))
  expect_identical(hit, c(TRUE, FALSE, FALSE))
  expect_identical(nrow(fr_data("5.3")), 15L)
})

test_that("francq2019-s52-s53: Section 5.2 (balanced) at the printed precision", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  fm <- fr_fit(fr_data("5.2"))
  expect_false(lme4::isSingular(fm))
  th <- .theta_cov(fm)
  expect_false(th$pi_df_fallback)
  expect_equal(signif(th$components[["run"]], 3), 0.000681)
  expect_equal(signif(th$components[["Residual"]], 4), 0.001253)
  W <- th$cov_theta
  expect_equal(signif(W[1, 1], 4), 5.123e-7)
  expect_equal(signif(W[2, 2], 4), 2.617e-7)
  expect_equal(signif(W[1, 2], 3), -8.72e-8)
  expect_equal(signif(th$sigmaT2, 4), 0.001934)
  expect_equal(signif(sum(W), 4), 5.995e-7)
  expect_equal(round(th$df_pi, 3), 12.484)

  ci <- ci_lmm(fm, level = 0.95)
  ck <- ci_lmm(fm, level = 0.95, control = lmm_interval_control(ci_df = "kenward-roger"))
  expect_equal(round(ci$estimate, 3), 0.981)
  expect_equal(round(sqrt(stats::vcov(fm)[1, 1]), 5), 0.01353)
  expect_equal(round(sqrt(pbkrtest::vcovAdj(fm)[1, 1]), 5), 0.01353)   # balanced: no KR inflation
  for (x in list(ci, ck)) {
    expect_equal(attr(x, "diagnostics")$df, 5, tolerance = 1e-6)
    expect_equal(round(c(x$lower, x$upper), 3), c(0.946, 1.016))
  }
  pi <- pi_lmm(fm, level = 0.95)
  expect_equal(round(attr(pi, "diagnostics")$df, 3), 12.484)
  expect_equal(round(c(pi$lower, pi$upper), 3), c(0.881, 1.081))
  ti <- ti_lmm(fm, P = 0.95, level = 0.90)
  te <- ti_lmm(fm, P = 0.95, level = 0.90, method = "ems-mls")
  expect_identical(ti$method, "anova-mls")      # default: Hoffman-Kringle MLS
  expect_identical(te$method, "ems-mls")        # Francq EMS MLS
  for (x in list(ti, te)) expect_equal(round(c(x$lower, x$upper), 3), c(0.845, 1.117))
})

test_that("francq2019-s52-s53: Section 5.3 (unbalanced) at the printed precision; mismatches diagnosed", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  d <- fr_data("5.3")
  fm <- fr_fit(d)
  expect_false(lme4::isSingular(fm))
  th <- .theta_cov(fm)
  expect_false(th$pi_df_fallback)
  expect_equal(signif(th$components[["run"]], 3), 0.000575)
  expect_equal(signif(th$components[["Residual"]], 4), 0.001594)
  W <- th$cov_theta
  expect_equal(signif(W[1, 1], 4), 8.605e-7)
  expect_equal(signif(W[2, 2], 4), 6.151e-7)
  expect_equal(signif(W[1, 2], 3), -3.22e-7)
  expect_equal(signif(th$sigmaT2, 4), 0.002169)
  ## EMS_A = sigma_e^2 + 2.48 sigma_a^2; H_A on I - 1 = 5, H_e on N - I = 9 df
  expect_equal(.design_of(fm)$n0, 2.48, tolerance = 1e-12)
  expect_identical(c(nlevels(droplevels(d$run)) - 1L, nrow(d) - nlevels(droplevels(d$run))), c(5L, 9L))

  ## Var(total), printed 8.316e-7: package 8.317e-7; the print sums the rounded entries
  expect_equal(signif(sum(W), 4), 8.317e-7)
  expect_equal(8.605e-7 + 6.151e-7 - 2 * 3.22e-7, 8.316e-7, tolerance = 1e-12)
  ## df_PI, printed 11.31448: package 11.31143 (2 dp: 11.31); the print is the
  ## Satterthwaite formula on the rounded total variance and the rounded Var(total)
  expect_equal(round(th$df_pi, 2), 11.31)
  expect_equal(round(th$df_pi, 5), 11.31143)
  expect_equal(round(2 * signif(th$sigmaT2, 4)^2 /
                       (signif(W[1, 1], 4) + signif(W[2, 2], 4) + 2 * signif(W[1, 2], 3)), 5),
               11.31448)

  ci <- ci_lmm(fm, level = 0.95)
  ck <- ci_lmm(fm, level = 0.95, control = lmm_interval_control(ci_df = "kenward-roger"))
  expect_equal(round(ci$estimate, 4), 0.9822)
  ## SE, printed 0.01441: the package's (model-based) SE is 0.01429; pbkrtest's
  ## KR-adjusted SE (expected information) 0.01439; KR with the observed
  ## information (the package's .theta_cov) reproduces the print
  expect_equal(round(sqrt(stats::vcov(fm)[1, 1]), 5), 0.01429)
  expect_equal(round(sqrt(pbkrtest::vcovAdj(fm)[1, 1]), 5), 0.01439)
  se_kr_obs <- fr_kr_se(fm, W)
  expect_equal(round(se_kr_obs, 5), 0.01441)
  ## CI, printed [0.943, 1.022]
  expect_equal(round(attr(ci, "diagnostics")$df, 2), 4.18)
  expect_equal(round(c(ci$lower, ci$upper), 3), c(0.943, 1.021))          # default: upper differs
  expect_equal(round(attr(ck, "diagnostics")$df, 2), 4.91)
  expect_equal(round(c(ck$lower, ck$upper), 3), c(0.945, 1.019))          # KR arm: both differ
  expect_equal(round(ci$estimate + c(-1, 1) * stats::qt(0.975, attr(ci, "diagnostics")$df) *
                       se_kr_obs, 3), c(0.943, 1.022))                     # the print, reconstructed

  pi <- pi_lmm(fm, level = 0.95)
  expect_equal(round(c(pi$lower, pi$upper), 3), c(0.875, 1.089))
  ti <- ti_lmm(fm, P = 0.95, level = 0.90)
  te <- ti_lmm(fm, P = 0.95, level = 0.90, method = "ems-mls")
  expect_identical(ti$method, "anova-mls")      # default route: Hoffman-Kringle, not the print
  expect_equal(round(c(ti$lower, ti$upper), 3), c(0.832, 1.136))
  expect_identical(te$method, "ems-mls")
  expect_equal(round(c(te$lower, te$upper), 3), c(0.838, 1.127))
})
