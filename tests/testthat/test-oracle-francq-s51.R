# oracle: francq2019-s51-esm (fixtures/oracles.csv)
#
# Francq BG, Lin D, Hoyer W (2019). Statistics in Medicine 38(30):5603-5622.
# doi:10.1002/sim.8386. Section 5.1, the orthopedic data (Supporting Information
# Table 1; fixtures/francq2019_esm/data.csv), against Table 2 as corrected in
# Statistics in Medicine 2020;39(14):2015. The uncorrected 2019 Table 2 has sign
# and digit misprints and is not used.
#
# Model (the authors' SAS PROC MIXED per procedure: MODEL ESM = slice / NOINT,
# RANDOM surgeon, DDFM=KR): esm ~ 0 + factor(slice) + (1 | surgeon), REML, fitted
# separately for each procedure. The REML likelihood factorises over procedures
# because no parameter is shared. Each fit is balanced (23 surgeons x 4 slices):
# H_surgeon on 22 df, H_residual on 92 - 4 - 22 = 66 df, k = 1/n0 = 1/4.
#
# Every printed value is reproduced at its printed precision except the
# navigated PI df, printed 84.93, package 84.909 (84.91). On this balanced fit
# the observed-information, expected-information and ANOVA (22 / 66 df) Var(total)
# all equal 0.38790, so the choice of information matrix is not the cause. The
# print corresponds to Var(total) = 2 x 4.058114^2 / 84.93 = 0.38781 (0.023%
# smaller), which rounds to the printed 0.388 too. The printed digits cannot
# separate the two, so the package value is asserted at its own 2 dp. A typo or
# a slightly different SAS optimum are both possible.

esm_data <- function() {
  d <- utils::read.csv(test_path("fixtures", "francq2019_esm", "data.csv"), comment.char = "#")
  d$slice <- factor(d$slice); d$surgeon <- factor(d$surgeon)
  d
}
esm_fit <- function(d, proc) suppressMessages(
  lme4::lmer(esm ~ 0 + slice + (1 | surgeon), data = d[d$procedure == proc, ], REML = TRUE))

## corrected Table 2 (2 dp), rows slice 1-4, per procedure (1 = free hand, 2 = navigated)
esm_table2 <- list(
  `1` = data.frame(est = c(0.11, 3.53, -2.49, -1.91),
                   ci_l = c(-1.98, 1.43, -4.59, -4.01), ci_u = c(2.21, 5.63, -0.40, 0.18),
                   pi_l = c(-10.16, -6.74, -12.77, -12.19), pi_u = c(10.39, 13.80, 7.78, 8.36),
                   ti_l = c(-10.89, -7.47, -13.50, -12.92), ti_u = c(11.11, 14.53, 8.51, 9.09)),
  `2` = data.frame(est = c(1.29, 1.28, -0.89, 0.19),
                   ci_l = c(0.45, 0.45, -1.72, -0.64), ci_u = c(2.12, 2.12, -0.05, 1.03),
                   pi_l = c(-2.80, -2.81, -4.98, -3.90), pi_u = c(5.38, 5.37, 3.20, 4.28),
                   ti_l = c(-3.06, -3.07, -5.24, -4.16), ti_u = c(5.64, 5.63, 3.46, 4.54)))
## section 5.1 text
esm_text <- list(
  `1` = list(vc = c(7.349, 18.074), total = 25.423, var_total = 18.372, pi_df = 70.36,
             kr_df = 70.4, se1 = 1.051, est1 = 0.114),
  `2` = list(vc = c(0.447, 3.611), total = 4.058, var_total = 0.388, pi_df = 84.93,
             kr_df = 84.9))

test_that("francq2019-s51-esm: the data are 23 surgeons x 2 procedures x 4 slices", {
  d <- esm_data()
  expect_identical(dim(d), c(184L, 5L))
  expect_identical(nlevels(d$surgeon), 23L)
  expect_true(all(table(d$surgeon, d$procedure, d$slice) == 1L))
  expect_true(all(tapply(d$seniority, d$surgeon, function(x) length(unique(x))) == 1L))
})

test_that("francq2019-s51-esm: section 5.1 components, Var(total), df and the corrected Table 2", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  d <- esm_data()
  nd <- data.frame(slice = factor(1:4))
  for (p in c("1", "2")) {
    fm <- esm_fit(d, as.integer(p))
    tx <- esm_text[[p]]; tb <- esm_table2[[p]]
    expect_false(lme4::isSingular(fm))
    expect_identical(.design_of(fm)$type, "oneway")
    expect_equal(.design_of(fm)$n0, 4)                                   # k = 1/4
    expect_identical(c(nlevels(droplevels(fm@flist$surgeon)) - 1L,
                       nrow(fm@frame) - 4L - 22L), c(22L, 66L))           # H df
    th <- .theta_cov(fm)
    expect_false(th$pi_df_fallback)
    expect_equal(round(unname(th$components), 3), tx$vc, info = p)
    expect_equal(round(th$sigmaT2, 3), tx$total, info = p)
    expect_equal(round(sum(th$cov_theta), 3), tx$var_total, info = p)
    if (p == "1") {
      expect_equal(round(th$df_pi, 2), tx$pi_df)
    } else {
      ## printed 84.93; package 84.91 (see the header). On this balanced fit the
      ## ANOVA closed form (22 / 66 df) gives the same Var(total) as the Hessian.
      expect_equal(round(th$df_pi, 2), 84.91)
      s <- d[d$procedure == 2, ]
      ym <- tapply(s$esm, s$surgeon, mean)
      mse <- sum(stats::residuals(stats::lm(esm ~ slice + surgeon, s))^2) / 66
      msa <- 4 * sum((ym - mean(ym))^2) / 22
      expect_equal(sum(th$cov_theta), 2 * (msa / 4)^2 / 22 + 2 * (3 / 4 * mse)^2 / 66,
                   tolerance = 1e-6)
      expect_equal(round(2 * th$sigmaT2^2 / 84.93, 5), 0.38781)           # what the print implies
    }

    ck <- ci_lmm(fm, nd, level = 0.95, control = lmm_interval_control(ci_df = "kenward-roger"))
    cs <- ci_lmm(fm, nd, level = 0.95)
    expect_equal(round(attr(ck, "diagnostics")$df, 1), rep(tx$kr_df, 4), info = p)
    expect_equal(attr(cs, "diagnostics")$df, attr(ck, "diagnostics")$df, tolerance = 1e-6) # balanced
    if (p == "1") {
      expect_equal(round(ck$estimate[1], 3), tx$est1)
      expect_equal(round(sqrt(as.matrix(pbkrtest::vcovAdj(fm))[1, 1]), 3), tx$se1)
    }
    pi <- pi_lmm(fm, nd, level = 0.95)
    ti <- ti_lmm(fm, nd, P = 0.95, level = 0.80)
    expect_identical(ti$method, rep("ems-mls", 4))                       # default route
    for (x in list(ck, cs)) {
      expect_equal(round(x$estimate, 2), tb$est, info = p)
      expect_equal(round(x$lower, 2), tb$ci_l, info = p)
      expect_equal(round(x$upper, 2), tb$ci_u, info = p)
    }
    expect_equal(round(pi$lower, 2), tb$pi_l, info = p)
    expect_equal(round(pi$upper, 2), tb$pi_u, info = p)
    expect_equal(round(ti$lower, 2), tb$ti_l, info = p)
    expect_equal(round(ti$upper, 2), tb$ti_u, info = p)
  }
})
