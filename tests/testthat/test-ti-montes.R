# The new M1 ANOVA/MLS engine. Its whole reason to exist is reproducing the
# Montes-Burdick-Leblond (2019) Appendix A oracle EXACTLY -- the published TI and
# its printed intermediates -- on unbalanced data with single-observation lots.

test_that("ti_anova_raw reproduces the MBL Appendix A oracle to the published decimals", {
  d <- mbl_appendix_a
  out <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, P = 0.99, level = 0.95)

  expect_equal(out$lower, 57.25, tolerance = 0.01)   # published two-sided TI
  expect_equal(out$upper, 75.41, tolerance = 0.01)
  # printed intermediates (the published worked example). The paper prints the slope as
  # a magnitude (0.457); the genuine regression slope is negative (degradation),
  # so the engine returns the signed value -0.457.
  expect_equal(attr(out, "diagnostics")$beta_hat, -0.457, tolerance = 1e-3)
  expect_equal(abs(attr(out, "diagnostics")$beta_hat), 0.457, tolerance = 1e-3)
  expect_equal(unname(attr(out, "diagnostics")$n_E), 6.729, tolerance = 0.01)
  expect_equal(attr(out, "diagnostics")$U,        10.820, tolerance = 0.01)
})

test_that("ti_anova_raw returns the documented one-row-per-t0 schema", {
  d <- mbl_appendix_a
  out <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = c(0, 6, 12))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_equal(out$time, c(0, 6, 12))
  expect_true(all(out$type == "ti"))
  expect_true(all(out$design == "random_intercept_fixed_slope"))
  expect_true(all(out$method == "anova-mls"))
  expect_true(all(attr(out, "diagnostics")$df == min(11 - 1, 26 - 11 - 1)))   # min(s, r)
  expect_length(attr(out, "diagnostics")$n_E, 3L)                    # one effective n per t0
})

test_that("no ANOVA-estimate truncation: the mean squares are SS/df, non-negative by construction", {
  # S2_L = sum(adj^2)/(I-1) and S2_E = within-SS/r are both sums of squares over
  # df, so Vhat_Y and U are >= 0 with no truncation step. Confirm on a design
  # whose lot means are nearly identical (tiny between-lot signal): no NaN/<0.
  set.seed(1)
  lots <- rep(LETTERS[1:6], each = 4)
  time <- rep(c(0, 3, 6, 9), 6)
  y    <- 100 - 0.2 * time + stats::rnorm(24, 0, 0.3)    # tight, near-zero S2_L
  out  <- ti_anova_raw(y, time, lots, t0 = 9)
  expect_gte(attr(out, "diagnostics")$Vhat_Y, 0)
  expect_gte(attr(out, "diagnostics")$U, 0)
  expect_true(is.finite(out$lower) && is.finite(out$upper))
})

test_that("one-sided intervals open the unbounded side", {
  d <- mbl_appendix_a
  up <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, sides = "upper")
  lo <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, sides = "lower")
  expect_identical(up$lower, -Inf)
  expect_true(is.finite(up$upper))
  expect_identical(lo$upper, Inf)
  expect_true(is.finite(lo$lower))
})

test_that("input guards: P/level in (0,1), finite t0, >= 2 lots, within-lot time variation", {
  d <- mbl_appendix_a
  expect_error(ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, P = 1), "between 0 and 1")
  expect_error(ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, level = 0), "(0, 1)")
  expect_error(ti_anova_raw(d$Y, d$Month, d$Lot, t0 = Inf), "finite")
  expect_error(ti_anova_raw(d$Y, d$Month, d$Lot, t0 = numeric(0)), "non-empty")
  expect_error(ti_anova_raw(1:3, 0:2, rep("A", 3), t0 = 1), ">= 2 lots")
  # all observations at a single time point -> no within-lot slope is estimable
  expect_error(
    ti_anova_raw(c(1, 2, 3, 4), c(0, 0, 0, 0), c("A", "A", "B", "B"), t0 = 0),
    "within-lot time variation"
  )
})

test_that("one-sided limits are Hoffman's (2010) H: Montes et al. (2019) Tables I-II on mbl_appendix_a", {
  d <- mbl_appendix_a; r2 <- function(x) round(x, 2)
  ## intercept-only (the slope-free / Hoffman-Kringle branch), P = 0.99, conf 0.95
  two <- ti_anova_raw(d$Y, time = NULL, lot = d$Lot, P = 0.99, level = 0.95)
  lo  <- ti_anova_raw(d$Y, time = NULL, lot = d$Lot, P = 0.99, level = 0.95, sides = "lower")
  up  <- ti_anova_raw(d$Y, time = NULL, lot = d$Lot, P = 0.99, level = 0.95, sides = "upper")
  expect_equal(r2(c(two$lower, two$upper)), c(60.44, 81.97))
  expect_equal(r2(c(lo$lower, up$upper)), c(60.63, 81.78))
  g <- attr(two, "diagnostics")
  yb <- tapply(d$Y, d$Lot, mean)
  expect_equal(round(sum((yb - mean(yb))^2) / 10, 3), 5.917)                     # S2_L
  expect_equal(round(sum((d$Y - yb[d$Lot])^2) / (nrow(d) - 11), 3), 3.636)       # S2_E
  expect_equal(round(g$n_E, 3), 13.110)
  expect_equal(round(g$U, 3), 16.232)
  expect_equal(round(g$C1, 3), 1.538)                                            # H1
  expect_equal(round(g$C2, 3), 1.066)                                            # H2
  ## fixed slope at t0 = 12
  lo12 <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, P = 0.99, level = 0.95, sides = "lower")
  up12 <- ti_anova_raw(d$Y, d$Month, d$Lot, t0 = 12, P = 0.99, level = 0.95, sides = "upper")
  expect_equal(r2(c(lo12$lower, up12$upper)), c(57.36, 75.31))
  ## not the superseded "HK1" form (two-sided formula with z_P): 61.48 / 80.92, 58.13 / 74.53
  expect_false(isTRUE(all.equal(c(lo$lower, up$upper), c(61.48, 80.92), tolerance = 1e-3)))
  expect_false(isTRUE(all.equal(c(lo12$lower, up12$upper), c(58.13, 74.53), tolerance = 1e-3)))
  expect_gt(abs(lo$lower - 61.48), 0.5); expect_gt(abs(lo12$lower - 58.13), 0.5)
})
