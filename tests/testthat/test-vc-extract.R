# fit_52() is provided by helper-fixtures.R.

test_that("vc_extract reproduces Section 5.2 components and df_PI FROM the Hessian", {
  skip_if_not_installed("lme4")
  comp <- vc_extract(fit_52())
  expect_s3_class(comp, "re_components")
  expect_equal(unname(comp$components[["run"]]), 0.000681, tolerance = 1e-4)
  expect_equal(unname(comp$components[["Residual"]]), 0.001253, tolerance = 1e-4)
  # df_PI from the observed-information Hessian. Locked tightly to the ported
  # route's value 12.486 (paper prints 12.484). This is
  # the no-drift gate: full-precision REML estimates, so no rounding excuse.
  expect_equal(comp$dfs[["pi"]], 12.486, tolerance = 1e-3)
  expect_identical(comp$design, "oneway")
  expect_false(is.null(comp$ems))
  expect_false(attr(comp, "singular"))
  expect_false(attr(comp, "pi_df_fallback"))
})

test_that("vc_extract end-to-end reproduces Section 5.2 CI / PI / TI", {
  skip_if_not_installed("lme4")
  comp <- vc_extract(fit_52())
  ci <- ci_mean(comp, conf = 0.95)
  pi <- pi_newobs(comp, level = 0.95)
  ti <- ti_francq(comp, content = 0.95, conf = 0.90)
  expect_equal(c(ci$lower, ci$upper), c(0.946, 1.016), tolerance = 2e-3)
  expect_equal(c(pi$lower, pi$upper), c(0.881, 1.081), tolerance = 2e-3)
  expect_equal(c(ti$lower, ti$upper), c(0.845, 1.117), tolerance = 2e-3)
})

test_that("vc_extract CI ddf equals the balanced between-groups df (A - 1 = 5)", {
  skip_if_not_installed("lme4")
  comp <- vc_extract(fit_52())
  expect_equal(comp$dfs[["ci"]], 5, tolerance = 1e-2)
})

test_that("balanced REML components equal the closed-form ANOVA estimates", {
  skip_if_not_installed("lme4")
  dat <- build_52_data()
  comp <- vc_extract(lme4::lmer(y ~ 1 + (1 | run), data = dat, REML = TRUE))
  grand <- mean(dat$y)
  runmeans <- tapply(dat$y, dat$run, mean)
  msb <- 3 * sum((runmeans - grand)^2) / (6 - 1)
  msw <- sum((dat$y - runmeans[dat$run])^2) / (18 - 6)
  expect_equal(unname(comp$components[["Residual"]]), msw, tolerance = 1e-6)
  expect_equal(unname(comp$components[["run"]]), (msb - msw) / 3, tolerance = 1e-6)
})

test_that("vc_extract rejects random-slope models", {
  skip_if_not_installed("lme4")
  dat <- build_52_data()
  dat$t <- rep(c(-1, 0, 1), 6)
  m <- suppressWarnings(lme4::lmer(y ~ t + (t | run), data = dat, REML = TRUE))
  expect_error(vc_extract(m), "random-intercept")
})

test_that("Section 5.3 unbalanced Hessian df_PI passes the sanity gate", {
  skip_if_not_installed("lme4")
  # Exact 11.31 needs genuine raw data (see deferred test below). This only
  # gates the unbalanced Hessian path: finite, positive, below the balanced df,
  # and roughly in range -- so the path has coverage, not a precision claim.
  dat <- build_52_data()
  d53 <- dat[!((dat$run %in% c(1, 2, 6)) & dat$rep == 1), ]
  comp52 <- vc_extract(fit_52())
  comp53 <- vc_extract(lme4::lmer(y ~ 1 + (1 | run), data = d53, REML = TRUE))
  df53 <- comp53$dfs[["pi"]]
  expect_true(is.finite(df53) && df53 > 0)
  expect_lt(df53, comp52$dfs[["pi"]])
  expect_gt(df53, 10)
  expect_lt(df53, 12)
  expect_false(attr(comp53, "pi_df_fallback"))     # interior fit -> Hessian route
})

test_that("Section 5.3 df_PI == 11.31 from a real fit (DEFERRED)", {
  skip(paste("Deferred: exact 5.3 df_PI (11.31) needs the genuine",
             "Hoffman-Kringle raw assay values; constructed unbalanced data",
             "reaches ~11.13. Sanity-gated in the test above."))
})
