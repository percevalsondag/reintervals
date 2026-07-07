# Uses the shared `fixture_52()` (Francq Section 5.2) from helper-fixtures.R.

test_that("pi_newobs reproduces Francq Section 5.2: 95% PI = [0.881, 1.081]", {
  pi <- pi_newobs(fixture_52(), level = 0.95)
  expect_s3_class(pi, "re_interval")
  expect_equal(round(pi$lower, 3), 0.881)
  expect_equal(round(pi$upper, 3), 1.081)
})

test_that("pi_newobs equals the closed-form Eq. 21 half-width", {
  comp <- fixture_52()
  pi <- pi_newobs(comp, level = 0.95)
  sigma2_t <- 0.000681 + 0.001253                      # = 0.001934
  hw <- stats::qt(0.975, df = 12.484) * sqrt(comp$var_mean + sigma2_t)
  expect_equal(pi$estimate, 0.981)
  expect_equal(pi$lower, 0.981 - hw)
  expect_equal(pi$upper, 0.981 + hw)
  expect_equal(pi$df, 12.484)
  expect_identical(pi$type, "PI")
  expect_identical(pi$level, 0.95)
  expect_identical(pi$sides, "two")
  expect_identical(pi$method, "francq")
  expect_identical(pi$target, "observable")
})

test_that("pi_newobs total variance is sum(coefs * components)", {
  # The PI must be wider than the CI: it adds sigma^2_T to var_mean.
  comp <- fixture_52()
  ci <- ci_mean(comp, conf = 0.95)
  pi <- pi_newobs(comp, level = 0.95)
  expect_gt(pi$upper - pi$lower, ci$upper - ci$lower)
})

test_that("the true_value target narrows the PI and is recorded", {
  # Observable: sigma^2_T = run + residual. True value: residual excluded.
  obs <- fixture_52()
  tv <- re_components(
    components = c(run = 0.000681, Residual = 0.001253),
    dfs        = c(ci = 5, pi = 12.484),
    mean       = 0.981,
    var_mean   = 0.01353^2,
    coefs      = c(run = 1, Residual = 0),
    target     = "true_value",
    design     = "oneway"
  )
  pi_obs <- pi_newobs(obs)
  pi_tv  <- pi_newobs(tv)
  expect_identical(pi_tv$target, "true_value")
  expect_lt(pi_tv$upper - pi_tv$lower, pi_obs$upper - pi_obs$lower)
})

test_that("pi_newobs is symmetric and widens with the prediction level", {
  comp <- fixture_52()
  p90 <- pi_newobs(comp, level = 0.90)
  p99 <- pi_newobs(comp, level = 0.99)
  expect_equal(p90$estimate - p90$lower, p90$upper - p90$estimate)
  expect_gt(p99$upper - p99$lower, p90$upper - p90$lower)
})

test_that("infinite PI df falls back to the normal quantile", {
  comp <- re_components(
    components = c(Residual = 3),
    dfs        = c(ci = Inf, pi = Inf),
    mean = 10, var_mean = 1
  )
  pi <- pi_newobs(comp, level = 0.95)
  expect_equal(pi$upper, 10 + stats::qnorm(0.975) * sqrt(1 + 3))
})

test_that("pi_newobs validates its arguments", {
  expect_error(pi_newobs(list(mean = 1)), "re_components")
  expect_error(pi_newobs(fixture_52(), level = 0), "between 0 and 1")
  expect_error(pi_newobs(fixture_52(), level = NA), "between 0 and 1")

  no_pi <- re_components(
    components = c(run = 0.000681, Residual = 0.001253),
    dfs        = c(ci = 5),                          # no "pi"
    mean = 0.981, var_mean = 0.01353^2
  )
  expect_error(pi_newobs(no_pi), "\"pi\" entry")
})
