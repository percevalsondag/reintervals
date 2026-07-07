# Uses the shared `fixture_52()` (Francq Section 5.2) from helper-fixtures.R.

test_that("ci_mean reproduces Francq Section 5.2: 95% CI = [0.946, 1.016]", {
  ci <- ci_mean(fixture_52(), conf = 0.95)
  expect_s3_class(ci, "re_interval")
  expect_equal(round(ci$lower, 3), 0.946)
  expect_equal(round(ci$upper, 3), 1.016)
})

test_that("ci_mean equals the closed-form Eq. 12 half-width", {
  comp <- fixture_52()
  ci <- ci_mean(comp, conf = 0.95)
  hw <- stats::qt(0.975, df = 5) * sqrt(comp$var_mean)
  expect_equal(ci$estimate, 0.981)
  expect_equal(ci$lower, 0.981 - hw)
  expect_equal(ci$upper, 0.981 + hw)
  expect_equal(ci$df, 5)
  expect_identical(ci$type, "CI")
  expect_identical(ci$conf, 0.95)
  expect_identical(ci$sides, "two")
  expect_identical(ci$method, "francq")
})

test_that("ci_mean is symmetric about the estimate and widens with confidence", {
  comp <- fixture_52()
  ci90 <- ci_mean(comp, conf = 0.90)
  ci99 <- ci_mean(comp, conf = 0.99)
  expect_equal(ci90$estimate - ci90$lower, ci90$upper - ci90$estimate)
  expect_lt(ci99$upper - ci99$lower, Inf)
  expect_gt(ci99$upper - ci99$lower, ci90$upper - ci90$lower)
})

test_that("infinite CI df falls back to the normal quantile", {
  comp <- re_components(
    components = c(Residual = 1),
    dfs        = c(Residual = Inf, ci = Inf, pi = Inf),
    mean = 10, var_mean = 4
  )
  ci <- ci_mean(comp, conf = 0.95)
  expect_equal(ci$upper, 10 + stats::qnorm(0.975) * 2)
})

test_that("as.data.frame.re_interval yields a one-row tidy frame", {
  df <- as.data.frame(ci_mean(fixture_52()))
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), 1L)
  expect_true(all(c("type", "estimate", "lower", "upper", "df", "conf") %in%
                    names(df)))
  expect_identical(df$type, "CI")
})

test_that("ci_mean validates its arguments", {
  expect_error(ci_mean(list(mean = 1)), "re_components")
  expect_error(ci_mean(fixture_52(), conf = 1.2), "between 0 and 1")
  expect_error(ci_mean(fixture_52(), conf = NA), "between 0 and 1")

  no_ci <- re_components(
    components = c(run = 0.00046, Residual = 0.00147),
    dfs        = c(run = 5, Residual = 12, pi = 12.484),  # no "ci"
    mean = 0.981, var_mean = 0.01353^2
  )
  expect_error(ci_mean(no_ci), "\"ci\" entry")
})
