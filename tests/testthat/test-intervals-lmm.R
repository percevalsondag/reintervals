# fit_52() / build_52_data() are provided by helper-fixtures.R.

test_that("intervals_lmm reproduces Section 5.2 CI and PI from a real fit", {
  skip_if_not_installed("lme4")
  out <- intervals_lmm(fit_52())
  ci <- out[out$type == "CI", ]
  pi <- out[out$type == "PI", ]
  expect_equal(c(ci$lower, ci$upper), c(0.946, 1.016), tolerance = 2e-3)
  expect_equal(c(pi$lower, pi$upper), c(0.881, 1.081), tolerance = 2e-3)
})

test_that("intervals_lmm preserves the data.frame columns and attributes (GxP contract)", {
  skip_if_not_installed("lme4")
  out <- intervals_lmm(fit_52())
  expect_s3_class(out, "data.frame")
  expect_true(all(c("type", "estimate", "lower", "upper", "df", "level") %in%
                    names(out)))
  expect_setequal(out$type, c("CI", "PI"))
  # the GxP-read attributes
  expect_false(is.null(attr(out, "components")))
  expect_named(attr(out, "components"), c("run", "Residual"))
  expect_false(attr(out, "pi_df_fallback"))
  expect_false(attr(out, "singular"))
  # df_PI on the PI row comes from the Hessian (12.486)
  expect_equal(out$df[out$type == "PI"], 12.486, tolerance = 1e-3)
})

test_that("which = 'CI' computes no Hessian and sets no PI attributes", {
  skip_if_not_installed("lme4")
  out <- intervals_lmm(fit_52(), which = "CI")
  expect_setequal(out$type, "CI")
  expect_null(attr(out, "components"))
  expect_null(attr(out, "pi_df_fallback"))
  expect_null(attr(out, "singular"))
})

test_that("intervals_lmm handles multiple newdata rows and keeps predictor columns", {
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  out <- intervals_lmm(fm, newdata = data.frame(Days = c(0, 9)))
  expect_equal(nrow(out), 4L)                # 2 rows x (CI + PI)
  expect_true("Days" %in% names(out))
  expect_equal(sum(out$type == "PI"), 2L)
})

test_that("intervals_lmm rejects random-slope models", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(
    lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  )
  expect_error(intervals_lmm(fm), "random-intercept")
})
