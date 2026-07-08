# fit_52() is provided by helper-fixtures.R.

test_that("ti_lmm reproduces Section 5.2 TI (content = 0.95, conf = 0.90)", {
  skip_if_not_installed("lme4")
  out <- .ti_ems(fit_52(), level = 0.95, conf = 0.90)
  expect_equal(c(out$lower, out$upper), c(0.845, 1.117), tolerance = 2e-3)
  expect_identical(attr(out, "design"), "oneway")
  expect_true(attr(out, "balanced"))
  expect_false(attr(out, "singular"))
})

test_that("ti_lmm preserves the data.frame columns and attributes", {
  skip_if_not_installed("lme4")
  out <- .ti_ems(fit_52(), level = 0.95, conf = 0.90)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("type", "estimate", "lower", "upper", "level", "conf") %in%
                    names(out)))
  expect_identical(unique(out$type), "TI")
})

test_that("ti_lmm rejects the true_value target for the fixed-slope design", {
  skip_if_not_installed("lme4")
  expect_error(.ti_ems(fit_52(), target = "true_value"), "not supported")
})

test_that("ti_lmm errors for a genuinely off-catalog design (crossed, no interaction)", {
  skip_if_not_installed("lme4")
  # two crossed random mains WITHOUT an interaction term -> outside the catalog
  # (unbalanced nested/crossed WITH interaction are now supported in M2). No
  # closed-form variance decomposition exists, so ti_lmm signals a clean error
  # (previously it returned NA bounds + a warning).
  dat <- expand.grid(a = factor(1:4), b = factor(1:4))
  dat$y <- as.numeric(dat$a) + as.numeric(dat$b)
  m <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ 1 + (1 | a) + (1 | b), data = dat,
               control = lme4::lmerControl(check.conv.singular = "ignore"))
  ))
  expect_error(.ti_ems(m), "not available")
})

test_that("ti_lmm rejects random-slope models", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(
    lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  )
  expect_error(.ti_ems(fm), "random-intercept")
})
