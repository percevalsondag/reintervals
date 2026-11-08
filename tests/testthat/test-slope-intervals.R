# Closed-form checks for the three slope intervals (CI / PI / new_group_mean).
# These pin the exact formula per branch so the suite stays green on CI:
#   CI             : mu +/- t(nu_F) * sqrt(V_F)
#   PI             : mu +/- t(nu_T) * sqrt(V_F + V_T)
#   new_group_mean : mu +/- t(nu_G) * sqrt(V_F + V_G)     (PI minus residual)
# vfun_extract's df machinery is tested separately (test-vfun-extract.R); these
# tests pin that .vfun_interval applies the right (df, se) per kind.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)

slope_fit <- function() {
  suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
                              data = lme4::sleepstudy, control = lc))
}

test_that(".vfun_interval matches the independent base-R recompute for CI / PI / new_group_mean", {
  skip_if_not_installed("lme4")
  fit  <- slope_fit()
  comp <- vfun_extract(fit, t0 = 5)
  mu   <- comp$mean
  tt   <- function(nu) stats::qt(0.975, nu)            # two-sided 95%

  ci  <- .vfun_interval(comp, kind = "CI",             level = 0.95)
  pi  <- .vfun_interval(comp, kind = "PI",             level = 0.95)
  ngm <- .vfun_interval(comp, kind = "new_group_mean", level = 0.95)

  # CI: df = nu_F, se = sqrt(V_F)
  expect_equal(ci$df, comp$dfs[["F"]])
  expect_equal(ci$se, sqrt(comp$var_mean))
  expect_equal(ci$lower, mu - tt(comp$dfs[["F"]]) * sqrt(comp$var_mean))
  expect_equal(ci$upper, mu + tt(comp$dfs[["F"]]) * sqrt(comp$var_mean))

  # PI: df = nu_T, se = sqrt(V_F + V_T)
  expect_equal(pi$df, comp$dfs[["T"]])
  expect_equal(pi$se, sqrt(comp$var_mean + comp$V_T))
  expect_equal(pi$lower, mu - tt(comp$dfs[["T"]]) * sqrt(comp$var_mean + comp$V_T))

  # new_group_mean: df = nu_G, se = sqrt(V_F + V_G)  (PI minus residual)
  expect_equal(ngm$df, comp$dfs[["G"]])
  expect_equal(ngm$se, sqrt(comp$var_mean + comp$V_G))
  expect_equal(ngm$lower, mu - tt(comp$dfs[["G"]]) * sqrt(comp$var_mean + comp$V_G))
})

test_that("the three intervals nest correctly: CI inside new_group_mean inside PI", {
  skip_if_not_installed("lme4")
  comp <- vfun_extract(slope_fit(), t0 = 5)
  ci  <- .vfun_interval(comp, "CI")
  pi  <- .vfun_interval(comp, "PI")
  ngm <- .vfun_interval(comp, "new_group_mean")
  # se ordering: V_F < V_F+V_G < V_F+V_T  (V_T = V_G + s2e >= V_G)
  expect_lt(ci$se, ngm$se)
  expect_lt(ngm$se, pi$se)
  # new_group_mean is PI with the residual s2e removed from the content variance
  expect_equal(pi$se^2 - ngm$se^2, comp$V_T - comp$V_G)   # = s2e
})

test_that("one-sided intervals open the correct side", {
  skip_if_not_installed("lme4")
  comp <- vfun_extract(slope_fit(), t0 = 5)
  for (k in c("CI", "PI", "new_group_mean")) {
    up <- .vfun_interval(comp, k, sides = "upper")
    lo <- .vfun_interval(comp, k, sides = "lower")
    expect_identical(up$lower, -Inf); expect_true(is.finite(up$upper))
    expect_identical(lo$upper,  Inf); expect_true(is.finite(lo$lower))
  }
})

test_that(".vfun_interval input guards", {
  skip_if_not_installed("lme4")
  comp <- vfun_extract(slope_fit(), t0 = 5)
  expect_error(.vfun_interval(comp, "CI", level = 1), "between 0 and 1")
  expect_error(.vfun_interval(list(), "CI"), "vfun_components")
  expect_error(.vfun_interval(comp, "bogus"), "should be one of")
})

# ---- dispatch (.interval_lmm; the public verbs wrap this) ----
test_that(".interval_lmm routes random-slope models to the slope ports", {
  skip_if_not_installed("lme4")
  fit <- slope_fit()
  for (k in c("CI", "PI", "new_group_mean")) {
    out  <- .interval_lmm(fit, newdata = data.frame(Days = c(0, 5)), kind = k)
    comp <- vfun_extract(fit, t0 = 5)
    expect_equal(nrow(out), 2L)
    expect_true(all(out$type == k))
    expect_equal(out$design[1], "M2c")
    # second row (Days = 5) equals the direct slope-constructor call
    direct <- .vfun_interval(comp, kind = k)
    expect_equal(out$lower[2], direct$lower)
    expect_equal(out$upper[2], direct$upper)
  }
})

test_that(".interval_lmm routes random-intercept CI/PI to the EMS engine", {
  skip_if_not_installed("lme4")
  m <- fit_52()                                          # y ~ 1 + (1 | run)
  ci <- .interval_lmm(m, kind = "CI")
  pi <- .interval_lmm(m, kind = "PI")
  ems_ci <- suppressWarnings(intervals_lmm(m, which = "CI"))
  ems_pi <- suppressWarnings(intervals_lmm(m, which = "PI"))
  expect_equal(ci$lower, ems_ci$lower)                   # same EMS numbers
  expect_equal(pi$upper, ems_pi$upper)
})

test_that("random-intercept new_group_mean (true mean of a new group) is narrower than the PI", {
  skip_if_not_installed("lme4")
  m   <- fit_52()
  ngm <- .ngm_ri(m, level = 0.95)                # native value-frame (no list)
  pi  <- suppressWarnings(intervals_lmm(m, which = "PI", level = 0.95))
  expect_s3_class(ngm, "data.frame")             # S5: not an re_interval list
  # between-level-only content -> strictly narrower than the future-observation PI
  expect_lt(ngm$upper - ngm$lower, pi$upper - pi$lower)
  # and the verb surfaces it with type = "new_group_mean" in the unified schema
  expect_identical(new_group_mean_lmm(m)$type, "new_group_mean")
})

# The closed-form random-slope tolerance interval (ti_vfun, "reml-mls") and its
# two tests (Graybill-Wang recompute; GPQ cross-method anchor) were removed in
# 1.0.0 with the method itself. The random-slope TI is GPQ (test-gpq.R,
# test-dispatch.R).
