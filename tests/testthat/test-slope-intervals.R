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
  comp <- vfun_extract(fit, t0 = 5, target = "observable")
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
    comp <- vfun_extract(fit, t0 = 5, target = "observable")
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

test_that("EMS new_group_mean (true_value path) drops the residual vs the observable PI", {
  skip_if_not_installed("lme4")
  m   <- fit_52()
  ngm <- .ngm_ems(m, level = 0.95)               # native value-frame (no list)
  pi  <- suppressWarnings(intervals_lmm(m, which = "PI", level = 0.95))
  expect_s3_class(ngm, "data.frame")             # S5: not an re_interval list
  # between-level-only content -> strictly narrower than the future-observation PI
  expect_lt(ngm$upper - ngm$lower, pi$upper - pi$lower)
  # and the verb surfaces it with type = "new_group_mean" in the unified schema
  expect_identical(new_group_mean_lmm(m)$type, "new_group_mean")
})

# ---- Independent anchors for the reml-mls slope TI (ti_vfun) ------------------
# Two checks validate the slope TI without relying on a stored fixture:
#   (1) an exact recompute of the published Graybill-Wang MLS formula, reading
#       only the separately-anchored component fields (catches implementation
#       drift); and
#   (2) a cross-method agreement check against the GPQ engine (Oliva-Hauser, a
#       different derivation) -- the non-circular anchor (catches a wrong formula,
#       which (1) alone could not, since it recomputes the same formula).
#
# Why not a reduce-to-Montes equality anchor: on balanced fixed-slope data the
# reml-mls and anova-mls (Montes) tolerance intervals differ by ~0.094. They
# share the variance-component decomposition exactly (same V_T, same n_E) but
# apply different tolerance factors (Graybill-Wang two-component vs the Montes
# factor), so the two engines are close but not equal by construction. An
# equality anchor would therefore be wrong; hence (1) an independent formula
# recompute + (2) the GPQ cross-method check instead.

.m2i_fit <- function() {
  di <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18, 24),
                 s2_0 = 2, s2_1 = 0.01, s01 = 0, s2e = 0.7, seed = 5)
  suppressWarnings(lme4::lmer(y ~ time + (1 | batch) + (0 + time | batch),
                              data = di, control = lc))
}

test_that("ti_vfun reproduces the published Graybill-Wang MLS TI formula", {
  skip_if_not_installed("lme4")
  P <- 0.99; g <- 0.95
  cases <- list(M2c = list(fit = slope_fit(), t0 = 5),
                M2i = list(fit = .m2i_fit(),  t0 = 12))
  for (nm in names(cases)) {
    comp <- vfun_extract(cases[[nm]]$fit, t0 = cases[[nm]]$t0, target = "observable")
    got  <- ti_vfun(comp, content = P, conf = g)

    # INDEPENDENT re-derivation of the two-component MLS tolerance interval from
    # the PUBLISHED Graybill-Wang formula, using ONLY base-R stats and the
    # component fields (which are anchored elsewhere: extraction vs the published
    # Montes oracle, test-vfun-extract.R). It calls NO package internal and does
    # NOT read ti_vfun's output -- k is rebuilt from scratch.
    mu  <- comp$mean; V_G <- comp$V_G; V_T <- comp$V_T; s2e <- comp$phi[["s2e"]]
    nuG <- comp$dfs[["G"]]; nue <- comp$dfs[["e"]]; nE <- comp$n_E
    H_G <- nuG / stats::qchisq(1 - g, nuG) - 1          # between-component tail factor
    H_e <- nue / stats::qchisq(1 - g, nue) - 1          # residual-component tail factor
    U   <- V_T + sqrt((H_G * V_G)^2 + (H_e * s2e)^2)     # content-variance upper bound
    hw  <- stats::qnorm((1 + P) / 2) * sqrt(1 + 1 / nE) * sqrt(U)

    expect_equal(got$lower, mu - hw, tolerance = 1e-9, info = paste(nm, "lower"))
    expect_equal(got$upper, mu + hw, tolerance = 1e-9, info = paste(nm, "upper"))
    expect_equal(got$estimate, mu,   tolerance = 1e-9, info = paste(nm, "estimate"))
  }
})

test_that("ti_vfun (reml-mls closed form) agrees with the GPQ engine on the observable TI (cross-method anchor)", {
  skip_if_not_installed("lme4")
  P <- 0.99; g <- 0.95; t0 <- 5
  fit  <- slope_fit()                                   # M2c (sleepstudy)
  reml <- ti_vfun(vfun_extract(fit, t0 = t0, target = "observable"), content = P, conf = g)
  # GPQ is a fully independent method (Monte-Carlo generalized pivotal quantities,
  # Oliva-Aviles & Hauser) computed from the RAW data -- not the closed form
  # c_resid = 1 is the future-observation (observable) target ti_vfun uses.
  gpq  <- ti_gpq_raw(lme4::sleepstudy$Reaction, lme4::sleepstudy$Days,
                     lme4::sleepstudy$Subject, t0 = t0, P = P, level = g, c_resid = 1,
                     control = lmm_interval_control(seed = 1, M = 20000))
  width <- reml$upper - reml$lower
  gap   <- max(abs(c(reml$lower - gpq$lower, reml$upper - gpq$upper)))
  # the two methods agree to well under 1% of the interval width (observed ~0.2%);
  # 2% is a robust cross-method-agreement bound (MC + closed-form-vs-pivotal gap).
  expect_lt(gap, 0.02 * width)
})
