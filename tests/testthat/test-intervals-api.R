# public surface: the four verbs, the unified door, control, the shared
# arg-resolution helpers, the route-independent provenance, and the V_G floor.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
slope_fit <- function() {
  suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
                              data = lme4::sleepstudy, control = lc))
}

test_that("lmm_interval_control validates its tuning knobs", {
  ctl <- lmm_interval_control(seed = 1, M = 5000)
  expect_s3_class(ctl, "lmm_interval_control")
  expect_identical(ctl$M, 5000L)
  expect_error(lmm_interval_control(M = 0), "positive")
  expect_error(lmm_interval_control(seed = "x"), "seed")
})

test_that(".resolve_alpha implements the mutually-exclusive alpha/level pair", {
  expect_equal(.resolve_alpha(NULL, NULL), 0.05)        # neither -> default
  expect_equal(.resolve_alpha(0.10, NULL), 0.10)        # alpha given
  expect_equal(.resolve_alpha(NULL, 0.90), 0.10)        # level -> 1 - level
  expect_error(.resolve_alpha(0.05, 0.95), "only one")  # both -> error
  expect_error(.resolve_alpha(NULL, 1.2), "\\(0, 1\\)")
})

test_that("defensive class check fires at the front door of every verb", {
  for (f in list(ci_lmm, pi_lmm, new_group_mean_lmm, ti_lmm, lmm_predict)) {
    expect_error(f(1L), "lme4 model|lmerMod")
  }
})

test_that("P lives on ti_lmm only; the coverage verbs have no P argument", {
  m <- fit_52()
  expect_error(ci_lmm(m, P = 0.99), "unused argument")
  expect_error(pi_lmm(m, P = 0.99), "unused argument")
  expect_error(new_group_mean_lmm(m, P = 0.99), "unused argument")
  expect_silent(suppressWarnings(ti_lmm(m, P = 0.99)))  # ti_lmm accepts P
})

test_that(".model_class is route-independent and matches the design tags", {
  skip_if_not_installed("lme4")
  expect_identical(.model_class(fit_52()), "random_intercept")
  dm1 <- vfun_sim(B = 8, times = c(0, 3, 6, 9, 12), s2_0 = 2, s2_1 = 0, s2e = 0.7, seed = 11)
  fm1 <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch), data = dm1, control = lc))
  expect_identical(.model_class(fm1), "random_intercept_fixed_slope")
  expect_identical(.model_class(slope_fit()), "random_slope_correlated")
  di <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2, s2_1 = 0.01,
                 s01 = 0, s2e = 0.7, seed = 5)
  f2i <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch) + (0 + time | batch),
                                     data = di, control = lc))
  expect_identical(.model_class(f2i), "random_slope_independent")
  # grouping attribute: oneway for intercept models, NA for slope models
  expect_identical(.grouping_of(fit_52()), "oneway")
  expect_true(is.na(.grouping_of(slope_fit())))
})

test_that("lmm_predict dispatches to the verbs and gates P/over to tolerance", {
  skip_if_not_installed("lme4")
  m <- fit_52()
  expect_equal(lmm_predict(m, "confidence")$lower, ci_lmm(m)$lower)
  expect_equal(lmm_predict(m, "prediction")$upper, pi_lmm(m)$upper)
  # P / over only valid for interval = "tolerance"
  expect_error(lmm_predict(m, "confidence", P = 0.99), "tolerance")
  expect_error(lmm_predict(m, "prediction", over = "group_mean"), "tolerance")
  expect_silent(suppressWarnings(lmm_predict(m, "tolerance", P = 0.99)))
})

test_that("the four verbs carry harmonised provenance (method/design/grouping)", {
  skip_if_not_installed("lme4")
  fs <- slope_fit()
  ci  <- ci_lmm(fs, newdata = data.frame(Days = 5))
  expect_identical(ci$method, "reml-mls")
  expect_identical(ci$design, "random_slope_correlated")
  expect_true(is.na(attr(ci, "grouping")))
  # EMS verb: ems-mls + random_intercept + grouping
  ti <- suppressWarnings(ti_lmm(fit_52(), P = 0.99))
  expect_identical(ti$method, "ems-mls")
  expect_identical(ti$design, "random_intercept")
  expect_null(attr(ti, "design"))                       # legacy attr dropped (now a column)
})

# ---- V_G floor-and-flag (ratified): floor at 0 + warn(t0) + flag --------------
test_that(".floor_vg floors a model-implied negative V_G, flags it, and warns naming t0", {
  expect_warning(r <- .floor_vg(-0.5, t0 = 24), "V_G\\(t0 = 24\\)")
  expect_equal(r$value, 0); expect_true(r$floored)
  # non-negative passes through silently, unflagged
  expect_silent(r2 <- .floor_vg(3.2, t0 = 5))
  expect_equal(r2$value, 3.2); expect_false(r2$floored)
  # numerical round-off is absorbed (max(0,.)) WITHOUT a spurious warning/flag
  expect_silent(r3 <- .floor_vg(-1e-15, t0 = 5))
  expect_equal(r3$value, 0); expect_false(r3$floored)
})

test_that("slope interval rows surface the vg_floored flag (FALSE in the normal V_G>=0 regime)", {
  skip_if_not_installed("lme4")
  comp <- vfun_extract(slope_fit(), t0 = 5)
  expect_false(isTRUE(attr(comp, "vg_floored")))        # PD Sigma -> V_G >= 0
  expect_false(ti_vfun(comp)$vg_floored)
  expect_false(.vfun_interval(comp, "new_group_mean")$vg_floored)
})

# ---- the unified output contract ---------------------------------------------
test_that("every verb on every engine returns the identical core column contract (no list)", {
  skip_if_not_installed("lme4")
  core <- c("type","estimate","lower","upper","level","P","sides",
            "method","design","singular","boundary","vg_floored")
  m  <- fit_52()                                   # EMS
  fs <- slope_fit()                                # slope
  verbs_ri <- list(ci_lmm(m), pi_lmm(m), suppressWarnings(ti_lmm(m, P = 0.99)),
                   new_group_mean_lmm(m))
  verbs_sl <- list(ci_lmm(fs, newdata = data.frame(Days = 5)),
                   pi_lmm(fs, newdata = data.frame(Days = 5)),
                   ti_lmm(fs, newdata = data.frame(Days = 5), P = 0.99),
                   new_group_mean_lmm(fs, newdata = data.frame(Days = 5)))
  for (x in c(verbs_ri, verbs_sl)) {
    expect_s3_class(x, "lmm_interval")
    expect_s3_class(x, "data.frame")              # never a list
    expect_true(all(core %in% names(x)))          # full core contract
    # the trailing 12 columns are exactly the core, in order (eval cols precede)
    expect_identical(tail(names(x), length(core)), core)
  }
  # type tags are the lower-case kinds
  expect_identical(vapply(verbs_sl, function(z) z$type[1], ""),
                   c("ci", "pi", "ti", "new_group_mean"))
  # P is NA for coverage intervals, the value for TI
  expect_true(is.na(ci_lmm(m)$P))
  expect_equal(suppressWarnings(ti_lmm(m, P = 0.99))$P, 0.99)
})

test_that("engine diagnostics live in attr(.,'diagnostics'), off the core frame", {
  skip_if_not_installed("lme4")
  g <- ti_lmm(slope_fit(), newdata = data.frame(Days = 5), over = "group_mean",
              control = lmm_interval_control(seed = 1, M = 500))
  d <- attr(g, "diagnostics")
  expect_type(d, "list")
  expect_true(all(c("tau2_med", "s2e_med", "frac_boundary") %in% names(d)))  # GPQ audit
  expect_false(any(c("tau2_med", "frac_boundary") %in% names(g)))            # not columns
})

# ---- A3: multi-covariate EMS (random intercept, several fixed effects) --------
test_that("A3: all four verbs compute on a multi-fixed-covariate random-intercept model", {
  skip_if_not_installed("lme4")
  set.seed(7)
  n <- 120; g <- factor(rep(1:12, each = 10))
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  b  <- stats::rnorm(12, 0, 2)
  y  <- 50 + 2 * x1 - 1.5 * x2 + b[as.integer(g)] + stats::rnorm(n)
  fm <- suppressWarnings(lme4::lmer(y ~ x1 + x2 + (1 | g), control = lc))
  nd <- data.frame(x1 = 0.5, x2 = -0.5)
  core <- c("type","estimate","lower","upper","level","P","sides",
            "method","design","singular","boundary","vg_floored")

  ci  <- ci_lmm(fm, newdata = nd)
  pi  <- pi_lmm(fm, newdata = nd)
  ti  <- suppressWarnings(ti_lmm(fm, newdata = nd, P = 0.99))
  ngm <- new_group_mean_lmm(fm, newdata = nd)
  for (x in list(ci, pi, ti, ngm)) {
    expect_s3_class(x, "lmm_interval")
    expect_true(all(core %in% names(x)))
    expect_identical(x$method, "ems-mls")
    expect_identical(x$design, "random_intercept")
    expect_true(all(c("x1","x2") %in% names(x)))   # both eval covariates carried
    expect_true(is.finite(x$lower) && x$upper > x$lower)
  }
  # nesting sanity: CI inside new_group_mean inside PI at the same point
  expect_lt(ci$upper - ci$lower, ngm$upper - ngm$lower)
  expect_lt(ngm$upper - ngm$lower, pi$upper - pi$lower)
})

# ---- A4: EMS new_group_mean with MULTI-ROW newdata ----------------------------
test_that("A4: new_group_mean_lmm handles multi-row newdata (each row correct, schema holds)", {
  skip_if_not_installed("lme4")
  fm  <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a,
                                     control = lc))
  nd  <- data.frame(Month = c(0, 6, 12))
  # MBL is single-obs fixed-slope -> the C1 anti-conservative warning fires (tested
  # in the C1 block); suppress it here since this test is about multi-row mechanics.
  out <- suppressWarnings(new_group_mean_lmm(fm, newdata = nd, level = 0.95))
  expect_s3_class(out, "lmm_interval")
  expect_equal(nrow(out), 3L)
  expect_equal(out$Month, c(0, 6, 12))
  expect_true(all(out$type == "new_group_mean"))
  expect_true(all(out$method == "ems-mls"))
  expect_true(all(is.finite(out$lower)) && all(out$upper > out$lower))
  # each row equals the single-row call at that Month (the loop is row-correct)
  for (i in seq_len(3)) {
    one <- suppressWarnings(new_group_mean_lmm(fm, newdata = nd[i, , drop = FALSE],
                                               level = 0.95))
    expect_equal(out$lower[i], one$lower)
    expect_equal(out$upper[i], one$upper)
  }
})

# ---- B5: print.lmm_interval flag note ----------------------------------------
test_that("B5: print.lmm_interval shows the flag note iff a row is flagged", {
  skip_if_not_installed("lme4")
  x <- ci_lmm(fit_52())
  expect_false(any(grepl("^note:", capture.output(print(x)))))  # nothing flagged
  x$singular <- TRUE                                             # force a flag
  expect_true(any(grepl("note:.*flag", capture.output(print(x)))))
})

# ---- B8: internal-class S3 methods dispatch (registered via S3method, not export)
test_that("B8: print/as.data.frame dispatch for the internal classes", {
  skip_if_not_installed("lme4")
  # NAMESPACE registers these via S3method(...) (no bare export); verify dispatch
  rc <- fixture_52()                                   # re_components
  expect_s3_class(rc, "re_components")
  expect_gt(length(capture.output(print(rc))), 0L)     # print.re_components fires
  ri <- ci_mean(rc)                                    # re_interval
  expect_s3_class(ri, "re_interval")
  expect_gt(length(capture.output(print(ri))), 0L)     # print.re_interval fires
  expect_s3_class(as.data.frame(ri), "data.frame")     # as.data.frame.re_interval fires
  vc <- vfun_extract(slope_fit(), t0 = 5)              # vfun_components
  expect_s3_class(vc, "vfun_components")
  expect_true(any(grepl("vfun_components", capture.output(print(vc)))))
})

# ---- C1: anti-conservative warning on single-obs / unbalanced fixed-slope -----
test_that("C1 FIRES: ci/pi/new_group_mean warn on single-obs (MBL) fixed-slope data", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a,
                                    control = lc))
  nd <- data.frame(Month = 12)
  expect_warning(ci_lmm(fm, newdata = nd),             "anti-conservative")
  expect_warning(pi_lmm(fm, newdata = nd),             "anti-conservative")
  expect_warning(new_group_mean_lmm(fm, newdata = nd), "anti-conservative")
})

test_that("C1 DOES NOT FIRE: clean balanced M1, plain random-intercept, slope, or ti_lmm", {
  skip_if_not_installed("lme4")
  # clean balanced fixed-slope: EMS numbers are exact -> NO warning
  db <- vfun_sim(B = 8, times = c(0, 3, 6, 9, 12), s2_0 = 2, s2_1 = 0, s2e = 0.7, seed = 11)
  fb <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch), data = db, control = lc))
  expect_warning(ci_lmm(fb, newdata = data.frame(time = 12)),             NA)
  expect_warning(pi_lmm(fb, newdata = data.frame(time = 12)),             NA)
  expect_warning(new_group_mean_lmm(fb, newdata = data.frame(time = 12)), NA)
  # plain random-intercept (no fixed slope): NO warning
  expect_warning(ci_lmm(fit_52()), NA)
  # random-slope model never goes through the EMS CI/PI path: NO warning
  expect_warning(ci_lmm(slope_fit(), newdata = data.frame(Days = 5)), NA)
  # ti_lmm on MBL: TI is re-routed to anova-mls, so NO anti-conservative warning
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a,
                                    control = lc))
  expect_warning(ti_lmm(fm, newdata = data.frame(Month = 12), P = 0.99), NA)
})

test_that("C1: the warning is ADDITIVE -- MBL interval numbers are unchanged", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a,
                                    control = lc))
  nd <- data.frame(Month = 12)
  ci <- suppressWarnings(ci_lmm(fm, newdata = nd))
  pi <- suppressWarnings(pi_lmm(fm, newdata = nd))
  # equal the underlying EMS engine output (the warning changes no computation)
  ems_ci <- suppressWarnings(intervals_lmm(fm, newdata = nd, which = "CI"))
  ems_pi <- suppressWarnings(intervals_lmm(fm, newdata = nd, which = "PI"))
  expect_equal(ci$lower, ems_ci$lower); expect_equal(ci$upper, ems_ci$upper)
  expect_equal(pi$lower, ems_pi$lower); expect_equal(pi$upper, ems_pi$upper)
})

# ---- S2: clear error on malformed newdata (names the missing predictor) -------
test_that("S2 FIRES: missing-predictor newdata stops naming the column, all verbs", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 | Subject),
                                    data = lme4::sleepstudy, control = lc))
  bad <- data.frame(Nonsense = 5)
  # the missing column NAME (Days) must appear in the message, for every verb
  expect_error(ci_lmm(fm, newdata = bad),             "Days")
  expect_error(pi_lmm(fm, newdata = bad),             "Days")
  expect_error(new_group_mean_lmm(fm, newdata = bad), "Days")
  expect_error(suppressWarnings(ti_lmm(fm, newdata = bad, P = 0.99)), "Days")
  expect_error(ci_lmm(fm, newdata = bad), "missing required predictor")
  # bare-vector newdata (N1) -> clear message, not a base-R data.frame error
  expect_error(ci_lmm(fm, newdata = 5), "must be a data frame")
  # multi-covariate: each missing predictor is named
  set.seed(1)
  d <- data.frame(y = stats::rnorm(60), x1 = stats::rnorm(60),
                  x2 = stats::rnorm(60), g = factor(rep(1:6, 10)))
  fm2 <- suppressWarnings(lme4::lmer(y ~ x1 + x2 + (1 | g), data = d, control = lc))
  expect_error(ci_lmm(fm2, newdata = data.frame(z = 1)), "x1, x2")
  expect_error(ci_lmm(fm2, newdata = data.frame(x1 = 1)), "x2")
})

test_that("S2 DOES NOT FIRE on valid newdata; the guard moves no interval number", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 | Subject),
                                    data = lme4::sleepstudy, control = lc))
  nd <- data.frame(Days = 5)
  # valid call is silent and equals the underlying EMS engine output (guard is
  # a pure input check -- it changes no computation)
  ci  <- expect_silent(ci_lmm(fm, newdata = nd))
  ems <- intervals_lmm(fm, newdata = nd, which = "CI")
  expect_equal(ci$lower, ems$lower); expect_equal(ci$upper, ems$upper)
  # NULL newdata and list newdata still work (unchanged behavior)
  expect_silent(ci_lmm(fm))
  expect_equal(ci_lmm(fm, newdata = list(Days = 5))$lower, ci$lower)
})
