# dispatch through the public verb ti_lmm.
# Checks: (1) EMS numeric core unchanged through ti_lmm; (2) the
# vc-extract random-slope stop() still fires direct; (3) every random-slope
# tolerance interval is GPQ (1.0.0: the closed form reml-mls is defunct); (4)
# M1-overlap routing
# (balanced -> EMS, single-obs -> anova-mls oracle); (5) the EMS and slope df fallbacks stay distinct; plus the
# harmonized provenance (method/design/grouping) and over/P/method args.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)

slope_fit <- function() {
  suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
                              data = lme4::sleepstudy, control = lc))
}

test_that(".has_random_slope classifies fits by structure (intercept FALSE, slope TRUE)", {
  skip_if_not_installed("lme4")
  expect_false(.has_random_slope(fit_52()))                       # y ~ 1 + (1|run)
  expect_true(.has_random_slope(slope_fit()))                     # M2c
  di <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24),
                 s2_0 = 2, s2_1 = 0.01, s01 = 0, s2e = 0.7, seed = 5)
  fm2i <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch) + (0 + time | batch),
                                      data = di, control = lc))
  expect_true(.has_random_slope(fm2i))                            # M2i
})

# ---- (1) EMS numeric core unchanged through the renamed verb -------------------
test_that("random-intercept ti_lmm matches the frozen EMS engine numerically", {
  skip_if_not_installed("lme4")
  m <- fit_52()
  direct <- suppressWarnings(.ti_ems(m, level = 0.95, conf = 0.90))   # frozen EMS
  # fit_52() is intercept-only one-way, whose default TI is now Hoffman-Kringle
  # (test-hk.R); request the EMS engine to check the verb against its core.
  via    <- suppressWarnings(ti_lmm(m, P = 0.95, level = 0.90, method = "ems-mls"))
  expect_equal(via$estimate, direct$estimate)
  expect_equal(via$lower, direct$lower)
  expect_equal(via$upper, direct$upper)
  expect_equal(via$P, direct$level)      # content P (unified) == EMS `level` arg
  expect_equal(via$level, direct$conf)   # confidence (unified) == EMS `conf`
  # harmonized provenance
  expect_identical(via$method, "ems-mls")
  expect_identical(via$design, "random_intercept")
  expect_identical(attr(via, "grouping"), "oneway")
  expect_identical(ti_lmm(m, P = 0.95, level = 0.90)$method, "anova-mls")   # default
})

test_that("the vc-extract random-slope stop() still EXISTS and fires direct; dispatch routes around it", {
  skip_if_not_installed("lme4")
  fs <- slope_fit()
  expect_error(vc_extract(fs), "random-intercept")
  expect_error(.ti_ems(fs), "random-intercept")
  expect_error(ti_lmm(fs, newdata = data.frame(Days = 0),          # does not error
                      control = lmm_interval_control(seed = 1, M = 2000)), NA)  # small M: not an oracle
})

# ---- (3) random slope -> GPQ for every target / side / fit --------------------
# These are routing tests, not oracles: M = 2000 keeps them fast (the validated
# default is M = 20000).
test_that("every random-slope TI is GPQ: both targets x all sides, non-singular and singular", {
  skip_if_not_installed("lme4")
  ctl <- lmm_interval_control(seed = 3, M = 2000L)        # routing test, not an oracle
  set.seed(11); d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:6))
  d$y <- 100 - 2 * d$time + rnorm(nrow(d))                 # no between-batch variation
  fz <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | batch), data = d)))
  fits <- list(nonsingular = slope_fit(), singular = fz)
  nds  <- list(nonsingular = data.frame(Days = c(0, 5)), singular = data.frame(time = c(0, 6)))
  for (f in names(fits)) for (ov in c("observation", "group_mean"))
    for (sd in c("two.sided", "lower", "upper")) {
      out <- ti_lmm(fits[[f]], newdata = nds[[f]], over = ov, sides = sd, control = ctl)
      info <- paste(f, ov, sd)
      expect_identical(out$method, c("gpq", "gpq"), info = info)
      expect_identical(attr(out, "diagnostics")$route_reason, c("auto", "auto"), info = info)
      expect_identical(out$sides, c(sd, sd), info = info)
      if (sd != "upper") expect_true(all(is.finite(out$lower)), info = info)
      if (sd != "lower") expect_true(all(is.finite(out$upper)), info = info)
    }
})

test_that("slope future-obs (non-singular) routes to GPQ, one row per t0", {
  skip_if_not_installed("lme4")
  fs <- slope_fit()
  out <- ti_lmm(fs, newdata = data.frame(Days = c(0, 5)), P = 0.99,
                control = lmm_interval_control(seed = 1, M = 2000L))   # not an oracle
  expect_equal(nrow(out), 2L)
  expect_equal(out$Days, c(0, 5))        # eval column named by the actual covariate
  expect_true(all(out$method == "gpq"))
  expect_true(all(out$design == "random_slope_correlated"))
  expect_true(all(is.finite(out$lower)) && all(out$upper > out$lower))
})

test_that("slope batch-mean (over = group_mean) routes to GPQ even when non-singular", {
  skip_if_not_installed("lme4")
  out <- ti_lmm(slope_fit(), newdata = data.frame(Days = 5), over = "group_mean",
                control = lmm_interval_control(seed = 11, M = 2000L))
  expect_equal(out$method, "gpq")
  expect_equal(out$design, "random_slope_correlated")
  expect_true(is.finite(out$lower) && out$upper > out$lower)
})

test_that("method override: gpq is forceable; reml-mls is defunct; EMS methods are rejected", {
  skip_if_not_installed("lme4")
  fs <- slope_fit()
  ctl <- lmm_interval_control(seed = 7, M = 2000L)          # not an oracle
  forced <- ti_lmm(fs, newdata = data.frame(Days = 5), method = "gpq", control = ctl)
  auto   <- ti_lmm(fs, newdata = data.frame(Days = 5), control = ctl)
  expect_equal(forced$method, "gpq")
  expect_identical(attr(forced, "diagnostics")$route_reason, "method")
  expect_equal(c(forced$lower, forced$upper), c(auto$lower, auto$upper))  # same engine, same seed
  for (ov in c("observation", "group_mean")) {
    expect_error(ti_lmm(fs, newdata = data.frame(Days = 5), over = ov, method = "reml-mls"),
                 class = "reintervals_error_defunct")
  }
  expect_error(ti_lmm(fit_52(), method = "reml-mls"), class = "reintervals_error_defunct")
  # an EMS-only algorithm on a slope model is rejected
  expect_error(ti_lmm(fs, newdata = data.frame(Days = 5), method = "ems-mls"),
               "not available")
})

test_that("GPQ draws are reproducible through ti_lmm (control seed honoured)", {
  skip_if_not_installed("lme4")
  fs <- slope_fit()
  a <- ti_lmm(fs, newdata = data.frame(Days = 5), method = "gpq",
              control = lmm_interval_control(seed = 99, M = 2000L))
  b <- ti_lmm(fs, newdata = data.frame(Days = 5), method = "gpq",
              control = lmm_interval_control(seed = 99, M = 2000L))
  expect_equal(a$lower, b$lower); expect_equal(a$upper, b$upper)
})

# ---- (4) M1-overlap routing ---------------------------------------------------
test_that("clean balanced M1 routes to EMS; numeric core matches the EMS engine", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 8, times = c(0, 3, 6, 9, 12), beta0 = 100, beta1 = -0.25,
                s2_0 = 2.0, s2_1 = 0, s01 = 0, s2e = 0.7, seed = 11)
  fb <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch), data = d,
                                    REML = TRUE, control = lc))
  expect_true(.m1_is_clean_balanced(fb, "time"))
  via    <- suppressWarnings(ti_lmm(fb, newdata = data.frame(time = 12), P = 0.99, level = 0.95))
  direct <- suppressWarnings(.ti_ems(fb, newdata = data.frame(time = 12), level = 0.99, conf = 0.95))
  expect_equal(via$lower, direct$lower); expect_equal(via$upper, direct$upper)
  expect_identical(via$method, "ems-mls")
  expect_identical(via$design, "random_intercept_fixed_slope")     # route-independent
  expect_identical(attr(via, "grouping"), "oneway")
})

test_that("M1 single-obs (MBL Appendix A) routes to the anova-mls oracle, not the EMS/REML drift", {
  skip_if_not_installed("lme4")
  d  <- mbl_appendix_a
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = d,
                                    REML = TRUE, control = lc))
  expect_false(.m1_is_clean_balanced(fm, "Month"))
  out <- ti_lmm(fm, newdata = data.frame(Month = 12), P = 0.99, level = 0.95)
  expect_equal(out$method, "anova-mls")
  expect_equal(out$design, "random_intercept_fixed_slope")
  expect_equal(out$lower, 57.25, tolerance = 0.01)   # published Montes oracle
  expect_equal(out$upper, 75.41, tolerance = 0.01)
  expect_false(isTRUE(all.equal(out$lower, 57.39, tolerance = 1e-3)))  # not EMS drift
  expect_false(isTRUE(all.equal(out$lower, 56.92, tolerance = 1e-3)))  # not REML
})

test_that("batch-mean tolerance (over = group_mean) is DNE for M1; error names the alternative", {
  skip_if_not_installed("lme4")
  d  <- mbl_appendix_a
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = d, control = lc))
  expect_error(
    ti_lmm(fm, newdata = data.frame(Month = 12), over = "group_mean"),
    "not implemented|new_group_mean_lmm"
  )
})

# ---- (5) the two df fallbacks are NOT unified ---------------------------------
test_that("EMS and slope df fallbacks are distinct functions (not shared)", {
  expect_true(is.function(.ddf_containment))         # EMS fallback
  expect_true(is.function(.vfun_df_containment))     # slope fallback
  expect_false(identical(.ddf_containment, .vfun_df_containment))
  parsed <- list(B = 11, k_re = 1, n = 26, p = 2)
  expect_identical(.vfun_df_containment(parsed, w_between = 1),   max(1, 11 - 1))
  expect_identical(.vfun_df_containment(parsed, w_between = 0), max(1, 26 - 2 - 11))
})

test_that("a non-singular random-slope fit with an unavailable Satterthwaite df: TI is GPQ, slope NGM reports containment", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # the fit sits near the variance-component boundary (the
                   # mechanism behind the df failure), so its singular /
                   # fallback status could flip with platform-level lme4 jitter
  # One of the 361 E3 fits (paper repo data/defects/e3_fallback_fits.rds): sigma1
  # = 0.10, 20 batches, rep 1084, seed = MASTER_SEED + 100000 * k + rep with k = 3
  # (the coverage study's DGP, scripts/coverage_reml.R). Non-singular by
  # lme4::isSingular(), yet the Hessian of Cov(phi) is unusable, so the
  # Satterthwaite df of nu_G falls back to containment.
  set.seed(20261026L + 100000L * 3L + 1084L)
  d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:20))
  b0 <- rnorm(20, 0, 3); b1 <- rnorm(20, 0, 0.10)
  d$y <- 100 - 2 * d$time + b0[d$batch] + b1[d$batch] * d$time + rnorm(nrow(d), 0, 2)
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | batch), data = d)))
  expect_false(lme4::isSingular(fm))
  expect_true(attr(vfun_extract(fm, t0 = 12), "df_fallback")[["G"]])

  nd <- data.frame(time = c(6, 12)); ctl <- lmm_interval_control(seed = 1, M = 2000)
  x <- ti_lmm(fm, nd, P = 0.99, control = ctl)              # not an oracle
  expect_identical(x$method, c("gpq", "gpq"))
  expect_identical(attr(x, "diagnostics")$route_reason, c("auto", "auto"))
  ngm <- new_group_mean_lmm(fm, nd)
  expect_identical(ngm$method, c("reml", "reml"))
  expect_identical(attr(ngm, "diagnostics")$df_rule, c("containment", "containment"))
})

test_that("route_reason is 'auto' for every dispatched random-slope TI and 'method' when gpq is forced", {
  skip_if_not_installed("lme4")
  fs <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  ctl <- lmm_interval_control(seed = 1, M = 2000)           # not an oracle
  x  <- ti_lmm(fs, newdata = data.frame(Days = 5), control = ctl)
  gm <- ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean", control = ctl)
  fo <- ti_lmm(fs, newdata = data.frame(Days = 5), method = "gpq", control = ctl)
  expect_identical(attr(x,  "diagnostics")$route_reason, "auto")
  expect_identical(attr(gm, "diagnostics")$route_reason, "auto")
  expect_identical(attr(fo, "diagnostics")$route_reason, "method")
})

test_that(".is_pd_cov accepts a PD covariance and rejects indefinite / singular / malformed ones", {
  expect_true(.is_pd_cov(matrix(c(2, 0.5, 0.5, 1), 2)))                  # PD
  expect_false(.is_pd_cov(matrix(c(1, 2, 2, 1), 2)))                     # eigen 3, -1
  expect_false(.is_pd_cov(matrix(c(1, 1, 1, 1), 2)))                     # singular
  expect_false(.is_pd_cov(matrix(c(1, 0, 0, 1e-12), 2)))                 # below 1e-8 * max
  expect_true(.is_pd_cov(matrix(c(1, 0, 0, 1e-6), 2)))                   # above it
  expect_false(.is_pd_cov(matrix(c(2, 0.5, 0.4, 1), 2)))                 # not symmetric
  expect_false(.is_pd_cov(NULL))
  expect_false(.is_pd_cov(matrix(c(1, NA, NA, 1), 2)))
})

test_that("an indefinite Cov(phi) is unavailable: containment df at every t0 (TI is GPQ regardless)", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # near-boundary fit (the mechanism): status could flip with lme4 jitter
  # One of the 7 E3 fits whose Satterthwaite df failed at one t0 only, because
  # the inverse Hessian is indefinite (smallest eigenvalue ~ -1e-9): sigma1 =
  # 0.05, 20 batches, rep 229, seed = MASTER_SEED + 100000 * k + rep with k = 2
  # (the coverage study's DGP, scripts/coverage_reml.R).
  set.seed(20261026L + 100000L * 2L + 229L)
  d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:20))
  b0 <- rnorm(20, 0, 3); b1 <- rnorm(20, 0, 0.05)
  d$y <- 100 - 2 * d$time + b0[d$batch] + b1[d$batch] * d$time + rnorm(nrow(d), 0, 2)
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | batch), data = d)))
  expect_false(lme4::isSingular(fm))
  p <- .vfun_parse(fm)
  H <- numDeriv::hessian(p$reml_nll, p$phi)
  expect_false(.is_pd_cov(solve(H)))                 # the raw inverse is indefinite
  for (t0 in c(6, 12)) {
    cp <- vfun_extract(fm, t0 = t0)
    expect_null(cp$cov)
    expect_true(attr(cp, "df_fallback")[["G"]], info = t0)
  }
  x <- ti_lmm(fm, data.frame(time = c(6, 12)), P = 0.99,
              control = lmm_interval_control(seed = 1, M = 2000))   # not an oracle
  expect_identical(x$method, c("gpq", "gpq"))
  expect_identical(attr(x, "diagnostics")$route_reason, c("auto", "auto"))
})
