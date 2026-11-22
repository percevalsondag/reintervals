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

test_that("a boundary-adjacent non-singular random-slope fit: the half-step stability guard rejects Cov(phi), TI is GPQ", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # the fit sits near the variance-component boundary (the
                   # mechanism behind the df failure), so its singular /
                   # fallback status could flip with platform-level lme4 jitter
  # One of the 361 E3 fits (paper repo data/defects/e3_fallback_fits.rds): sigma1
  # = 0.10, 20 batches, rep 1084, seed = MASTER_SEED + 100000 * k + rep with k = 3
  # (the coverage study's DGP, scripts/coverage_reml.R). Non-singular by
  # lme4::isSingular() (the package's singular criterion, kept as is; the fit is
  # boundary-adjacent, delta = -9.9e-6), and the default Hessian of Cov(phi) is
  # unusable (its first 10% step leaves the valid region). The d = 1e-3 retry
  # gives a PD Cov(phi), but s2_1 / s2e = 4.2e-6 is below numDeriv's default
  # zero.tol, so the half-step stability guard (.nd_unstable_psi) compares the
  # F, T, G dfs at the call's t0 from the d = 1e-3 and d = 5e-4 Hessians: at
  # t0 = 12, nu_F 19.13 vs 19.08 and nu_T 24.19 vs 24.12 (0.26% and 0.33%,
  # above the 1e-3 tolerance; nu_G at t0 = 12 is 12.16 at d = 0.05 and 11.87
  # at d = 2.5e-4), so Cov(phi) is unavailable ("unstable numerical Hessian")
  # and nu_G falls back to containment. (Before the scale-free derivatives the
  # retry was accepted here only because s2_1 = 1.83e-5 sat just above
  # zero.tol in these units; at y x 0.9 it was not.)
  set.seed(20261026L + 100000L * 3L + 1084L)
  d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:20))
  b0 <- rnorm(20, 0, 3); b1 <- rnorm(20, 0, 0.10)
  d$y <- 100 - 2 * d$time + b0[d$batch] + b1[d$batch] * d$time + rnorm(nrow(d), 0, 2)
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | batch), data = d)))
  expect_false(lme4::isSingular(fm))
  p <- .vfun_parse(fm)
  expect_true(.nd_unstable_psi(p$phi))
  nd <- data.frame(time = c(6, 12))
  cp <- vfun_extract(fm, t0 = 12, guard_t0 = nd$time)
  expect_identical(attr(cp, "vc_cov_route"), "unavailable")
  expect_identical(attr(cp, "vc_cov_reason"), "unstable numerical Hessian")
  expect_null(cp$cov)
  expect_true(attr(cp, "df_fallback")[["G"]])

  ctl <- lmm_interval_control(seed = 1, M = 2000)
  x <- ti_lmm(fm, nd, P = 0.99, control = ctl)              # not an oracle
  expect_identical(x$method, c("gpq", "gpq"))
  expect_identical(attr(x, "diagnostics")$route_reason, c("auto", "auto"))
  expect_warning(ngm <- new_group_mean_lmm(fm, nd), "unstable numerical Hessian",
                 class = "reintervals_warning_df_fallback")
  expect_identical(ngm$method, c("reml", "reml"))
  expect_identical(attr(ngm, "diagnostics")$df_rule, c("containment", "containment"))
  expect_identical(attr(ngm, "diagnostics")$vc_cov_route, c("unavailable", "unavailable"))
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
  ## judged on the correlation form, so scale-free: a tiny variance of a
  ## parameter in fine units is not "numerically singular"
  expect_true(.is_pd_cov(matrix(c(1, 0, 0, 1e-12), 2)))
  expect_true(.is_pd_cov(matrix(c(1, 0, 0, 1e-6), 2)))
  r <- 1 - 1e-10                                                         # corr eigen 1e-10
  expect_false(.is_pd_cov(matrix(c(1, r, r, 1), 2)))                     # below 1e-8 * max
  expect_false(.is_pd_cov(matrix(c(1, r, r, 1) * c(1, 1e-3, 1e-3, 1e-6), 2)))  # same, rescaled
  r <- 1 - 1e-6
  expect_true(.is_pd_cov(matrix(c(1, r, r, 1), 2)))                      # above it
  expect_false(.is_pd_cov(matrix(c(1, 0, 0, 0), 2)))                     # zero variance
  expect_false(.is_pd_cov(matrix(c(2, 0.5, 0.4, 1), 2)))                 # not symmetric
  expect_false(.is_pd_cov(NULL))
  expect_false(.is_pd_cov(matrix(c(1, NA, NA, 1), 2)))
})

test_that("an indefinite default Cov(phi) is recovered by the Hessian retry (TI is GPQ regardless)", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # near-boundary fit (the mechanism): status could flip with lme4 jitter
  # One of the 7 E3 fits whose Satterthwaite df failed at one t0 only, because
  # the inverse Hessian is indefinite (smallest eigenvalue ~ -1e-9): sigma1 =
  # 0.05, 20 batches, rep 229, seed = MASTER_SEED + 100000 * k + rep with k = 2
  # (the coverage study's DGP, scripts/coverage_reml.R). The fit is interior
  # (UN, rho = -0.95): the default Hessian is contaminated by the reml_nll
  # sentinel at its 10% first step; the d = 1e-3 retry recovers Cov(phi). The
  # broad-inference df is checked against the step-13 diagnostic (two
  # independent Hessians: reml_nll at d = 1e-3 and lme4's theta-scale Hessian by
  # the delta method, agreeing to 3 digits): 16.7 at t0 = 6, 11.3 at t0 = 12.
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
    expect_false(is.null(cp$cov), info = t0)
    expect_true(.is_pd_cov(cp$cov), info = t0)
    expect_identical(attr(cp, "vc_cov_route"), "retry-d1e-3", info = t0)
    expect_false(attr(cp, "df_fallback")[["G"]], info = t0)
  }
  bi <- suppressWarnings(new_group_mean_lmm(fm, data.frame(time = c(6, 12)),
          control = lmm_interval_control(slope_ngm_method = "broad-inference")))
  db <- attr(bi, "diagnostics")
  expect_identical(db$vc_cov_route, c("retry-d1e-3", "retry-d1e-3"))
  expect_false(any(db$ngm_df_failed))
  expect_equal(db$df, c(16.7, 11.3), tolerance = 1e-2)
  x <- ti_lmm(fm, data.frame(time = c(6, 12)), P = 0.99,
              control = lmm_interval_control(seed = 1, M = 2000))   # not an oracle
  expect_identical(x$method, c("gpq", "gpq"))
  expect_identical(attr(x, "diagnostics")$route_reason, c("auto", "auto"))
})

## ---- Hessian retry (.vfun_vc_cov): sentinel contamination -------------------
## reml_nll evaluated at phi, recording whether any evaluation of a numDeriv
## Hessian at method.args `ma` returned the 1e10 sentinel
sentinel_hit <- function(p, ma = list()) {
  hit <- FALSE
  f <- function(q) { v <- p$reml_nll(q); if (!is.finite(v) || v >= 1e10) hit <<- TRUE; v }
  H <- numDeriv::hessian(f, p$phi, method.args = ma)
  list(hit = hit, H = H)
}

test_that("a constructed UN random-slope fit with |rho| in 0.90-0.98: default Hessian hits the sentinel, the retry recovers Cov(phi)", {
  skip_if_not_installed("lme4")
  set.seed(2); K <- 12; tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- expand.grid(time = tt, lot = factor(1:K))
  z1 <- rnorm(K); z2 <- 0.95 * z1 + sqrt(1 - 0.95^2) * rnorm(K)
  d$y <- 100 + (2 * z1)[d$lot] + (-0.3 + (0.15 * z2)[d$lot]) * d$time + rnorm(nrow(d), 0, 1)
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | lot), data = d)))
  p <- .vfun_parse(fm); ph <- p$phi
  rho <- ph[["s01"]] / sqrt(ph[["s2_0"]] * ph[["s2_1"]])
  skip_if(isTRUE(lme4::isSingular(fm)) || abs(rho) < 0.90 || abs(rho) > 0.98,
          "the optimum differs on this platform")
  ## the default Hessian evaluates reml_nll outside the PD cone of Sigma
  def <- sentinel_hit(p)
  skip_if_not(def$hit, "the default Hessian does not hit the sentinel on this platform")
  expect_true(def$hit)
  expect_false(sentinel_hit(p, list(d = 1e-3))$hit)     # the retry stays inside
  vc <- .vfun_vc_cov(p)
  expect_identical(vc$vc_cov_route, "retry-d1e-3")
  expect_true(.is_pd_cov(vc$cov))
  nd <- data.frame(time = c(12, 24))
  bi <- suppressWarnings(new_group_mean_lmm(fm, nd,
          control = lmm_interval_control(slope_ngm_method = "broad-inference")))
  expect_true(all(is.finite(attr(bi, "diagnostics")$df)))
  expect_false(any(attr(bi, "diagnostics")$ngm_df_failed))
  expect_true(all(is.finite(c(bi$lower, bi$upper))))
  ## default reml arm: Satterthwaite nu_G (no containment) when nu_G >= 1
  ngm <- quiet_df_fallback(new_group_mean_lmm(fm, nd))
  dg <- attr(ngm, "diagnostics")
  expect_identical(dg$vc_cov_route, c("retry-d1e-3", "retry-d1e-3"))
  expect_identical(dg$df_rule[dg$df >= 1], rep("satterthwaite", sum(dg$df >= 1)))
  expect_true(all(dg$df >= 1))
  ## ci_lmm and pi_lmm report the same route
  expect_identical(attr(ci_lmm(fm, nd), "diagnostics")$vc_cov_route, c("retry-d1e-3", "retry-d1e-3"))
  expect_identical(attr(pi_lmm(fm, nd), "diagnostics")$vc_cov_route, c("retry-d1e-3", "retry-d1e-3"))
})

test_that("when the retry also hits the sentinel, Cov(phi) stays unavailable and the fallbacks apply", {
  skip_if_not_installed("lme4")
  ## |rho| > 0.9999 on a fit lme4 calls non-singular: the PD margin of Sigma
  ## (1 - rho^2) is below even the retry's 0.1% step on s01
  set.seed(13); K <- 12; tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- expand.grid(time = tt, lot = factor(1:K))
  z1 <- rnorm(K); z2 <- 0.9999 * z1 + sqrt(1 - 0.9999^2) * rnorm(K)
  d$y <- 100 + (2 * z1)[d$lot] + (-0.3 + (0.15 * z2)[d$lot]) * d$time + rnorm(nrow(d), 0, 1)
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (time | lot), data = d)))
  p <- .vfun_parse(fm); ph <- p$phi
  rho <- ph[["s01"]] / sqrt(ph[["s2_0"]] * ph[["s2_1"]])
  skip_if(isTRUE(lme4::isSingular(fm)) || abs(rho) < 0.9999,
          "the optimum differs on this platform")
  expect_true(sentinel_hit(p)$hit)
  expect_true(sentinel_hit(p, list(d = 1e-3))$hit)
  vc <- .vfun_vc_cov(p)
  expect_null(vc$cov)
  expect_identical(vc$vc_cov_route, "unavailable")
  nd <- data.frame(time = c(12, 24))
  ngm <- quiet_df_fallback(new_group_mean_lmm(fm, nd))
  expect_identical(attr(ngm, "diagnostics")$df_rule, c("containment", "containment"))
  expect_identical(attr(ngm, "diagnostics")$vc_cov_route, c("unavailable", "unavailable"))
  expect_identical(attr(ngm, "diagnostics")$df, rep(K - 2, 2))
  expect_true(all(is.finite(c(ngm$lower, ngm$upper))))
  bi <- suppressWarnings(new_group_mean_lmm(fm, nd,
          control = lmm_interval_control(slope_ngm_method = "broad-inference")))
  expect_true(all(attr(bi, "diagnostics")$ngm_df_failed))
  expect_true(all(is.na(bi$lower) & is.na(bi$upper)))
  x <- ti_lmm(fm, nd, P = 0.99, control = lmm_interval_control(seed = 1, M = 2000))  # not an oracle
  expect_identical(x$method, c("gpq", "gpq"))
})

test_that("a fit whose default Hessian is valid keeps vc_cov_route = 'default' (sleepstudy pins)", {
  skip_if_not_installed("lme4")
  fm <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy))
  expect_identical(.vfun_vc_cov(.vfun_parse(fm))$vc_cov_route, "default")
  expect_identical(attr(new_group_mean_lmm(fm, data.frame(Days = 5)), "diagnostics")$vc_cov_route,
                   "default")
})
