# Random-slope and GPQ diagnostics added in the 1.0.0 cleanup:
#   * df_fallback_reason on the random-slope CI rows (every slope_ci_df arm)
#     and new-group-mean rows (every slope_ngm_method arm), with the semantics
#     of the random-slope PI rows (test-slope-pi-method.R): why the default df
#     fell back to containment, NA where it did not and on arms with no
#     fallback;
#   * M (the number of Monte-Carlo draws) on every GPQ row: the random-slope
#     tolerance interval, the GPQ new-group-mean arm and ti_gpq_raw() (the
#     random-intercept GPQ route already reported it).
# Every other output was checked identical() to the previous commit when these
# fields were added.

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))
dg <- function(x) attr(x, "diagnostics")
## muffle the warnings these calls signal by design; any other surfaces
quiet_rs <- function(expr) {
  withCallingHandlers(expr,
    reintervals_warning_unvalidated   = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_limitation    = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"))
}
## the singular || fixture of test-slope-pi-method.R
sing_rs <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  lmer_q(y ~ time + (1 | batch) + (0 + time | batch), dn, control = lc)
}
sing_nd <- data.frame(time = c(0, 6, 12, 24))
ss_fit <- function() lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
ss_nd <- data.frame(Days = c(0, 4.5, 9))
ci_arm <- function(a) lmm_interval_control(slope_ci_df = a)
ngm_arm <- function(a) lmm_interval_control(slope_ngm_method = a, seed = 1, M = 2000)

test_that("random-slope CI rows carry df_fallback_reason on every slope_ci_df arm", {
  skip_if_not_installed("lme4")
  fs <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fs)), "fixture fit is not singular here")
  ## default: the singular fit's df falls back, with its reason
  r <- quiet_rs(ci_lmm(fs, sing_nd))
  expect_identical(dg(r)$df_rule, rep("containment", 4))
  expect_identical(dg(r)$df_fallback_reason, rep("singular fit", 4))
  ## the reason is the one the df-fallback warning reports
  w <- tryCatch(ci_lmm(fs, sing_nd), reintervals_warning_df_fallback = function(w) w)
  expect_identical(w$reason, "singular fit")
  ## the arms have no fallback: NA
  expect_identical(dg(quiet_rs(ci_lmm(fs, sing_nd, control = ci_arm("k-minus-1"))))$df_fallback_reason,
                   rep(NA_character_, 4))
  skip_if_not_installed("pbkrtest")
  expect_identical(dg(quiet_rs(ci_lmm(fs, sing_nd, control = ci_arm("kenward-roger"))))$df_fallback_reason,
                   rep(NA_character_, 4))
  ## interior fit: NA on every arm
  fm <- ss_fit()
  for (a in c("hessian-satterthwaite", "kenward-roger", "k-minus-1"))
    expect_identical(dg(quiet_rs(ci_lmm(fm, ss_nd, control = ci_arm(a))))$df_fallback_reason,
                     rep(NA_character_, 3), info = a)
})

test_that("random-slope new-group-mean rows carry df_fallback_reason on every slope_ngm_method arm", {
  skip_if_not_installed("lme4")
  fs <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fs)), "fixture fit is not singular here")
  r <- quiet_rs(new_group_mean_lmm(fs, sing_nd))
  expect_identical(dg(r)$df_rule, rep("containment", 4))
  expect_identical(dg(r)$df_fallback_reason, rep("singular fit", 4))
  for (a in c("broad-inference", "gpq", "fv-anova"))
    expect_identical(dg(quiet_rs(new_group_mean_lmm(fs, sing_nd, control = ngm_arm(a))))$df_fallback_reason,
                     rep(NA_character_, 4), info = a)
  fm <- ss_fit()
  for (a in c("reml", "broad-inference", "gpq", "fv-anova"))
    expect_identical(dg(quiet_rs(new_group_mean_lmm(fm, ss_nd, control = ngm_arm(a))))$df_fallback_reason,
                     rep(NA_character_, 3), info = a)
  ## vectorised: row-aligned (4 t0 x 2 sides x 2 levels), GPQ arm included
  v <- quiet_rs(new_group_mean_lmm(fs, sing_nd, level = c(0.9, 0.95), sides = c("two.sided", "lower")))
  expect_identical(dg(v)$df_fallback_reason, rep("singular fit", 16))
  g <- quiet_rs(new_group_mean_lmm(fs, sing_nd, level = c(0.9, 0.95), sides = c("two.sided", "lower"),
                                   control = ngm_arm("gpq")))
  expect_identical(dg(g)$df_fallback_reason, rep(NA_character_, 16))
})

test_that("random-slope CI and NGM df_fallback_reason: the other reasons", {
  skip_if_not_installed("lme4")
  fm <- ss_fit()
  reason <- function(v) dg(quiet_rs(v(fm, ss_nd)))$df_fallback_reason
  ## a non-singular fit without a covariance: the covariance reason
  local({
    local_mocked_bindings(.is_pd_cov = function(M, tol = 1e-8) FALSE)
    expect_identical(reason(ci_lmm), rep("Hessian not PD or not invertible", 3))
    expect_identical(reason(new_group_mean_lmm), rep("Hessian not PD or not invertible", 3))
  })
  ## a failed singularity check, read as singular
  local({
    local_mocked_bindings(.is_singular = function(model) NA)
    expect_identical(reason(ci_lmm), rep("singularity check failed; treated as singular", 3))
    expect_identical(reason(new_group_mean_lmm),
                     rep("singularity check failed; treated as singular", 3))
  })
  ## a Satterthwaite df below 1, and a non-finite one
  local({
    local_mocked_bindings(.vfun_df_satterthwaite = function(value, grad_vec, vc_cov) 0.5)
    expect_identical(reason(ci_lmm), rep("df below 1", 3))
    expect_identical(reason(new_group_mean_lmm), rep("df below 1", 3))
  })
  local({
    local_mocked_bindings(.vfun_df_satterthwaite = function(value, grad_vec, vc_cov) Inf)
    expect_identical(reason(ci_lmm), rep("non-finite df", 3))
    expect_identical(reason(new_group_mean_lmm), rep("non-finite df", 3))
  })
})

test_that("M is in the diagnostics of every GPQ row", {
  skip_if_not_installed("lme4")
  fm <- ss_fit()
  ctl <- lmm_interval_control(seed = 1, M = 2000)
  for (ov in c("observation", "group_mean")) {
    t <- quiet_rs(ti_lmm(fm, ss_nd, over = ov, control = ctl))
    expect_identical(dg(t)$M, 2000L, info = ov)
  }
  expect_identical(dg(quiet_rs(ti_lmm(fm, ss_nd, control = lmm_interval_control(seed = 1, M = 3000))))$M,
                   3000L)
  g <- quiet_rs(new_group_mean_lmm(fm, ss_nd, control = ngm_arm("gpq")))
  expect_identical(dg(g)$M, 2000L)
  d <- lme4::sleepstudy
  r <- quiet_rs(ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = c(2, 5), control = ctl))
  expect_identical(dg(r)$M, 2000L)
  ## the random-intercept GPQ route (unchanged)
  fi <- lmer_q(Reaction ~ Days + (1 | Subject), d)
  ri <- quiet_rs(ti_lmm(fi, data.frame(Days = 0), method = "gpq", control = ctl))
  expect_identical(dg(ri)$M, 2000L)
  ## non-GPQ rows do not carry it
  expect_null(dg(quiet_rs(new_group_mean_lmm(fm, ss_nd)))$M)
})
