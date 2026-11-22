# Singular fits never get a variance-parameter covariance (both engines), and
# why a covariance / a df is unavailable is reported per row in
# attr(., "diagnostics"): vc_cov_reason (random slope, next to vc_cov_route) and
# df_fallback_reason (random intercept).

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))
quiet <- function(expr) suppressWarnings(expr)
dg <- function(x) attr(x, "diagnostics")
bi_ctl <- lmm_interval_control(slope_ngm_method = "broad-inference")
ri_bi_ctl <- lmm_interval_control(ngm_method = "broad-inference")

## the singular || fixture of test-vfun-extract.R / test-a5-warnings.R
sing_rs <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  lmer_q(y ~ time + (1 | batch) + (0 + time | batch), dn, control = lc)
}
rs_nd <- data.frame(time = c(6, 12))
## every random-slope verb / arm that carries vc_cov_route
rs_calls <- function(fm, nd) list(
  ci  = quiet(ci_lmm(fm, nd)),
  pi  = quiet(pi_lmm(fm, nd)),
  ngm = quiet(new_group_mean_lmm(fm, nd)),
  bi  = quiet(new_group_mean_lmm(fm, nd, control = bi_ctl)),
  lp_ci  = quiet(lmm_predict(fm, "confidence", nd)),
  lp_pi  = quiet(lmm_predict(fm, "prediction", nd)),
  lp_ngm = quiet(lmm_predict(fm, "new_group_mean", nd)))
expect_reason <- function(calls, route, reason) {
  for (v in names(calls)) {
    d <- dg(calls[[v]]); n <- nrow(calls[[v]])
    expect_identical(d$vc_cov_route, rep(route, n), info = v)
    expect_identical(d$vc_cov_reason, rep(reason, n), info = v)
  }
}

test_that("a singular || random-slope fit: no covariance, reason 'singular fit', containment / NA limits", {
  skip_if_not_installed("lme4")
  fm <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fm)), "fit did not land singular in this environment")
  cp <- vfun_extract(fm, t0 = 12)
  expect_null(cp$cov)
  expect_identical(attr(cp, "vc_cov_route"), "unavailable")
  expect_identical(attr(cp, "vc_cov_reason"), "singular fit")
  expect_identical(unname(attr(cp, "df_fallback")), rep(TRUE, 3))
  r <- rs_calls(fm, rs_nd)
  expect_reason(r, "unavailable", "singular fit")
  for (v in c("ci", "pi", "ngm", "lp_ci", "lp_pi", "lp_ngm")) {
    expect_identical(dg(r[[v]])$df_rule, rep("containment", 2), info = v)
    expect_true(all(is.finite(c(r[[v]]$lower, r[[v]]$upper))), info = v)
  }
  ## the containment dfs: B - 2 between-lot pieces, as .vfun_df_containment
  expect_identical(dg(r$ngm)$df, rep(10 - 2, 2))
  expect_identical(dg(r$ci)$df, rep(10 - 2, 2))
  ## broad-inference: no df, no fallback -> NA limits, ngm_df_failed
  expect_true(all(dg(r$bi)$ngm_df_failed))
  expect_true(all(is.na(r$bi$lower) & is.na(r$bi$upper)))
  expect_true(all(is.na(dg(r$bi)$df)))
  ## the df-fallback warning names the reason
  expect_warning(ci_lmm(fm, rs_nd), "(singular fit)", fixed = TRUE,
                 class = "reintervals_warning_df_fallback")
})

test_that("a singular random-slope fit evaluates no Hessian (numDeriv / .nd_deriv mocked to error)", {
  skip_if_not_installed("lme4")
  fm <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fm)), "fit did not land singular in this environment")
  ref <- rs_calls(fm, rs_nd)
  n_hess <- 0L
  local_mocked_bindings(hessian = function(...) {
    n_hess <<- n_hess + 1L; stop("numDeriv::hessian must not be called")
  }, .package = "numDeriv")
  got <- rs_calls(fm, rs_nd)
  expect_identical(n_hess, 0L)
  for (v in names(ref)) {
    expect_identical(got[[v]]$lower, ref[[v]]$lower, info = v)
    expect_identical(got[[v]]$upper, ref[[v]]$upper, info = v)
    expect_identical(dg(got[[v]])$df, dg(ref[[v]])$df, info = v)
    expect_identical(dg(got[[v]])$df_rule, dg(ref[[v]])$df_rule, info = v)
  }
  ## the mock does intercept the Hessian: a non-singular fit calls it
  fs <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  quiet(ci_lmm(fs, data.frame(Days = 5)))
  expect_gt(n_hess, 0L)
})

test_that("a singular random-slope fit succeeds when every numerical derivative errors", {
  skip_if_not_installed("lme4")
  fm <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fm)), "fit did not land singular in this environment")
  ref <- rs_calls(fm, rs_nd)
  local_mocked_bindings(.nd_deriv = function(...) stop("no derivatives"))
  got <- rs_calls(fm, rs_nd)
  expect_reason(got, "unavailable", "singular fit")
  for (v in names(ref)) {
    expect_identical(got[[v]]$lower, ref[[v]]$lower, info = v)
    expect_identical(got[[v]]$upper, ref[[v]]$upper, info = v)
  }
})

test_that("vc_cov_reason is NA when a covariance was obtained (sleepstudy, default route)", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  expect_reason(rs_calls(fs, data.frame(Days = c(2, 5))), "default", NA_character_)
})

test_that("vc_cov_reason 'singularity check failed; treated as singular' (isSingular mocked to error)", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  local_mocked_bindings(isSingular = function(x, ...) stop("boom"), .package = "lme4")
  cp <- vfun_extract(fs, t0 = 5)
  expect_null(cp$cov)
  expect_identical(attr(cp, "vc_cov_reason"), "singularity check failed; treated as singular")
  r <- rs_calls(fs, data.frame(Days = c(2, 5)))
  expect_reason(r, "unavailable", "singularity check failed; treated as singular")
  expect_identical(dg(r$ngm)$df_rule, rep("containment", 2))
  expect_true(all(dg(r$bi)$ngm_df_failed))
})

test_that("vc_cov_reason 'Hessian not PD or not invertible' (mocked: no PD inverse; Hessian error)", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(2, 5))
  local({
    local_mocked_bindings(.is_pd_cov = function(M, tol = 1e-8) FALSE)
    r <- rs_calls(fs, nd)
    expect_reason(r, "unavailable", "Hessian not PD or not invertible")
    expect_identical(dg(r$pi)$df_rule, rep("containment", 2))
  })
  local({
    real <- .nd_deriv
    local_mocked_bindings(.nd_deriv = function(type, ...)
      if (type == "hessian") stop("no Hessian") else real(type, ...))
    cp <- vfun_extract(fs, t0 = 5)
    expect_null(cp$cov)
    expect_identical(attr(cp, "vc_cov_reason"), "Hessian not PD or not invertible")
  })
})

test_that("vc_cov_reason 'sentinel': default and retried Hessians both leave the valid region", {
  skip_if_not_installed("lme4")
  ## the |rho| > 0.9999 fixture of test-dispatch.R
  set.seed(13); K <- 12; tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- expand.grid(time = tt, lot = factor(1:K))
  z1 <- rnorm(K); z2 <- 0.9999 * z1 + sqrt(1 - 0.9999^2) * rnorm(K)
  d$y <- 100 + (2 * z1)[d$lot] + (-0.3 + (0.15 * z2)[d$lot]) * d$time + rnorm(nrow(d), 0, 1)
  fm <- lmer_q(y ~ time + (time | lot), d)
  ph <- .vfun_parse(fm)$phi
  rho <- ph[["s01"]] / sqrt(ph[["s2_0"]] * ph[["s2_1"]])
  skip_if(isTRUE(lme4::isSingular(fm)) || abs(rho) < 0.9999,
          "the optimum differs on this platform")
  r <- rs_calls(fm, data.frame(time = c(12, 24)))
  expect_reason(r, "unavailable", "sentinel")
  expect_true(all(dg(r$bi)$ngm_df_failed))
})

test_that("vc_cov_reason 'unstable numerical Hessian' (rep-1084 boundary-adjacent fit)", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # boundary-adjacent: see test-dispatch.R
  set.seed(20261026L + 100000L * 3L + 1084L)
  d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:20))
  b0 <- rnorm(20, 0, 3); b1 <- rnorm(20, 0, 0.10)
  d$y <- 100 - 2 * d$time + b0[d$batch] + b1[d$batch] * d$time + rnorm(nrow(d), 0, 2)
  fm <- lmer_q(y ~ time + (time | batch), d)
  skip_if(isTRUE(lme4::isSingular(fm)), "the fit is singular on this platform")
  r <- rs_calls(fm, data.frame(time = c(6, 12)))
  expect_reason(r, "unavailable", "unstable numerical Hessian")
})

## ---- random intercept: df_fallback_reason per row --------------------------
## identical run means -> zero between-run variance (test-a5-warnings.R)
sing_ri <- function() {
  d <- sqrt(0.001253)
  y <- as.vector(sapply(rep(0.981, 6), function(m) m + c(-d, 0, d)))
  lmer_q(y ~ 1 + (1 | run), data.frame(y = y, run = factor(rep(1:6, each = 3))))
}
## the per-row column against the engine's own attr(., "df_fallback_reason")
expect_row_reason <- function(model, kind, newdata = NULL, expected) {
  native <- quiet(.interval_lmm(model, newdata = newdata, kind = kind))
  out <- quiet(switch(kind, CI = ci_lmm(model, newdata),
                      PI = pi_lmm(model, newdata),
                      new_group_mean = new_group_mean_lmm(model, newdata)))
  lp <- quiet(lmm_predict(model, switch(kind, CI = "confidence", PI = "prediction",
                                        new_group_mean = "new_group_mean"), newdata))
  r <- dg(out)$df_fallback_reason
  expect_identical(r, expected, info = kind)
  expect_identical(dg(lp)$df_fallback_reason, expected, info = kind)
  u <- unique(stats::na.omit(r))
  if (length(u)) expect_identical(as.vector(u), attr(native, "df_fallback_reason"), info = kind)
  else expect_null(attr(native, "df_fallback_reason"), info = kind)
}

test_that("random-intercept df_fallback_reason in diagnostics equals attr(., 'df_fallback_reason')", {
  skip_if_not_installed("lme4")
  m <- sing_ri()
  expect_true(lme4::isSingular(m))
  expect_row_reason(m, "CI", expected = "singular fit")
  expect_row_reason(m, "PI", expected = "singular fit")
  ## it is also the reason the df-fallback warning reports
  w <- tryCatch(ci_lmm(m), reintervals_warning_df_fallback = function(w) w)
  expect_identical(w$reason, dg(quiet(ci_lmm(m)))$df_fallback_reason)
  ## a healthy fit: NA on every row, no attribute
  fs <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  for (k in c("CI", "PI", "new_group_mean"))
    expect_row_reason(fs, k, nd, expected = rep(NA_character_, 2))
})

test_that("random-intercept df_fallback_reason: mocked 'df below 1' (CI, PI) and the crossed NGM", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  local({
    local_mocked_bindings(.ci_ddf_satterthwaite = function(model, Lrow, theta) 0.004)
    expect_row_reason(fs, "CI", nd, expected = rep("df below 1", 2))
  })
  local({
    local_mocked_bindings(.pi_df_hessian = function(sigmaT2, vS) 0.004)
    expect_row_reason(fs, "PI", nd, expected = rep("df below 1", 2))
  })
  ## crossed new-group mean (test-a5-warnings.R): k_A = k_B = 0 -> containment
  set.seed(1); A <- 3; B <- 4; n <- 2
  g <- expand.grid(r = 1:n, a = factor(1:A), b = factor(1:B)); g$ab <- interaction(g$a, g$b)
  g$y <- 10 + rnorm(A, 0, 2)[g$a] + rnorm(B, 0, 2)[g$b] + rnorm(A * B, 0, 1)[g$ab] + rnorm(nrow(g))
  fc <- lmer_q(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), g)
  expect_row_reason(fc, "new_group_mean", expected = NA_character_)
  real <- .design_components_fit
  local_mocked_bindings(.design_components_fit = function(...) {
    dc <- real(...); dc$ems$k[1:2] <- 0; dc
  })
  expect_row_reason(fc, "new_group_mean",
                    expected = "k_A EMS_A + k_B EMS_B is not positive")
})

test_that("random-intercept broad-inference arm: df_fallback_reason is why the df could not be formed", {
  skip_if_not_installed("lme4")
  m <- sing_ri()
  bi <- quiet(new_group_mean_lmm(m, control = ri_bi_ctl))
  expect_true(all(dg(bi)$ngm_df_failed))
  expect_true(all(is.na(bi$lower) & is.na(bi$upper)))
  expect_identical(dg(bi)$df_fallback_reason, "singular fit")
  ## formed df: NA reason
  fs <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  ok <- quiet(new_group_mean_lmm(fs, nd, control = ri_bi_ctl))
  expect_false(any(dg(ok)$ngm_df_failed))
  expect_identical(dg(ok)$df_fallback_reason, rep(NA_character_, 2))
  ## a failed gradient
  local_mocked_bindings(.contrast_var_grad = function(...) NULL)
  ng <- quiet(new_group_mean_lmm(fs, nd, control = ri_bi_ctl))
  expect_true(all(dg(ng)$ngm_df_failed))
  expect_identical(dg(ng)$df_fallback_reason, rep("numerical derivative failed", 2))
})
