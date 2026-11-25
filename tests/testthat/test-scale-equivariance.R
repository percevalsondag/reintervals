# Scale equivariance: for y -> c * y (c > 0) every interval and df the package
# returns must transform as the data do. Limits, estimates and SEs scale by c,
# variances (W, tau2_med, s2e_med) by c^2; dfs, rules, routes, methods and
# fallback flags are unchanged. Before the fix numDeriv's absolute step below
# zero.tol (~1.8e-5) made the derivatives (and the dfs built on them) depend on
# the measurement units once a variance component fell below ~1.8e-5; every
# derivative is now taken on psi = phi / s, s = the residual variance, with
# relative steps for every non-zero psi_j (.nd_deriv). Time units: for
# t -> c * t (newdata t0 rescaled the same way) every interval, SE, variance
# and df is unchanged.
#
# Tolerances (element-wise relative): 1e-8 for limits / estimates / SEs /
# variances on closed-form and Satterthwaite rows, 1e-10 on GPQ rows (same
# seed), 1e-6 for dfs and other derived numbers.
#
# The scaled fits are built at the exactly transformed theta with no
# optimization (lmerControl(optimizer = NULL)): for y -> c y the relative theta
# is unchanged; for t -> c t the slope rows of the relative Cholesky factor are
# divided by c. lme4's optimizer is not itself equivariant: refitting y x c
# moves the variance components by up to ~3e-8 relative (the fixed-slope
# fixture), and refitting t x 12 or t x 100 on sleepstudy (||) stops 2-16% away
# from the transformed optimum with convergence warnings. That is a property
# of the fit, not of the package's derived quantities, which are what is tested
# here. (The two regression tests below use genuine lme4 refits.)

y_scales <- c(1e-4, 1e-2, 1, 1e2, 1e4)
t_scales <- c(1e-2, 1 / 12, 1, 12, 1e2)

## columns (data or diagnostics) and their power of c under y -> c y
## (y / Reaction: the response column of the model-frame row used when
## newdata = NULL; U, Vhat_Y: the MLS / Graybill-Wang variance pieces)
sc_c1 <- c("estimate", "lower", "upper", "se", "mc_se_lower", "mc_se_upper",
           "y", "Reaction")
sc_c2 <- c("W", "tau2_med", "s2e_med", "U", "Vhat_Y")
## integer / bookkeeping numbers compared exactly
sc_exact <- c("level", "P", "gpq_seed", "mc_sections", "K", "df_resid", "n_na_draws")

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))

## element-wise relative agreement of x with r (matching Inf / NA count as equal)
rel_ok <- function(x, r, tol) {
  if (length(x) != length(r)) return(FALSE)
  same_inf <- is.infinite(x) & is.infinite(r) & sign(x) == sign(r)
  same_na <- is.na(x) & is.na(r)
  ok <- same_inf | same_na |
    (is.finite(x) & is.finite(r) & abs(x - r) <= tol * pmax(abs(r), .Machine$double.xmin))
  all(ok)
}

## compare two interval outputs: base (c = 1) and scaled (c). `kind` = "y"
## (y -> c y) or "time" (t -> c t; `time_col` is the time column).
## `tol`: list(lim, df) overriding the closed-form limit (1e-8) and df (1e-6)
## tolerances for one fixture (the retry fixture; see sc_fixtures()).
expect_equivariant <- function(x1, xc, c, info, kind = c("y", "time"), time_col = NULL,
                               tol = list(lim = 1e-8, df = 1e-6)) {
  kind <- match.arg(kind)
  gpq <- if ("method" %in% names(x1)) x1$method == "gpq" else rep(FALSE, nrow(x1))
  cmp <- function(a, b, where) {
    expect_identical(sort(names(b)), sort(names(a)), info = paste(info, where))
    for (nm in names(a)) {
      va <- a[[nm]]; vb <- b[[nm]]; lab <- paste(info, where, nm)
      ## the advisory near-boundary flag compares s2_1 (y^2 per time^2) with the
      ## variance total; it is documented as time-scale dependent
      ## (?vfun_components) and drives no route or df: not compared under t -> c t
      if (kind == "time" && nm == "boundary") next
      if (!is.numeric(va)) { expect_identical(vb, va, info = lab); next }
      tl <- if (nm %in% c("lower", "upper")) tol$lim else 1e-8
      tol_lim <- if (length(va) == length(gpq)) ifelse(gpq, 1e-10, tl) else
        if (any(gpq)) 1e-10 else tl
      if (!is.null(time_col) && nm == time_col) {
        expect_equal(vb, va * c, tolerance = 1e-14, info = lab)
      } else if (nm %in% sc_exact) {
        expect_identical(vb, va, info = lab)
      } else if (nm %in% c(sc_c1, sc_c2)) {
        k <- if (kind == "time") 0 else if (nm %in% sc_c1) 1 else 2
        expect_true(all(mapply(rel_ok, vb, va * c^k, tol_lim)), info = lab)
      } else {
        ## dfs and other derived numbers (phi1, frac_boundary, c_resid, ...)
        expect_true(rel_ok(vb, va, tol$df), info = lab)
      }
    }
  }
  cmp(as.list(as.data.frame(x1)), as.list(as.data.frame(xc)), "[data]")
  d1 <- attr(x1, "diagnostics"); dc <- attr(xc, "diagnostics")
  expect_identical(is.null(dc), is.null(d1), info = info)
  if (!is.null(d1)) cmp(d1, dc, "[diagnostics]")
  for (a in c("df_fallback_reason", "pi_df_fallback", "singular")) {
    expect_identical(attr(xc, a), attr(x1, a), info = paste(info, a))
  }
}

## call a verb, keeping its value and every warning class (each warning
## muffled); an error is kept as its class
run_verb <- function(expr) {
  .reset_once_warnings()      # once-per-session warnings fire on every call here
  w <- character(0)
  v <- tryCatch(withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, class(cnd)[1]); invokeRestart("muffleWarning")
  }), error = function(e) structure(list(), class = c("verb_error", class(e)[1])))
  list(value = v, warnings = sort(unique(w)))
}

## ---- fixtures (data at c = 1; each fit is refitted on c * y) ----------------
sc_fixtures <- function() {
  ss <- lme4::sleepstudy
  set.seed(11)
  g <- factor(rep(1:9, c(2, 4, 3, 5, 2, 6, 3, 4, 3)))   # unbalanced one-way
  ow <- data.frame(g = g, y = 0.12 + rnorm(9, 0, 0.003)[g] + rnorm(length(g), 0, 0.002))
  set.seed(1); A <- 5; B <- 3; n <- 2
  a <- factor(rep(1:A, each = B * n)); b <- factor(rep(1:(A * B), each = n))
  ne <- data.frame(a = a, b = b, y = 10 + rnorm(A, 0, 2)[a] + rnorm(A * B, 0, 1)[b] +
                     rnorm(A * B * n))
  set.seed(1); A <- 3; B <- 4; n <- 2
  cr <- expand.grid(r = 1:n, a = factor(1:A), b = factor(1:B))
  cr$ab <- interaction(cr$a, cr$b)
  cr$y <- 10 + rnorm(A, 0, 2)[cr$a] + rnorm(B, 0, 2)[cr$b] + rnorm(A * B, 0, 1)[cr$ab] +
    rnorm(nrow(cr))
  set.seed(5); K <- 6; tt <- c(0, 3, 6, 9, 12)
  fs <- expand.grid(time = tt, lot = factor(1:K))
  fs$y <- 100 - 0.4 * fs$time + rnorm(K, 0, 1.5)[fs$lot] + rnorm(nrow(fs), 0, 0.8)
  ## the K = 12, seed-2, rho ~ 0.95 Hessian-retry fixture of test-dispatch.R
  set.seed(2); K <- 12; tt <- c(0, 3, 6, 9, 12, 18, 24)
  rt <- expand.grid(time = tt, lot = factor(1:K))
  z1 <- rnorm(K); z2 <- 0.95 * z1 + sqrt(1 - 0.95^2) * rnorm(K)
  rt$y <- 100 + (2 * z1)[rt$lot] + (-0.3 + (0.15 * z2)[rt$lot]) * rt$time + rnorm(nrow(rt), 0, 1)
  list(
    rs_un  = list(f = Reaction ~ Days + (Days | Subject), d = ss, y = "Reaction",
                  nd = data.frame(Days = c(0, 4.5, 9)), slope = TRUE),
    rs_ind = list(f = Reaction ~ Days + (Days || Subject), d = ss, y = "Reaction",
                  nd = data.frame(Days = c(0, 4.5, 9)), slope = TRUE),
    oneway = list(f = y ~ 1 + (1 | g), d = ow, y = "y", nd = NULL, slope = FALSE),
    nested = list(f = y ~ 1 + (1 | a/b), d = ne, y = "y", nd = NULL, slope = FALSE),
    crossed = list(f = y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), d = cr, y = "y",
                   nd = NULL, slope = FALSE),
    fixed_slope = list(f = y ~ time + (1 | lot), d = fs, y = "y",
                       nd = data.frame(time = c(0, 12, 24)), slope = FALSE),
    ## Tolerance: its Cov(phi) comes from the 0.1%-step retry Hessian on a
    ## rho = 0.95 fit, whose second differences amplify the floating-point
    ## noise of reml_nll (its magnitude shifts with the units by
    ## (n - p) log c^2): measured over y x 1e-4 ... 1e4 and t x 1e-2 ... 1e2,
    ## dfs agree to 4.9e-5 and limits to 1.1e-6 relative (estimates, SEs and
    ## variances, which do not use the Hessian, to 1e-8; GPQ limits to
    ## 1.6e-15). Limits 1e-5, dfs 1e-4.
    retry = list(f = y ~ time + (time | lot), d = rt, y = "y",
                 nd = data.frame(time = c(12, 24)), slope = TRUE,
                 tol = list(lim = 1e-5, df = 1e-4))
  )
}

## y -> c y: the relative theta is unchanged; the scaled fit is built at the
## c = 1 theta with no optimization (see the header)
sc_fit <- function(fx, c) {
  f1 <- lmer_q(fx$f, fx$d)
  if (c == 1) return(f1)
  d <- fx$d; d[[fx$y]] <- d[[fx$y]] * c
  lmer_q(fx$f, d, start = list(theta = lme4::getME(f1, "theta")),
         control = lme4::lmerControl(optimizer = NULL))
}

## every verb / arm / side, as named thunks over a fit
sc_calls <- function(fm, fx) {
  nd <- fx$nd
  ctl <- function(...) lmm_interval_control(seed = 1, M = 2000, ...)
  arms <- if (fx$slope) c("reml", "broad-inference", "gpq", "fv-anova") else
    c("package", "broad-inference", "hts", "classical")
  ngm <- lapply(arms, function(a) {
    force(a)
    function() new_group_mean_lmm(fm, nd, control = if (fx$slope)
      ctl(slope_ngm_method = a) else ctl(ngm_method = a))
  })
  names(ngm) <- paste0("ngm:", arms)
  c(list(
    ci = function() ci_lmm(fm, nd),
    pi = function() pi_lmm(fm, nd),
    ngm_default = function() new_group_mean_lmm(fm, nd),
    ti_two = function() ti_lmm(fm, nd, P = 0.99, control = ctl()),
    ti_lower = function() ti_lmm(fm, nd, P = 0.99, sides = "lower", control = ctl()),
    ti_upper = function() ti_lmm(fm, nd, P = 0.99, sides = "upper", control = ctl()),
    pred_conf = function() lmm_predict(fm, "confidence", nd),
    pred_pred = function() lmm_predict(fm, "prediction", nd),
    pred_tol = function() lmm_predict(fm, "tolerance", nd, P = 0.99, control = ctl()),
    pred_ngm = function() lmm_predict(fm, "new_group_mean", nd)
  ), ngm)
}

## run every call on the c = 1 fit and on each scaled fit, and compare
check_equivariance <- function(fx, nm, scales, make_fit, kind, time_col = NULL,
                               scale_nd = function(nd, c) nd) {
  f1 <- make_fit(fx, 1)
  fx1 <- fx
  base <- lapply(sc_calls(f1, fx1), function(th) run_verb(th()))
  if (nm == "retry") {
    p <- .vfun_parse(f1); rho <- p$phi[["s01"]] / sqrt(p$phi[["s2_0"]] * p$phi[["s2_1"]])
    skip_if(isTRUE(lme4::isSingular(f1)) || abs(rho) < 0.90,
            "the retry fixture's optimum differs on this platform")
    expect_identical(attr(base$ngm_default$value, "diagnostics")$vc_cov_route,
                     c("retry-d1e-3", "retry-d1e-3"))
  }
  for (c in setdiff(scales, 1)) {
    fc <- make_fit(fx, c)
    fxc <- fx; fxc$nd <- scale_nd(fx$nd, c)
    sc <- lapply(sc_calls(fc, fxc), function(th) run_verb(th()))
    for (v in names(base)) {
      info <- sprintf("%s / %s / %s x %g", nm, v, kind, c)
      b <- base[[v]]; s <- sc[[v]]
      expect_identical(s$warnings, b$warnings, info = info)
      expect_identical(class(s$value), class(b$value), info = info)
      if (!inherits(b$value, "verb_error"))
        expect_equivariant(b$value, s$value, c, info, kind = kind, time_col = time_col,
                           tol = if (is.null(fx$tol)) list(lim = 1e-8, df = 1e-6) else fx$tol)
    }
  }
}

test_that("every interval and df is equivariant under y -> c y (c = 1e-4 ... 1e4)", {
  skip_if_not_installed("lme4")
  skip_on_cran()   # 35 fits and ~600 interval calls
  fxs <- sc_fixtures()
  for (nm in names(fxs)) check_equivariance(fxs[[nm]], nm, y_scales, sc_fit, "y")
})

## t -> c t at the transformed theta (no optimization; see the header)
sc_fit_time <- function(fx, c) {
  f1 <- sc_fit(fx, 1)
  if (c == 1) return(f1)
  d <- fx$d; d[[fx$time]] <- d[[fx$time]] * c
  lmer_q(fx$f, d, start = list(theta = lme4::getME(f1, "theta") * fx$theta_scale(c)),
         control = lme4::lmerControl(optimizer = NULL))
}

test_that("every interval and df is invariant under t -> c t (c = 1e-2, 1/12, 1, 12, 1e2)", {
  skip_if_not_installed("lme4")
  skip_on_cran()
  fxs <- sc_fixtures()
  ## relative-Cholesky theta under t -> c t: the slope row is divided by c
  fxs$rs_un$theta_scale  <- function(c) c(1, 1 / c, 1 / c)   # L11, L21, L22
  fxs$rs_ind$theta_scale <- function(c) c(1, 1 / c)          # intercept, slope
  fxs$retry$theta_scale  <- function(c) c(1, 1 / c, 1 / c)
  fxs$fixed_slope$theta_scale <- function(c) 1
  for (nm in c("rs_un", "rs_ind", "retry", "fixed_slope")) {
    fx <- fxs[[nm]]; fx$time <- names(fx$nd)[1]
    check_equivariance(fx, nm, t_scales, sc_fit_time, "time", time_col = fx$time,
                       scale_nd = function(nd, c) { nd[[1]] <- nd[[1]] * c; nd })
  }
})

test_that("a one-way fit in small units keeps its Satterthwaite PI df (no containment fallback)", {
  skip_if_not_installed("lme4")
  fx <- sc_fixtures()$oneway                         # genuine lme4 refits
  p1 <- run_verb(pi_lmm(lmer_q(fx$f, fx$d)))
  pk <- run_verb(pi_lmm(lmer_q(fx$f, transform(fx$d, y = y * 0.001))))
  d1 <- attr(p1$value, "diagnostics"); dk <- attr(pk$value, "diagnostics")
  ## the variance components at c = 1 are ~1e-5 (impurity in %, SD ~ 0.003) and
  ## ~1e-11 at c = 1e-3: below numDeriv's zero.tol at both
  vc <- as.data.frame(lme4::VarCorr(lmer_q(fx$f, fx$d)))$vcov
  expect_true(all(vc < 1.78e-5))
  expect_false(attr(p1$value, "diagnostics")$pi_df_fallback)
  expect_false(attr(pk$value, "diagnostics")$pi_df_fallback)
  expect_equal(dk$df, d1$df, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(d1$df, rep(.ddf_containment(lmer_q(fx$f, fx$d)), length(d1$df)))))
  expect_identical(pk$warnings, p1$warnings)
  expect_false("reintervals_warning_df_fallback" %in% pk$warnings)
})

test_that("sleepstudy x 1e-4: new_group_mean_lmm() runs and keeps the y x 1 df rule", {
  skip_if_not_installed("lme4")
  ss <- lme4::sleepstudy; nd <- data.frame(Days = c(0, 5, 9))
  f1 <- lmer_q(Reaction ~ Days + (Days | Subject), ss)
  ss$Reaction <- ss$Reaction * 1e-4
  fk <- lmer_q(Reaction ~ Days + (Days | Subject), ss)
  r1 <- new_group_mean_lmm(f1, nd)
  rk <- expect_no_error(new_group_mean_lmm(fk, nd))
  expect_identical(attr(rk, "diagnostics")$df_rule, attr(r1, "diagnostics")$df_rule)
  expect_equal(attr(rk, "diagnostics")$df, attr(r1, "diagnostics")$df, tolerance = 1e-6)
})

test_that("a failed numerical gradient of V_F gives the containment df, reason recorded, one warning", {
  skip_if_not_installed("lme4")
  fm <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  real <- .nd_deriv
  local_mocked_bindings(.nd_deriv = function(type, f, x, s, ...) {
    if (identical(type, "grad")) stop("leading minor not positive (mocked)")
    real(type, f, x, s, ...)
  })
  cp <- vfun_extract(fm, t0 = 5)
  expect_true(attr(cp, "df_fallback")[["F"]])
  expect_identical(attr(cp, "df_fallback_why")[["F"]], "numerical derivative failed")
  expect_false(attr(cp, "df_fallback")[["G"]])          # V_G's gradient is analytic
  ws <- list()
  x <- withCallingHandlers(ci_lmm(fm, data.frame(Days = 5)),
         warning = function(w) { ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning") })
  fb <- Filter(function(w) inherits(w, "reintervals_warning_df_fallback"), ws)
  expect_length(fb, 1L)
  expect_match(conditionMessage(fb[[1]]), "numerical derivative failed", fixed = TRUE)
  expect_identical(attr(x, "diagnostics")$df_rule, "containment")
  ## broad-inference: the df cannot be formed -> NA limits (its no-fallback rule)
  b <- suppressWarnings(new_group_mean_lmm(fm, data.frame(Days = 5),
         control = lmm_interval_control(slope_ngm_method = "broad-inference")))
  expect_true(attr(b, "diagnostics")$ngm_df_failed)
  expect_true(is.na(b$lower) && is.na(b$upper))
})

test_that(".nd_deriv back-transforms exactly: gradient / s, Hessian / s^2", {
  ## homogeneous cubic on parameters ~1e-6 (below numDeriv's default zero.tol)
  f <- function(p) 1e12 * (p[[1]]^2 * p[[2]] + p[[2]]^3)
  x <- c(a = 2e-6, b = 3e-6); s <- 3e-6
  g <- .nd_deriv("grad", f, x, s)
  H <- .nd_deriv("hessian", f, x, s)
  expect_equal(g, 1e12 * c(2 * x[[1]] * x[[2]], x[[1]]^2 + 3 * x[[2]]^2), tolerance = 1e-8)
  expect_equal(H, 1e12 * matrix(c(2 * x[[2]], 2 * x[[1]], 2 * x[[1]], 6 * x[[2]]), 2),
               tolerance = 1e-8)
  ## the same derivative at any scale s (scale-free), and relative steps
  expect_equal(.nd_deriv("hessian", f, x, 1), H, tolerance = 1e-8)
  expect_error(.nd_deriv("grad", f, x, 0))
  expect_error(.nd_deriv("grad", f, x, NA_real_))
  expect_identical(.nd_scale(c(s2_0 = 4, s01 = 1, s2e = 2)), 2)
  expect_identical(.nd_scale(c(run = 4, Residual = 3)), 3)
})

test_that("the half-step stability guard: a disagreeing half-step Hessian makes Cov unavailable, with its reason", {
  skip_if_not_installed("lme4")
  ## random slope: force the guard on (as if some |psi_j| < 1.78e-5) and the
  ## two Hessians' dfs to disagree
  fm <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  local_mocked_bindings(.nd_unstable_psi = function(x, s) TRUE)
  ok <- quiet_df_fallback(new_group_mean_lmm(fm, nd))        # guard passes: stable
  expect_identical(attr(ok, "diagnostics")$vc_cov_route, c("default", "default"))
  local_mocked_bindings(.df_stable = function(a, b, tol = 1e-3) FALSE)
  cp <- vfun_extract(fm, t0 = 0, guard_t0 = nd$Days)
  expect_null(cp$cov)
  expect_identical(attr(cp, "vc_cov_route"), "unavailable")
  expect_identical(attr(cp, "vc_cov_reason"), "unstable numerical Hessian")
  expect_warning(x <- new_group_mean_lmm(fm, nd), "unstable numerical Hessian",
                 class = "reintervals_warning_df_fallback")
  expect_identical(attr(x, "diagnostics")$df_rule, c("containment", "containment"))
  expect_identical(attr(x, "diagnostics")$vc_cov_route, c("unavailable", "unavailable"))
  ## random intercept: the PI df falls back with the same reason
  fr <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  th <- .theta_cov(fr, L = .fixed(fr, nd)$L)
  expect_null(th$cov_theta)
  expect_true(th$pi_df_fallback)
  expect_identical(th$fallback_reason, "unstable numerical Hessian")
  expect_warning(p <- pi_lmm(fr, nd), "unstable numerical Hessian",
                 class = "reintervals_warning_df_fallback")
  expect_true(isTRUE(attr(p, "diagnostics")$pi_df_fallback))
})

test_that(".nd_unstable_psi and .df_stable", {
  expect_false(.nd_unstable_psi(c(s2_0 = 1, s2_1 = 1e-3, s2e = 1)))
  expect_true(.nd_unstable_psi(c(s2_0 = 1, s2_1 = 1e-6, s2e = 1)))
  expect_true(.nd_unstable_psi(c(s2_0 = 1e-2, s2_1 = 1e-8, s2e = 1e-2)))  # ratio 1e-6
  expect_false(.nd_unstable_psi(c(s2_0 = 1e-8, s2_1 = 1e-9, s2e = 1e-8)))  # ratio 0.1
  expect_true(.nd_unstable_psi(c(s2_0 = 0, s2e = 1)))                      # singular
  expect_true(.df_stable(c(10, 20), c(10.005, 20)))
  expect_false(.df_stable(c(10, 20), c(10.02, 20)))
  expect_true(.df_stable(c(NA, 20), c(NA, 20)))
  expect_false(.df_stable(c(NA, 20), c(5, 20)))
  expect_true(.df_stable(numeric(0), numeric(0)))
})

## slope_pi_method = "fv-anova" (and "hybrid"): limits, estimates and SEs scale
## by c and W by c^2 under y -> c y; every number is unchanged under t -> c t;
## the df (Satterthwaite on S2_L and S2_E) and the effective df of S2_L
## (fv_df_lots) are invariant. rs_unbal (three lots on 0..24, three on 0..12)
## has unequal leverages, so fv_df_lots < K - 1 there.
test_that("slope_pi_method = 'fv-anova' is equivariant under y -> c y and invariant under t -> c t", {
  skip_if_not_installed("lme4")
  skip_on_cran()
  fxs <- sc_fixtures()
  fxs$rs_un$theta_scale  <- function(c) c(1, 1 / c, 1 / c)
  fxs$rs_ind$theta_scale <- function(c) c(1, 1 / c)
  fxs$retry$theta_scale  <- function(c) c(1, 1 / c, 1 / c)
  set.seed(3); full <- c(0, 3, 6, 9, 12, 18, 24); short <- c(0, 3, 6, 9, 12)
  ub <- do.call(rbind, lapply(1:6, function(k) {
    tt <- if (k <= 3) full else short
    data.frame(lot = factor(k), time = tt,
               y = 100 + rnorm(1) + (-0.3 + rnorm(1, 0, 0.05)) * tt + rnorm(length(tt)))
  }))
  fxs$rs_unbal <- list(f = y ~ time + (time | lot), d = ub, y = "y",
                       nd = data.frame(time = c(6, 24)), slope = TRUE,
                       theta_scale = function(c) c(1, 1 / c, 1 / c))
  pi_calls <- function(fm, nd) lapply(c(fv = "fv-anova", hybrid = "hybrid"), function(m)
    lapply(c("two.sided", "lower", "upper"), function(sd)
      run_verb(pi_lmm(fm, nd, sides = sd, control = lmm_interval_control(slope_pi_method = m)))))
  for (nm in c("rs_un", "rs_ind", "retry", "rs_unbal")) {
    fx <- fxs[[nm]]; fx$time <- names(fx$nd)[1]
    b <- pi_calls(sc_fit(fx, 1), fx$nd)
    expect_identical(b$fv[[1]]$value$method, rep("fv-anova", nrow(fx$nd)))
    dfl <- attr(b$fv[[1]]$value, "diagnostics")$fv_df_lots
    if (nm == "rs_unbal") expect_true(all(dfl < 6 - 1)) else
      expect_identical(dfl, rep(nlevels(factor(fx$d[[if (nm == "retry") "lot" else "Subject"]])) - 1,
                                nrow(fx$nd)))
    for (kind in c("y", "time")) {
      scales <- if (kind == "y") y_scales else t_scales
      for (c in setdiff(scales, 1)) {
        if (kind == "y") {
          fc <- sc_fit(fx, c); nd <- fx$nd
        } else {
          fc <- sc_fit_time(fx, c); nd <- fx$nd; nd[[1]] <- nd[[1]] * c
        }
        s <- pi_calls(fc, nd)
        for (m in names(b)) for (j in 1:3) {
          info <- sprintf("%s / pi:%s side %d / %s x %g", nm, m, j, kind, c)
          expect_identical(s[[m]][[j]]$warnings, b[[m]][[j]]$warnings, info = info)
          if (m == "fv")
            expect_true(rel_ok(attr(s[[m]][[j]]$value, "diagnostics")$fv_df_lots,
                               attr(b[[m]][[j]]$value, "diagnostics")$fv_df_lots, 1e-8),
                        info = paste(info, "fv_df_lots"))
          expect_equivariant(b[[m]][[j]]$value, s[[m]][[j]]$value, c, info, kind = kind,
                             time_col = if (kind == "time") fx$time,
                             tol = if (is.null(fx$tol)) list(lim = 1e-8, df = 1e-6) else fx$tol)
        }
      }
    }
  }
})

test_that("ti_lmm(method = 'gpq') on a random-intercept fit is equivariant under y -> c y and invariant under t -> c t (same seed)", {
  skip_if_not_installed("lme4")
  set.seed(3)
  nk <- c(5, 1, 4, 7, 1, 3, 6, 2); tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- data.frame(lot = factor(rep(seq_along(nk), nk)))
  d$time <- unlist(lapply(nk, function(k) tt[seq_len(k)]))
  d$y <- 100 - 0.3 * d$time + rnorm(length(nk), 0, 1.2)[d$lot] + rnorm(nrow(d), 0, 0.8)
  nd <- data.frame(time = c(0, 12, 30))
  ctl <- lmm_interval_control(seed = 7L, M = 2000L)
  ## the GPQ uses only (y, X, lot), never the REML fit, so genuine refits are used
  run <- function(dd, ndd, ov) run_verb(ti_lmm(lmer_q(y ~ time + (1 | lot), dd), ndd,
                                               P = c(0.9, 0.99), sides = c("two.sided", "lower"),
                                               over = ov, method = "gpq", control = ctl))
  for (ov in c("observation", "group_mean")) {
    b <- run(d, nd, ov)
    expect_identical(b$value$method, rep("gpq", nrow(b$value)))
    for (c in y_scales) {
      s <- run(transform(d, y = y * c), nd, ov)
      info <- sprintf("gpq-ri %s y x %g", ov, c)
      expect_identical(s$warnings, b$warnings, info = info)
      expect_equivariant(b$value, s$value, c, info, kind = "y")
    }
    for (c in t_scales) {
      s <- run(transform(d, time = time * c), transform(nd, time = time * c), ov)
      info <- sprintf("gpq-ri %s t x %g", ov, c)
      expect_identical(s$warnings, b$warnings, info = info)
      expect_equivariant(b$value, s$value, c, info, kind = "time", time_col = "time")
    }
  }
})
