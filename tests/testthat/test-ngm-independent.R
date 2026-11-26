# New-group-mean intervals against an INDEPENDENT IMPLEMENTATION written in this
# file from the published formulas. It is not a package call: it builds the
# marginal covariance V, the REML log-likelihood, C11 = (X' V^-1 X)^-1 and its
# gradient itself, and uses no reintervals internal. Oracles
# indep-ngm-ri-broad and indep-ngm-rs-reml-broad.
#
# Parameterisation (the one the package documents): the variance components on
# the variance scale, residual last -- (sigma2_1, ..., sigma2_e) for the
# random-intercept models; phi = (s2_0, s01, s2_1, s2e) (correlated) or
# (s2_0, s2_1, s2e) (independent) for the random-slope models. A = the inverse
# of minus the numerical Hessian (numDeriv, default Richardson settings) of the
# REML log-likelihood in these parameters at the REML estimates; the gradient
# of l C11 l' is analytic, d C11 / d p_j = C11 X' V^-1 (dV / d p_j) V^-1 X C11.
#
# (a) Random-intercept broad-inference prediction of a new group mean (McLean,
#     Sanders & Stroup 1991): W = l C11 l' + sum of the between-group
#     components, df = 2 W^2 / (g' A g), g the gradient of W
#     (ngm_method = "broad-inference").
# (b) Random-slope mean mu(t0) = x0' beta, W = x0' C11 x0 + V_G(t0),
#     V_G(t0) = z0' Sigma z0, z0 = (1, t0): "reml" (the default) takes the
#     Satterthwaite df of V_G(t0) alone, 2 V_G^2 / (gG' A gG); "broad-inference"
#     the Satterthwaite df of the whole W.
# Limits: estimate +/- t_{0.975, df} sqrt(W).
#
# Fits (all non-singular): (a) the one-way unbalanced, two-fold nested and
# staggered fixed-slope fits of fixtures/sim_designs/; (b) lme4::sleepstudy with
# a correlated random slope, the balanced stability design with independent
# random intercept and slope (stab_bal.csv), and the staggered stability design
# with a correlated random slope (stab_unbal_v2_a.csv).
#
# Measured on 2026-11-26: (a) max relative difference df 1.0e-10, SE 2e-16,
# limits 6e-13; (b) df 7.7e-9, SE 3e-15, limits 2.4e-10. Both sides take a
# numerical REML Hessian with numDeriv's Richardson scheme (in different code
# and, for the package, in standardised parameters), and the package's
# gradient of V_F is numerical where this one is analytic. Tolerance: df 1e-6
# relative, limits and SE 1e-8 relative (two orders of magnitude above the
# measured agreement, for the platform dependence of the numerical derivatives).

## REML log-likelihood (up to a constant) at variance-scale parameters `par`,
## Vfun(par) the marginal covariance of y; beta at its GLS value
reml_ll <- function(Vfun, X, y) function(par) {
  V <- Vfun(par)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(NA_real_)
  Vi <- chol2inv(R)
  XtViX <- crossprod(X, Vi %*% X)
  b <- solve(XtViX, crossprod(X, Vi %*% y))
  r <- y - X %*% b
  -0.5 * (2 * sum(log(diag(R))) + as.numeric(determinant(XtViX)$modulus) +
            as.numeric(crossprod(r, Vi %*% r)))
}
## C11 = (X' V^-1 X)^-1 and its analytic derivative in parameter j:
## dC11/dp_j = C11 X' V^-1 (dV/dp_j) V^-1 X C11
c11_and_grad <- function(Vfun, dV, X, par, l) {
  Vi <- solve(Vfun(par)); C <- solve(crossprod(X, Vi %*% X))
  P <- C %*% crossprod(X, Vi)
  g <- vapply(dV, function(D) as.numeric(t(l) %*% P %*% D %*% t(P) %*% l), numeric(1))
  list(value = as.numeric(t(l) %*% C %*% l), grad = g)
}

## (a) random-intercept broad inference:
## W = l C11 l' + sum_j sigma_j^2 (between-group components), df = 2 W^2 / g'Ag
ri_broad_ref <- function(fm, L, level = 0.95) {
  y <- lme4::getME(fm, "y"); X <- as.matrix(lme4::getME(fm, "X")); N <- length(y)
  Zl <- lapply(lme4::getME(fm, "Ztlist"), function(z) as.matrix(Matrix::crossprod(z)))
  vc <- as.data.frame(lme4::VarCorr(fm))
  gn <- sub("\\.\\(Intercept\\)$", "", names(Zl))
  par <- c(vc$vcov[match(gn, vc$grp)], vc$vcov[vc$grp == "Residual"])
  k <- length(par)
  Vfun <- function(p) Reduce(`+`, Map(`*`, Zl, p[-k])) + diag(p[k], N)
  dV <- c(Zl, list(diag(N)))
  A <- solve(-numDeriv::hessian(reml_ll(Vfun, X, y), par))
  do.call(rbind, lapply(seq_len(nrow(L)), function(i) {
    l <- L[i, ]; cg <- c11_and_grad(Vfun, dV, X, par, l)
    W <- cg$value + sum(par[-k]); g <- cg$grad + c(rep(1, k - 1), 0)
    df <- 2 * W^2 / as.numeric(t(g) %*% A %*% g)
    est <- sum(l * lme4::fixef(fm)); hw <- qt(1 - (1 - level) / 2, df) * sqrt(W)
    data.frame(estimate = est, se = sqrt(W), df = df, lower = est - hw, upper = est + hw)
  }))
}

## (b) random slope y ~ t + (1 + t | g) [cor] or (1 | g) + (0 + t | g) [ind]:
## phi = (s2_0, s01, s2_1, s2e) or (s2_0, s2_1, s2e); V_G(t0) = z0' Sigma z0,
## W = x0' C11 x0 + V_G(t0); "reml": df = 2 V_G^2 / gG'A gG; "broad": df on W.
rs_ref <- function(fm, time, group, t0s, correlated, level = 0.95) {
  mf <- stats::model.frame(fm)
  y <- lme4::getME(fm, "y"); X <- cbind(1, mf[[time]]); N <- length(y)
  tt <- mf[[time]]; G <- outer(mf[[group]], mf[[group]], "==") * 1
  D00 <- G; D01 <- G * (outer(tt, rep(1, N)) + outer(rep(1, N), tt)); D11 <- G * outer(tt, tt)
  vc <- as.data.frame(lme4::VarCorr(fm))
  s2e <- vc$vcov[vc$grp == "Residual"]
  if (correlated) {
    s2_0 <- vc$vcov[which(vc$var1 == "(Intercept)" & is.na(vc$var2))]
    s2_1 <- vc$vcov[which(vc$var1 == time & is.na(vc$var2))]
    s01 <- vc$vcov[which(!is.na(vc$var2))]
    par <- c(s2_0 = s2_0, s01 = s01, s2_1 = s2_1, s2e = s2e)
    Vfun <- function(p) p[1] * D00 + p[2] * D01 + p[3] * D11 + diag(p[4], N)
    dV <- list(D00, D01, D11, diag(N))
    gG <- function(t0) c(1, 2 * t0, t0^2, 0)
    VG <- function(t0) par[1] + 2 * t0 * par[2] + t0^2 * par[3]
  } else {
    s2_0 <- vc$vcov[which(vc$var1 == "(Intercept)")]; s2_1 <- vc$vcov[which(vc$var1 == time)]
    par <- c(s2_0 = s2_0, s2_1 = s2_1, s2e = s2e)
    Vfun <- function(p) p[1] * D00 + p[2] * D11 + diag(p[3], N)
    dV <- list(D00, D11, diag(N))
    gG <- function(t0) c(1, t0^2, 0)
    VG <- function(t0) par[1] + t0^2 * par[2]
  }
  A <- solve(-numDeriv::hessian(reml_ll(Vfun, X, y), par))
  do.call(rbind, lapply(t0s, function(t0) {
    x0 <- c(1, t0); cg <- c11_and_grad(Vfun, dV, X, par, x0)
    vg <- unname(VG(t0)); W <- cg$value + vg
    df_reml <- 2 * vg^2 / as.numeric(t(gG(t0)) %*% A %*% gG(t0))
    g <- cg$grad + gG(t0)
    df_broad <- 2 * W^2 / as.numeric(t(g) %*% A %*% g)
    est <- sum(x0 * lme4::fixef(fm))
    q <- function(df) qt(1 - (1 - level) / 2, df) * sqrt(W)
    data.frame(time = t0, estimate = est, se = sqrt(W), df_reml = df_reml, df_broad = df_broad,
               lower_reml = est - q(df_reml), upper_reml = est + q(df_reml),
               lower_broad = est - q(df_broad), upper_broad = est + q(df_broad))
  }))
}

sd_read <- function(f) utils::read.csv(test_path("fixtures", "sim_designs", f),
                                       stringsAsFactors = TRUE)
## package values in a common form: estimate, SE (from the two-sided width), df, limits
pkg_ngm <- function(fm, nd, control) {
  p <- quiet_extrapolation(new_group_mean_lmm(fm, newdata = nd, control = control))
  df <- attr(p, "diagnostics")$df
  list(estimate = p$estimate, se = (p$upper - p$lower) / 2 / stats::qt(0.975, df),
       df = df, lower = p$lower, upper = p$upper,
       df_rule = attr(p, "diagnostics")$df_rule)
}

# oracle: indep-ngm-ri-broad
test_that("random-intercept broad-inference NGM matches the independent implementation", {
  skip_if_not_installed("lme4")
  mo <- c(0, 12, 24, 36, 48)
  fits <- list(
    OW_UNBAL = list(y ~ 1 + (1 | grp), sd_read("oneway_unbal.csv"), NULL),
    NEST_BAL = list(y ~ 1 + (1 | a) + (1 | a:b), sd_read("nested_bal.csv"), NULL),
    FS_UNBAL = list(y ~ month + (1 | lot), sd_read("stab_unbal.csv"), data.frame(month = mo)))
  for (tag in names(fits)) {
    f <- fits[[tag]]
    fm <- lme4::lmer(f[[1]], data = f[[2]], REML = TRUE)
    expect_false(lme4::isSingular(fm), info = tag)
    L <- if (is.null(f[[3]])) matrix(1, 1, 1) else cbind(1, f[[3]][[1]])
    p <- pkg_ngm(fm, f[[3]], lmm_interval_control(ngm_method = "broad-inference"))
    r <- ri_broad_ref(fm, L)
    expect_identical(p$df_rule, rep("satterthwaite-broad", nrow(L)), info = tag)
    expect_equal(p$estimate, r$estimate, tolerance = 1e-8, info = tag)
    expect_equal(p$se, r$se, tolerance = 1e-8, info = tag)
    expect_equal(p$df, r$df, tolerance = 1e-6, info = tag)
    expect_equal(p$lower, r$lower, tolerance = 1e-8, info = tag)
    expect_equal(p$upper, r$upper, tolerance = 1e-8, info = tag)
  }
})

# oracle: indep-ngm-rs-reml-broad
test_that("random-slope reml (default) and broad-inference NGM match the independent implementation", {
  skip_if_not_installed("lme4")
  mo <- c(0, 12, 24, 36)
  fits <- list(
    SLEEP_UN  = list(Reaction ~ Days + (Days | Subject), lme4::sleepstudy,
                     "Days", "Subject", TRUE, c(0, 4.5, 9)),
    RS_VC_BAL = list(y ~ month + (1 | lot) + (0 + month | lot), sd_read("stab_bal.csv"),
                     "month", "lot", FALSE, mo),
    RS_UN_V2  = list(y ~ month + (1 + month | lot), sd_read("stab_unbal_v2_a.csv"),
                     "month", "lot", TRUE, mo))
  for (tag in names(fits)) {
    f <- fits[[tag]]
    fm <- lme4::lmer(f[[1]], data = f[[2]], REML = TRUE)
    expect_false(lme4::isSingular(fm), info = tag)
    nd <- stats::setNames(data.frame(f[[6]]), f[[3]])
    r <- rs_ref(fm, f[[3]], f[[4]], f[[6]], correlated = f[[5]])
    for (arm in c("reml", "broad-inference")) {
      p <- pkg_ngm(fm, nd, lmm_interval_control(slope_ngm_method = arm))
      sfx <- if (arm == "reml") "reml" else "broad"
      info <- paste(tag, arm)
      expect_identical(p$df_rule, rep(if (arm == "reml") "satterthwaite" else "satterthwaite-broad",
                                      length(f[[6]])), info = info)
      expect_equal(p$estimate, r$estimate, tolerance = 1e-8, info = info)
      expect_equal(p$se, r$se, tolerance = 1e-8, info = info)
      expect_equal(p$df, r[[paste0("df_", sfx)]], tolerance = 1e-6, info = info)
      expect_equal(p$lower, r[[paste0("lower_", sfx)]], tolerance = 1e-8, info = info)
      expect_equal(p$upper, r[[paste0("upper_", sfx)]], tolerance = 1e-8, info = info)
    }
  }
})
