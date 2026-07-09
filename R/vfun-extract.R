# =============================================================================
# Slope-engine extraction layer: the path that builds a populated
# vfun_components from a fitted random-slope (or fixed-slope) lmer model.
#
# This is the SLOPE-side counterpart to vc-extract.R: an extraction layer, so it
# may touch lme4 / Matrix / numDeriv (not a method-pure constructor). It is a
# SEPARATE path -- it does not call, modify, or route through vc_extract(), and
# the EMS random-slope stop() there is left untouched.
#
# Ported from the Montes/GPQ source (parse / vccov / df / intervals),
# adapted to populate vfun_components. The slope engine keeps its OWN
# singular/boundary handling and its OWN Satterthwaite->containment df fallback;
# these are NOT unified with the EMS-side guards in vc-extract.R.
#
# Safeguards preserved verbatim in behavior:
#   (1) the reml_nll non-PD-Sigma / failed-Cholesky sentinel (returns 1e10 so
#       the numerical Hessian stays well-defined);
#   (2) isSingular + near-boundary detection;
#   (3) the Satterthwaite -> containment df fallback on degenerate df.
# =============================================================================

## ---- V_G(t0) value and gradient w.r.t. phi (port of .value_VG / .grad_VG) ----
.vfun_value_VG <- function(type, phi, t0) {
  switch(type,
    M1  = phi[["s2_0"]],
    M2c = phi[["s2_0"]] + 2 * t0 * phi[["s01"]] + t0^2 * phi[["s2_1"]],
    M2i = phi[["s2_0"]] + t0^2 * phi[["s2_1"]])
}
.vfun_grad_VG <- function(type, t0) {
  switch(type,
    M1  = c(s2_0 = 1, s2e = 0),
    M2c = c(s2_0 = 1, s01 = 2 * t0, s2_1 = t0^2, s2e = 0),
    M2i = c(s2_0 = 1, s2_1 = t0^2, s2e = 0))
}
.vfun_grad_VT <- function(type, t0) {
  g <- .vfun_grad_VG(type, t0); g[["s2e"]] <- 1; g
}

## ---- df helpers (port of df.R) ----------------------------------------------
## generalized Satterthwaite df for a scalar variance v = f(phi)
.vfun_df_satterthwaite <- function(value, grad_vec, vc_cov) {
  if (is.null(vc_cov)) return(NA_real_)
  varv <- as.numeric(t(grad_vec) %*% vc_cov %*% grad_vec)
  if (!is.finite(varv) || varv <= 0) return(NA_real_)
  nu <- 2 * value^2 / varv
  if (!is.finite(nu) || nu <= 0) return(NA_real_)
  nu
}
.vfun_df_between  <- function(parsed) max(1, parsed$B - parsed$k_re)
.vfun_df_residual <- function(parsed) max(1, parsed$n - parsed$p - parsed$B * parsed$k_re)
.vfun_df_containment <- function(parsed, w_between) {
  if (w_between >= 0.5) .vfun_df_between(parsed) else .vfun_df_residual(parsed)
}

## Near-boundary predicate: a random-effect VARIANCE (s2_0, s2_1) is
## small-but-positive relative to total variance. Pure function of (phi,
## boundary_tol) so it is testable without a fit. 
.vfun_is_boundary <- function(phi, boundary_tol) {
  total <- sum(phi[grep("^s2", names(phi))])
  re_var_names <- setdiff(grep("^s2", names(phi), value = TRUE), "s2e")
  any(phi[re_var_names] < boundary_tol * total)
}

## ---- model parse: classify + Sigma/V/REML-NLL closures  --
.vfun_parse <- function(model) {
  if (!inherits(model, "merMod"))
    stop("vfun_extract(): `model` must be an lme4 merMod fit.", call. = FALSE)

  cnms <- lme4::getME(model, "cnms")
  gnames <- names(cnms)
  if (length(unique(gnames)) != 1L)
    stop("vfun_extract(): exactly one grouping factor is supported (got: ",
         paste(unique(gnames), collapse = ", "), ").", call. = FALSE)
  group_name <- unique(gnames)

  re_terms <- unlist(cnms, use.names = FALSE)
  has_slope <- length(re_terms) >= 2L
  n_terms   <- length(cnms)

  slope_name <- setdiff(re_terms, "(Intercept)")
  if (has_slope && length(slope_name) != 1L)
    stop("vfun_extract(): expected a single slope term; got: ",
         paste(slope_name, collapse = ", "), call. = FALSE)

  if (!has_slope) {
    type <- "M1"; time_name <- NA_character_
  } else if (n_terms == 1L) {
    type <- "M2c"; time_name <- slope_name
  } else {
    type <- "M2i"; time_name <- slope_name
  }

  fe <- lme4::fixef(model)
  fe_names <- names(fe)
  if (!("(Intercept)" %in% fe_names))
    stop("vfun_extract(): fixed part must contain an intercept.", call. = FALSE)
  if (has_slope && !(time_name %in% fe_names))
    stop("vfun_extract(): the random-slope variable '", time_name,
         "' must also be a fixed effect.", call. = FALSE)
  fe_time <- setdiff(fe_names, "(Intercept)")
  if (length(fe_time) != 1L)
    stop("vfun_extract(): fixed part must be exactly ~ 1 + time (found: ",
         paste(fe_names, collapse = ", "), ").", call. = FALSE)
  if (is.na(time_name)) time_name <- fe_time

  mf <- stats::model.frame(model)
  g  <- factor(mf[[group_name]])
  tt <- as.numeric(mf[[time_name]])
  y  <- as.numeric(lme4::getME(model, "y"))
  X  <- lme4::getME(model, "X")
  X  <- X[, c("(Intercept)", fe_time), drop = FALSE]
  n  <- length(y); p <- ncol(X)
  B  <- nlevels(g); lev <- as.integer(g)
  k_re <- if (has_slope) 2L else 1L

  ## rebuild Z with ordering [lvl1:(int,slope), lvl2:(int,slope), ...]
  if (!has_slope) {
    Z <- Matrix::sparseMatrix(i = seq_len(n), j = lev, x = 1, dims = c(n, B))
  } else {
    Zi <- Matrix::sparseMatrix(i = seq_len(n), j = 2L * lev - 1L, x = 1,  dims = c(n, 2L * B))
    Zs <- Matrix::sparseMatrix(i = seq_len(n), j = 2L * lev,      x = tt, dims = c(n, 2L * B))
    Z  <- Zi + Zs
  }

  ## fitted variance components from theta + sigma
  theta <- unname(lme4::getME(model, "theta"))
  s2e   <- unname(stats::sigma(model)^2)
  if (type == "M1") {
    s2_0 <- s2e * theta[1]^2
    phi  <- c(s2_0 = s2_0, s2e = s2e)
  } else if (type == "M2c") {
    th <- theta[1:3]
    L  <- matrix(c(th[1], th[2], 0, th[3]), 2, 2)
    Sig <- s2e * (L %*% t(L))
    phi  <- c(s2_0 = Sig[1, 1], s01 = Sig[1, 2], s2_1 = Sig[2, 2], s2e = s2e)
  } else { # M2i
    th <- theta[1:2]
    phi  <- c(s2_0 = s2e * th[1]^2, s2_1 = s2e * th[2]^2, s2e = s2e)
  }

  make_Sigma <- function(phi) {
    switch(type,
      M1  = matrix(phi[["s2_0"]], 1, 1),
      M2c = matrix(c(phi[["s2_0"]], phi[["s01"]], phi[["s01"]], phi[["s2_1"]]), 2, 2),
      M2i = matrix(c(phi[["s2_0"]], 0, 0, phi[["s2_1"]]), 2, 2))
  }

  In <- diag(n)
  IB <- Matrix::Diagonal(B)
  make_V <- function(phi) {
    Sig <- make_Sigma(phi)
    s2e <- phi[["s2e"]]
    if (type == "M1") {
      Vr <- phi[["s2_0"]] * (Z %*% Matrix::t(Z))
    } else {
      Gfull <- Matrix::kronecker(IB, Matrix::Matrix(Sig))
      Vr <- Z %*% Gfull %*% Matrix::t(Z)
    }
    as.matrix(Vr) + s2e * In
  }

  ## REML negative log-likelihood as a function of phi.
  ## SAFEGUARD (1): returns the 1e10 sentinel on non-positive s2e, non-PD Sigma,
  ## or any failed Cholesky, so numDeriv::hessian() stays well-defined.
  reml_nll <- function(phi) {
    Sig <- make_Sigma(phi)
    if (phi[["s2e"]] <= 0) return(1e10)
    if (length(Sig) > 1 && min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      return(1e10)
    V <- tryCatch(make_V(phi), error = function(e) NULL)
    if (is.null(V)) return(1e10)
    Vc <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(Vc)) return(1e10)
    Vinv     <- chol2inv(Vc)
    XtVinvX  <- t(X) %*% Vinv %*% X
    Xc       <- tryCatch(chol(XtVinvX), error = function(e) NULL)
    if (is.null(Xc)) return(1e10)
    beta <- solve(XtVinvX, t(X) %*% (Vinv %*% y))
    r    <- y - X %*% beta
    quad <- as.numeric(t(r) %*% Vinv %*% r)
    logdetV <- 2 * sum(log(diag(Vc)))
    logdetP <- 2 * sum(log(diag(Xc)))
    0.5 * ((n - p) * log(2 * pi) + logdetV + logdetP + quad)
  }

  ## value-of-V_F as a function of phi (for the mean-line Satterthwaite df)
  vF_fun <- function(phi, x0) {
    V <- make_V(phi); Vc <- chol(V); Vinv <- chol2inv(Vc)
    Vbeta <- solve(t(X) %*% Vinv %*% X)
    as.numeric(t(x0) %*% Vbeta %*% x0)
  }

  list(
    model = model, type = type,
    group_name = group_name, time_name = time_name, fe_time = fe_time,
    y = y, X = X, Z = Z, n = n, p = p, B = B, k_re = k_re, has_slope = has_slope,
    fixef = fe, vcov_beta = as.matrix(stats::vcov(model)),
    phi = phi, s2e = s2e,
    make_Sigma = make_Sigma, make_V = make_V, reml_nll = reml_nll, vF_fun = vF_fun
  )
}

## ---- Cov(phi) via inverse observed information  ---------
## SAFEGUARD (2): isSingular + near-boundary detection. cov is NULL when the
## Hessian fails or is not PD; `boundary` flags a small-but-positive component.
.vfun_vc_cov <- function(parsed, boundary_tol = 1e-4) {
  phi <- parsed$phi
  boundary <- .vfun_is_boundary(phi, boundary_tol)  # variances-only (see helper)

  H <- tryCatch(numDeriv::hessian(parsed$reml_nll, phi), error = function(e) NULL)
  ok <- !is.null(H) && all(is.finite(H))
  cov <- NULL; singular <- TRUE
  if (ok) {
    cov <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(cov) && all(is.finite(cov))) {
      ev <- tryCatch(eigen(cov, symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) NA_real_)
      singular <- any(!is.finite(ev)) || any(ev <= 0)
    } else {
      cov <- NULL
    }
  }
  if (is.null(cov)) singular <- TRUE

  list(cov = cov, phi = phi, boundary = boundary, singular = singular,
       ok = !is.null(cov))
}

#' Build a `vfun_components` from a fitted slope (or fixed-slope) lmer model
#'
#' The slope-engine extraction path: it parses a fitted `merMod` stability model
#' (random intercept + fixed slope, correlated random slope, or independent
#' random slope; see [reintervals-models]), reconstructs `Cov(phi)` from the
#' observed information (the numerical Hessian of the profiled REML negative
#' log-likelihood), evaluates the time-quadratic between-batch variance
#' `V_G(t0) = z0' Sigma z0` and the Satterthwaite/containment degrees of freedom
#' at a single time point `t0`, and returns a populated, validated
#' `vfun_components`.
#'
#' This is a separate path from `vc_extract()` (the EMS engine): it does not
#' route through it and does not relax its random-slope `stop()`. The slope
#' engine keeps its own singular handling (advisory `singular`/`boundary` flags)
#' and its own Satterthwaite -> containment df fallback; neither is shared with
#' the EMS path.
#'
#' @param model A fitted `lme4::lmer` stability model with one grouping factor
#'   and fixed part `~ 1 + time`: `y ~ time + (1 | g)` (random intercept + fixed
#'   slope), `y ~ time + (1 + time | g)` (correlated random slope), or
#'   `y ~ time + (1 | g) + (0 + time | g)` (independent random slope).
#' @param t0 A single time point at which to evaluate the payload.
#' @param target `"observable"` (default; content variance `V_T = V_G + s2e`,
#'   the future-observation target) or `"true_value"` (content variance `V_G`,
#'   the batch-mean target). Selects which content variance the effective sample
#'   size `n_E` is formed from; recorded as `attr(., "target")`.
#' @param df_method `"satterthwaite"` (default; generalized Satterthwaite with a
#'   containment fallback on degenerate df) or `"containment"` (integer df).
#' @param boundary_tol A component is flagged near-boundary if below this
#'   fraction of total variance. Default `1e-4`.
#'
#' @return A `vfun_components` object at `t0`, with `attr(., "target")` set.
#' @seealso `vfun_components()`, `vc_extract()` (the EMS path).
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'                  data = lme4::sleepstudy)
#' vfun_extract(fm, t0 = 5)
#' @noRd
vfun_extract <- function(model, t0,
                         target = c("observable", "true_value"),
                         df_method = c("satterthwaite", "containment"),
                         boundary_tol = 1e-4) {
  target    <- match.arg(target)
  df_method <- match.arg(df_method)
  if (!is.numeric(t0) || length(t0) != 1L || !is.finite(t0))
    stop("`t0` must be a single finite numeric value (the container holds one ",
         "time point).", call. = FALSE)

  parsed <- .vfun_parse(model)
  vc     <- .vfun_vc_cov(parsed, boundary_tol = boundary_tol)
  sing   <- isTRUE(tryCatch(lme4::isSingular(model), error = function(e) NA))

  phi   <- parsed$phi; type <- parsed$type; s2e <- parsed$s2e
  beta  <- parsed$fixef; Vbeta <- parsed$vcov_beta; fe_time <- parsed$fe_time
  cov   <- vc$cov

  ## fixed-effect design row in (Intercept, time) order
  x0 <- stats::setNames(numeric(length(beta)), names(beta))
  x0[["(Intercept)"]] <- 1; x0[[fe_time]] <- t0
  x0 <- x0[c("(Intercept)", fe_time)]

  mu  <- as.numeric(crossprod(x0, beta))
  V_F <- as.numeric(t(x0) %*% Vbeta %*% x0)
  ## V_G(t0) = z0' Sigma z0 is a quadratic form; for a PD Sigma it is >= 0 at
  ## every t0, so this normally returns the raw value untouched. The floor +
  ## flag is the defensive case: a model-implied
  ## NEGATIVE V_G at this t0 (the signed-s01 quadratic dipping below 0 -- NOT
  ## singularity, which is a fit-stage condition handled upstream). Flooring to
  ## 0 returns a finite (degenerate, possibly anti-conservative) interval where
  ## carrying it raw would give NaN; the flag/warning make that t0 auditable.
  vg <- .floor_vg(.vfun_value_VG(type, phi, t0), t0)
  V_G <- vg$value
  V_T <- V_G + s2e
  content_var <- if (target == "observable") V_T else V_G
  n_E <- if (V_F > 0) content_var / V_F else Inf

  ## df pieces, aligned to phi ordering, with the Satterthwaite->containment
  ## fallback (SAFEGUARD 3).
  gVG <- .vfun_grad_VG(type, t0)[names(phi)]
  gVT <- .vfun_grad_VT(type, t0)[names(phi)]
  gVe <- stats::setNames(as.numeric(names(phi) == "s2e"), names(phi))
  grad_vF <- numDeriv::grad(function(p) parsed$vF_fun(p, x0), phi)

  resolve_df <- function(value, grad_vec, w_between) {
    if (df_method == "containment") return(.vfun_df_containment(parsed, w_between))
    d <- .vfun_df_satterthwaite(value, grad_vec, cov)
    if (is.finite(d) && d > 0) return(d)
    .vfun_df_containment(parsed, w_between)        # fallback on degenerate df
  }
  nu_F <- resolve_df(V_F, grad_vF, w_between = 1)
  nu_T <- resolve_df(V_T, gVT, w_between = if (V_T > 0) V_G / V_T else 0)
  nu_G <- resolve_df(V_G, gVG, w_between = 1)
  nu_e <- resolve_df(s2e, gVe, w_between = 0)

  out <- vfun_components(
    type = type, phi = phi, cov = cov,
    mean = mu, var_mean = V_F,
    V_G = V_G, V_T = V_T, n_E = n_E,
    dfs = c(F = nu_F, T = nu_T, G = nu_G, e = nu_e),
    t0 = t0, singular = sing, boundary = isTRUE(vc$boundary)
  )
  attr(out, "target") <- target
  attr(out, "vg_floored") <- vg$floored
  out
}

## V_G floor-and-flag. raw >= -tol: pass through
## (max(0, .) absorbs round-off without flagging). raw < -tol: floor to 0 and
## flag `floored = TRUE` so the caller can warn (naming t0) and surface it on the
## result. Pure -> unit-testable without a fit. tol guards against flagging
## numerical round-off (~ -1e-15) as a genuine negative group variance.
.floor_vg <- function(raw_vg, t0, tol = 1e-8) {
  if (is.finite(raw_vg) && raw_vg < -tol) {
    warning(sprintf(
      "V_G(t0 = %g) is model-implied negative (%.3g); floored to 0. The interval ",
      t0, raw_vg),
      "at this t0 is degenerate and may be anti-conservative.", call. = FALSE)
    return(list(value = 0, floored = TRUE))
  }
  list(value = max(0, raw_vg), floored = FALSE)
}
