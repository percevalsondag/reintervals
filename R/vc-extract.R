# =============================================================================
# Extraction layer: the ONLY place reintervals touches lme4 / Matrix / numDeriv.
# Constructors (ci-mean.R, pi-newobs.R, ti-francq.R) and design adapters
# (design-*.R) stay method-pure. The wrappers (intervals_lmm / ti_lmm) reach
# lme4 ONLY through vc_extract() and the internal .fixed() helper here.
#
# .theta_cov() computes the covariance of the variance-component estimates,
# following the Francq, Lin & Hoyer (2019) observed-information construction. Do not rewrite the Hessian
# construction: it must reproduce df_PI to full precision (12.486 on the
# Section 5.2 balanced fixture).
# =============================================================================

## ---- scope guard: random intercepts only -----------------------------------
.assert_random_intercepts <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  if (!all(is.na(vc$var2))) {
    stop("Only random-intercept models are supported (no random slopes or ",
         "correlated random effects).", call. = FALSE)
  }
  vc
}

.level_counts <- function(flist) vapply(flist, nlevels, numeric(1))

## Cheap variance components from VarCorr (no Hessian) -- same naming/order as
## .theta_cov. Used by the wrappers when the prediction-interval df is not
## needed (CI-only intervals_lmm; the normal-based ti_lmm).
.var_components <- function(model) {
  vc <- as.data.frame(lme4::VarCorr(model))
  zl <- lme4::getME(model, "Ztlist")
  gn <- sub("\\.\\(Intercept\\)$", "", names(zl))
  stats::setNames(vc$vcov[match(c(gn, "Residual"), vc$grp)], c(gn, "Residual"))
}

.design_of <- function(model) {
  classify_design(lme4::getME(model, "flist"), length(lme4::getME(model, "y")))
}

## Cell-count table for the unbalanced EMS synthesis (plain data for the pure
## adapters). Crossed: an A-by-B integer matrix n_ij with rows = the grp_A factor
## and cols = grp_B (orientation is load-bearing -- the synthesis is asymmetric).
## Nested: a list whose i-th element is the replicate counts of the fine levels
## within coarse level i.
.cell_counts <- function(model, desc) {
  fl <- lme4::getME(model, "flist")
  if (desc$type == "crossed") {
    tab <- table(fl[[desc$grp_A]], fl[[desc$grp_B]])
    matrix(as.integer(tab), nrow = nrow(tab))
  } else {
    tab <- table(fl[[desc$grp_alpha]], fl[[desc$grp_beta]])
    lapply(seq_len(nrow(tab)), function(i) {
      r <- as.integer(tab[i, ])
      r[r > 0]
    })
  }
}

## counts needed only for unbalanced nested/crossed (synthesis); NULL otherwise
.counts_if_unbalanced <- function(model, desc) {
  if (desc$type %in% c("nested", "crossed") && !isTRUE(desc$balanced)) {
    .cell_counts(model, desc)
  } else {
    NULL
  }
}

.is_singular <- function(model) isTRUE(lme4::isSingular(model))

## ---- bounded Satterthwaite df for sigma_T^2 (singular/boundary fallback) -----
## Bounded generalized-Satterthwaite df for the total variance. Used only when
## the observed-information route is unreliable; keeps df_PI finite as a component -> 0.
.satterthwaite_df <- function(model, sigmaT2) {
  des <- tryCatch(
    classify_design(lme4::getME(model, "flist"), length(lme4::getME(model, "y"))),
    error = function(e) list(type = "unknown")
  )
  vc <- as.data.frame(lme4::VarCorr(model))
  s2e <- vc$vcov[vc$grp == "Residual"]
  comp_df <- function(ms, df) if (df > 0 && ms > 0) ms^2 / df else 0

  terms <- NULL
  if (des$type == "oneway") {
    A <- des$A; N <- des$N; n0 <- des$n0
    s2a <- vc$vcov[vc$grp == des$grp]
    emsa <- s2e + n0 * s2a
    terms <- c(comp_df((1 / n0) * emsa, A - 1),
               comp_df((1 - 1 / n0) * s2e, N - A))
  } else if (des$type == "nested" && isTRUE(des$balanced)) {
    A <- des$A; B <- des$B; n <- des$n
    s2a <- vc$vcov[vc$grp == des$grp_alpha]
    s2b <- vc$vcov[vc$grp == des$grp_beta]
    emsb <- s2e + n * s2b
    emsa <- s2e + n * s2b + n * B * s2a
    terms <- c(comp_df((1 / (n * B)) * emsa, A - 1),
               comp_df((1 / n - 1 / (n * B)) * emsb, A * (B - 1)),
               comp_df((1 - 1 / n) * s2e, A * B * (n - 1)))
  } else if (des$type == "crossed" && isTRUE(des$balanced)) {
    A <- des$A; B <- des$B; n <- des$n
    v1 <- vc$vcov[vc$grp == des$grp_A]
    v2 <- vc$vcov[vc$grp == des$grp_B]
    s2ab <- vc$vcov[vc$grp == des$grp_AB]
    emsab <- s2e + n * s2ab
    emsa <- s2e + n * s2ab + n * B * v1
    emsb <- s2e + n * s2ab + n * A * v2
    terms <- c(comp_df((1 / (n * B)) * emsa, A - 1),
               comp_df((1 / (n * A)) * emsb, B - 1),
               comp_df((1 / n - 1 / (n * A) - 1 / (n * B)) * emsab, (A - 1) * (B - 1)),
               comp_df((1 - 1 / n) * s2e, A * B * (n - 1)))
  } else {
    N <- length(lme4::getME(model, "y"))
    p <- length(lme4::fixef(model))
    q <- sum(vc$grp != "Residual")
    return(max(1, N - p - q))
  }
  denom <- sum(terms)
  if (!is.finite(denom) || denom <= 0) return(Inf)
  sigmaT2^2 / denom
}

## ---- total variance, its variance, df_PI (observed Fisher information) -------
## Cov(theta_hat) from the observed information: a numerical Hessian
## of the profiled REML log-likelihood (== PROC MIXED ASYCOV). 
.theta_cov <- function(model) {
  vc <- .assert_random_intercepts(model)
  Zl <- lme4::getME(model, "Ztlist")
  gn <- sub("\\.\\(Intercept\\)$", "", names(Zl))          # grouping name per block
  Al <- lapply(Zl, function(zt) Matrix::tcrossprod(Matrix::t(zt)))  # A_g = Z_g Z_g'
  th <- vc$vcov[match(c(gn, "Residual"), vc$grp)]          # block order; residual last
  k <- length(th)
  y <- lme4::getME(model, "y")
  X <- lme4::getME(model, "X")
  N <- length(y)
  ll <- function(t) {
    V <- as.matrix(Reduce(`+`, Map(function(a, s) s * a, Al, t[-k])) +
                     Matrix::Diagonal(N, t[k]))
    lc <- chol(V)
    vi <- chol2inv(lc)
    xtvi <- crossprod(X, vi)
    xtvix <- xtvi %*% X
    b <- solve(xtvix, xtvi %*% y)
    r <- y - X %*% b
    -0.5 * (2 * sum(log(diag(lc))) +
              as.numeric(determinant(xtvix, logarithm = TRUE)$modulus) +
              as.numeric(crossprod(r, vi %*% r)))
  }
  sigmaT2 <- sum(th)
  singular <- isTRUE(lme4::isSingular(model))
  hmat <- if (singular) NULL else tryCatch(numDeriv::hessian(ll, th),
                                           error = function(e) NULL)
  cov_theta <- if (is.null(hmat)) NULL else tryCatch(solve(-hmat),
                                                     error = function(e) NULL)
  vS <- if (is.null(cov_theta)) NA_real_ else sum(cov_theta)
  use_hessian <- isTRUE(is.finite(vS) && vS > 0) && !singular
  df_pi <- if (use_hessian) 2 * sigmaT2^2 / vS else .satterthwaite_df(model, sigmaT2)
  list(sigmaT2 = sigmaT2,
       df_pi = df_pi,
       pi_df_fallback = !use_hessian,
       singular = singular,
       components = stats::setNames(th, c(gn, "Residual")),
       cov_theta = if (use_hessian) cov_theta else NULL,
       parts = list(Al = Al, th = th, X = X, N = N, k = k))
}

## ---- fixed-effect prediction l'beta and Var(l'beta) = l C11 l' --------------
## Port of `.lmm_fixed`.
## Validate user-supplied `newdata` before it reaches model.matrix(), which
## otherwise fails with an opaque "object '<var>' not found". Guards MALFORMED
## input only; valid `newdata` (and NULL) pass through unchanged, so no
## valid-input interval is affected.
.check_newdata <- function(model, newdata) {
  if (is.null(newdata)) return(invisible(NULL))
  if (!is.data.frame(newdata) && !is.list(newdata)) {
    stop("`newdata` must be a data frame with one column per model predictor; ",
         "got ", class(newdata)[1L], ". ",
         "Supply e.g. data.frame(", paste(all.vars(stats::delete.response(
           stats::terms(model, fixed.only = TRUE))), collapse = " = , "),
         " = ).", call. = FALSE)
  }
  need <- all.vars(stats::delete.response(stats::terms(model, fixed.only = TRUE)))
  miss <- setdiff(need, names(newdata))
  if (length(miss)) {
    stop("`newdata` is missing required predictor(s): ",
         paste(miss, collapse = ", "),
         ". `newdata` must supply a column for each model predictor (",
         paste(need, collapse = ", "), ").", call. = FALSE)
  }
  invisible(NULL)
}

.fixed <- function(model, newdata = NULL) {
  .check_newdata(model, newdata)
  if (is.null(newdata)) newdata <- stats::model.frame(model)[1, , drop = FALSE]
  terms <- stats::delete.response(stats::terms(model, fixed.only = TRUE))
  l <- stats::model.matrix(terms, newdata)
  c11 <- as.matrix(stats::vcov(model))
  list(newdata = newdata, L = l,
       fit = as.vector(l %*% lme4::fixef(model)),
       var_fix = diag(l %*% c11 %*% t(l)))
}

## ---- CI denominator df: Kenward-Roger, else Satterthwaite (never NA) --------
.ci_ddf <- function(model, Lrow, theta = NULL, ci_df = NULL) {
  if (!is.null(ci_df)) return(as.numeric(ci_df))
  if (requireNamespace("pbkrtest", quietly = TRUE)) {
    ddf <- tryCatch(pbkrtest::get_Lb_ddf(model, as.numeric(Lrow)),
                    error = function(e) NA_real_)
    if (isTRUE(is.finite(ddf) && ddf > 0)) return(ddf)
  }
  .ci_ddf_satterthwaite(model, Lrow, theta)
}

## Satterthwaite df for the contrast variance l C11 l', reusing the observed-
## information covariance from .theta_cov(). Falls back to the residual df when
## that covariance is unavailable (e.g. singular fit). Never returns NA.
.ci_ddf_satterthwaite <- function(model, Lrow, theta) {
  l <- as.numeric(Lrow)
  resid_df <- function() {
    max(1, length(lme4::getME(model, "y")) - length(lme4::fixef(model)))
  }
  if (is.null(theta) || is.null(theta$cov_theta)) return(resid_df())
  Al <- theta$parts$Al; th <- theta$parts$th
  X <- theta$parts$X; N <- theta$parts$N; k <- theta$parts$k
  cvar <- function(t) {
    V <- as.matrix(Reduce(`+`, Map(function(a, s) s * a, Al, t[-k])) +
                     Matrix::Diagonal(N, t[k]))
    c11 <- solve(crossprod(X, solve(V, X)))
    as.numeric(t(l) %*% c11 %*% l)
  }
  g <- tryCatch(numDeriv::grad(cvar, th), error = function(e) NULL)
  if (is.null(g)) return(resid_df())
  var_c <- as.numeric(t(g) %*% theta$cov_theta %*% g)
  if (!is.finite(var_c) || var_c <= 0) return(resid_df())
  2 * cvar(th)^2 / var_c
}

#' Extract variance components and degrees of freedom from a fitted model
#'
#' The single bridge from a fitted `lmerMod` to the method-pure interval
#' constructors: it is the only function in the package that touches `lme4`,
#' `Matrix`, or `numDeriv`. It rejects random slopes, detects the design,
#' estimates the variance components and the observed-information
#' prediction-interval degrees of freedom (the numerical Hessian of the profiled
#' REML log-likelihood, Francq 2019 Eq. 23-24), builds the design's
#' expected-mean-square decomposition for the tolerance interval, and returns a
#' canonical `re_components` object evaluated at one prediction point.
#'
#' The returned object carries two attributes: `singular`
#' (`lme4::isSingular()`) and `pi_df_fallback` (`TRUE` when the prediction df
#' came from the bounded-Satterthwaite guard rather than the Hessian --- the
#' singular-fit width-explosion fix). When the design has no closed-form EMS
#' decomposition (unbalanced nested/crossed in v1, or off-catalog), `ems` is
#' `NULL` and a `note` attribute explains why; the CI and PI are still returned.
#'
#' @param model A fitted `lmerMod` (random intercepts only).
#' @param design Optional pre-computed design descriptor (as from
#'   `classify_design()`); `NULL` (default) detects it from the fit.
#' @param target `"observable"` (default; future-value variance = all
#'   components) or `"true_value"` (between-level components only).
#' @param newdata Optional one-row data frame giving the fixed-effect
#'   combination `l` to evaluate the mean at; `NULL` uses the first model-frame
#'   row.
#' @return An `re_components` object with `singular`/`pi_df_fallback` attributes.
#' @seealso `ci_mean()`, `pi_newobs()`, `ti_francq()`
#' @noRd
vc_extract <- function(model, design = NULL,
                       target = c("observable", "true_value"),
                       newdata = NULL) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted `lmerMod` (from lme4::lmer).", call. = FALSE)
  }
  target <- match.arg(target)
  .assert_random_intercepts(model)
  flist <- lme4::getME(model, "flist")
  n_obs <- length(lme4::getME(model, "y"))
  desc <- if (is.null(design)) classify_design(flist, n_obs) else design

  theta <- .theta_cov(model)
  dc <- .design_components(desc, theta$components, target,
                          .counts_if_unbalanced(model, desc))
  fx <- .fixed(model, newdata)
  ci_ddf <- .ci_ddf(model, fx$L[1, ], theta = theta)

  comp <- re_components(
    components = theta$components,
    dfs        = c(ci = ci_ddf, pi = theta$df_pi),
    mean       = fx$fit[1],
    var_mean   = fx$var_fix[1],
    coefs      = dc$coefs,
    ems        = dc$ems,
    target     = target,
    n_levels   = .level_counts(flist),
    estimator  = "reml",
    design     = desc$type
  )
  attr(comp, "singular") <- theta$singular
  attr(comp, "pi_df_fallback") <- theta$pi_df_fallback
  if (!is.null(dc$note)) attr(comp, "note") <- dc$note
  comp
}
