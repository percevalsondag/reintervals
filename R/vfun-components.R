#' Variance-function container for the slope (Montes / GPQ) engine
#'
#' The parallel container to the frozen `re_components` object. Where
#' `re_components` carries the expected-mean-square (EMS) representation for the
#' random-intercept Francq engine, `vfun_components` carries the
#' *variance-function* representation for the random-slope stability models
#' (all three slope-model classes): the variance-parameter vector `phi`, its REML-Hessian
#' covariance `cov`, and the quadratic-form payload evaluated at a single time
#' point `t0`.
#'
#' This is a **replace**-architecture container, not a generalization of
#' `re_components`, for one load-bearing reason: the slope representation
#' includes the intercept--slope **covariance** `s01`, which can legitimately be
#' **negative**. `re_components` truncates every component at zero; applying that
#' here would corrupt `s01`. The `vfun_components` validator therefore enforces
#' sign rules *by parameter role* --- variances non-negative, the covariance
#' unrestricted --- rather than the blanket truncation.
#'
#' The slot layout mirrors the Montes/GPQ representation:
#' \describe{
#'   \item{`type`}{Model code (see [reintervals-models] for the code-to-model
#'     mapping): random intercept, correlated random slope, or independent
#'     random slope.}
#'   \item{`phi`}{The variance-parameter vector on the natural variance scale,
#'     named in canonical order per `type` (`s2_0` intercept variance, `s01`
#'     intercept--slope covariance, `s2_1` slope variance, `s2e` residual).
#'     `s2_0`, `s2_1`, `s2e` are variances (>= 0); `s01` **may be negative**.}
#'   \item{`cov`}{`Cov(phi)`, the inverse-observed-information (REML-Hessian)
#'     reconstruction; a symmetric `length(phi) x length(phi)` matrix, or `NULL`
#'     when the fit is singular and the Hessian is unavailable.}
#'   \item{`mean`, `var_mean`}{The per-`t0` mean-line scalars `mu(t0) = x0'beta`
#'     and `V_F(t0) = x0' Vbeta x0`.}
#'   \item{`V_G`, `V_T`, `n_E`}{The quadratic-form payload at `t0`:
#'     `V_G(t0) = z0' Sigma z0` (between-batch variance, residual excluded),
#'     `V_T = V_G + s2e` (future-observation variance), and the effective sample
#'     size `n_E = V_T / V_F` (`Inf` when `V_F = 0`).}
#'   \item{`dfs`}{Named numeric of the Satterthwaite / containment degrees of
#'     freedom (e.g. the mean-line, future-obs, between, and residual df). `Inf`
#'     allowed; values positive.}
#'   \item{`t0`}{The single time point at which the payload is evaluated.}
#'   \item{`singular`}{`lme4::isSingular()` flag --- the fit is *exactly*
#'     degenerate (a variance component on the boundary at zero).}
#'   \item{`boundary`}{Near-boundary flag --- a variance component is
#'     small-but-positive (below `boundary_tol` of total variance). Distinct
#'     from `singular`: this is the regime where the Satterthwaite df collapses
#'     (Karl, Rushing, Burdick & Hofer 2026) and the batch-mean closed form
#'     breaks, so the dispatch layer must be able to see it.}
#' }
#'
#' @param type Character model code (see [reintervals-models]).
#' @param phi Named numeric variance-parameter vector (see layout above).
#' @param mean Numeric scalar `mu(t0)`.
#' @param var_mean Numeric scalar `V_F(t0)` (>= 0).
#' @param V_G Numeric scalar `V_G(t0)`. Finite; normally `>= 0`, but carried
#'   raw (it can be negative at a non-PD / boundary `Sigma`, because it is the
#'   raw quadratic form --- PSD-ness is guarded upstream at the fit).
#' @param V_T Numeric scalar `V_T(t0) = V_G + s2e`. Finite; carried raw.
#' @param n_E Numeric scalar effective sample size (> 0; `Inf` allowed).
#' @param dfs Named numeric vector of degrees of freedom (positive; `Inf`
#'   allowed; no `NA`).
#' @param t0 Numeric scalar time point.
#' @param cov `Cov(phi)` matrix, or `NULL` (singular fit). Default `NULL`.
#' @param singular Logical scalar; the fit is exactly degenerate. Default `FALSE`.
#' @param boundary Logical scalar; a variance component is small-but-positive
#'   (near-boundary). Distinct from `singular`. Default `FALSE`.
#'
#' @return An object of class `vfun_components`.
#'
#' @seealso `re_components()` for the frozen EMS container.
#' @examples
#' # Hand-built correlated-random-slope payload at t0 = 24 with a NEGATIVE intercept-slope
#' # covariance s01 -- the validator accepts it (it would be truncated by
#' # re_components, which is exactly why this is a separate container).
#' vc <- vfun_components(
#'   type = "M2c",
#'   phi  = c(s2_0 = 2.0, s01 = -0.05, s2_1 = 0.01, s2e = 0.7),
#'   mean = 98.5, var_mean = 0.12,
#'   V_G = 2.0 + 2 * 24 * -0.05 + 24^2 * 0.01,   # z0' Sigma z0 at t0 = 24
#'   V_T = (2.0 + 2 * 24 * -0.05 + 24^2 * 0.01) + 0.7,
#'   n_E = 30,
#'   dfs = c(F = 8, T = 6.4, G = 5.1, e = 40),
#'   t0 = 24
#' )
#' vc
#' @noRd
vfun_components <- function(type, phi, mean, var_mean, V_G, V_T, n_E, dfs,
                            t0, cov = NULL, singular = FALSE, boundary = FALSE) {
  validate_vfun_components(
    type = type, phi = phi, cov = cov, mean = mean, var_mean = var_mean,
    V_G = V_G, V_T = V_T, n_E = n_E, dfs = dfs, t0 = t0, singular = singular,
    boundary = boundary
  )
}

## canonical phi names per model type (matches the Montes/GPQ source parse)
.vfun_phi_names <- function(type) {
  switch(type,
    M1  = c("s2_0", "s2e"),
    M2c = c("s2_0", "s01", "s2_1", "s2e"),
    M2i = c("s2_0", "s2_1", "s2e"),
    stop("`type` must be a recognized model code (see ?reintervals-models).", call. = FALSE)
  )
}

## Validate and assemble a vfun_components object (internal constructor). The
## sign rules are the whole point (see the class doc):
## variances (s2_*) must be non-negative, the covariance s01 is unrestricted.
validate_vfun_components <- function(type, phi, cov, mean, var_mean,
                                     V_G, V_T, n_E, dfs, t0, singular,
                                     boundary = FALSE) {
  ## --- type --------------------------------------------------------------
  if (!is.character(type) || length(type) != 1L ||
        !type %in% c("M1", "M2c", "M2i")) {
    stop("`type` must be a single recognized model code (see ?reintervals-models).",
         call. = FALSE)
  }
  expected <- .vfun_phi_names(type)

  ## --- phi: names, order, finiteness, and ROLE-AWARE sign rules ----------
  if (!is.numeric(phi) || is.null(names(phi))) {
    stop("`phi` must be a named numeric vector.", call. = FALSE)
  }
  if (!identical(names(phi), expected)) {
    stop(sprintf(paste("`phi` names must be exactly c(%s) in that order for this",
                       "model code (see ?reintervals-models)."),
                 paste0("\"", expected, "\"", collapse = ", ")),
         call. = FALSE)
  }
  if (any(!is.finite(phi))) {
    stop("`phi` entries must be finite.", call. = FALSE)
  }
  ## variances (everything except the covariance s01) must be non-negative;
  ## s01 is a covariance and MAY be negative -- do NOT truncate it.
  var_names <- setdiff(names(phi), "s01")
  neg_var <- var_names[phi[var_names] < 0]
  if (length(neg_var)) {
    stop("Variance component(s) must be non-negative: ",
         paste(neg_var, collapse = ", "),
         ". (The intercept-slope covariance `s01` may be negative; ",
         "variances may not.)", call. = FALSE)
  }

  ## --- cov: NULL (singular) or a symmetric length(phi)-square matrix ------
  if (!is.null(cov)) {
    if (!is.matrix(cov) || !is.numeric(cov)) {
      stop("`cov` must be a numeric matrix or NULL.", call. = FALSE)
    }
    k <- length(phi)
    if (nrow(cov) != k || ncol(cov) != k) {
      stop(sprintf("`cov` must be %d x %d to match `phi`; got %d x %d.",
                   k, k, nrow(cov), ncol(cov)), call. = FALSE)
    }
    if (any(!is.finite(cov))) {
      stop("`cov` entries must be finite (use `cov = NULL` for a singular fit).",
           call. = FALSE)
    }
    if (!isSymmetric(unname(cov), tol = 1e-8)) {
      stop("`cov` must be symmetric (it is a covariance matrix).", call. = FALSE)
    }
  }

  ## --- scalar payload -----------------------------------------------------
  .scalar <- function(x, nm, nonneg = FALSE, positive = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be a single non-NA numeric value.", nm),
           call. = FALSE)
    }
    if (positive && !(x > 0)) {
      stop(sprintf("`%s` must be positive (Inf allowed).", nm), call. = FALSE)
    }
    if (nonneg && (is.finite(x) && x < 0)) {
      stop(sprintf("`%s` must be non-negative.", nm), call. = FALSE)
    }
    if (!positive && nm != "n_E" && !is.finite(x)) {
      stop(sprintf("`%s` must be finite.", nm), call. = FALSE)
    }
    unname(x)
  }
  mean     <- .scalar(mean, "mean")
  var_mean <- .scalar(var_mean, "var_mean", nonneg = TRUE)  # V_F = x0'Vbeta x0, PSD by construction
  ## V_G(t0) = z0' Sigma z0 and V_T = V_G + s2e are carried RAW, not forced >= 0.
  ## PSD-ness (which normally keeps V_G >= 0) is guarded upstream at the fit
  ## (the eigen sentinel in vfun_extract), not re-imposed here: a non-PD or
  ## boundary Sigma can give a negative V_G, and the container carries it as-is.
  ## Forcing V_G >= 0 here while permitting negative s01 would be inconsistent
  ## (permit the cause, forbid the effect), so don't add a floor here.
  V_G      <- .scalar(V_G, "V_G")                          # finite; may be < 0
  V_T      <- .scalar(V_T, "V_T")                          # finite; may be < 0
  n_E      <- .scalar(n_E, "n_E", positive = TRUE)        # Inf allowed
  t0       <- .scalar(t0, "t0")

  ## --- dfs: same contract as re_components$dfs (positive, Inf ok, no NA) --
  if (!is.numeric(dfs) || length(dfs) < 1L ||
        is.null(names(dfs)) || any(names(dfs) == "")) {
    stop("`dfs` must be a non-empty named numeric vector.", call. = FALSE)
  }
  if (any(is.na(dfs)) || any(dfs <= 0)) {
    stop("`dfs` must be positive (Inf allowed); NA is not permitted.",
         call. = FALSE)
  }

  ## --- singular / boundary flags ------------------------------------------
  if (!is.logical(singular) || length(singular) != 1L || is.na(singular)) {
    stop("`singular` must be a single TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(boundary) || length(boundary) != 1L || is.na(boundary)) {
    stop("`boundary` must be a single TRUE/FALSE.", call. = FALSE)
  }

  structure(
    list(type = type, phi = phi, cov = cov,
         mean = mean, var_mean = var_mean,
         V_G = V_G, V_T = V_T, n_E = n_E,
         dfs = dfs, t0 = t0, singular = singular, boundary = boundary),
    class = "vfun_components"
  )
}

#' @export
print.vfun_components <- function(x, ...) {
  flags <- c(if (isTRUE(x$singular)) "singular",
             if (isTRUE(x$boundary)) "near-boundary")
  cat("<vfun_components>", paste0("[", x$type, "]"),
      if (length(flags)) paste0("(", paste(flags, collapse = ", "), ")"), "\n")
  cat("  phi      :",
      paste(names(x$phi), signif(x$phi, 4), sep = "=", collapse = ", "), "\n")
  cat("  cov      :",
      if (is.null(x$cov)) "<NULL> (singular / Hessian unavailable)"
      else sprintf("%d x %d", nrow(x$cov), ncol(x$cov)), "\n")
  cat("  @ t0     :", format(x$t0), "\n")
  cat("  mean     :", format(x$mean), " var_mean:", format(x$var_mean), "\n")
  cat("  V_G      :", format(x$V_G), " V_T:", format(x$V_T),
      " n_E:", format(x$n_E), "\n")
  cat("  dfs      :",
      paste(names(x$dfs), signif(x$dfs, 4), sep = "=", collapse = ", "), "\n")
  invisible(x)
}
