#' Canonical variance-component container
#'
#' Builds the `re_components` object that the method-pure interval constructors
#' (`ci_mean()`, `pi_newobs()`, `ti_francq()`) consume. It carries everything an
#' interval needs as plain numbers, so the constructors never touch `lme4`: the
#' fitted point `lb-hat` and its variance, the estimated variance components, the
#' degrees of freedom for each interval, and -- for the tolerance interval -- the
#' expected-mean-square (EMS) linear combination. 
#'
#' This slot layout is the shared contract across all three constructors and is
#' frozen:
#' \describe{
#'   \item{`components`, `coefs`}{the prediction interval reads
#'     `sigma^2_T = sum(coefs * components)` (`pi_newobs()`, Eq. 21).}
#'   \item{`dfs["ci"]`, `dfs["pi"]`}{the confidence- and prediction-interval
#'     degrees of freedom (Eq. 12; Eq. 23-24).}
#'   \item{`ems`}{the tolerance interval reads the EMS terms `ms`, their
#'     coefficients `k`, and their ANOVA degrees of freedom `df`
#'     (`ti_francq()`, Eq. 26).}
#' }
#'
#' @param components Named numeric vector of variance-component estimates
#'   `sigma^2_j` (residual included, e.g. `c(run = 0.000681, Residual = 0.001253)`).
#'   Negative ANOVA estimates are truncated at 0 with a warning.
#' @param dfs Named numeric vector of interval degrees of freedom. The reserved
#'   entries are `"ci"` (the CI denominator df, Eq. 12) and `"pi"` (the
#'   prediction-interval df for the total variance, Eq. 23-24). `Inf` is allowed
#'   (the residual-only df limit); values must be positive. Extra named entries
#'   are permitted and ignored.
#' @param mean Numeric scalar: the fitted linear combination `lb-hat` predicted at
#'   the requested `newdata` row.
#' @param var_mean Numeric scalar: `Var(lb-hat) = l C11 l'`, the variance of the
#'   fitted mean. Must be non-negative.
#' @param coefs Named numeric target selector `a_i` aligned to `components`: `1`
#'   includes a component in the total target variance, `0` excludes it. The
#'   "observable" target is all ones; the "true value" target zeroes the
#'   residual. Defaults to all ones.
#' @param ems Tolerance-interval EMS decomposition, or `NULL` when no TI applies
#'   (CI/PI-only fixtures, off-catalog designs). When supplied it is a list of
#'   three aligned named numeric vectors: `ms` (the EMS values `EMS_j`), `k`
#'   (their linear-combination coefficients `k_j`), and `df` (their ANOVA degrees
#'   of freedom `r_j`, finite and positive). The TI total variance is
#'   `sum(k * ms)`.
#' @param target Character scalar, `"observable"` (default) or `"true_value"` --
#'   the future-value distribution the PI/TI address. The extraction layer sets
#'   `coefs`/`ems` to match; recorded on the returned interval.
#' @param n_levels Optional named numeric of per-factor level counts (`NULL` if
#'   not supplied).
#' @param estimator Character scalar, `"reml"` (default) or `"anova"`.
#' @param design Character scalar naming the detected random structure
#'   (e.g. `"oneway"`, `"nested"`, `"crossed"`, or `""` when unknown).
#'
#' @return An object of class `re_components`: a list with the validated fields
#'   above.
#'
#' @examples
#' # Hand-built from Francq, Lin & Hoyer (2019) Section 5.2 (balanced one-way,
#' # 6 runs x 3 reps): run 0.000681, residual 0.001253, sigma^2_T 0.001934,
#' # df_PI 12.484, intercept 0.981 (SE 0.01353), CI df 5.
#' comp <- re_components(
#'   components = c(run = 0.000681, Residual = 0.001253),
#'   dfs        = c(ci = 5, pi = 12.484),
#'   mean       = 0.981,
#'   var_mean   = 0.01353^2,
#'   ems = list(
#'     ms = c(run = 0.001253 + 3 * 0.000681, Residual = 0.001253),
#'     k  = c(run = 1 / 3,                    Residual = 2 / 3),
#'     df = c(run = 5,                        Residual = 12)
#'   ),
#'   design = "oneway"
#' )
#' comp
#' @noRd
re_components <- function(components, dfs, mean, var_mean, coefs = NULL,
                          ems = NULL, target = "observable", n_levels = NULL,
                          estimator = "reml", design = "") {
  ## --- components ----------------------------------------------------------
  if (!is.numeric(components) || length(components) < 1L ||
        is.null(names(components)) || any(names(components) == "")) {
    stop("`components` must be a non-empty named numeric vector.", call. = FALSE)
  }
  if (any(!is.finite(components))) {
    stop("`components` must be finite.", call. = FALSE)
  }
  if (any(components < 0)) {
    neg <- names(components)[components < 0]
    warning("Negative variance-component estimate(s) truncated at 0: ",
            paste(neg, collapse = ", "), ".", call. = FALSE)
    components[components < 0] <- 0
  }

  ## --- dfs ------------------------------------------------------------------
  if (!is.numeric(dfs) || length(dfs) < 1L ||
        is.null(names(dfs)) || any(names(dfs) == "")) {
    stop("`dfs` must be a non-empty named numeric vector.", call. = FALSE)
  }
  if (any(is.na(dfs)) || any(dfs <= 0)) {
    stop("`dfs` must be positive (Inf allowed); NA is not permitted.",
         call. = FALSE)
  }

  ## --- coefs (target selector) ---------------------------------------------
  if (is.null(coefs)) {
    coefs <- stats::setNames(rep(1, length(components)), names(components))
  } else {
    if (!is.numeric(coefs) || is.null(names(coefs))) {
      stop("`coefs` must be a named numeric vector or NULL.", call. = FALSE)
    }
    if (!setequal(names(coefs), names(components))) {
      stop("`coefs` names must match `components` names.", call. = FALSE)
    }
    coefs <- coefs[names(components)]
  }

  ## --- mean / var_mean ------------------------------------------------------
  if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean)) {
    stop("`mean` must be a single finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(var_mean) || length(var_mean) != 1L ||
        !is.finite(var_mean) || var_mean < 0) {
    stop("`var_mean` must be a single finite non-negative numeric value.",
         call. = FALSE)
  }
  mean <- unname(mean)            # keep interval bounds free of stray names
  var_mean <- unname(var_mean)

  ## --- ems (TI EMS decomposition) ------------------------------------------
  if (!is.null(ems)) {
    ems <- validate_ems(ems)
  }

  ## --- target / estimator / design -----------------------------------------
  target    <- match.arg(target, c("observable", "true_value"))
  estimator <- match.arg(estimator, c("reml", "anova"))
  if (!is.character(design) || length(design) != 1L) {
    stop("`design` must be a character scalar.", call. = FALSE)
  }
  if (!is.null(n_levels) && !is.numeric(n_levels)) {
    stop("`n_levels` must be NULL or a numeric vector.", call. = FALSE)
  }

  structure(
    list(components = components, dfs = dfs, coefs = coefs,
         mean = mean, var_mean = var_mean, ems = ems, target = target,
         n_levels = n_levels, estimator = estimator, design = design),
    class = "re_components"
  )
}

## Validate and normalize the `ems` sub-list. Internal.
validate_ems <- function(ems) {
  if (!is.list(ems) || !setequal(names(ems), c("ms", "k", "df"))) {
    stop('`ems` must be NULL or a list with elements "ms", "k", "df".',
         call. = FALSE)
  }
  for (nm in c("ms", "k", "df")) {
    v <- ems[[nm]]
    if (!is.numeric(v) || length(v) < 1L || is.null(names(v)) ||
          any(names(v) == "")) {
      stop(sprintf("`ems$%s` must be a non-empty named numeric vector.", nm),
           call. = FALSE)
    }
    if (any(!is.finite(v))) {
      stop(sprintf("`ems$%s` must be finite.", nm), call. = FALSE)
    }
  }
  if (!setequal(names(ems$ms), names(ems$k)) ||
        !setequal(names(ems$ms), names(ems$df))) {
    stop("`ems$ms`, `ems$k`, and `ems$df` must share the same names.",
         call. = FALSE)
  }
  if (any(ems$df <= 0)) {
    stop("`ems$df` (ANOVA degrees of freedom) must be positive.", call. = FALSE)
  }
  if (any(ems$ms < 0)) {
    stop("`ems$ms` (expected mean squares) must be non-negative.", call. = FALSE)
  }
  ord <- names(ems$ms)
  list(ms = ems$ms, k = ems$k[ord], df = ems$df[ord])
}

#' @export
print.re_components <- function(x, ...) {
  cat("<re_components>", if (nzchar(x$design)) paste0("[", x$design, "]"), "\n")
  cat("  estimator :", x$estimator, " target:", x$target, "\n")
  cat("  mean      :", format(x$mean), " var_mean:", format(x$var_mean), "\n")
  cat("  components:",
      paste(names(x$components), signif(x$components, 4), sep = "=",
            collapse = ", "), "\n")
  cat("  dfs       :",
      paste(names(x$dfs), signif(x$dfs, 4), sep = "=", collapse = ", "), "\n")
  if (is.null(x$ems)) {
    cat("  ems       : <none> (no tolerance interval)\n")
  } else {
    cat("  ems$ms    :",
        paste(names(x$ems$ms), signif(x$ems$ms, 4), sep = "=",
              collapse = ", "), "\n")
  }
  invisible(x)
}
