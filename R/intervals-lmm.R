#' Confidence and prediction intervals for a linear mixed model
#'
#' Confidence interval (Francq, Lin & Hoyer 2019, Eq. 12) and prediction interval
#' for a future observation (Eq. 21). Both are general across random-intercept
#' designs, balanced or unbalanced: the prediction-interval degrees of freedom
#' come from the observed Fisher information (the numerical Hessian of the
#' REML log-likelihood on the variance-component scale, with beta at its GLS
#' value, Eq. 23-24), not from a balanced-design formula.
#'
#' This is a thin wrapper: it reaches `lme4` only through `vc_extract()`'s
#' internals and delegates the arithmetic to the method-pure constructors
#' `ci_mean()` and `pi_newobs()`.
#'
#' @param model An `lmerMod` fit (random intercepts only; random slopes are
#'   rejected).
#' @param newdata Row(s) giving the fixed-effect combination(s) to predict at.
#'   Defaults to the first model-frame row.
#' @param level Prediction level (`1 - psi`) for the PI. Default 0.95.
#' @param conf Confidence level for the CI. Default 0.95.
#' @param which Character subset of `c("CI", "PI")`. Default both.
#' @param ci_df Optional numeric overriding the CI denominator df (internal).
#'   `NULL` (default) uses `ci_method`.
#' @param ci_method CI df rule: `"satterthwaite"` (default; Satterthwaite with
#'   the observed-information covariance of the variance components, falling
#'   back to the SAS-style containment df when it is unavailable) or
#'   `"kenward-roger"` (requires \pkg{pbkrtest}; errors otherwise).
#'
#' @return A data frame with the `newdata` predictor columns plus `type`,
#'   `estimate`, `lower`, `upper`, `df`, `df_method` (the CI df source:
#'   `"satterthwaite"`, `"containment"`, or `"kenward-roger"`; `NA` on PI rows),
#'   `level`. When the PI is computed the
#'   variance components are attached as `attr(., "components")` and
#'   `attr(., "pi_df_fallback")` flags the singular-fit df guard (see
#'   `vc_extract()`). `attr(., "singular")` (`lme4::isSingular()`) is always
#'   set, for the CI and the PI alike.
#'
#' @seealso [ti_lmm()] for tolerance intervals; `ci_mean()`, `pi_newobs()`.
#' @references Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and
#'   tolerance in linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' intervals_lmm(fm, newdata = data.frame(Days = c(0, 9)))
#' @noRd
intervals_lmm <- function(model, newdata = NULL, level = 0.95, conf = 0.95,
                          which = c("CI", "PI"), ci_df = NULL,
                          ci_method = c("satterthwaite", "kenward-roger")) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted `lmerMod` (an lmer fit).", call. = FALSE)
  }
  which <- match.arg(which, c("CI", "PI"), several.ok = TRUE)
  ci_method <- match.arg(ci_method)
  fx <- .fixed(model, newdata)
  ## theta (observed information) is needed for the PI df and for the default
  ## Satterthwaite CI df; a Kenward-Roger or overridden CI-only call skips it.
  need_theta <- "PI" %in% which ||
    ("CI" %in% which && is.null(ci_df) && ci_method == "satterthwaite")
  theta <- if (need_theta) .theta_cov(model) else NULL
  comps <- if (!is.null(theta)) theta$components else .var_components(model)
  nr <- nrow(fx$L)
  base <- fx$newdata

  ## one container per prediction row, reusing the model-level components
  make_comp <- function(i, ci_ddf) {
    dfs <- if (!is.null(theta)) c(ci = ci_ddf, pi = theta$df_pi) else c(ci = ci_ddf)
    re_components(components = comps, dfs = dfs, mean = fx$fit[i],
                  var_mean = fx$var_fix[i], target = "observable")
  }

  blocks <- list()
  if ("CI" %in% which) {
    dd <- lapply(seq_len(nr), function(i)
      .ci_ddf(model, fx$L[i, ], theta = theta, ci_df = ci_df, df_method = ci_method))
    ddf <- vapply(dd, function(x) x$df, numeric(1))
    cis <- lapply(seq_len(nr), function(i) ci_mean(make_comp(i, ddf[i]), conf = conf))
    blocks$CI <- cbind(base, data.frame(
      type = "CI",
      estimate = vapply(cis, function(x) x$estimate, numeric(1)),
      lower = vapply(cis, function(x) x$lower, numeric(1)),
      upper = vapply(cis, function(x) x$upper, numeric(1)),
      df = ddf, df_method = vapply(dd, function(x) x$method, character(1)),
      level = conf, stringsAsFactors = FALSE
    ), row.names = NULL)
  }
  if ("PI" %in% which) {
    ## ci df is irrelevant to the PI; the container carries a placeholder
    pis <- lapply(seq_len(nr), function(i) pi_newobs(make_comp(i, Inf), level = level))
    blocks$PI <- cbind(base, data.frame(
      type = "PI",
      estimate = vapply(pis, function(x) x$estimate, numeric(1)),
      lower = vapply(pis, function(x) x$lower, numeric(1)),
      upper = vapply(pis, function(x) x$upper, numeric(1)),
      df = theta$df_pi, df_method = NA_character_, level = level,
      stringsAsFactors = FALSE
    ), row.names = NULL)
  }

  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  if ("PI" %in% which) {
    attr(out, "components") <- theta$components
    attr(out, "pi_df_fallback") <- isTRUE(theta$pi_df_fallback)
  }
  ## set whether or not theta was computed (a CI-only call skips .theta_cov), so
  ## the CI reports lme4::isSingular() like the PI / new-group-mean / TI.
  attr(out, "singular") <- .is_singular(model)
  out
}
