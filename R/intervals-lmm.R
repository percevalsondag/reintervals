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
#'   back to the SAS-style containment df when it is unavailable or below 1) or
#'   `"kenward-roger"` (requires \pkg{pbkrtest}; errors otherwise): the KR
#'   (1997) df with the KR-adjusted covariance `pbkrtest::vcovAdj()` for the
#'   standard error.
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
                          ci_method = c("satterthwaite", "kenward-roger",
                                        "containment", "capped-kr"),
                          pi_method = c("francq-satterthwaite", "containment",
                                        "capped", "ratio-guard"),
                          pi_df_ratio = 1.5) {
  if (!inherits(model, "lmerMod")) {
    .abort("`model` must be a fitted `lmerMod` (an lmer fit).",
           class = "reintervals_error_input")
  }
  which <- match.arg(which, c("CI", "PI"), several.ok = TRUE)
  ci_method <- match.arg(ci_method)
  pi_method <- match.arg(pi_method)
  fx <- .fixed(model, newdata)
  ## theta (observed information) is needed for the PI df and for the default
  ## Satterthwaite CI df; a Kenward-Roger or overridden CI-only call skips it.
  need_theta <- "PI" %in% which ||
    ("CI" %in% which && is.null(ci_df) && ci_method == "satterthwaite")
  theta <- if (need_theta) .theta_cov(model, L = fx$L) else NULL
  comps <- if (!is.null(theta)) theta$components else .var_components(model)
  nr <- nrow(fx$L)
  base <- fx$newdata

  ## one container per prediction row, reusing the model-level components.
  ## `var_mean` overrides the unadjusted l C11 l' (Kenward-Roger CI only).
  make_comp <- function(i, ci_ddf, var_mean = fx$var_fix[i]) {
    dfs <- if (!is.null(theta)) c(ci = ci_ddf, pi = theta$df_pi) else c(ci = ci_ddf)
    re_components(components = comps, dfs = dfs, mean = fx$fit[i],
                  var_mean = var_mean, target = "observable")
  }

  ## the containment df (the rank rule), reported on every CI / PI row as
  ## df_C; a diagnostic only on the arms that do not use it
  df_C <- tryCatch(.ddf_containment(model), error = function(e) NA_real_)
  blocks <- list()
  if ("CI" %in% which) {
    ## ci_df = "capped-kr" (step 13b): the Kenward-Roger SE and df, with the df
    ## capped at the containment df, df = min(df_KR, df_C)
    capped <- is.null(ci_df) && ci_method == "capped-kr"
    kr_like <- is.null(ci_df) && ci_method %in% c("kenward-roger", "capped-kr")
    ## Kenward-Roger: one vcovAdj() per call, shared by every newdata row
    kr <- if (kr_like) .kr_parts(model) else NULL
    dd <- lapply(seq_len(nr), function(i)
      .ci_ddf(model, fx$L[i, ], theta = theta, ci_df = ci_df,
              df_method = if (capped) "kenward-roger" else ci_method, kr = kr))
    ddf <- vapply(dd, function(x) x$df, numeric(1))
    df_KR <- if (kr_like) ddf else rep(NA_real_, nr)
    cap_binds <- rep(NA, nr)
    if (capped) {
      if (is.na(df_C)) df_C <- .ddf_containment(model)   # its own error, if any
      cap_binds <- df_C < df_KR
      ddf <- pmin(df_KR, df_C)
      for (i in which(cap_binds)) dd[[i]]$method <- "containment"
    }
    ## per-row reason a default Satterthwaite CI df fell back (NA otherwise):
    ## the .theta_cov reason when there is no covariance, else .ci_ddf's
    ci_fb <- is.null(ci_df) && ci_method == "satterthwaite"
    ci_reason <- vapply(dd, function(x) {
      if (!ci_fb || !identical(x$method, "containment")) return(NA_character_)
      if (is.null(theta$cov_theta)) theta$fallback_reason
      else if (is.null(x$reason)) NA_character_ else x$reason
    }, character(1))
    cis <- lapply(seq_len(nr), function(i) {
      vm <- if (is.null(dd[[i]]$var_mean)) fx$var_fix[i] else dd[[i]]$var_mean
      ci_mean(make_comp(i, ddf[i], var_mean = vm), conf = conf)
    })
    blocks$CI <- cbind(base, data.frame(
      type = "CI",
      estimate = vapply(cis, function(x) x$estimate, numeric(1)),
      lower = vapply(cis, function(x) x$lower, numeric(1)),
      upper = vapply(cis, function(x) x$upper, numeric(1)),
      df = ddf, df_method = vapply(dd, function(x) x$method, character(1)),
      df_fallback_reason = ci_reason,
      level = conf, df_KR = df_KR, df_C = df_C, ci_cap_binds = cap_binds,
      stringsAsFactors = FALSE
    ), row.names = NULL)
  }
  if ("PI" %in% which) {
    ## ci df is irrelevant to the PI; the container carries a placeholder
    ## pi_df = "containment": the same variance (V_F + sigma2_T), with the
    ## containment df (.ddf_containment, the rank rule) in place of the Francq
    ## generalized-Satterthwaite df and its fallback.
    ## pi_df = "capped" / "ratio-guard" (step 13b): the Francq centre and
    ## variance with df = min(df_F, df_C), always ("capped") or when
    ## df_F / df_C >= pi_df_ratio ("ratio-guard"); df_F is the Francq df after
    ## its own fallback (df_F = df_C on a fallback, so the guard never fires)
    df_F <- theta$df_pi
    if (pi_method %in% c("capped", "ratio-guard") && is.na(df_C))
      df_C <- .ddf_containment(model)                     # its own error, if any
    guard <- switch(pi_method, capped = TRUE,
                    `ratio-guard` = df_F / df_C >= pi_df_ratio, NA)
    if (pi_method == "containment") theta$df_pi <- .ddf_containment(model)
    if (isTRUE(guard)) theta$df_pi <- min(df_F, df_C)
    ## the PI counterpart of ci_cap_binds: df_C < df_F, on the two guarded arms
    cap_binds <- if (pi_method %in% c("capped", "ratio-guard")) df_C < df_F else NA
    pis <- lapply(seq_len(nr), function(i) pi_newobs(make_comp(i, Inf), level = level))
    blocks$PI <- cbind(base, data.frame(
      type = "PI",
      estimate = vapply(pis, function(x) x$estimate, numeric(1)),
      lower = vapply(pis, function(x) x$lower, numeric(1)),
      upper = vapply(pis, function(x) x$upper, numeric(1)),
      df = theta$df_pi, df_method = NA_character_,
      ## the Francq df fell back (default pi_df only): .theta_cov's reason
      df_fallback_reason = if (pi_method != "containment" &&
                                 isTRUE(theta$pi_df_fallback))
        theta$fallback_reason else NA_character_,
      level = level, df_F = df_F, df_C = df_C, pi_guard_fired = guard,
      pi_cap_binds = cap_binds, stringsAsFactors = FALSE
    ), row.names = NULL)
  }

  ## the CI and PI blocks carry their own diagnostics columns: fill the other
  ## block's with NA before binding (a which = c("CI", "PI") call)
  if (length(blocks) > 1L) {
    cols <- unique(unlist(lapply(blocks, names)))
    blocks <- lapply(blocks, function(b) {
      for (cn in setdiff(cols, names(b))) b[[cn]] <- NA
      b[cols]
    })
  }
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  if ("PI" %in% which) {
    attr(out, "components") <- theta$components
    ## the singular-fit df guard is a property of the default (Francq) df only
    attr(out, "pi_df_fallback") <- if (pi_method == "containment") NA else
      isTRUE(theta$pi_df_fallback)
    if (isTRUE(attr(out, "pi_df_fallback")))
      attr(out, "df_fallback_reason") <- theta$fallback_reason
  }
  ## reason for a default Satterthwaite CI df that fell back to containment
  ## (the union of the per-row column df_fallback_reason)
  if ("CI" %in% which && is.null(ci_df) && ci_method == "satterthwaite" &&
        any(blocks$CI$df_method == "containment")) {
    attr(out, "df_fallback_reason") <- c(attr(out, "df_fallback_reason"),
      unique(stats::na.omit(ci_reason)))
  }
  ## set whether or not theta was computed (a CI-only call skips .theta_cov), so
  ## the CI reports lme4::isSingular() like the PI / new-group-mean / TI.
  attr(out, "singular") <- .is_singular(model)
  out
}
