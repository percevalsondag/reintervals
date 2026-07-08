# =============================================================================
# Public surface: the four interval verbs (ci_lmm / pi_lmm / ti_lmm /
# new_group_mean_lmm), the predict-style wrapper (lmm_predict), and the control
# constructor (lmm_interval_control). These wrappers translate the uniform
# argument contract into the engine calls and stamp the harmonized
# provenance (method / design / grouping) onto the result; the numbers come from
# the engines unchanged.
# =============================================================================

#' Tuning controls for interval Monte-Carlo (GPQ)
#'
#' Tuning knobs only --- the RNG seed and Monte-Carlo size used when the
#' GPQ algorithm runs. Algorithm *selection* is the `method` argument, not this.
#'
#' @param seed Optional integer seed for the GPQ draws (reproducibility / GxP
#'   traceability). `NULL` uses the ambient RNG state.
#' @param M Number of GPQ Monte-Carlo realizations. Default 10000.
#' @return An `lmm_interval_control` list (passed to a verb's `control`
#'   argument); a list with elements `seed` and `M`.
#' @seealso [ci_lmm()], [ti_lmm()], [ti_gpq_raw()].
#' @examples
#' lmm_interval_control(seed = 1, M = 5000)
#' @export
lmm_interval_control <- function(seed = NULL, M = 10000L) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  }
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || M < 1) {
    stop("`M` must be a single positive integer.", call. = FALSE)
  }
  structure(list(seed = seed, M = as.integer(M)), class = "lmm_interval_control")
}

## alpha / level as a mutually-exclusive pair (the rgamma rate/scale pattern).
## Returns alpha. both -> error; neither -> 0.05; level -> 1 - level.
.resolve_alpha <- function(alpha, level) {
  if (!is.null(alpha) && !is.null(level)) {
    stop("Supply only one of `alpha` or `level` (they are the same quantity).",
         call. = FALSE)
  }
  if (is.null(alpha) && is.null(level)) return(0.05)
  if (!is.null(level)) {
    if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
          level <= 0 || level >= 1) {
      stop("`level` must be a single number in (0, 1).", call. = FALSE)
    }
    return(1 - level)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  alpha
}

## defensive front-door class check (clear error at the door).
.assert_lmm <- function(model) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted lme4 model (an `lmerMod` from lmer()).",
         call. = FALSE)
  }
  invisible(TRUE)
}

## route-independent model-class tag (the `design` column vocabulary). Keyed on
## the random/fixed STRUCTURE, not the engine that computed the interval, so a
## given model gets the same design tag whichever route ran (kills oneway-vs-fixed-slope).
.model_class <- function(model) {
  cnms <- lme4::getME(model, "cnms")
  has_slope <- any(vapply(cnms, function(t) any(t != "(Intercept)"), logical(1)))
  if (has_slope) {
    return(if (length(cnms) == 1L) "random_slope_correlated"
           else "random_slope_independent")
  }
  fe <- lme4::fixef(model)
  if (length(setdiff(names(fe), "(Intercept)")) == 1L) {
    return("random_intercept_fixed_slope")
  }
  "random_intercept"
}

## the random-INTERCEPT topology (oneway/nested/crossed), moved off `design`
## into its own attribute. Defined only for models with NO random slope; NA
## otherwise.
.grouping_of <- function(model) {
  cnms <- lme4::getME(model, "cnms")
  has_slope <- any(vapply(cnms, function(t) any(t != "(Intercept)"), logical(1)))
  if (has_slope) return(NA_character_)
  tryCatch(.design_of(model)$type, error = function(e) NA_character_)
}

## interval-kind tag for the `type` column (always lower-case / snake).
.type_tag <- function(kind) {
  switch(kind, CI = , ci = "ci", PI = , pi = "pi", TI = , ti = "ti",
         new_group_mean = "new_group_mean", tolower(kind))
}

## Value columns an engine may emit; everything else in the native frame is an
## eval-point (predictor) column. "time" is the slope/raw engines' generic
## placeholder -- renamed to the actual covariate when a model is in hand.
.value_cols <- c("time", "type", "estimate", "lower", "upper", "level", "conf",
                 "content", "P", "sides", "df", "se", "method", "design",
                 "singular", "boundary", "vg_floored", "tau2_med", "s2e_med",
                 "frac_boundary", "K", "df_resid", "phi1", "c_resid", "n_E")

## The output unifier. Maps any engine's native frame to the single tidy
## contract; relabels only, never recomputes. Universal columns, in order:
##   <eval cols...>, type, estimate, lower, upper, level, P, sides, method,
##   design, singular, boundary, vg_floored
## Engine-specific diagnostics (df, se, beta_hat, GPQ audit, ...) go to
## attr(.,"diagnostics"). `model = NULL` serves the raw engines (design taken
## from the native row, grouping NA, the generic "time" eval column kept).
.unify_intervals <- function(native, model, kind, P_value = NA_real_,
                             method_tag = NULL) {
  n <- nrow(native)
  ## eval-point columns (predictors); rename the generic "time" to the covariate
  eval_df <- native[, setdiff(names(native), .value_cols), drop = FALSE]
  if ("time" %in% names(native)) {
    tv <- if (!is.null(model)) {
      setdiff(names(lme4::fixef(model)), "(Intercept)")[1]
    } else "time"
    eval_df[[tv]] <- native$time
  }

  level_out  <- if (kind %in% c("ti", "TI")) native$conf else native$level
  method_out <- if ("method" %in% names(native)) native$method else method_tag
  design_out <- if (!is.null(model)) .model_class(model)
                else if ("design" %in% names(native)) native$design else NA_character_
  pick <- function(col) if (col %in% names(native)) native[[col]] else rep(NA, n)
  singular_out <- if ("singular" %in% names(native)) native$singular
                  else rep(isTRUE(attr(native, "singular")), n)

  core <- data.frame(
    type = .type_tag(kind), estimate = native$estimate,
    lower = native$lower, upper = native$upper,
    level = level_out, P = rep(if (kind %in% c("ti", "TI")) P_value else NA_real_, n),
    sides = if ("sides" %in% names(native)) native$sides else rep("two.sided", n),
    method = method_out, design = design_out,
    singular = singular_out, boundary = pick("boundary"),
    vg_floored = pick("vg_floored"), stringsAsFactors = FALSE
  )
  out <- cbind(eval_df, core); rownames(out) <- NULL

  ## diagnostics -> attribute (off the core shape, still retrievable)
  diag <- list()
  for (cn in c("df", "se", "tau2_med", "s2e_med", "frac_boundary", "K",
               "df_resid", "phi1", "c_resid")) {
    if (cn %in% names(native)) diag[[cn]] <- native[[cn]]
  }
  for (a in c("beta_hat", "n_E", "U", "Vhat_Y")) {
    v <- attr(native, a); if (!is.null(v)) diag[[a]] <- v
  }
  attr(out, "diagnostics") <- diag
  attr(out, "grouping") <- if (!is.null(model)) .grouping_of(model) else NA_character_
  class(out) <- c("lmm_interval", "data.frame")
  out
}

#' @export
print.lmm_interval <- function(x, ...) {
  kind <- if (length(unique(x$type)) == 1L) unique(x$type) else "mixed"
  cat(sprintf("<lmm_interval> %s  (method = %s, design = %s)\n",
              kind, paste(unique(x$method), collapse = "/"),
              paste(unique(x$design), collapse = "/")))
  core <- as.data.frame(x)[, !names(x) %in% c("design", "method"), drop = FALSE]
  print(core, row.names = FALSE, ...)
  if (any(stats::na.omit(x$singular)) || any(stats::na.omit(x$boundary)) ||
        any(stats::na.omit(x$vg_floored))) {
    cat("note: some rows carry a singular/boundary/vg_floored flag ",
        "(degenerate fit at that point); see the flag columns.\n", sep = "")
  }
  d <- attr(x, "diagnostics")
  if (length(d)) cat(sprintf("diagnostics: %s  (see attr(., \"diagnostics\"))\n",
                             paste(names(d), collapse = ", ")))
  invisible(x)
}

#' Confidence interval for the mean of a linear mixed model
#'
#' One of the four interval verbs ([ci_lmm()], [pi_lmm()], [ti_lmm()],
#' [new_group_mean_lmm()]). `ci_lmm()` gives a confidence interval for the mean
#' response at `newdata`. It dispatches by random structure: random-intercept
#' designs (including a fixed slope) use the expected-mean-square engine;
#' random-slope designs use the REML closed form. See [reintervals-models] for
#' the model-class to engine map and [lmm_predict()] for the unified entry point.
#'
#' @param model A fitted `lmerMod` (from [lme4::lmer()]).
#' @param newdata Row(s) giving the fixed-effect combination(s) to evaluate at.
#'   `NULL` uses the first model-frame row.
#' @param alpha,level The interval confidence, as a mutually-exclusive pair: give
#'   at most one. `level` is `1 - alpha`; the default is `alpha = 0.05` (95%).
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`.
#' @param method `"auto"` (default) chooses the algorithm by structure; one of
#'   `"ems-mls"`, `"reml-mls"`, `"anova-mls"`, `"gpq"` forces that engine and
#'   errors if it is incompatible with the model.
#' @param control An [lmm_interval_control()] object (the GPQ `seed` and `M`).
#'
#' @return An `lmm_interval` object (a data frame), one row per evaluation point,
#'   with columns, in order:
#'   \itemize{
#'     \item the eval-point column(s), named by the actual model predictor(s)
#'       (e.g. the time covariate for a slope model, or each `newdata` column);
#'     \item `type` (`"ci"`, `"pi"`, `"ti"`, or `"new_group_mean"`), `estimate`,
#'       `lower`, `upper`;
#'     \item `level` (the confidence, `1 - alpha`); `P` (content proportion ---
#'       the value for `ti`, `NA` for the coverage intervals); `sides`;
#'     \item `method` (the engine: `"ems-mls"`, `"reml-mls"`, `"anova-mls"`, or
#'       `"gpq"`); `design` (the route-independent model class, see
#'       [reintervals-models]);
#'     \item the per-row flags `singular`, `boundary`, `vg_floored`.
#'   }
#'   `attr(., "diagnostics")` is a list of engine-specific extras (degrees of
#'   freedom `df`, standard error `se`, and any GPQ / ANOVA intermediates);
#'   `attr(., "grouping")` is the random-intercept topology (`"oneway"`,
#'   `"nested"`, `"crossed"`, or `NA` for slope models).
#'
#' @section Anti-conservative drift on single-observation fixed-slope data:
#' On a fixed-slope model fitted to single-observation / unbalanced lots (the
#' release-only stability structure), the EMS confidence, prediction, and
#' new-group-mean intervals can be anti-conservative (too narrow), because
#' REML precision-weighting down-weights the single-observation lots.
#' [ci_lmm()], [pi_lmm()], and [new_group_mean_lmm()] emit a warning when
#' that structure is detected (the interval is still returned, unchanged). The
#' *tolerance* interval ([ti_lmm()]) is instead routed to the validated Montes
#' ANOVA closed form on this structure; a validated closed form for the
#' confidence / prediction / new-group-mean case is planned for a future version.
#'
#' @seealso [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' ci_lmm(fm, newdata = data.frame(Days = 0))
#' @export
ci_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "CI", level = lvl, sides = sides,
                 method = method, control = control)
}

#' Prediction interval for a future observation of a linear mixed model
#'
#' One of the four interval verbs ([ci_lmm()], [pi_lmm()], [ti_lmm()],
#' [new_group_mean_lmm()]). `pi_lmm()` gives a prediction interval for one future
#' observation at `newdata`. Dispatch and arguments are as for [ci_lmm()].
#'
#' @inheritParams ci_lmm
#' @inherit ci_lmm return
#' @inheritSection ci_lmm Anti-conservative drift on single-observation fixed-slope data
#' @seealso [ci_lmm()], [ti_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' pi_lmm(fm, newdata = data.frame(Days = 0))
#' @export
pi_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "PI", level = lvl, sides = sides,
                 method = method, control = control)
}

#' New-group-mean (CInew / new-batch-mean) interval
#'
#' One of the four interval verbs: a prediction interval for the mean of one new,
#' as-yet-unobserved group/batch --- the CInew / new-batch-mean interval. It
#' is the prediction interval with the residual component removed,
#' \eqn{\mu(t_0) \pm t_{\nu_G}\sqrt{V_F + V_G}}; coverage is `level`/`alpha`, so
#' it carries no content `P` (unlike [ti_lmm()]). Dispatch and arguments are
#' as for [ci_lmm()].
#'
#' @inheritParams ci_lmm
#' @inherit ci_lmm return
#' @inheritSection ci_lmm Anti-conservative drift on single-observation fixed-slope data
#' @seealso [ci_lmm()], [pi_lmm()], [ti_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'                  data = lme4::sleepstudy)
#' new_group_mean_lmm(fm, newdata = data.frame(Days = 5))
#' @export
new_group_mean_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                               sides = c("two.sided", "lower", "upper"),
                               method = "auto",
                               control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "new_group_mean", level = lvl,
                 sides = sides, method = method, control = control)
}

## shared body for the three coverage verbs: dispatch via the engine entry
## (.interval_lmm) then stamp harmonized provenance. method tags: slope ->
## "reml-mls" (already on the row); EMS -> "ems-mls".
.coverage_verb <- function(model, newdata, kind, level, sides, method, control) {
  if (!identical(method, "auto")) {
    ok <- if (.has_random_slope(model)) "reml-mls" else "ems-mls"
    if (!method %in% ok) {
      stop(sprintf(
        "method = \"%s\" is not available for a %s model on this interval; ",
        method, .model_class(model)),
        sprintf("use \"auto\" or \"%s\".", ok), call. = FALSE)
    }
  }
  ## Anti-conservative-drift warning, scoped to EXACTLY the drifting cells:
  ## fixed-slope model with single-observation / unbalanced lots. On that
  ## structure the EMS CI / PI / new-group-mean drift narrow (the same REML
  ## precision-weighting that pushes the EMS TI ~2.4% narrow vs the Montes
  ## oracle). We make the drift VISIBLE -- we do NOT re-route (a Montes/GPQ
  ## CI/PI/new-group-mean is possible new methodology, recorded as v2). Reuses
  ## the same .is_single_obs_m1() detector ti_lmm uses to re-route the TI, so it
  ## cannot fire on clean-balanced M1 (exact), plain random-intercept, or slope
  ## models, and it does not touch ti_lmm (whose TI is already re-routed).
  if (.is_single_obs_m1(model)) {
    warning("This is a fixed-slope model with single-observation / unbalanced ",
            "lots; the EMS ", kind, " interval may be anti-conservative (too ",
            "narrow) on such data, because REML precision-weighting down-weights ",
            "the single-observation lots. The tolerance interval (ti_lmm) is ",
            "routed to the validated Montes ANOVA closed form on this structure; ",
            "no validated closed form exists yet for the confidence / prediction ",
            "/ new-group-mean interval here (planned v2).", call. = FALSE)
  }
  out <- .interval_lmm(model, newdata = newdata, kind = kind, level = level,
                       sides = sides)
  tag <- if (.has_random_slope(model)) NULL else "ems-mls"   # slope row already tagged
  .unify_intervals(out, model, kind = kind, method_tag = tag)
}

#' Unified interface to the four interval verbs
#'
#' A predict-style wrapper over the four interval verbs ([ci_lmm()], [pi_lmm()],
#' [ti_lmm()], [new_group_mean_lmm()]), selected by `interval`. This is a plain
#' function, not an S3 method --- it does not override `lme4`'s
#' `predict.merMod`.
#'
#' @inheritParams ti_lmm
#' @param interval Which interval: `"confidence"`, `"prediction"`, `"tolerance"`,
#'   or `"new_group_mean"`.
#' @param P,over Content proportion and tolerance target; honored only for
#'   `interval = "tolerance"` (forwarded to [ti_lmm()]). Supplying a non-default
#'   `P`/`over` with any other `interval` is an error.
#' @inherit ci_lmm return
#' @seealso [ci_lmm()], [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' lmm_predict(fm, interval = "confidence", newdata = data.frame(Days = 0))
#' @export
lmm_predict <- function(model,
                        interval = c("confidence", "prediction", "tolerance",
                                     "new_group_mean"),
                        newdata = NULL, P = 0.95,
                        over = c("observation", "group_mean"),
                        alpha = NULL, level = NULL,
                        sides = c("two.sided", "lower", "upper"),
                        method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  interval <- match.arg(interval)
  sides    <- match.arg(sides)
  over     <- match.arg(over)
  if (interval != "tolerance" && (!isTRUE(all.equal(P, 0.95)) ||
        over != "observation")) {
    stop("`P` and `over` apply to interval = \"tolerance\" only.", call. = FALSE)
  }
  switch(interval,
    confidence     = ci_lmm(model, newdata, alpha, level, sides, method, control),
    prediction     = pi_lmm(model, newdata, alpha, level, sides, method, control),
    new_group_mean = new_group_mean_lmm(model, newdata, alpha, level, sides,
                                        method, control),
    tolerance      = ti_lmm(model, newdata, P = P, over = over, alpha = alpha,
                            level = level, sides = sides, method = method,
                            control = control)
  )
}
