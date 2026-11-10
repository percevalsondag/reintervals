# =============================================================================
# Classed conditions. Every error the package raises inherits from
# `reintervals_error`, and every warning from `reintervals_warning`, so callers
# can catch or muffle them by class (tryCatch / withCallingHandlers /
# testthat::expect_error(class = )). Base R only (errorCondition /
# warningCondition); no rlang.
#
# Hierarchy (each class also carries all of its ancestors):
#
#   reintervals_error
#   +- reintervals_error_input          bad arguments, newdata, alpha/level/P, control
#   +- reintervals_error_unsupported    model / target / interval combination not available
#   |  +- reintervals_error_design      no closed-form decomposition for this design
#   |  +- reintervals_error_method      forced `method` incompatible with the model/interval
#   +- reintervals_error_data           the data cannot support the interval
#   |  +- reintervals_error_groups      fewer than 2 groups / lots / batches
#   |  +- reintervals_error_time_points a lot / batch without >= 2 distinct time points
#   |  +- reintervals_error_residual_df non-positive residual df
#   +- reintervals_error_numerical      a quantity is undefined (W <= 0, sum(k EMS) <= 0,
#   |                                   singular synthesis, too many NA Monte-Carlo draws)
#   +- reintervals_error_dependency     a Suggests package (pbkrtest) is missing
#   +- reintervals_error_internal       an internal invariant failed (a bug)
#   +- reintervals_error_defunct        reserved: defunct method names
#
#   reintervals_warning
#   +- reintervals_warning_fallback     a rule fell back (df, variance)
#   |  +- reintervals_warning_truncation  a negative estimate was truncated / floored at 0
#   |  +- reintervals_warning_df_fallback a DEFAULT df rule fell back to containment
#   +- reintervals_warning_singular     singular (boundary) fit
#   +- reintervals_warning_limitation   a documented known limitation applies; the
#   |                                   interval is still returned unchanged
#   +- reintervals_warning_deprecated   reserved: deprecated method names
#   +- reintervals_warning_mc           reserved: Monte-Carlo error above threshold
#   +- reintervals_warning_unvalidated  a non-default validation arm (lmm_interval_control)
#   |                                   was applied; its coverage is under evaluation
#   +- reintervals_warning_extrapolation  a slope-model t0 lies outside the observed
#                                       time range (the interval extrapolates the model)
# =============================================================================

.reintervals_condition_parents <- c(
  reintervals_error              = NA,
  reintervals_error_input        = "reintervals_error",
  reintervals_error_unsupported  = "reintervals_error",
  reintervals_error_design       = "reintervals_error_unsupported",
  reintervals_error_method       = "reintervals_error_unsupported",
  reintervals_error_data         = "reintervals_error",
  reintervals_error_groups       = "reintervals_error_data",
  reintervals_error_time_points  = "reintervals_error_data",
  reintervals_error_residual_df  = "reintervals_error_data",
  reintervals_error_numerical    = "reintervals_error",
  reintervals_error_dependency   = "reintervals_error",
  reintervals_error_internal     = "reintervals_error",
  reintervals_error_defunct      = "reintervals_error",
  reintervals_warning            = NA,
  reintervals_warning_fallback   = "reintervals_warning",
  reintervals_warning_truncation = "reintervals_warning_fallback",
  reintervals_warning_singular   = "reintervals_warning",
  reintervals_warning_limitation = "reintervals_warning",
  reintervals_warning_deprecated = "reintervals_warning",
  reintervals_warning_mc         = "reintervals_warning",
  reintervals_warning_unvalidated = "reintervals_warning",
  reintervals_warning_extrapolation = "reintervals_warning",
  reintervals_warning_df_fallback = "reintervals_warning_fallback"
)

## Package-level state for once-per-session warnings (flags keyed by name).
## Reset with .reset_once_warnings() (used by the tests).
.reintervals_state <- new.env(parent = emptyenv())

## TRUE the first time it is called for `key` in this R session, FALSE after.
.once <- function(key) {
  flag <- paste0("warned_", key)
  if (isTRUE(get0(flag, envir = .reintervals_state, inherits = FALSE))) return(FALSE)
  assign(flag, TRUE, envir = .reintervals_state)
  TRUE
}

## Clear every once-per-session flag (internal; for tests).
.reset_once_warnings <- function() {
  rm(list = ls(.reintervals_state, all.names = TRUE), envir = .reintervals_state)
  invisible(NULL)
}

## Leaf class -> c(leaf, parent, ..., root). Errors on an unknown class so a
## typo at a call site fails loudly in the tests.
.condition_chain <- function(class) {
  if (!class %in% names(.reintervals_condition_parents)) {
    stop(errorCondition(sprintf("Internal: unknown condition class '%s'.", class),
                        class = c("reintervals_error_internal", "reintervals_error"),
                        call = NULL))
  }
  out <- character(0)
  while (!is.na(class)) {
    out <- c(out, class)
    class <- .reintervals_condition_parents[[class]]
  }
  out
}

## Signal a classed error. `...` is pasted like stop()'s message arguments;
## `data` adds fields to the condition object. No call is recorded (the
## equivalent of stop(call. = FALSE)).
.abort <- function(..., class = "reintervals_error", data = list()) {
  cnd <- do.call(errorCondition,
                 c(list(.makeMessage(...), class = .condition_chain(class),
                        call = NULL), data))
  stop(cnd)
}

## Signal a classed warning (warning(call. = FALSE) equivalent).
.warn <- function(..., class = "reintervals_warning", data = list()) {
  cnd <- do.call(warningCondition,
                 c(list(.makeMessage(...), class = .condition_chain(class),
                        call = NULL), data))
  warning(cnd)
}

## ---- slope-model newdata guard and extrapolation warning ----------------------
## The time variable of a slope model: list(var, random) for a random-slope model
## (the one non-intercept fixed effect) or a fixed-slope model (.m1_time_name),
## NULL when the model has no time covariate (random-intercept designs; also a
## random-slope fit with several fixed effects, which the engines reject).
.slope_time_info <- function(model) {
  if (.has_random_slope(model)) {
    tv <- setdiff(names(lme4::fixef(model)), "(Intercept)")
    if (length(tv) != 1L) return(NULL)
    return(list(var = tv, random = TRUE))
  }
  tv <- .m1_time_name(model)
  if (is.na(tv)) return(NULL)
  list(var = tv, random = FALSE)
}

## Slope-model `newdata`: it must carry the time variable as a numeric, finite,
## non-missing column (extra columns are ignored). Returns the requested t0
## values (NULL when newdata is NULL: the engines then use the first
## model-frame row, which is inside the observed range).
.check_time_newdata <- function(model, newdata, info) {
  if (is.null(newdata)) return(NULL)
  tv <- info$var
  if (!is.data.frame(newdata) && !is.list(newdata)) .check_newdata(model, newdata)
  if (!tv %in% names(newdata)) {
    if (!info$random) .check_newdata(model, newdata)   # the existing message
    .abort("Could not determine the evaluation time point(s) for this interval; ",
           "supply `newdata` with a numeric '", tv,
           "' column giving the time(s) to evaluate at.", class = "reintervals_error_input")
  }
  v <- newdata[[tv]]
  bad <- if (!is.numeric(v) || NCOL(v) != 1L) {
    sprintf("is of class '%s', not numeric", class(v)[1L])
  } else if (!length(v)) {
    "has no rows"
  } else if (anyNA(v)) {
    "contains missing values (NA)"
  } else if (any(!is.finite(v))) {
    "contains non-finite values (Inf)"
  }
  if (!is.null(bad)) {
    .abort("The time variable '", tv, "' in `newdata` ", bad, ". `newdata$", tv,
           "` must be a numeric column of finite, non-missing time point(s).",
           class = "reintervals_error_input", data = list(variable = tv))
  }
  as.numeric(v)
}

## Warn ONCE (per call) when any requested t0 lies outside the observed time
## range [min, max] of the fitted time variable.
.warn_extrapolation <- function(t0, time, random, var = "time") {
  if (is.null(t0) || is.null(time)) return(invisible(NULL))
  t0 <- suppressWarnings(as.numeric(t0)); time <- as.numeric(time)
  rng <- range(time[is.finite(time)])
  out <- t0[is.finite(t0) & (t0 < rng[1] | t0 > rng[2])]
  if (!length(out) || any(!is.finite(rng))) return(invisible(NULL))
  model_txt <- if (random) {
    paste0("the random-slope model, whose between-group variance V_G(t0) is ",
           "quadratic in t0, holds beyond the data")
  } else {
    "the linear (fixed-slope) model holds beyond the data"
  }
  out_u <- unique(out)
  .warn(sprintf("Extrapolation: requested %s = %s %s outside the observed time range [%s, %s]. ",
                var, paste(format(out_u, trim = TRUE), collapse = ", "),
                if (length(out_u) == 1L) "is" else "are",
                format(rng[1]), format(rng[2])),
        "The interval there assumes that ", model_txt, ".",
        class = "reintervals_warning_extrapolation",
        data = list(t0 = out_u, range = rng, variable = var))
}

## Model-based guard used by every verb: validate newdata's time column and warn
## on extrapolation. A no-op for models without a time covariate.
.slope_newdata_guard <- function(model, newdata) {
  info <- .slope_time_info(model)
  if (is.null(info)) return(invisible(NULL))
  t0 <- .check_time_newdata(model, newdata, info)
  .warn_extrapolation(t0, stats::model.frame(model)[[info$var]],
                      random = info$random, var = info$var)
  invisible(t0)
}

## ---- default-df fallback warning ---------------------------------------------
## Warn ONCE (per call) when a DEFAULT df rule fell back to containment on any
## row of a coverage-interval result `out` (from .interval_lmm). Validation arms
## never fall back, and a user-chosen containment arm is not a fallback, so both
## are skipped. Reads the per-engine fallback markers:
##   random slope (CI / PI / NGM, default arms): df_rule == "containment"
##     (vfun_extract's Hessian-Satterthwaite -> .vfun_df_containment);
##   random-intercept CI (ci_df = "satterthwaite"): df_method == "containment";
##   random-intercept PI (pi_df = "francq-satterthwaite"): attr pi_df_fallback;
##   random-intercept NGM (ngm_method = "package", crossed Gv):
##     df_rule == "containment".
.warn_df_fallback <- function(out, model, kind, applied) {
  verb <- switch(kind, CI = "ci_lmm()", PI = "pi_lmm()",
                 new_group_mean = "new_group_mean_lmm()")
  arm <- function(opt, default) {
    if (opt %in% names(applied)) applied[[opt]] else default
  }
  rows <- logical(0); rule <- NULL
  if (.has_random_slope(model)) {
    default <- switch(kind,
      CI = arm("slope_ci_df", "hessian-satterthwaite") == "hessian-satterthwaite",
      PI = TRUE,
      new_group_mean = arm("slope_ngm_method", "reml") == "reml")
    if (default && "df_rule" %in% names(out)) rows <- out$df_rule %in% "containment"
    rule <- "Hessian-Satterthwaite"
  } else if (kind == "CI") {
    if (arm("ci_df", "satterthwaite") == "satterthwaite" && "df_method" %in% names(out))
      rows <- out$df_method %in% "containment"
    rule <- "Satterthwaite"
  } else if (kind == "PI") {
    if (arm("pi_df", "francq-satterthwaite") == "francq-satterthwaite" &&
          isTRUE(attr(out, "pi_df_fallback"))) rows <- rep(TRUE, nrow(out))
    rule <- "Francq generalized-Satterthwaite"
  } else {
    if (arm("ngm_method", "package") == "package" && "df_rule" %in% names(out))
      rows <- out$df_rule %in% "containment"
    rule <- "Satterthwaite (crossed Gv)"
  }
  if (!any(rows)) return(invisible(NULL))
  dfs <- unique(out$df[rows])
  reason <- unique(stats::na.omit(attr(out, "df_fallback_reason")))
  .warn(verb, ": the default ", rule, " df could not be used",
        if (length(reason)) paste0(" (", paste(reason, collapse = "; "), ")"),
        "; the containment df was used instead (rule: containment, df = ",
        paste(format(dfs, digits = 6, trim = TRUE), collapse = ", "), ").",
        class = "reintervals_warning_df_fallback",
        data = list(verb = verb, rule = "containment", from = rule, df = dfs,
                    reason = if (length(reason)) reason else NA_character_))
}
