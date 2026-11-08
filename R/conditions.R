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
#   +- reintervals_warning_singular     singular (boundary) fit
#   +- reintervals_warning_limitation   a documented known limitation applies; the
#   |                                   interval is still returned unchanged
#   +- reintervals_warning_deprecated   reserved: deprecated method names
#   +- reintervals_warning_mc           reserved: Monte-Carlo error above threshold
#   +- reintervals_warning_unvalidated  a non-default validation arm (lmm_interval_control)
#                                       was applied; its coverage is under evaluation
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
  reintervals_warning_unvalidated = "reintervals_warning"
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
