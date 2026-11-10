# Muffle ONLY the documented known-limitation warning (class
# reintervals_warning_limitation) at call sites where it is not what the test is
# about, e.g. pi_lmm() on a design whose top random factor has <= 5 levels. It
# also muffles the two A5 warnings (reintervals_warning_extrapolation,
# reintervals_warning_df_fallback), for call sites that deliberately extrapolate
# or use a singular fit. Any other warning still surfaces. That these warnings
# are raised is tested in test-conditions.R and test-a5-warnings.R.
quiet_limitation <- function(expr) {
  withCallingHandlers(expr,
    reintervals_warning_limitation    = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))
}

# Muffle ONLY the extrapolation warning: a call site that deliberately asks for a
# t0 outside the observed time range (e.g. an oracle evaluated at month 48).
quiet_extrapolation <- function(expr) {
  withCallingHandlers(expr, reintervals_warning_extrapolation = function(w)
    invokeRestart("muffleWarning"))
}

# Muffle ONLY the default-df fallback warning: a call site that deliberately uses
# a singular (or otherwise degenerate) fit whose default df falls back.
quiet_df_fallback <- function(expr) {
  withCallingHandlers(expr, reintervals_warning_df_fallback = function(w)
    invokeRestart("muffleWarning"))
}
