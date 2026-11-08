# Muffle ONLY the documented known-limitation warning (class
# reintervals_warning_limitation) at call sites where it is not what the test is
# about, e.g. pi_lmm() on a design whose top random factor has <= 5 levels. Any
# other warning still surfaces. That the warning is raised is tested in
# test-conditions.R.
quiet_limitation <- function(expr) {
  withCallingHandlers(expr, reintervals_warning_limitation = function(w)
    invokeRestart("muffleWarning"))
}
