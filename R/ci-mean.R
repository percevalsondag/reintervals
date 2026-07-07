#' Confidence interval for the mean of a linear mixed model
#'
#' The two-sided confidence interval for the fixed-effect linear combination
#' `lb` of Francq, Lin & Hoyer (2019), Equation 12:
#' \deqn{l\hat\beta \pm t_{1 - \alpha/2,\, \nu}\, \sqrt{l \hat C_{11} l'}.}
#' The denominator degrees of freedom `nu` are supplied by the extraction layer
#' (Kenward-Roger via \pkg{pbkrtest}, or a Satterthwaite fallback) and read from
#' `comp$dfs["ci"]`. This constructor is method-pure: it consumes the numbers in
#' an `re_components` object and never touches `lme4`, so it can be checked
#' directly against a paper's printed values.
#'
#' @param comp An `re_components` object. Uses `comp$mean` (`l\hat\beta`),
#'   `comp$var_mean` (`l \hat C_{11} l'`), and `comp$dfs["ci"]` (the CI
#'   denominator df).
#' @param conf Confidence level `1 - alpha`. Default 0.95.
#'
#' @return An [re_interval] of `type = "CI"`.
#'
#' @references Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and
#'   tolerance in linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#'
#' @examples
#' # Francq, Lin & Hoyer (2019) Section 5.2 balanced one-way:
#' # intercept 0.981 (SE 0.01353), CI df = 5 -> 95% CI [0.946, 1.016].
#' comp <- re_components(
#'   components = c(run = 0.000681, Residual = 0.001253),
#'   dfs        = c(ci = 5, pi = 12.484),
#'   mean       = 0.981,
#'   var_mean   = 0.01353^2,
#'   design     = "oneway"
#' )
#' ci_mean(comp, conf = 0.95)
#' @noRd
ci_mean <- function(comp, conf = 0.95) {
  if (!inherits(comp, "re_components")) {
    stop("`comp` must be an `re_components` object (see ?re_components).",
         call. = FALSE)
  }
  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) ||
        conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!"ci" %in% names(comp$dfs)) {
    stop("`comp$dfs` must contain a \"ci\" entry (the CI denominator df). ",
         "The extraction layer supplies it via Kenward-Roger or a ",
         "Satterthwaite fallback.", call. = FALSE)
  }

  df <- unname(comp$dfs[["ci"]])
  hw <- stats::qt(1 - (1 - conf) / 2, df = df) * sqrt(comp$var_mean)

  new_re_interval(
    estimate = comp$mean,
    lower    = comp$mean - hw,
    upper    = comp$mean + hw,
    type     = "CI",
    conf     = conf,
    df       = df,
    sides    = "two",
    design   = comp$design,
    call     = match.call()
  )
}
