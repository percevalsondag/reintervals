#' Prediction interval for a future observation
#'
#' The two-sided prediction interval of Francq, Lin & Hoyer (2019),
#' Equation 21:
#' \deqn{l\hat\beta \pm t_{1 - \psi/2,\, \nu_{PI}}\,
#'        \sqrt{l \hat C_{11} l' + \hat\sigma^2_T},}
#' where the total variance of the predicted future value is
#' `sigma^2_T = sum(coefs * components)` -- all variance components for an
#' observable future value, or the between-level components only for an
#' unobservable level "true value" (selected by `comp$coefs`/`comp$target`).
#' The degrees of freedom `nu_PI` come from the observed Fisher information
#' (Eq. 23-24) and are read from `comp$dfs["pi"]`; computing them is the
#' extraction layer's job, so this constructor stays method-pure and never
#' touches `lme4`. Unlike the tolerance interval, the prediction interval is
#' **general** across random-intercept designs (balanced or unbalanced).
#'
#' @param comp An `re_components` object. Uses `comp$mean`, `comp$var_mean`,
#'   `comp$components`, `comp$coefs`, and `comp$dfs["pi"]`.
#' @param level Prediction level `1 - psi`. Default 0.95.
#'
#' @return An [re_interval] of `type = "PI"`.
#'
#' @references Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and
#'   tolerance in linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#'
#' @examples
#' # Francq, Lin & Hoyer (2019) Section 5.2 balanced one-way:
#' # sigma^2_T = 0.001934, df_PI = 12.484 -> 95% PI [0.881, 1.081].
#' comp <- re_components(
#'   components = c(run = 0.000681, Residual = 0.001253),
#'   dfs        = c(ci = 5, pi = 12.484),
#'   mean       = 0.981,
#'   var_mean   = 0.01353^2,
#'   design     = "oneway"
#' )
#' pi_newobs(comp, level = 0.95)
#' @noRd
pi_newobs <- function(comp, level = 0.95) {
  if (!inherits(comp, "re_components")) {
    stop("`comp` must be an `re_components` object (see ?re_components).",
         call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!"pi" %in% names(comp$dfs)) {
    stop("`comp$dfs` must contain a \"pi\" entry (the prediction-interval df ",
         "from the observed Fisher information). The extraction layer supplies ",
         "it.", call. = FALSE)
  }

  sigma2_t <- sum(comp$coefs * comp$components)
  df <- unname(comp$dfs[["pi"]])
  hw <- stats::qt(1 - (1 - level) / 2, df = df) * sqrt(comp$var_mean + sigma2_t)

  new_re_interval(
    estimate = comp$mean,
    lower    = comp$mean - hw,
    upper    = comp$mean + hw,
    type     = "PI",
    level    = level,
    df       = df,
    sides    = "two",
    target   = comp$target,
    design   = comp$design,
    call     = match.call()
  )
}
