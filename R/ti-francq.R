#' Tolerance interval for a linear mixed model
#'
#' The two-sided beta-content / gamma-confidence tolerance interval of Francq,
#' Lin & Hoyer (2019), Equation 26, via the modified large-sample (MLS) method on
#' the expected-mean-square (EMS) linear combination:
#' \deqn{l\hat\beta \pm z_{(1+\beta)/2}\, \sqrt{l \hat C_{11} l' + \hat\sigma^2_T}\,
#'        \sqrt{1 + \frac{1}{\hat\sigma^2_T}\sqrt{\textstyle\sum_j H_j^2 k_j^2 EMS_j^2}}, }
#' with `sigma^2_T = sum(k_j * EMS_j)` and the chi-square tail factors
#' `H_j = r_j / qchisq(1 - gamma, r_j) - 1` on the ANOVA degrees of freedom
#' `r_j`. The EMS terms (`EMS_j`, `k_j`, `r_j`) are supplied design-by-design in
#' `comp$ems` by the extraction layer; this constructor is method-pure and never
#' touches `lme4`.
#'
#' Unlike the prediction interval, the tolerance interval is design-restricted:
#' it is defined only for the designs with a closed-form EMS decomposition
#' (one random factor; two nested or two crossed random factors). When
#' `comp$ems` is `NULL` (off-catalog design, or a CI/PI-only container) the
#' tolerance interval is unavailable and this errors.
#'
#' Only the **observable** target is supported in v1. An `re_components` carrying
#' `target = "true_value"` is rejected with an error rather than silently
#' returning the observable interval: the between-level-only tolerance interval
#' needs the difference-of-EMS (sign-unrestricted MLS) variant, which has no
#' worked oracle in the paper and is deferred to a focused follow-up. The
#' prediction interval (`pi_newobs()`) does support `true_value`.
#'
#' @param comp An `re_components` object with a non-`NULL` `ems` slot. Uses
#'   `comp$mean`, `comp$var_mean`, and `comp$ems`.
#' @param content Content fraction `beta` (proportion of the population covered).
#'   Default 0.95.
#' @param conf Confidence level `gamma`. Default 0.95.
#' @param sides Character: `"two"` (the default and only option in v1). One-sided
#'   tolerance limits are not implemented.
#'
#' @return An [re_interval] of `type = "TI"`.
#'
#' @references Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and
#'   tolerance in linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#'
#' @examples
#' # Francq, Lin & Hoyer (2019) Section 5.2 balanced one-way EMS decomposition:
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
#' ti_francq(comp, content = 0.90, conf = 0.95)
#' @noRd
ti_francq <- function(comp, content = 0.95, conf = 0.95,
                      sides = c("two", "lower", "upper")) {
  if (!inherits(comp, "re_components")) {
    stop("`comp` must be an `re_components` object (see ?re_components).",
         call. = FALSE)
  }
  sides <- match.arg(sides)
  if (sides != "two") {
    stop("Only two-sided tolerance intervals are implemented in v1 ",
         "(`sides = \"two\"`).", call. = FALSE)
  }
  for (nm in c("content", "conf")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0 || v >= 1) {
      stop(sprintf("`%s` must be a single number strictly between 0 and 1.", nm),
           call. = FALSE)
    }
  }
  if (identical(comp$target, "true_value")) {
    stop("The true_value (between-level only) tolerance interval is not ",
         "supported for the fixed-slope design; only `target = \"observable\"` is implemented. ",
         "(`pi_newobs()` does support true_value.) The true_value TI requires ",
         "the difference-of-EMS / sign-unrestricted MLS variant and is deferred ",
         "to a focused follow-up with its own derivation review and validation.",
         call. = FALSE)
  }
  if (is.null(comp$ems)) {
    stop("`comp$ems` is NULL: no expected-mean-square decomposition is ",
         "available for this design/target, so the Francq tolerance interval ",
         "(Eq. 26) cannot be computed. CI and PI remain available.",
         call. = FALSE)
  }

  ems <- comp$ems
  sigma2_t <- sum(ems$k * ems$ms)
  if (sigma2_t <= 0) {
    stop("Total variance `sum(k * EMS)` is not positive; the tolerance ",
         "interval is undefined (all variance components at zero).",
         call. = FALSE)
  }

  a <- 1 - conf                                   # chi-square lower-tail prob
  h <- ems$df / stats::qchisq(a, df = ems$df) - 1
  inner <- sum(h^2 * ems$k^2 * ems$ms^2)
  z <- stats::qnorm((1 + content) / 2)
  hw <- z * sqrt(comp$var_mean + sigma2_t) * sqrt(1 + sqrt(inner) / sigma2_t)

  new_re_interval(
    estimate = comp$mean,
    lower    = comp$mean - hw,
    upper    = comp$mean + hw,
    type     = "TI",
    content  = content,
    conf     = conf,
    df       = NA_real_,                          # normal-based; no t df
    sides    = "two",
    target   = comp$target,
    design   = comp$design,
    call     = match.call()
  )
}
