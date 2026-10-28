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
#' @param sides Character: `"two"` (default), `"lower"`, or `"upper"`. One-sided
#'   limits are the package's extension of Hoffman's (2010) one-sided limit to
#'   this construction, validated by simulation.
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
  for (nm in c("content", "conf")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0 || v >= 1) {
      stop(sprintf("`%s` must be a single number strictly between 0 and 1.", nm),
           call. = FALSE)
    }
  }
  if (identical(comp$target, "true_value")) {
    stop("The true_value (between-level only) tolerance interval is not ",
         "supported for any design; only `target = \"observable\"` is implemented. ",
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
  ## One-sided: the package's extension of Hoffman's (2010) one-sided limit to
  ## the Francq construction (validated by simulation):
  ##   mean -/+ ( z_P sqrt(U) + z_{1-alpha} sqrt(V_F) ),
  ## with U = sigma2_T + sqrt(sum H_j^2 k_j^2 EMS_j^2), Francq's MLS upper bound
  ## on sigma2_T. Hoffman D (2010), Technometrics 52(3):303-312.
  hw <- if (sides == "two") {
    stats::qnorm((1 + content) / 2) * sqrt(comp$var_mean + sigma2_t) *
      sqrt(1 + sqrt(inner) / sigma2_t)
  } else {
    stats::qnorm(content) * sqrt(sigma2_t + sqrt(inner)) +
      stats::qnorm(conf) * sqrt(comp$var_mean)
  }

  new_re_interval(
    estimate = comp$mean,
    lower    = if (sides == "upper") -Inf else comp$mean - hw,
    upper    = if (sides == "lower")  Inf else comp$mean + hw,
    type     = "TI",
    content  = content,
    conf     = conf,
    df       = NA_real_,                          # normal-based; no t df
    sides    = sides,
    target   = comp$target,
    design   = comp$design,
    call     = match.call()
  )
}
