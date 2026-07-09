#' Closed-form tolerance interval from a slope (vfun) container
#'
#' The method-pure slope-engine counterpart to `ti_francq()`: it consumes a
#' `vfun_components` object and computes the modified large-sample (MLS)
#' tolerance interval for the random-slope stability model, evaluated at the
#' container's time point `t0`. It is the closed-form (deterministic) slope TI. 
#' Like `ti_francq()` it is method-pure --- it reads the numbers off the container 
#' and never touches `lme4`.
#'
#' This constructor serves the **observable** (future-observation) target, content
#' variance `V_T = V_G + s2e`. The batch-mean (`true_value`) target is not served
#' here; the batch-mean interval is the GPQ engine's job, via the dispatcher.
#'
#' The interval is
#' \deqn{\mu(t_0) \pm z_{(1+P)/2}\,\sqrt{1 + 1/n_E}\,\sqrt{U},}
#' with the Graybill--Wang content-variance bound
#' `U = V_T + sqrt((H_G V_G)^2 + (H_e s2e)^2)`, chi-square tail factors
#' `H_j = nu_j / qchisq(1 - gamma, nu_j) - 1`, and effective sample size
#' `n_E = V_T / V_F` (carried on the container).
#'
#' @param comp A `vfun_components` object (from `vfun_extract()`). Uses `mean`,
#'   `var_mean`, `V_G`, `V_T`, `n_E`, `phi[["s2e"]]`, and `dfs[c("G","e")]`.
#' @param content Content proportion `P` (e.g. 0.99).
#' @param conf Confidence `gamma` (e.g. 0.95).
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`.
#' @return A one-row data frame: `time`, `type = "TI"`, `estimate`, `lower`,
#'   `upper`, `content`, `conf`, `df`, `sides`, `design` (the model type),
#'   `method = "reml-mls"`, and the `singular` / `boundary` flags.
#' @seealso `ti_francq()` (the EMS counterpart), `vfun_extract()`.
#' @references Graybill FA, Wang CM (1980). Confidence intervals on nonnegative
#'   linear combinations of variances. Journal of the American Statistical
#'   Association 75(372):869-873.
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'                  data = lme4::sleepstudy)
#' ti_vfun(vfun_extract(fm, t0 = 5), content = 0.99, conf = 0.95)
#' @noRd
ti_vfun <- function(comp, content = 0.99, conf = 0.95,
                    sides = c("two.sided", "lower", "upper")) {
  if (!inherits(comp, "vfun_components")) {
    stop("`comp` must be a `vfun_components` object (see ?vfun_extract).",
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
  tgt <- attr(comp, "target")
  if (!is.null(tgt) && !identical(tgt, "observable")) {
    stop("ti_vfun() implements the observable (future-observation) target only ",
         "in this stage; the batch-mean (true_value) target is the GPQ ",
         "engine's job.", call. = FALSE)
  }
  if (!all(c("G", "e") %in% names(comp$dfs))) {
    stop("`comp$dfs` must carry the \"G\" (between) and \"e\" (residual) df ",
         "entries; build the container with vfun_extract().", call. = FALSE)
  }

  mu  <- comp$mean
  V_F <- comp$var_mean
  V_G <- comp$V_G
  V_T <- comp$V_T
  s2e <- comp$phi[["s2e"]]
  nu_G <- comp$dfs[["G"]]
  nu_e <- comp$dfs[["e"]]
  nE   <- comp$n_E

  ## Graybill-Wang two-component content-variance bound (Graybill & Wang 1980)
  H_G <- nu_G / stats::qchisq(1 - conf, nu_G) - 1
  H_e <- nu_e / stats::qchisq(1 - conf, nu_e) - 1
  U   <- V_T + sqrt((H_G * V_G)^2 + (H_e * s2e)^2)
  zc  <- if (sides == "two.sided") stats::qnorm((1 + content) / 2) else stats::qnorm(content)
  hw  <- zc * sqrt(1 + 1 / nE) * sqrt(U)

  data.frame(
    time = comp$t0, type = "TI", estimate = mu,
    lower = if (sides == "upper") -Inf else mu - hw,
    upper = if (sides == "lower")  Inf else mu + hw,
    content = content, conf = conf, df = min(nu_G, nu_e), sides = sides,
    design = comp$type, method = "reml-mls",
    singular = isTRUE(comp$singular), boundary = isTRUE(comp$boundary),
    vg_floored = isTRUE(attr(comp, "vg_floored")),
    stringsAsFactors = FALSE
  )
}
