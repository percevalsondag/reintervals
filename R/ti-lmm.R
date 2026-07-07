#' Tolerance interval for a linear mixed model (dispatched by design)
#'
#' The beta-content / gamma-confidence tolerance interval of Francq, Lin & Hoyer
#' (2019, Eq. 26), dispatched by the detected random structure (one random
#' factor; two nested or two crossed random factors). Designs outside that
#' closed-form catalog --- and, in v1, unbalanced nested/crossed --- return `NA`
#' bounds with an explanatory `note` attribute; use `intervals_lmm()` for the
#' prediction interval there.
#'
#' This is a thin wrapper: it reaches `lme4` only through `vc_extract()`'s
#' internals and delegates the arithmetic to the method-pure `ti_francq()`.
#'
#' @param model An `lmerMod` fit (random intercepts only; random slopes are
#'   rejected).
#' @param newdata Row(s) giving the fixed-effect combination(s). Defaults to the
#'   first model-frame row.
#' @param level Content level (`beta`, the proportion of the population covered).
#'   Default 0.95.
#' @param conf Confidence level (`gamma`). Default 0.95.
#' @param target `"observable"` (default; future-value variance = all
#'   components). `"true_value"` (between-level only) is **not supported for the
#'   tolerance interval for the fixed-slope design** and is rejected; the prediction interval does
#'   support it.
#'
#' @return A data frame with the `newdata` predictor columns plus `type = "TI"`,
#'   `estimate`, `lower`, `upper`, `level`, `conf`. `attr(., "design")` is the
#'   detected structure, `attr(., "balanced")` the balance flag,
#'   `attr(., "singular")` whether `lme4::isSingular()` flags the fit, and
#'   `attr(., "note")` any reason the TI was not computed.
#'
#' @seealso `intervals_lmm()` for CI/PI; `ti_francq()`.
#' @references Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and
#'   tolerance in linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' ti_lmm(fm, newdata = data.frame(Days = 0), P = 0.99)
#' @keywords internal
#' @noRd
.ti_ems <- function(model, newdata = NULL, level = 0.95, conf = 0.95,
                    target = c("observable", "true_value")) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted `lmerMod` (an lmer fit).", call. = FALSE)
  }
  target <- match.arg(target)
  if (target == "true_value") {
    stop("true_value tolerance intervals are not supported for the random-intercept fixed-slope design (only ",
         "`target = \"observable\"`). The prediction interval (intervals_lmm / ",
         "pi_newobs) does support true_value.", call. = FALSE)
  }
  .assert_random_intercepts(model)
  desc <- .design_of(model)
  comps <- .var_components(model)
  dc <- .design_components(desc, comps, target, .counts_if_unbalanced(model, desc))
  fx <- .fixed(model, newdata)
  singular <- .is_singular(model)
  nr <- nrow(fx$L)

  bounds <- function(i) {
    if (is.null(dc$ems)) return(c(NA_real_, NA_real_))
    comp <- re_components(
      components = comps, dfs = c(pi = Inf),     # df unused by ti_francq (normal-based)
      mean = fx$fit[i], var_mean = fx$var_fix[i],
      coefs = dc$coefs, ems = dc$ems, target = target, design = desc$type
    )
    ti <- ti_francq(comp, content = level, conf = conf)
    c(ti$lower, ti$upper)
  }
  bm <- vapply(seq_len(nr), bounds, numeric(2))

  out <- cbind(fx$newdata, data.frame(
    type = "TI", estimate = fx$fit, lower = bm[1, ], upper = bm[2, ],
    level = level, conf = conf
  ), row.names = NULL)
  attr(out, "design") <- desc$type
  attr(out, "balanced") <- isTRUE(desc$balanced)
  attr(out, "singular") <- singular
  if (!is.null(dc$note)) attr(out, "note") <- dc$note

    if (is.null(dc$ems)) {
    warning(dc$note, call. = FALSE)
  } else if (singular) {
    warning("Singular fit: a variance component is at the boundary (zero); the ",
            "tolerance interval reflects a near-degenerate model.", call. = FALSE)
  }
  out
}
