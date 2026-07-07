#' EMS decomposition for the one-random-factor design (pure)
#'
#' Builds the expected-mean-square linear combination that `ti_francq()` (Eq. 26)
#' consumes, for `y ~ FE + (1|a)`. Valid for **balanced and unbalanced** data:
#' the effective replicate `n0` (Burdick-Graybill) absorbs the imbalance
#' (Francq 2019, Section 5.3). Pure --- no `lme4`; testable from hand-entered
#' counts and components.
#'
#' @param desc A one-way design descriptor from `classify_design()` (uses `grp`,
#'   `n0`, `A`, `N`).
#' @param components Named numeric variance components, including `Residual`.
#' @param target `"observable"` or `"true_value"` (selects `coefs`).
#' @return `list(ems = list(ms, k, df), coefs)`.
#' @noRd
ems_oneway <- function(desc, components, target = "observable") {
  grp <- desc$grp
  s2a <- components[[grp]]
  s2e <- components[["Residual"]]
  n0 <- desc$n0
  nm <- c(grp, "Residual")
  ems <- list(
    ms = stats::setNames(c(s2e + n0 * s2a, s2e), nm),
    k  = stats::setNames(c(1 / n0, 1 - 1 / n0), nm),
    df = stats::setNames(c(desc$A - 1, desc$N - desc$A), nm)
  )
  list(ems = ems, coefs = .target_coefs(components, target))
}
