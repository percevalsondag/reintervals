#' EMS decomposition for the two-fold nested design (pure, balanced)
#'
#' Builds the expected-mean-square linear combination for `y ~ FE + (1|a/b)`
#' (Francq 2019, Eq. 29-31), balanced data. Pure --- no `lme4`. The unbalanced
#' nested EMS coefficients are deferred to v2.
#'
#' @param desc A nested design descriptor from `classify_design()` (uses
#'   `grp_alpha`, `grp_beta`, `A`, `B`, `n`).
#' @param components Named numeric variance components, including `Residual`.
#' @param target `"observable"` or `"true_value"` (selects `coefs`).
#' @return `list(ems = list(ms, k, df), coefs)`.
#' @noRd
ems_nested <- function(desc, components, target = "observable") {
  ga <- desc$grp_alpha
  gb <- desc$grp_beta
  s2a <- components[[ga]]
  s2b <- components[[gb]]
  s2e <- components[["Residual"]]
  A <- desc$A
  B <- desc$B
  n <- desc$n
  emsb <- s2e + n * s2b
  emsa <- s2e + n * s2b + n * B * s2a
  nm <- c(ga, gb, "Residual")
  ems <- list(
    ms = stats::setNames(c(emsa, emsb, s2e), nm),
    k  = stats::setNames(c(1 / (n * B), 1 / n - 1 / (n * B), 1 - 1 / n), nm),
    df = stats::setNames(c(A - 1, A * (B - 1), A * B * (n - 1)), nm)
  )
  list(ems = ems, coefs = .target_coefs(components, target))
}
