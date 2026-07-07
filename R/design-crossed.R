#' EMS decomposition for the two-crossed-with-interaction design (pure, balanced)
#'
#' Builds the expected-mean-square linear combination for
#' `y ~ FE + (1|a) + (1|b) + (1|a:b)` (Francq 2019, Eq. 32-34), balanced data.
#' Pure --- no `lme4`. The unbalanced crossed EMS coefficients are deferred to
#' v2. 
#'
#' @param desc A crossed design descriptor from `classify_design()` (uses
#'   `grp_A`, `grp_B`, `grp_AB`, `A`, `B`, `n`).
#' @param components Named numeric variance components, including `Residual`.
#' @param target `"observable"` or `"true_value"` (selects `coefs`).
#' @return `list(ems = list(ms, k, df), coefs)`.
#' @noRd
ems_crossed <- function(desc, components, target = "observable") {
  ga <- desc$grp_A
  gb <- desc$grp_B
  gab <- desc$grp_AB
  v1 <- components[[ga]]
  v2 <- components[[gb]]
  s2ab <- components[[gab]]
  s2e <- components[["Residual"]]
  A <- desc$A
  B <- desc$B
  n <- desc$n
  emsab <- s2e + n * s2ab
  emsa <- s2e + n * s2ab + n * B * v1
  emsb <- s2e + n * s2ab + n * A * v2
  nm <- c(ga, gb, gab, "Residual")
  ems <- list(
    ms = stats::setNames(c(emsa, emsb, emsab, s2e), nm),
    k  = stats::setNames(
      c(1 / (n * B), 1 / (n * A), 1 / n - 1 / (n * A) - 1 / (n * B), 1 - 1 / n),
      nm
    ),
    df = stats::setNames(
      c(A - 1, B - 1, (A - 1) * (B - 1), A * B * (n - 1)), nm
    )
  )
  list(ems = ems, coefs = .target_coefs(components, target))
}
