#' Montes-Burdick-Leblond (2019) Appendix A stability data
#'
#' The worked-example stability dataset from Montes, Burdick & Leblond (2019),
#' Appendix A: 11 lots measured over time, mixing genuine multi-timepoint
#' stability lots (A-D) with release-only single-observation lots (E-K). It is
#' the published oracle for the random-intercept + fixed-slope tolerance
#' interval: a two-sided P = 0.99 / gamma = 0.95 interval at `Month = 12`
#' reproduces the published `(57.25, 75.41)` (via [ti_anova_raw()] or
#' [ti_lmm()]); intermediates `beta = 0.457`, `n_E = 6.729`, `U = 10.820`.
#'
#' @format A data frame with 26 rows and 3 columns:
#' \describe{
#'   \item{Lot}{lot identifier (character), 11 lots `A`-`K`.}
#'   \item{Month}{time since release (numeric, months).}
#'   \item{Y}{measured response.}
#' }
#' @source Montes RO, Burdick RK, Leblond DJ (2019). Simple approach to calculate
#'   random effects model tolerance intervals to set release and shelf-life
#'   specification limits of pharmaceutical products. *PDA Journal of Pharmaceutical Science and Technology*. \doi{10.5731/pdajpst.2018.008839}.
#'   Appendix A (reproduced from their Table II).
#' @seealso [ti_anova_raw()], [ti_lmm()], [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a)
#' ti_lmm(fm, newdata = data.frame(Month = 12), P = 0.99)
"mbl_appendix_a"

#' Oliva-Aviles & Hauser (2025) Section 5 DGP specification (GPQ reference)
#'
#' The published data-generating specification for the random-intercept +
#' random-slope single-temperature stability model of Oliva-Aviles & Hauser
#' (2025), Section 5, at variance ratio `kappa = 1`. The paper reports a
#' simulation / coverage study rather than a single raw dataset, so the
#' specification (not a dataset) is stored: it makes the closed-form content
#' variance exactly reproducible and lets reference datasets be simulated
#' deterministically. The GPQ tolerance interval targets the content variance
#' `tau^2(t0) = z' V z`, `z = (1, t0)`; at `t0 = 20` this is `8.439388`
#' (the paper's `8.439`), reproduced exactly by the engine's `V_G` machinery.
#'
#' @format A named list:
#' \describe{
#'   \item{beta0, beta1}{fixed intercept and slope (100, -0.03).}
#'   \item{s2_0, s01, s2_1}{random intercept variance, intercept-slope
#'     covariance, slope variance of `V` (at `kappa = 1`).}
#'   \item{s2e}{residual variance (1).}
#'   \item{kappa}{the variance ratio (1).}
#'   \item{times}{the staggered measurement times.}
#'   \item{t0}{the evaluation time (20).}
#'   \item{tau2_at_t0}{the closed-form content variance at `t0` (8.439388).}
#'   \item{citation}{the source reference.}
#' }
#' @source Oliva-Aviles C, Hauser P (2025). *Technometrics* 67(2):193-202 \doi{10.1080/00401706.2024.2407324},
#'   Section 5.
#' @seealso [ti_gpq_raw()], [reintervals-models].
#' @examples
#' # the published content-variance functional, reproduced exactly:
#' oliva_hauser_dgp$tau2_at_t0
"oliva_hauser_dgp"

#' Hoffman & Kringle (2007) Table III assay data
#'
#' The one-way random-effects worked example of Hoffman & Kringle (2007),
#' Table III (their Example 1): an assay concentration measured in 6 analytical
#' runs with 3 replicates each. It is the published oracle for the
#' Hoffman-Kringle (2005) two-sided tolerance interval: with content
#' `P = 0.95` and confidence `0.90`, [ti_lmm()] (default `"anova-mls"`) or
#' [ti_anova_raw()] with `time = NULL` reproduces the published
#' `(0.846, 1.114)`; with `P = 0.667` it reproduces `(0.914, 1.046)`.
#'
#' The same data appear in Hoffman & Kringle (2005), Table 6, with run 5 /
#' replicate 3 misprinted as `1.020`. The correct value is `0.998`: only it is
#' consistent with the published mean `0.9798` and mean squares (between
#' `0.003251`, within `0.001167`). Francq, Lin & Hoyer (2019), Sections 5.2-5.3,
#' used the misprinted value, which is why their intervals for this example
#' differ slightly from Hoffman & Kringle's.
#'
#' @format A data frame with 18 rows and 3 columns:
#' \describe{
#'   \item{run}{analytical run (factor, 6 levels).}
#'   \item{replicate}{replicate number within run (integer, 1-3).}
#'   \item{conc}{measured concentration.}
#' }
#' @source Hoffman D, Kringle R (2007). A total error approach for the
#'   validation of quantitative analytical methods. *Pharmaceutical Research*
#'   24(6):1157-1164. \doi{10.1007/s11095-007-9242-3}. Table III.
#' @references Hoffman D, Kringle R (2005). Two-sided tolerance intervals for
#'   balanced and unbalanced random effects models. *Journal of
#'   Biopharmaceutical Statistics* 15(2):283-293. \doi{10.1081/BIP-200048826}
#'
#'   Francq BG, Lin D, Hoyer W (2019). Confidence, prediction, and tolerance in
#'   linear mixed models. *Statistics in Medicine* 38(30):5603-5622.
#'   \doi{10.1002/sim.8386}
#' @seealso [ti_lmm()], [ti_anova_raw()], [reintervals-models].
#' @examples
#' fm <- lme4::lmer(conc ~ 1 + (1 | run), data = hk_assay)
#' ti_lmm(fm, P = 0.95, level = 0.90)   # Hoffman-Kringle: (0.846, 1.114)
"hk_assay"
