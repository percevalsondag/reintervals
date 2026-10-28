#' reintervals: intervals for linear mixed models
#'
#' Confidence, prediction, tolerance, and new-group-mean (CInew / new-batch-mean)
#' intervals for linear mixed models fitted with [lme4::lmer()], returned through
#' one uniform tidy interface. The right interval engine is selected
#' automatically for the fitted model; but it can be changed.
#'
#' @section The four interval verbs:
#' \describe{
#'   \item{[ci_lmm()]}{confidence interval for the mean.}
#'   \item{[pi_lmm()]}{prediction interval for a future observation.}
#'   \item{[ti_lmm()]}{tolerance interval (content `P`); `over = "group_mean"`
#'     gives the batch-mean tolerance interval.}
#'   \item{[new_group_mean_lmm()]}{the CInew / new-batch-mean interval (an
#'     interval for the true mean of one new group: a new innermost cell, or a
#'     new lot at `t0`).}
#' }
#' [lmm_predict()] is a `predict()`-style wrapper over the four, selected by
#' `interval =`. [lmm_interval_control()] supplies the GPQ tuning (seed, `M`)
#' and the confidence-interval df rule (`ci_df`).
#' Every verb returns the same `lmm_interval` data frame (see [ci_lmm()] for the
#' column/attribute contract).
#'
#' @section Engine and scope map:
#' The engine is chosen from the fitted model's structure and can differ
#' between the verbs (full mapping in `vignette("choosing-an-interval")`;
#' design tags in [reintervals-models]). Engines: `"ems-mls"`
#' (expected-mean-square / modified large-sample), `"anova-mls"` (unweighted
#' lot-mean ANOVA: Montes et al. with a fixed slope, Hoffman & Kringle without
#' one), `"reml-mls"` (REML closed form), `"gpq"` (generalized pivotal quantity,
#' Monte Carlo).
#' \tabular{lllll}{
#'   **model / situation** \tab **CI** \tab **PI** \tab **TI** \tab
#'     **new-group mean** \cr
#'   random intercept, one-way, intercept only (`y ~ 1 + (1 | g)`) \tab
#'     `ems-mls` \tab `ems-mls` \tab `anova-mls` (Hoffman & Kringle) \tab
#'     `anova-mls` \cr
#'   random intercept, one-way with covariates \tab `ems-mls` \tab `ems-mls`
#'     \tab `ems-mls` \tab `ems-mls` (Gm) \cr
#'   random intercept, nested \tab `ems-mls` \tab `ems-mls` \tab `ems-mls`
#'     \tab `ems-mls` (Gm) \cr
#'   random intercept, crossed with interaction \tab `ems-mls` \tab `ems-mls`
#'     \tab `ems-mls` \tab `ems-mls` (Gv) \cr
#'   random intercept, crossed without interaction or 4+ grouping factors \tab
#'     `ems-mls` \tab `ems-mls` \tab error \tab error \cr
#'   random intercept + fixed slope, clean balanced \tab `ems-mls` \tab
#'     `ems-mls` \tab `ems-mls` \tab `anova-mls` \cr
#'   random intercept + fixed slope, single-observation / staggered \tab
#'     `ems-mls` (warns with single-observation lots) \tab `ems-mls` (warns
#'     likewise) \tab `anova-mls` \tab `anova-mls` \cr
#'   random slope, future-observation, non-singular \tab `reml-mls` \tab
#'     `reml-mls` \tab `reml-mls` \tab `reml-mls` \cr
#'   random slope, singular fit or TI Satterthwaite df unavailable \tab
#'     `reml-mls` \tab `reml-mls` \tab `gpq` \tab `reml-mls` \cr
#'   random slope, batch-mean (`over = "group_mean"`, TI only) \tab --- \tab
#'     --- \tab `gpq` \tab --- \cr
#' }
#' The random-slope tolerance interval switches to GPQ when the fit is
#' singular or, at a given `t0`, when the closed form's Satterthwaite df for the
#' between-group or residual variance is unavailable (including when the
#' variance-component covariance is not positive definite); the reason is in
#' `attr(., "diagnostics")$route_reason`. The random-intercept confidence
#' interval uses a Satterthwaite df, falling back to the SAS-style containment
#' df (`A - 1` one-way or nested, `min(A, B) - 1` crossed) when that df is
#' unavailable; Kenward-Roger df only on request
#' (`lmm_interval_control(ci_df = "kenward-roger")`).
#'
#' **One-sided intervals** (`sides = "lower"` / `"upper"`) are available from
#' every verb. CI, PI, and new-group mean keep the two-sided center, variance,
#' and df with `t_{1 - alpha}`. One-sided tolerance limits: `anova-mls` uses
#' Hoffman's (2010) one-sided limit; `ems-mls` and random-slope models
#' (`reml-mls`, the default on every fit) use the package's extension of
#' Hoffman's limit to the Francq and REML closed-form bounds, validated by
#' simulation. See [reintervals-models].
#'
#' **Scope.** General fixed effects are supported for random-intercept
#' models; random-slope models are supported for a single continuous covariate (the
#' stability design). The two tier-2 direct engines [ti_anova_raw()] and
#' [ti_gpq_raw()] compute from raw `(y, time, group)` vectors without a fit.
#'
#' @section Validation:
#' Every shipped interval is checked against a published or external-source
#' oracle on CI: the Francq, Lin & Hoyer (2019) worked numbers; the
#' Montes-Burdick-Leblond (2019) Appendix A tolerance interval `(57.25, 75.41)`;
#' the random-slope intervals against committed reference values to ~1e-9; and the
#' GPQ content variance against the Oliva-Aviles & Hauser (2025) closed form
#' (`tau^2 = 8.439`).
#'
#' @keywords internal
"_PACKAGE"
