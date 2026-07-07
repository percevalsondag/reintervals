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
#'   \item{[new_group_mean_lmm()]}{the CInew / new-batch-mean interval (a
#'     prediction interval for one new group's mean).}
#' }
#' [lmm_predict()] is a `predict()`-style wrapper over the four, selected by
#' `interval =`. [lmm_interval_control()] supplies the GPQ tuning (seed, `M`).
#' Every verb returns the same `lmm_interval` data frame (see [ci_lmm()] for the
#' column/attribute contract).
#'
#' @section Engine and scope map:
#' The engine is chosen from the fitted model's structure (see
#' [reintervals-models] for the full mapping):
#' \tabular{ll}{
#'   **model class** \tab **engine (`method`)** \cr
#'   random intercept (any fixed effects; one-way / nested / crossed) \tab
#'     expected-mean-square / modified large-sample (`"ems-mls"`) \cr
#'   random intercept + fixed slope, clean balanced \tab `"ems-mls"` \cr
#'   random intercept + fixed slope, single-observation / unbalanced \tab
#'     unweighted lot-mean ANOVA closed form (`"anova-mls"`) \cr
#'   random slope, future-observation \tab REML closed form (`"reml-mls"`);
#'     GPQ (`"gpq"`) when the fit is singular \cr
#'   random slope, batch-mean \tab GPQ Monte Carlo (`"gpq"`) \cr
#' }
#' **v1 scope.** General fixed effects are supported for random-intercept
#' models; random-slope models are supported for a single time covariate (the
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
