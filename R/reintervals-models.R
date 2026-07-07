#' Mixed-model classes reintervals computes intervals for
#'
#' `reintervals` covers four linear-mixed-model classes for a single grouping
#' factor with fixed part `~ 1 + time`. The interval verbs ([ci_lmm()],
#' [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()]) detect the class from the
#' fitted model and route to the engine that is exact/validated for it; the class
#' is reported back, route-independently, in the `design` column of the result.
#' This page records the mapping from each model definition to its `design` tag.
#'
#' @section Model classes:
#'
#' \tabular{lll}{
#'   **lme4 formula** \tab **`design` tag** \tab **description** \cr
#'   `y ~ 1 + (1 | g)` (and nested/crossed) \tab `random_intercept` \tab random
#'     intercept only \cr
#'   `y ~ time + (1 | g)` \tab `random_intercept_fixed_slope` \tab random
#'     intercept with a fixed (population) slope \cr
#'   `y ~ time + (1 + time | g)` \tab `random_slope_correlated` \tab random
#'     intercept and random slope, correlated \cr
#'   `y ~ time + (1 | g) + (0 + time | g)` \tab `random_slope_independent` \tab
#'     random intercept and random slope, independent \cr
#' }
#'
#' @section Random-intercept topology (the `grouping` attribute):
#' For random-intercept models (no random slope) the one-way / nested / crossed
#' topology is reported separately in `attr(result, "grouping")`
#' (`"oneway"`, `"nested"`, or `"crossed"`); it is `NA` for random-slope models.
#'
#' @section Extrapolation in time (slope models):
#' For the random-slope and fixed-slope models the between-group variance
#' `V_G(t0) = z0' Sigma z0` is **quadratic in the time covariate**. Evaluating an
#' interval at a time point outside the observed range therefore extrapolates
#' that quadratic; the package computes it without a guard or warning, so
#' out-of-range `t0` should be requested knowingly.
#'
#' @section new_group_mean (CInew / new-batch-mean):
#' [new_group_mean_lmm()] is the **CInew / new-batch-mean** interval --- a
#' prediction interval for the mean of one new, as-yet-unobserved group/batch
#' (the prediction interval with the residual component removed). It is listed
#' here under both names so it is discoverable as "CInew" or "new batch mean".
#'
#' @seealso [ci_lmm()], [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()],
#'   [lmm_predict()], [ti_anova_raw()], [ti_gpq_raw()].
#' @name reintervals-models
NULL
