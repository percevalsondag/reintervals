# =============================================================================
# The public verbs ci_lmm, pi_lmm, ti_lmm, new_group_mean_lmm. Plus the
# lmm_predict wrapper and lmm_interval_control. These are thin wrappers: they
# take the shared arguments, dispatch to the right engine, and tag the result
# with its method, design, and grouping. The numbers are the engine's own; the
# wrappers don't change them.
# =============================================================================

#' Tuning controls for the interval engines
#'
#' Tuning knobs only --- the RNG seed and Monte-Carlo size used when the
#' GPQ algorithm runs, the degrees-of-freedom rules of the confidence and
#' prediction intervals, and the validation arms of the new-group-mean
#' interval. Algorithm *selection* is the `method` argument, not this.
#'
#' @section Validation arms:
#' The non-default values of `ci_df`, `slope_ci_df`, `pi_df`, `ngm_method` and
#' `slope_ngm_method` are validation
#' arms: their coverage is under evaluation in the 1.0.0 validation study, and
#' the defaults are the validated rules. Applying a non-default arm signals a
#' warning of class `reintervals_warning_unvalidated` on every call (muffle it
#' by class), except `ci_df = "kenward-roger"` on random-intercept models,
#' supported since 0.2.0. Each option applies to one model structure and verb:
#' `ci_df` to [ci_lmm()] on random-intercept models (including the fixed-slope
#' design), `slope_ci_df` to [ci_lmm()] on random-slope models, `pi_df` to
#' [pi_lmm()] on random-intercept models, `ngm_method` to
#' [new_group_mean_lmm()] on random-intercept models and `slope_ngm_method` to
#' [new_group_mean_lmm()] on random-slope models. Elsewhere an option is a
#' no-op, never an error, so one control object can be reused across models
#' and verbs; so is a new-group-mean arm that is not defined for the design
#' (`"hts"` on anything but an intercept-only one-way model with at least three
#' groups, `"classical"` on unbalanced groups): the default is computed. The
#' diagnostics of every confidence, prediction and new-group-mean interval
#' record each option as `ci_df_rule`, `slope_ci_df_rule`, `pi_df_rule`,
#' `ngm_arm` and `slope_ngm_arm`: the value applied, or `"not_applicable"`;
#' `df` is the df actually used. The `method` column keeps the default route's
#' tag; the arm is read from these fields.
#'
#' @param seed Optional integer seed for the GPQ draws (reproducibility / GxP
#'   traceability). With the default `NULL`, one integer seed is drawn from the
#'   caller's RNG stream (`sample.int(.Machine$integer.max, 1L)`), so a
#'   `set.seed()` before the call still controls the result; this advances the
#'   caller's stream by exactly one draw. The GPQ draws then run under that seed
#'   (or under the supplied one) and the caller's RNG state is restored
#'   afterwards. Either way the seed actually used is reported as
#'   `attr(., "diagnostics")$gpq_seed`, and the RNG kind in effect
#'   (the [RNGkind()] triple) as `$rng_kind`. Reproducibility contract:
#'   `lmm_interval_control(seed = <gpq_seed>)` with the same `M`, under the same
#'   `rng_kind`, returns an identical interval.
#' @param M Number of GPQ Monte-Carlo realizations, an integer `>= 2000`.
#'   Default 20000, the size validated by simulation for the random-slope
#'   tolerance interval (one- and two-sided, both targets). Smaller values run
#'   faster but are not validated. Values below 2000 are an error: the
#'   Monte-Carlo standard error of each limit (`mc_se_lower`, `mc_se_upper` in
#'   the diagnostics) is computed from 20 MC sections of consecutive draws, and
#'   below 2000 a section has fewer than 100 draws, too few for a 0.99 tail
#'   quantile. The MC SE is the run-to-run jitter of the Monte-Carlo
#'   computation, not statistical uncertainty: a limit within approximately 2
#'   MC SE of a specification limit has a seed-dependent pass/fail. At
#'   `M = 2000` the sectioning estimate ran about 8% below the observed
#'   run-to-run SD in a sleepstudy check (sections of 100 draws); the gap is
#'   expected to shrink at the default `M = 20000`. When `M` is not divisible by 20, the last `M %% 20` draws are
#'   left out of the sections only; the limit itself uses all `M` draws.
#' @param ci_df Denominator-df rule for the confidence interval ([ci_lmm()]) on
#'   random-intercept models: `"satterthwaite"` (default; Satterthwaite with the
#'   observed-information covariance of the variance components, falling back to
#'   the SAS-style containment df --- e.g. `A - 1` for a one-way or nested
#'   design, `min(A, B) - 1` for a crossed one --- when it is unavailable, as on
#'   a singular fit) or `"kenward-roger"` (Kenward & Roger 1997 as implemented
#'   in \pkg{pbkrtest}: the standard error uses the KR-adjusted covariance
#'   `pbkrtest::vcovAdj()` and the df is the KR df, `pbkrtest::Lb_ddf()`;
#'   requires \pkg{pbkrtest} and errors if it is not installed --- never a
#'   silent fallback). SAS PROC MIXED `DDFM=KR` can differ on unbalanced data:
#'   in the package's checks it reported the Satterthwaite df with an adjusted
#'   standard error. The df source is reported in
#'   `attr(., "diagnostics")$df_method`. A validation arm (see below):
#'   `"containment"`, the containment df alone (the rank rule above: the
#'   smallest random-stratum df, each fixed-effect column that is constant
#'   within a factor's levels costing one df there --- `A - 1` one-way and
#'   nested, `A - 2` with one between-group covariate, `min(A, B) - 1` crossed,
#'   `K - 1` for the fixed-slope design), with the unadjusted standard error.
#'   Not applied to random-slope models (see `slope_ci_df`).
#' @param slope_ci_df Denominator-df rule for the confidence interval
#'   ([ci_lmm()]) on random-slope models: `"hessian-satterthwaite"` (default;
#'   the generalized Satterthwaite df from the observed-information covariance
#'   of the variance components, falling back to `K - 2`, `K` lots, when it is
#'   unavailable; reported as `df_rule = "satterthwaite"` / `"containment"`).
#'   Validation arms: `"kenward-roger"` (\pkg{pbkrtest} `vcovAdj()` +
#'   `Lb_ddf()` on the random-slope fit; the standard error uses the
#'   KR-adjusted covariance; requires \pkg{pbkrtest}) and `"k-minus-1"`
#'   (`K - 1` with the unadjusted standard error: the exact df for a balanced
#'   random-coefficient model, Chow & Shao 1991). `df_rule` reports the arm.
#'   Not applied to random-intercept models (see `ci_df`).
#' @param pi_df Denominator-df rule for the prediction interval ([pi_lmm()]) on
#'   random-intercept models: `"francq-satterthwaite"` (default; the Francq et
#'   al. 2019 generalized Satterthwaite df, with the containment df as its
#'   singular-fit fallback) or the validation arm `"containment"` (the same
#'   variance with the containment df of `ci_df = "containment"`). There is no
#'   Kenward-Roger prediction interval: KR adjusts fixed-effect contrasts only,
#'   and the prediction variance includes estimated variance components. Not
#'   applied to random-slope models, whose prediction interval is unchanged.
#' @param ngm_method New-group-mean construction ([new_group_mean_lmm()]) on
#'   random-intercept models: `"package"` (default; the 0.2.0 constructions,
#'   see [new_group_mean_lmm()]) or a validation arm:
#'   `"broad-inference"` (the SAS PROC MIXED broad-inference prediction of a new
#'   group, McLean, Sanders & Stroup 1991: `W = V_F` plus the REML estimates of
#'   every between-group variance component, with the Satterthwaite df
#'   `2 W^2 / (g' A g)`, `g` the gradient of `W` in the variance components and
#'   `A` their observed-information covariance, as DDFM=SATTERTH; when that df is
#'   not finite or below 1 the limits are `NA` and `ngm_df_failed = TRUE`, with
#'   no fallback; every random-intercept design the default accepts),
#'   `"hts"` (Higgins, Thompson & Spiegelhalter 2009 prediction interval,
#'   `mu +/- t_{I-2} sqrt(sigma2_L + V_F)`, with the REML between-group variance
#'   in place of their DerSimonian-Laird `tau^2`; intercept-only one-way models
#'   with `I >= 3` groups) or `"classical"` (Hahn & Meeker 1991 prediction
#'   interval on the group means, `Ybar +/- t_{I-1} S_L sqrt(1 + 1/I)`;
#'   intercept-only one-way models with balanced groups, where it equals the
#'   default).
#' @param slope_ngm_method New-group-mean construction on random-slope models:
#'   `"reml"` (default; `mu(t0) +/- t(nu_G) sqrt(V_F + V_G)`) or a validation
#'   arm: `"broad-inference"` (the SAS PROC MIXED broad-inference prediction of
#'   a new lot, McLean, Sanders & Stroup 1991: the same center and variance
#'   `W = V_F + V_G` as `"reml"`, with the Satterthwaite df of the whole `W`,
#'   `2 W^2 / (g' A g)`, `g` the gradient of `V_F + V_G` in the
#'   variance-covariance parameters and `A` their observed-information
#'   covariance, as DDFM=SATTERTH; `"reml"` uses the df of `V_G` alone; when the
#'   df is not finite or below 1 the limits are `NA` and `ngm_df_failed = TRUE`,
#'   with no fallback; it accepts lots with a single time point, as the
#'   default does), `"gpq"` (a generalized prediction interval on the tolerance
#'   interval's pivots, Oliva-Aviles & Hauser 2025 with the generalized
#'   prediction construction of Krishnamoorthy & Mathew 2009:
#'   `R_mu(t0) + Z_new sqrt(R_VG(t0))`, limits its quantiles; the GPQ `seed` and
#'   `M` apply, with the same seed contract, common random numbers across `t0`
#'   and Monte-Carlo SE as [ti_lmm()]) or `"fv-anova"` (package-original:
#'   `Ybar_fit +/- t_{K-1} S_L sqrt(1 + 1/K)` on the per-lot OLS fitted values
#'   at `t0`, `K` lots; in a balanced design the exact t prediction interval
#'   for a new lot's fitted value at `t0`, whose variance
#'   `V_G + sigma^2 h(t0)` exceeds that of the new lot's true mean, so it is
#'   conservative for the new-group mean, by `sigma^2 h(t0) (1 + 1/K)`, and
#'   cannot truncate). `"gpq"` and `"fv-anova"` need two distinct time points
#'   per lot and use the per-lot OLS fits, not the fitted random-effects
#'   structure (an independent `||` fit gives the same interval as the
#'   correlated fit).
#' @return An `lmm_interval_control` list (passed to a verb's `control`
#'   argument); a list with elements `seed`, `M`, `ci_df`, `slope_ci_df`,
#'   `pi_df`, `ngm_method` and `slope_ngm_method`.
#' @seealso [ci_lmm()], [ti_lmm()], [ti_gpq_raw()].
#' @examples
#' lmm_interval_control(seed = 1, M = 5000)
#' lmm_interval_control(ci_df = "kenward-roger")
#' lmm_interval_control(slope_ci_df = "k-minus-1", pi_df = "containment")
#' lmm_interval_control(ngm_method = "broad-inference",
#'                      slope_ngm_method = "broad-inference")
#' @export
lmm_interval_control <- function(seed = NULL, M = 20000L,
                                 ci_df = c("satterthwaite", "kenward-roger",
                                           "containment"),
                                 slope_ci_df = c("hessian-satterthwaite",
                                                 "kenward-roger", "k-minus-1"),
                                 pi_df = c("francq-satterthwaite", "containment"),
                                 ngm_method = c("package", "broad-inference", "hts",
                                                "classical"),
                                 slope_ngm_method = c("reml", "broad-inference", "gpq",
                                                      "fv-anova")) {
  ci_df       <- match.arg(ci_df)
  slope_ci_df <- match.arg(slope_ci_df)
  pi_df       <- match.arg(pi_df)
  ngm_method  <- match.arg(ngm_method)
  slope_ngm_method <- match.arg(slope_ngm_method)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
                           abs(seed) > .Machine$integer.max)) {
    .abort("`seed` must be NULL or a single number in the integer range.",
           class = "reintervals_error_input")
  }
  ## M >= 2000: the MC standard error splits the draws into 20 MC sections, and
  ## below 2000 a section has fewer than 100 draws, too few for a 0.99 tail quantile
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || !is.finite(M) ||
        M != round(M) || M < 2000) {
    .abort("`M` must be a single integer >= 2000 (the Monte-Carlo standard error ",
           "uses 20 MC sections of at least 100 draws each).",
           class = "reintervals_error_input")
  }
  structure(list(seed = seed, M = as.integer(M), ci_df = ci_df,
                 slope_ci_df = slope_ci_df, pi_df = pi_df,
                 ngm_method = ngm_method, slope_ngm_method = slope_ngm_method),
            class = "lmm_interval_control")
}

## ---- validation arms (step-13 options) ------------------------------------
## The df arms of lmm_interval_control(), with their validated defaults. Each
## arm applies to one (model structure, verb) cell; anywhere else it is a no-op,
## recorded as "not_applicable" in the diagnostics field `<option>_rule`, so one
## control object can be reused across a simulation grid. The new-group-mean
## arms (B2) record the arm in `ngm_arm` / `slope_ngm_arm`.
.df_arms <- list(
  ci_df       = list(default = "satterthwaite", rule = "ci_df_rule",
                     ## the 0.2.0 options, which do not warn
                     no_warn = c("satterthwaite", "kenward-roger")),
  slope_ci_df = list(default = "hessian-satterthwaite", rule = "slope_ci_df_rule",
                     no_warn = "hessian-satterthwaite"),
  pi_df       = list(default = "francq-satterthwaite", rule = "pi_df_rule",
                     no_warn = "francq-satterthwaite"),
  ngm_method  = list(default = "package", rule = "ngm_arm", no_warn = "package"),
  slope_ngm_method = list(default = "reml", rule = "slope_ngm_arm", no_warn = "reml")
)

## The arm value of `opt` from a control object, or its default when the
## control predates the option (a hand-built list).
.arm_value <- function(control, opt) {
  v <- control[[opt]]
  if (is.null(v)) .df_arms[[opt]]$default else v
}

## Which arms apply to this model and verb: named character, option -> value
## for the applied arms only. A new-group-mean arm that is not defined for the
## design (e.g. ngm_method = "hts" on a nested fit, "classical" on unbalanced
## groups) is not applied: the default runs and the option is recorded as
## "not_applicable".
.applied_arms <- function(model, kind, control) {
  rs <- .has_random_slope(model)
  opts <- if (kind == "CI") (if (rs) "slope_ci_df" else "ci_df")
          else if (kind == "PI" && !rs) "pi_df"
          else if (kind == "new_group_mean") (if (rs) "slope_ngm_method" else "ngm_method")
          else character(0)
  out <- stats::setNames(vapply(opts, function(o) .arm_value(control, o), ""), opts)
  if ("ngm_method" %in% names(out) && !.ngm_arm_applies(model, out[["ngm_method"]])) {
    out <- out[names(out) != "ngm_method"]
  }
  out
}

## Diagnostics fields for every arm: the applied value or "not_applicable".
.arm_rules <- function(applied) {
  out <- lapply(names(.df_arms), function(o)
    if (o %in% names(applied)) applied[[o]] else "not_applicable")
  stats::setNames(out, vapply(.df_arms, function(a) a$rule, ""))
}

## One reintervals_warning_unvalidated per call for each applied non-default arm
## (not once per session: simulations muffle it by class).
.warn_unvalidated <- function(applied) {
  for (o in names(applied)) {
    a <- .df_arms[[o]]
    if (applied[[o]] %in% a$no_warn) next
    .warn(sprintf(paste0("%s = \"%s\" is a validation arm: its coverage is being ",
                         "evaluated for reintervals 1.0.0. The validated default is ",
                         "%s = \"%s\"."), o, applied[[o]], o, a$default),
          class = "reintervals_warning_unvalidated",
          data = list(option = o, value = applied[[o]]))
  }
  invisible(NULL)
}

## alpha / level as a mutually-exclusive pair.
## Returns alpha. both -> error; neither -> 0.05; level -> 1 - level.
.resolve_alpha <- function(alpha, level) {
  if (!is.null(alpha) && !is.null(level)) {
    .abort("Supply only one of `alpha` or `level` (they are the same quantity).",
           class = "reintervals_error_input")
  }
  if (is.null(alpha) && is.null(level)) return(0.05)
  if (!is.null(level)) {
    if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
          level <= 0 || level >= 1) {
      .abort("`level` must be a single number in (0, 1).",
             class = "reintervals_error_input")
    }
    return(1 - level)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
        alpha <= 0 || alpha >= 1) {
    .abort("`alpha` must be a single number in (0, 1).",
           class = "reintervals_error_input")
  }
  alpha
}

## defensive front-door class check (clear error at the door).
.assert_lmm <- function(model) {
  if (!inherits(model, "lmerMod")) {
    .abort("`model` must be a fitted lme4 model (an `lmerMod` from lmer()).",
           class = "reintervals_error_input")
  }
  invisible(TRUE)
}

## route-independent model-class tag (the `design` column vocabulary). Keyed on
## the random/fixed STRUCTURE, not the engine that computed the interval, so a
## given model gets the same design tag whichever route ran (kills oneway-vs-fixed-slope).
.model_class <- function(model) {
  cnms <- lme4::getME(model, "cnms")
  has_slope <- any(vapply(cnms, function(t) any(t != "(Intercept)"), logical(1)))
  if (has_slope) {
    return(if (length(cnms) == 1L) "random_slope_correlated"
           else "random_slope_independent")
  }
  fe <- lme4::fixef(model)
  if (length(setdiff(names(fe), "(Intercept)")) == 1L) {
    return("random_intercept_fixed_slope")
  }
  "random_intercept"
}

## the random-INTERCEPT topology (oneway/nested/crossed), moved off `design`
## into its own attribute. Defined only for models with NO random slope; NA
## otherwise.
.grouping_of <- function(model) {
  cnms <- lme4::getME(model, "cnms")
  has_slope <- any(vapply(cnms, function(t) any(t != "(Intercept)"), logical(1)))
  if (has_slope) return(NA_character_)
  ## a failed classification is NA with its reason kept (attr "reason"), which
  ## .unify_intervals() records as diagnostics$grouping_reason
  tryCatch(.design_of(model)$type, error = function(e)
    structure(NA_character_, reason = paste0("classification failed: ",
                                             conditionMessage(e))))
}

## interval-kind tag for the `type` column (always lower-case / snake).
.type_tag <- function(kind) {
  switch(kind, CI = , ci = "ci", PI = , pi = "pi", TI = , ti = "ti",
         new_group_mean = "new_group_mean", tolower(kind))
}

## Value columns an engine may emit; everything else in the native frame is an
## eval-point (predictor) column. "time" is the slope/raw engines' generic
## placeholder -- renamed to the actual covariate when a model is in hand.
.value_cols <- c("time", "type", "estimate", "lower", "upper", "level", "conf",
                 "content", "P", "sides", "df", "se", "method", "design",
                 "singular", "boundary", "vg_floored", "tau2_med", "s2e_med",
                 "frac_boundary", "K", "df_resid", "phi1", "c_resid", "n_E",
                 "df_method", "route_reason", "df_rule", "W", "mls_bound",
                 "n_na_draws", "mc_se_lower", "mc_se_upper", "ci_df_rule",
                 "slope_ci_df_rule", "pi_df_rule", "ngm_arm", "slope_ngm_arm",
                 "ngm_df_failed")

## The output unifier. Maps any engine's native frame to the single tidy
## contract; relabels only, never recomputes. Universal columns, in order:
##   <eval cols...>, type, estimate, lower, upper, level, P, sides, method,
##   design, singular, boundary, vg_floored
## Engine-specific diagnostics (df, se, beta_hat, GPQ audit, ...) go to
## attr(.,"diagnostics"). `model = NULL` serves the raw engines (design taken
## from the native row, grouping NA, the generic "time" eval column kept).
.unify_intervals <- function(native, model, kind, P_value = NA_real_,
                             method_tag = NULL, verb = NULL) {
  n <- nrow(native)
  ## NaN guard: an engine must never hand back a NaN estimate or limit (NA is
  ## allowed: e.g. a validation arm whose df failed reports NA limits; +/-Inf is
  ## a one-sided limit).
  for (col in c("estimate", "lower", "upper")) {
    bad <- which(is.nan(native[[col]]))
    if (length(bad)) {
      if (is.null(verb)) verb <- switch(tolower(kind), ci = "ci_lmm()", pi = "pi_lmm()",
                                        new_group_mean = "new_group_mean_lmm()",
                                        ti = "ti_lmm()", kind)
      .abort(verb, ": the engine returned NaN for `", col, "` in row(s) ",
             paste(bad, collapse = ", "), "; the interval is not available for this ",
             "fit.", class = "reintervals_error_numerical",
             data = list(verb = verb, column = col, rows = bad))
    }
  }
  ## eval-point columns (predictors); rename the generic "time" to the covariate
  eval_df <- native[, setdiff(names(native), .value_cols), drop = FALSE]
  if ("time" %in% names(native)) {
    tv <- if (!is.null(model)) {
      setdiff(names(lme4::fixef(model)), "(Intercept)")[1]
    } else "time"
    eval_df[[tv]] <- native$time
  }

  level_out  <- if (kind %in% c("ti", "TI")) native$conf else native$level
  method_out <- if ("method" %in% names(native)) native$method else method_tag
  design_out <- if (!is.null(model)) .model_class(model)
                else if ("design" %in% names(native)) native$design else NA_character_
  pick <- function(col) if (col %in% names(native)) native[[col]] else rep(NA, n)
  singular_out <- if ("singular" %in% names(native)) native$singular
                  else rep(if (is.null(attr(native, "singular"))) FALSE
                           else as.logical(attr(native, "singular"))[1], n)
  ## a failed lme4::isSingular() is reported as singular = NA, never FALSE
  sing_check <- if (!is.null(model)) .singular_status(model)$check else NA_character_
  if (!is.na(sing_check)) singular_out <- rep(NA, n)

  core <- data.frame(
    type = .type_tag(kind), estimate = native$estimate,
    lower = native$lower, upper = native$upper,
    level = level_out, P = rep(if (kind %in% c("ti", "TI")) P_value else NA_real_, n),
    sides = if ("sides" %in% names(native)) native$sides else rep("two.sided", n),
    method = method_out, design = design_out,
    singular = singular_out, boundary = pick("boundary"),
    vg_floored = pick("vg_floored"), stringsAsFactors = FALSE
  )
  out <- cbind(eval_df, core); rownames(out) <- NULL

  ## diagnostics -> attribute (off the core shape, still retrievable)
  diag <- list()
  for (cn in c("df", "df_method", "df_rule", "ci_df_rule", "slope_ci_df_rule",
               "pi_df_rule", "ngm_arm", "slope_ngm_arm", "ngm_df_failed", "W", "se", "tau2_med", "s2e_med",
               "frac_boundary", "K", "df_resid", "phi1", "c_resid",
               "route_reason", "mls_bound", "n_na_draws", "mc_se_lower",
               "mc_se_upper")) {
    if (cn %in% names(native)) diag[[cn]] <- native[[cn]]
  }
  for (a in c("beta_hat", "n_E", "U", "Vhat_Y", "C1", "C2", "pi_df_fallback",
              "gpq_seed", "rng_kind", "mc_sections")) {
    v <- attr(native, a); if (!is.null(v)) diag[[a]] <- v
  }
  if (!is.na(sing_check)) diag$singular_check <- sing_check
  grp <- if (!is.null(model)) .grouping_of(model) else NA_character_
  if (!is.null(attr(grp, "reason"))) diag$grouping_reason <- attr(grp, "reason")
  attr(out, "diagnostics") <- diag
  attr(out, "grouping") <- as.vector(grp)
  class(out) <- c("lmm_interval", "data.frame")
  out
}

#' @export
print.lmm_interval <- function(x, ...) {
  kind <- if (length(unique(x$type)) == 1L) unique(x$type) else "mixed"
  cat(sprintf("<lmm_interval> %s  (method = %s, design = %s)\n",
              kind, paste(unique(x$method), collapse = "/"),
              paste(unique(x$design), collapse = "/")))
  core <- as.data.frame(x)[, !names(x) %in% c("design", "method"), drop = FALSE]
  print(core, row.names = FALSE, ...)
  if (any(stats::na.omit(x$singular)) || any(stats::na.omit(x$boundary)) ||
        any(stats::na.omit(x$vg_floored))) {
    cat("note: some rows carry a singular/boundary/vg_floored flag ",
        "(degenerate fit at that point); see the flag columns.\n", sep = "")
  }
  d <- attr(x, "diagnostics")
  if (length(d)) cat(sprintf("diagnostics: %s  (see attr(., \"diagnostics\"))\n",
                             paste(names(d), collapse = ", ")))
  invisible(x)
}

#' Confidence interval for the mean of a population, from a linear mixed model
#'
#' One of the four interval verbs ([ci_lmm()], [pi_lmm()], [ti_lmm()],
#' [new_group_mean_lmm()]). `ci_lmm()` gives a confidence interval for the mean
#' response at `newdata`. It dispatches by random structure: random-intercept
#' designs (including a fixed slope) use the expected-mean-square engine;
#' random-slope designs use a REML-based t-interval (`"reml"`). See [reintervals] and
#' `vignette("choosing-an-interval")` for the engine map, [reintervals-models]
#' for the design tags, and [lmm_predict()] for the unified entry point.
#'
#' @param model A fitted `lmerMod` (from [lme4::lmer()]).
#' @param newdata Row(s) giving the fixed-effect combination(s) to evaluate at.
#'   `NULL` uses the first model-frame row.
#' @param alpha,level The interval confidence, as a mutually-exclusive pair: give
#'   at most one. `level` is `1 - alpha`; the default is `alpha = 0.05` (95%).
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`. A one-sided
#'   interval keeps the two-sided interval's center, variance, and df and uses
#'   `t_{1 - alpha}` in place of `t_{1 - alpha/2}`.
#' @param method `"auto"` (default) chooses the algorithm by structure; one of
#'   `"ems-mls"`, `"anova-mls"` (new-group mean only), `"reml"` (random-slope
#'   models) forces that engine and errors if it is incompatible with the model.
#' @param control An [lmm_interval_control()] object (the GPQ `seed` and `M`,
#'   the df rules `ci_df`, `slope_ci_df` and `pi_df`, and the new-group-mean
#'   arms `ngm_method` and `slope_ngm_method`; an option that does not apply to
#'   the model or verb is ignored and recorded as `"not_applicable"` in the
#'   diagnostics).
#'
#' @return An `lmm_interval` object (a data frame), one row per evaluation point,
#'   with columns, in order:
#'   \itemize{
#'     \item the eval-point column(s), named by the actual model predictor(s)
#'       (e.g. the time covariate for a slope model, or each `newdata` column);
#'     \item `type` (`"ci"`, `"pi"`, `"ti"`, or `"new_group_mean"`), `estimate`,
#'       `lower`, `upper`;
#'     \item `level` (the confidence, `1 - alpha`); `P` (content proportion ---
#'       the value for `ti`, `NA` for the coverage intervals); `sides`;
#'     \item `method` (the engine: `"ems-mls"`, `"anova-mls"`, `"reml"`, or
#'       `"gpq"`); `design` (the route-independent model class, see
#'       [reintervals-models]);
#'     \item the per-row flags `singular`, `boundary`, `vg_floored`.
#'   }
#'   `attr(., "diagnostics")` is a list of engine-specific extras (degrees of
#'   freedom `df`; for the random-intercept confidence interval its source
#'   `df_method` --- `"satterthwaite"`, `"containment"`, or `"kenward-roger"`;
#'   for the random-slope tolerance interval the GPQ `route_reason`; standard
#'   error `se`; and any GPQ / ANOVA intermediates);
#'   `attr(., "grouping")` is the random-intercept topology (`"oneway"`,
#'   `"nested"`, `"crossed"`, or `NA` for slope models).
#'
#' @section Known limitation: 
#' On a fixed-slope model with single-observation lots (the release-only
#' stability structure), the EMS confidence and prediction intervals can be
#' anti-conservative (too narrow), because REML precision-weighting
#' down-weights the single-observation lots. [ci_lmm()] and [pi_lmm()] emit a
#' warning when at least one lot has a single observation (the interval is
#' still returned, unchanged); staggered designs whose lots all have two or more
#' observations do not warn. The tolerance interval ([ti_lmm()]) and the
#' new-group-mean interval ([new_group_mean_lmm()]) use the unweighted lot-mean
#' (Montes) closed forms on single-observation and staggered fixed-slope
#' designs and do not warn; no validated alternative is currently available for
#' the confidence or prediction interval.
#'
#' @seealso [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' ci_lmm(fm, newdata = data.frame(Days = 0))
#' @export
ci_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "CI", level = lvl, sides = sides,
                 method = method, control = control)
}

#' Prediction interval for a future observation of a linear mixed model
#'
#' One of the four interval verbs ([ci_lmm()], [pi_lmm()], [ti_lmm()],
#' [new_group_mean_lmm()]). `pi_lmm()` gives a prediction interval for one future
#' observation at `newdata`. Dispatch and arguments are as for [ci_lmm()].
#'
#' @inheritParams ci_lmm
#' @inherit ci_lmm return
#' @inheritSection ci_lmm Known limitation
#' @section Known limitation (few top-level groups):
#' With few levels of the top random factor (five or fewer) and a dominant
#' between-group variance (an intraclass correlation near 0.9), the Francq et
#' al. (2019) prediction interval used for random-intercept models can
#' under-cover: the worst coverage observed in the package's validation was
#' about 0.93 for a nominal 0.95, consistent with Francq et al. (2019, Fig. 2).
#' When the top random factor of a random-intercept model has five or fewer
#' levels, `pi_lmm()` signals a warning of class
#' `reintervals_warning_limitation` (the interval itself is unchanged). The
#' warning is given at most once per R session.
#' @seealso [ci_lmm()], [ti_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' pi_lmm(fm, newdata = data.frame(Days = 0))
#' @export
pi_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "PI", level = lvl, sides = sides,
                 method = method, control = control)
}

#' New-group-mean (CInew / new-batch-mean) interval
#'
#' One of the four interval verbs: an interval for the mean of one new,
#' as-yet-unobserved group/batch (the CInew / new-batch-mean interval). Its
#' coverage is `level`/`alpha`, so it carries no content `P` (unlike
#' [ti_lmm()]).
#'
#' @section Target:
#' The target is the **true** mean of a new innermost cell: a new group in a
#' one-way design, a new B-within-A cell in a nested design, a new (a, b) cell
#' in a crossed design; for the fixed-slope and random-slope designs, the true
#' mean of a new lot at `t0`.
#'
#' @section Constructions:
#' The interval is `center +/- t_df * sqrt(W)`, with `W`, `df` and the rule
#' reported in `attr(., "diagnostics")` (`W`, `df`, `df_rule`):
#' \itemize{
#'   \item fixed slope (one grouping factor, any balance): the unweighted
#'     lot-mean (Montes) center `Yhat(t0)` with
#'     `W = S2_L (1 + 1/I) + S2_E (t0 - tbar)^2 / S_tt` and a Satterthwaite df
#'     on the two mean squares (`df_rule = "anova-satterthwaite"`, method
#'     `"anova-mls"`);
#'   \item intercept-only one-way: the unweighted mean of group means with
#'     `W = S2_L (1 + 1/I)` and `df = I - 1` (`"anova-single"`,
#'     `"anova-mls"`);
#'   \item nested, and one-way with other fixed effects:
#'     `W = V_F + sum_j k_j EMS_j` over the random terms (REML-implied expected
#'     mean squares, synthesis coefficients when unbalanced) with the top
#'     factor's df `A - 1` (`"top-factor"`, `"ems-mls"`);
#'   \item crossed with interaction: the same `W`, with a Satterthwaite df in
#'     which `V_F` is attached to the two main-effect terms
#'     (`"satterthwaite-vf"`, `"ems-mls"`);
#'   \item random slope: the REML-based t-interval
#'     \eqn{\mu(t_0) \pm t_{\nu_G}\sqrt{V_F + V_G}} (`"reml"`), with
#'     \eqn{\nu_G} the Satterthwaite df of \eqn{V_G} alone (SAS DDFM=SATTERTH
#'     uses the df of \eqn{V_F + V_G}, which is larger; that df is the
#'     validation arm `slope_ngm_method = "broad-inference"`).
#' }
#' A design without a closed-form variance decomposition (two crossed factors
#' without their interaction term, or four or more grouping factors) is an
#' error.
#' One-sided intervals (`sides = "lower"` / `"upper"`) use the same center, `W`,
#' and df with `t_{1 - alpha}`.
#'
#' @section Validation arms:
#' [lmm_interval_control()] selects alternative constructions for the
#' validation study: `ngm_method = "broad-inference"`, `"hts"` or
#' `"classical"` on random-intercept models and `slope_ngm_method =
#' "broad-inference"`, `"gpq"` or `"fv-anova"` on random-slope models (see there for the formulas). They warn
#' (`reintervals_warning_unvalidated`), and the arm used is recorded as
#' `attr(., "diagnostics")$ngm_arm` / `$slope_ngm_arm`. The `"gpq"` arm reports
#' `df = NA`, the Monte-Carlo SE of each limit (`mc_se_lower`, `mc_se_upper`)
#' and `gpq_seed`; `"broad-inference"` (both models) reports `ngm_df_failed`.
#'
#' @references
#' McLean RA, Sanders WL, Stroup WW (1991). A unified approach to mixed linear
#' models. *The American Statistician* 45(1):54-64.
#'
#' Higgins JPT, Thompson SG, Spiegelhalter DJ (2009). A re-evaluation of
#' random-effects meta-analysis. *Journal of the Royal Statistical Society A*
#' 172(1):137-159.
#'
#' Hahn GJ, Meeker WQ (1991). *Statistical Intervals: A Guide for
#' Practitioners*. Wiley.
#'
#' Krishnamoorthy K, Mathew T (2009). *Statistical Tolerance Regions: Theory,
#' Applications, and Computation*. Wiley.
#'
#' Oliva-Aviles C, Hauser P (2025). *Technometrics* 67(2):193-202.
#' \doi{10.1080/00401706.2024.2407324}
#'
#' @inheritParams ci_lmm
#' @inherit ci_lmm return
#' @seealso [ci_lmm()], [pi_lmm()], [ti_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
#'                  data = lme4::sleepstudy)
#' new_group_mean_lmm(fm, newdata = data.frame(Days = 5))
#' @export
new_group_mean_lmm <- function(model, newdata = NULL, alpha = NULL, level = NULL,
                               sides = c("two.sided", "lower", "upper"),
                               method = "auto",
                               control = lmm_interval_control()) {
  .assert_lmm(model)
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal coverage `level` (1 - alpha).
  lvl <- 1 - .resolve_alpha(alpha, level)
  .coverage_verb(model, newdata, kind = "new_group_mean", level = lvl,
                 sides = sides, method = method, control = control)
}

## shared body for the three coverage verbs: dispatch via the engine entry
## (.interval_lmm) then stamp harmonized provenance. method tags: slope ->
## "reml" (already on the row; "reml-mls" until 1.0.0); EMS -> "ems-mls".
.coverage_verb <- function(model, newdata, kind, level, sides, method, control) {
  if (identical(method, "reml-mls") && .has_random_slope(model)) {
    .abort("method = \"reml-mls\" is no longer a method of the random-slope ",
           "confidence, prediction and new-group-mean intervals: since 1.0.0 their ",
           "method is \"reml\" (\"reml-mls\" named the removed closed-form ",
           "tolerance interval). Use \"auto\" or \"reml\".",
           class = "reintervals_error_method")
  }
  if (!identical(method, "auto")) {
    ok <- if (.has_random_slope(model)) "reml"
          else if (kind == "new_group_mean" && .ngm_uses_anova(model)) "anova-mls"
          else "ems-mls"
    if (!method %in% ok) {
      .abort(sprintf(
        "method = \"%s\" is not available for a %s model on this interval; ",
        method, .model_class(model)),
             sprintf("use \"auto\" or \"%s\".", ok), class = "reintervals_error_method")
    }
  }
  ## slope models: newdata must carry a usable time column; warn (once) when a
  ## requested t0 lies outside the observed time range
  .slope_newdata_guard(model, newdata)
  ## Anti-conservative-drift warning, scoped to the drifting cells: the EMS CI
  ## and PI on a fixed-slope model with at least one single-observation lot
  ## (.m1_has_single_obs_lot; REML precision-weighting down-weights those lots,
  ## the same effect that pushes the EMS TI ~2.4% narrow vs the Montes oracle).
  ## Staggered designs whose lots all have >= 2 observations showed no shortfall
  ## in validation and do not warn. We make the drift VISIBLE -- we do NOT
  ## re-route. The new-group-mean interval on this structure uses the unweighted
  ## lot-mean construction (Fa, .ngm_ri) and does not warn.
  if (kind %in% c("CI", "PI") && .m1_has_single_obs_lot(model)) {
    .warn("This is a fixed-slope model with single-observation lots; the ",
          "EMS ", kind, " interval may be anti-conservative (too ",
          "narrow) on such data, because REML precision-weighting down-weights ",
          "the single-observation lots. The tolerance and new-group-mean ",
          "intervals use the unweighted lot-mean (Montes) closed forms on this ",
          "structure; no validated alternative is currently available for the ",
          "confidence or prediction interval here.",
          class = "reintervals_warning_limitation")
  }
  ## Known limitation (few top-level groups), documented in ?pi_lmm since
  ## 0.2.0: the Francq PI on a random-intercept model whose top (coarsest)
  ## random factor has five or fewer levels. The documented condition also
  ## names a dominant between-group variance (ICC near 0.9), but no ICC
  ## threshold is documented, so only the level count is checked. Emitted at
  ## most once per R session (.once, conditions.R).
  if (kind == "PI" && !.has_random_slope(model)) {
    n_top <- min(vapply(lme4::getME(model, "flist"), nlevels, numeric(1)))
    if (n_top <= 5L && .once("few_groups_pi")) {
      .warn(sprintf("The top random factor has %d levels (five or fewer). ", n_top),
            "With few levels of the top random factor (five or fewer) and a ",
            "dominant between-group variance (an intraclass correlation near ",
            "0.9), the Francq et al. (2019) prediction interval used for ",
            "random-intercept models can under-cover: the worst coverage observed ",
            "in the package's validation was about 0.93 for a nominal 0.95, ",
            "consistent with Francq et al. (2019, Fig. 2).",
            class = "reintervals_warning_limitation",
            data = list(n_top_levels = n_top))
    }
  }
  applied <- .applied_arms(model, kind, control)
  .warn_unvalidated(applied)
  out <- .interval_lmm(model, newdata = newdata, kind = kind, level = level,
                       sides = sides, control = control, applied = applied)
  .warn_df_fallback(out, model, kind, applied)
  ## the df arm of every option: the value applied here, or "not_applicable"
  for (r in names(rules <- .arm_rules(applied))) out[[r]] <- rules[[r]]
  tag <- if (.has_random_slope(model)) NULL else "ems-mls"   # slope row already tagged
  .unify_intervals(out, model, kind = kind, method_tag = tag)
}

#' Unified interface to the four interval verbs
#'
#' A predict-style wrapper over the four interval verbs ([ci_lmm()], [pi_lmm()],
#' [ti_lmm()], [new_group_mean_lmm()]), selected by `interval`. This is a plain
#' function, not an S3 method --- it does not override `lme4`'s
#' `predict.merMod`.
#'
#' @inheritParams ti_lmm
#' @param interval Which interval: `"confidence"`, `"prediction"`, `"tolerance"`,
#'   or `"new_group_mean"`.
#' @param P,over Content proportion and tolerance target; honored only for
#'   `interval = "tolerance"` (forwarded to [ti_lmm()]). Supplying a non-default
#'   `P`/`over` with any other `interval` is an error.
#' @inherit ci_lmm return
#' @seealso [ci_lmm()], [pi_lmm()], [ti_lmm()], [new_group_mean_lmm()],
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' lmm_predict(fm, interval = "confidence", newdata = data.frame(Days = 0))
#' @export
lmm_predict <- function(model,
                        interval = c("confidence", "prediction", "tolerance",
                                     "new_group_mean"),
                        newdata = NULL, P = 0.95,
                        over = c("observation", "group_mean"),
                        alpha = NULL, level = NULL,
                        sides = c("two.sided", "lower", "upper"),
                        method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  interval <- match.arg(interval)
  sides    <- match.arg(sides)
  over     <- match.arg(over)
  if (interval != "tolerance" && (!isTRUE(all.equal(P, 0.95)) ||
        over != "observation")) {
    .abort("`P` and `over` apply to interval = \"tolerance\" only.",
           class = "reintervals_error_input")
  }
  switch(interval,
    confidence     = ci_lmm(model, newdata, alpha, level, sides, method, control),
    prediction     = pi_lmm(model, newdata, alpha, level, sides, method, control),
    new_group_mean = new_group_mean_lmm(model, newdata, alpha, level, sides,
                                        method, control),
    tolerance      = ti_lmm(model, newdata, P = P, over = over, alpha = alpha,
                            level = level, sides = sides, method = method,
                            control = control)
  )
}
