# =============================================================================
# Top-level tolerance-interval dispatch. The classifier routes from the fitted
# model's random structure, its fixed-effect structure, the requested target,
# and the flag:
#
#   * any random slope present -> GPQ (Oliva-Aviles & Hauser 2025), for both
#     targets (future observation, batch mean), every side, singular or not.
#     The REML closed-form MLS interval (`reml-mls`) was removed in 1.0.0: it
#     under-covers (one-sided 0.916-0.934 when the residual variance dominates
#     and t0 is extrapolated; two-sided 0.84-0.93 where forced), and the
#     phase-0 study validated one-sided GPQ at M = 20000 in the regime that had
#     motivated the closed-form default (outcome A, simstudy d20b954).
#
#   * random intercept + a single fixed slope:
#       - balanced -> EMS engine: exact, paper-anchored.
#       - single-observation / release-only lots, or staggered / unbalanced
#         entry -> the Montes ANOVA/MLS closed form. The EMS path drifts ~2.4%
#         anti-conservative against the Montes oracle on single-obs lots (REML
#         precision-weighting down-weights them); the output therefore changes
#         discontinuously at this boundary by design.
#
#   * intercept-only fixed part with exactly ONE random-intercept factor (one-way)
#     -> the Hoffman-Kringle (2005) one-way interval, i.e. the slope-free branch
#     of the anova-mls core, balanced or not. method = "ems-mls" (Francq 2019)
#     stays available; the two coincide on balanced non-singular data.
#
#   * any other random-intercept model, no fixed slope (one-way with covariates /
#     nested / crossed) -> EMS engine.
#
# Each engine keeps its own singular guard (they do different jobs; not unified):
#   * EMS path  -> reintervals' Hessian -> stratum containment df guard (the
#     PI-width guard pinned by test-singular-guard.R).
#   * slope path -> the tolerance interval is GPQ, which has no df to fall back
#     on; the slope CI / PI / new-group mean keep the guard carried by
#     vfun_extract (its own Satterthwaite -> containment df fallback +
#     singular/boundary flags).
#
# The df fallback is also not shared: the EMS fallback (.ddf_containment,
# vc-extract.R) and the slope fallback (.vfun_df_containment, vfun-extract.R)
# are distinct functions validated against distinct oracles.
#
# The EMS random-slope stop() in vc-extract.R (.assert_random_intercepts) is the
# guard for anyone calling vc_extract / the EMS engine directly on a random-slope
# model; this dispatcher routes random-slope models to the slope path before they
# reach it.
# =============================================================================

## TRUE iff any random-effect term is a non-intercept (a random slope). Uses the
## fitted model's component-name structure, not a guess. Random-intercept
## designs (one-way / nested / crossed intercepts) are all FALSE -> EMS.
.has_random_slope <- function(model) {
  cnms <- lme4::getME(model, "cnms")
  any(vapply(cnms, function(terms) any(terms != "(Intercept)"), logical(1)))
}

## Slope-design tag from the random structure (one random-effect term ->
## correlated random slope; two terms -> independent). Labels slope-path rows.
.slope_design <- function(model) {
  if (length(lme4::getME(model, "cnms")) == 1L) "random_slope_correlated"
  else "random_slope_independent"
}

## Raw (y, time, batch) pulled off the fit, for the engines that take vectors
## (.gpq_core, .montes_core) rather than a re-parsed fit.
.model_yts <- function(model, time_name) {
  group_name <- names(lme4::getME(model, "cnms"))[1]
  mf <- stats::model.frame(model)
  list(
    y     = as.numeric(lme4::getME(model, "y")),
    time  = as.numeric(mf[[time_name]]),
    batch = factor(mf[[group_name]])
  )
}

## Fixed-slope "clean balanced" predicate. A design is clean
## balanced iff every lot has >= 2 observations AND all lots share the same
## multiset of time points. Single-observation (release-only) lots, or staggered
## / unbalanced entry, make it FALSE -> route to .montes_core.
.m1_is_clean_balanced <- function(model, time_name) {
  group_name <- names(lme4::getME(model, "cnms"))[1]
  mf <- stats::model.frame(model)
  g  <- factor(mf[[group_name]]); tt <- as.numeric(mf[[time_name]])
  if (any(table(g) < 2L)) return(FALSE)                  # single-obs / release-only
  tsets <- tapply(tt, g, function(v) paste(sort(v), collapse = ","))
  length(unique(tsets)) == 1L                            # identical time grids
}

## The fixed-slope (M1) time column, or NA if the fixed part is not a single
## numeric time covariate. Maps the one non-intercept fixed TERM to its
## model-frame column via the terms object (not the coefficient name, which for
## a factor / logical is e.g. "grpB" / "flagTRUE" and names no column), and
## requires that column to be numeric. A 2-level factor or a logical covariate
## also yields one non-intercept coefficient but is NOT a time design -> NA,
## and the model follows the EMS path. The column must also VARY WITHIN at
## least one level of a grouping factor: a covariate constant within every
## level (a between-group covariate, e.g. y ~ z + (1 | g)) carries no
## within-lot slope information, is not a time design -> NA, EMS path.
.m1_time_name <- function(model) {
  covars <- setdiff(names(lme4::fixef(model)), "(Intercept)")
  if (length(covars) != 1L) return(NA_character_)
  labs <- attr(stats::terms(model, fixed.only = TRUE), "term.labels")
  if (length(labs) != 1L) return(NA_character_)
  col <- stats::model.frame(model)[[labs]]
  if (is.null(col) || !is.numeric(col) || NCOL(col) != 1L) return(NA_character_)
  varies_within <- any(vapply(lme4::getME(model, "flist"), function(g)
    any(tapply(col, g, function(v) length(unique(v)) > 1L), na.rm = TRUE), logical(1)))
  if (!varies_within) return(NA_character_)
  labs
}

## TI ROUTING predicate: a fixed-slope model (no random slope, one numeric
## covariate varying within lots; see .m1_time_name) that is NOT clean balanced
## (single-observation lots, or staggered / unequal time grids). ti_lmm routes
## its tolerance interval to the Montes ANOVA closed form (anova-mls).
.m1_not_clean_balanced <- function(model) {
  if (.has_random_slope(model)) return(FALSE)
  time_name <- .m1_time_name(model)
  if (is.na(time_name)) return(FALSE)
  !.m1_is_clean_balanced(model, time_name)
}

## CI / PI WARNING predicate: a fixed-slope model with AT LEAST ONE
## single-observation lot. Narrower than the routing predicate: the validation
## study found the EMS CI / PI mildly anti-conservative on single-observation
## designs (E4) but not on staggered designs whose lots all have >= 2
## observations (V2, worst cell 0.943). Separate from the routing predicate on
## purpose: the TI stays on anova-mls for every non-clean-balanced design.
.m1_has_single_obs_lot <- function(model) {
  if (.has_random_slope(model)) return(FALSE)
  if (is.na(.m1_time_name(model))) return(FALSE)
  any(table(lme4::getME(model, "flist")[[1]]) < 2L)
}

## The Hoffman-Kringle one-way detector: no random slope, an intercept-only
## fixed part (no covariates), and exactly one random-intercept grouping factor.
.is_oneway_intercept_only <- function(model) {
  if (.has_random_slope(model)) return(FALSE)
  identical(names(lme4::fixef(model)), "(Intercept)") &&
    length(lme4::getME(model, "cnms")) == 1L
}

## Normalize .gpq_core output to the slope-path row schema (the cross-engine
## output contract is handled at the unify step). `boundary` is NA: the GPQ
## method has no fitted-boundary notion (frac_boundary is an internal MC
## diagnostic, reported in the diagnostics).
.gpq_rows <- function(g, design, sides, singular) {
  out <- data.frame(
    time = g$time, type = "TI", estimate = g$estimate,
    lower = g$lower, upper = g$upper,
    content = g$content, conf = g$conf, df = g$df_resid, sides = sides,
    design = design, method = "gpq",
    singular = singular, boundary = NA,
    ## GPQ audit fields carried through (-> attr(.,"diagnostics") at unify)
    tau2_med = g$tau2_med, s2e_med = g$s2e_med, frac_boundary = g$frac_boundary,
    n_na_draws = g$n_na_draws,
    mc_se_lower = g$mc_se_lower, mc_se_upper = g$mc_se_upper,
    K = g$K, df_resid = g$df_resid, phi1 = g$phi1, c_resid = g$c_resid,
    stringsAsFactors = FALSE
  )
  ## per-call (not per-row) Monte-Carlo provenance
  for (a in c("gpq_seed", "rng_kind", "mc_sections")) attr(out, a) <- attr(g, a)
  out
}

## Random-slope tolerance interval: GPQ at every requested time point, for
## both targets and all sides (one .gpq_core call over all t0). `forced` is TRUE
## when the user passed method = "gpq"; it only sets route_reason.
##
## route_reason now has two values: "auto" (the dispatcher chose GPQ, which it
## always does for a random slope) and "method" (method = "gpq" was passed).
## The 0.2.0 values ("none", "singular", "df_unavailable", "group_mean_target")
## described why the closed form was or was not used; with the closed form
## gone there is no per-t0 or per-fit choice left to report. The field is kept
## so code reading attr(., "diagnostics")$route_reason still finds it.
.ti_slope_path <- function(model, newdata, content, conf, target, sides,
                           forced, seed, M) {
  fe <- lme4::fixef(model)
  time_var <- setdiff(names(fe), "(Intercept)")
  if (length(time_var) != 1L) {
    .abort("random-slope intervals require a single continuous covariate (e.g. ~ 1 + time); ",
           "multiple fixed effects with a random slope are not currently supported.",
           class = "reintervals_error_unsupported")
  }
  t0s <- if (is.null(newdata)) stats::model.frame(model)[[time_var]][1] else
    newdata[[time_var]]
  if (is.null(t0s) || !length(t0s)) {
    .abort("Could not determine the evaluation time point(s) for this interval; ",
           "supply `newdata` with a numeric '", time_var,
           "' column giving the time(s) to evaluate at.", class = "reintervals_error_input")
  }

  singular <- isTRUE(tryCatch(lme4::isSingular(model), error = function(e) NA))
  yts <- .model_yts(model, time_var)
  c_resid <- if (target == "observable") 1 else 0       # future-obs vs batch-mean
  g <- .gpq_core(yts$y, yts$time, yts$batch, t0 = t0s,
                 content = content, conf = conf, sides = sides,
                 c_resid = c_resid, M = M, seed = seed)
  out <- .gpq_rows(g, design = .slope_design(model), sides = sides,
                   singular = singular)
  rownames(out) <- NULL
  out$route_reason <- if (forced) "method" else "auto"
  out
}

## method = "reml-mls" (the random-slope REML closed-form MLS tolerance
## interval) was removed in 1.0.0. The coverage figures are from the decision
## memo (v4), random-slope closed form, observation target.
.reml_mls_defunct <- function() {
  .abort("method = \"reml-mls\" was removed in reintervals 1.0.0: the REML ",
         "closed-form MLS tolerance interval for random-slope models under-covers ",
         "(one-sided 0.916-0.934 when the residual variance dominates and t0 is ",
         "extrapolated; two-sided 0.84-0.93 where forced). Use method = \"gpq\" ",
         "(the default; Oliva-Aviles & Hauser 2025).",
         class = "reintervals_error_defunct")
}

## Fixed-slope single-obs / staggered route -> the Montes ANOVA/MLS closed form.
.ti_m1_montes_path <- function(model, newdata, content, conf, sides, time_name) {
  yts <- .model_yts(model, time_name)
  t0s <- if (is.null(newdata)) yts$time[1] else newdata[[time_name]]
  if (is.null(t0s) || !length(t0s)) {
    .abort("Could not determine the evaluation time point(s) for this interval; ",
           "supply `newdata` with a numeric '", time_name,
           "' column giving the time(s) to evaluate at.", class = "reintervals_error_input")
  }
  .montes_core(yts$y, yts$time, yts$batch, t0 = t0s,
               content = content, conf = conf, sides = sides)
}

#' Tolerance interval for a linear mixed model
#'
#' One of the four interval verbs, and the dispatched tolerance interval. It
#' classifies the fitted model and routes by random/fixed structure and the
#' tolerance target (`over`) to the engine that is exact/validated for that
#' case --- all transparently:
#'
#' * **Random slope.** The generalized-pivotal-quantity (GPQ) interval of
#'   Oliva-Aviles & Hauser (2025), for both targets (`over = "observation"` and
#'   `"group_mean"`), two-sided and one-sided, on every fit (singular included).
#'   It is a Monte-Carlo interval: set `lmm_interval_control(seed = )` for
#'   reproducibility (see Details); the default `M = 20000` is the size
#'   validated by simulation. `attr(., "diagnostics")$route_reason` is `"auto"`, or
#'   `"method"` when `method = "gpq"` was passed. The REML closed-form interval
#'   (`method = "reml-mls"`) was removed in 1.0.0 because it under-covers; it
#'   now signals an error of class `reintervals_error_defunct`.
#' * **Random intercept + fixed slope.** Clean balanced -> EMS engine; single-
#'   observation / staggered / unbalanced -> the Montes ANOVA/MLS closed form.
#' * **Intercept-only one-way** (`y ~ 1 + (1 | g)`, one grouping factor, no
#'   covariates), balanced or not -> the Hoffman & Kringle (2005) closed form
#'   (`"anova-mls"`). `method = "ems-mls"` gives the Francq et
#'   al. (2019) interval instead; the two coincide on balanced data with a
#'   non-singular fit and differ on singular or unbalanced data.
#' * **Any other random-intercept model, no fixed slope** (one-way with
#'   covariates, nested, crossed) -> EMS engine.
#'
#' For slope models the between-group variance `V_G(t0)` is quadratic in the
#' time covariate, so evaluating at a `t0` outside the observed time range is
#' an extrapolation of that quadratic and should be done knowingly (the interval
#' is computed without a guard or warning). The same applies to the other slope
#' verbs ([ci_lmm()], [pi_lmm()], [new_group_mean_lmm()]).
#'
#' @details
#' **Monte-Carlo (GPQ) intervals.** Within one call, the GPQ pivot draws are
#' generated once and reused for every time point in `newdata` and for both
#' limits (common random numbers): for a given seed and `M`, the limits at a
#' time point are identical whatever other time points are requested, and in
#' any order. The seed used is in `attr(., "diagnostics")$gpq_seed` and the
#' RNG kind in `$rng_kind`; passing
#' `control = lmm_interval_control(seed = <gpq_seed>, M = <same M>)` under the
#' same RNG kind reproduces the interval exactly (with `seed = NULL`, one seed
#' is drawn from the caller's RNG stream, so `set.seed()` before the call also
#' reproduces it; see [lmm_interval_control()]). The diagnostics also carry the
#' Monte-Carlo standard error of each limit, `mc_se_lower` and `mc_se_upper`
#' (`NA` on the infinite side of a one-sided interval), computed by MC
#' sections (`mc_sections = 20` consecutive groups of draws; the standard
#' error is the standard deviation of the 20 section limits over
#' `sqrt(20)`). The MC SE is the run-to-run jitter of the Monte-Carlo
#' computation, not statistical uncertainty: a limit within approximately 2 MC
#' SE of a specification limit has a seed-dependent pass/fail. In a sleepstudy
#' check at `M = 2000` the sectioning estimate ran about 8% below the observed
#' run-to-run SD (sections of 100 draws); the gap is expected to shrink at the
#' default `M = 20000`. The GPQ route requires
#' every lot to have at least two distinct time points and otherwise signals
#' an error of class `reintervals_error_time_points`.
#'
#' @section Availability:
#' Tolerance intervals require a closed-form variance decomposition. They are
#' available for one-way, nested, and crossed-with-interaction random-intercept
#' designs, the fixed-slope stability design, and random-slope models. For a
#' design without such a decomposition (two crossed factors with no interaction
#' term, or four or more grouping factors), [ti_lmm()] and
#' [new_group_mean_lmm()] signal an error; [ci_lmm()] and [pi_lmm()] remain
#' available for those models.
#'
#' @param model A fitted `lmerMod`.
#' @param newdata Row(s) giving the fixed-effect combination(s); the time column
#'   supplies the evaluation point(s). `NULL` uses the first model-frame row.
#' @param P Content proportion (the population fraction covered). Default 0.95
#'   (pharma stability work typically passes `P = 0.99`).
#' @param over What the interval covers: `"observation"` (default; a proportion
#'   `P` of future observations --- the ordinary tolerance interval) or
#'   `"group_mean"` (a proportion `P` of new group means --- the batch-mean
#'   tolerance interval).
#' @param alpha,level The confidence, as a mutually-exclusive pair: give at most
#'   one. `level` is `1 - alpha`; default `alpha = 0.05`.
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`. One-sided
#'   tolerance limits: `anova-mls` uses Hoffman's (2010) one-sided limit;
#'   `ems-mls` uses the package's extension of it, validated by simulation;
#'   random-slope models use the one-sided GPQ limit. See [reintervals-models].
#' @param method `"auto"` (default) chooses the algorithm by structure; one of
#'   `"ems-mls"`, `"anova-mls"`, `"gpq"` forces it and errors if incompatible
#'   with the model. `"anova-mls"` is available for fixed-slope models and for
#'   intercept-only one-way models (Hoffman-Kringle); `"gpq"` for random-slope
#'   models. `"reml-mls"` is defunct (removed in 1.0.0) and errors.
#' @param control [lmm_interval_control()] tuning (GPQ seed / M).
#'
#' @inherit ci_lmm return
#' @seealso [ci_lmm()], [pi_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @references Oliva-Aviles C, Hauser P (2025). *Technometrics* 67(2):193-202.
#'   \doi{10.1080/00401706.2024.2407324}
#'
#'   Hoffman D (2010). One-sided tolerance limits for balanced and
#'   unbalanced random effects models. *Technometrics* 52(3):303-312.
#'   \doi{10.1198/TECH.2010.09171}
#'
#'   Hoffman D, Kringle R (2005). Two-sided tolerance intervals for
#'   balanced and unbalanced random effects models. *Journal of
#'   Biopharmaceutical Statistics* 15(2):283-293. \doi{10.1081/BIP-200048826}
#' @examples
#' # intercept-only one-way -> Hoffman-Kringle (anova-mls)
#' fh <- lme4::lmer(conc ~ 1 + (1 | run), data = hk_assay)
#' ti_lmm(fh, P = 0.95, level = 0.90)
#' # random intercept + fixed slope (clean balanced) -> EMS engine
#' fi <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' ti_lmm(fi, newdata = data.frame(Days = 0))
#' # random slope -> GPQ (both targets). M = 2000 keeps the example fast; the
#' # default M = 20000 is the validated size.
#' fs <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy)
#' ti_lmm(fs, newdata = data.frame(Days = 5),
#'        control = lmm_interval_control(seed = 1, M = 2000))
#' ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean",
#'        control = lmm_interval_control(seed = 1, M = 2000))
#' @export
ti_lmm <- function(model, newdata = NULL, P = 0.95,
                   over = c("observation", "group_mean"),
                   alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  over   <- match.arg(over)
  sides  <- match.arg(sides)
  ## "reml-mls" stays in the list only so the defunct error below can fire
  method <- match.arg(method, c("auto", "ems-mls", "reml-mls", "anova-mls", "gpq"))
  if (method == "reml-mls") .reml_mls_defunct()
  ## boundary: public alpha/level -> internal `conf` (gamma); public `P` is the
  ## content proportion passed straight to the engines as `content`.
  conf   <- 1 - .resolve_alpha(alpha, level)
  if (!is.numeric(P) || length(P) != 1L || is.na(P) || P <= 0 || P >= 1) {
    .abort("`P` must be a single number strictly between 0 and 1.",
           class = "reintervals_error_input")
  }
  seed <- control$seed
  M    <- control$M
  target <- if (over == "observation") "observable" else "true_value"

  ## ---- random slope present -> GPQ (both targets, all sides) -----------------
  if (.has_random_slope(model)) {
    if (!method %in% c("auto", "gpq")) {
      .abort(sprintf("method = \"%s\" is not available for a %s model; use ",
                     method, .model_class(model)),
             "\"auto\" or \"gpq\".", class = "reintervals_error_method")
    }
    out <- .ti_slope_path(model, newdata = newdata, content = P, conf = conf,
                          target = target, sides = sides,
                          forced = method == "gpq", seed = seed, M = M)
    return(.unify_intervals(out, model, kind = "ti", P_value = P))
  }

  ## ---- random intercept only -------------------------------------------------
  time_name <- .m1_time_name(model)
  is_m1 <- !is.na(time_name)

  ## batch-mean tolerance is GPQ-only (random slope); DNE for random intercept / fixed-slope
  if (over == "group_mean") {
    avail <- if (is_m1) {
      "ti_lmm(over = \"observation\") or new_group_mean_lmm() (the prediction interval)"
    } else {
      "ti_lmm(over = \"observation\")"
    }
    .abort(sprintf(paste0("The batch-mean tolerance interval (over = ",
                          "\"group_mean\") is not implemented for a %s model; available here: %s."),
                   .model_class(model), avail), class = "reintervals_error_unsupported")
  }

  if (method == "gpq") {
    .abort(sprintf("method = \"%s\" is not available for a %s model; use ",
                   method, .model_class(model)),
           "\"auto\", \"ems-mls\", or \"anova-mls\".", class = "reintervals_error_method")
  }

  ## Intercept-only one-way -> Hoffman-Kringle (the slope-free branch of the
  ## anova-mls core), unless method = "ems-mls" asks for Francq.
  if (.is_oneway_intercept_only(model) && method %in% c("auto", "anova-mls")) {
    ## one-sided: Hoffman's (2010) published limit (see .montes_core)
    out <- .montes_core(as.numeric(lme4::getME(model, "y")), time = NULL,
                        lot = lme4::getME(model, "flist")[[1]], t0 = NULL,
                        content = P, conf = conf, sides = sides)
    out$singular <- .is_singular(model)   # report the fit's flag, like the other verbs
    return(.unify_intervals(out, model, kind = "ti", P_value = P))  # carries anova-mls
  }

  ## Fixed-slope single-obs / staggered -> Montes ANOVA; else EMS. `method` can
  ## force. (.m1_not_clean_balanced is the routing predicate; the coverage-verb
  ## warning uses the narrower .m1_has_single_obs_lot.)
  use_montes <- .m1_not_clean_balanced(model)
  if (method == "anova-mls") {
    if (!is_m1) {
      .abort("method = \"anova-mls\" needs a fixed-slope (random intercept + ",
             "fixed slope) model or an intercept-only one-way model.",
             class = "reintervals_error_method")
    }
    use_montes <- TRUE
  }
  if (method == "ems-mls") use_montes <- FALSE

  if (use_montes) {
    out <- .ti_m1_montes_path(model, newdata = newdata, content = P, conf = conf,
                              sides = sides, time_name = time_name)
    return(.unify_intervals(out, model, kind = "ti", P_value = P))  # carries anova-mls
  }

  ## one-sided EMS TI: the package's extension of Hoffman (2010), validated by
  ## simulation (see ti_francq)
  out <- .ti_ems(model, newdata = newdata, level = P, conf = conf, target = target,
                 sides = sides)
  .unify_intervals(out, model, kind = "ti", P_value = P, method_tag = "ems-mls")
}
