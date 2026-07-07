# =============================================================================
# Top-level tolerance-interval dispatch. The classifier routes from the fitted
# model's random structure, its fixed-effect structure, the requested target,
# and the singular flag:
#
#   * any random slope present -> SLOPE path, target-aware:
#       - future-observation: closed-form MLS (ti_vfun) when the fit is
#         non-singular; GPQ when singular.
#       - batch-mean: GPQ always. The closed form has a verified non-singular
#         breakdown for this target (~0.89 coverage at small tau^2), so it is
#         never used here, singular or not.
#     The `method` argument can pin the algorithm; the default is the
#     target-aware auto above.
#
#   * random intercept + a single fixed slope:
#       - clean balanced -> EMS engine: exact, paper-anchored.
#       - single-observation / release-only lots, or staggered / unbalanced
#         entry -> the Montes ANOVA/MLS closed form. The EMS path drifts ~2.4%
#         anti-conservative against the Montes oracle on single-obs lots (REML
#         precision-weighting down-weights them); the output therefore changes
#         discontinuously at this boundary by design.
#
#   * random intercept, no fixed slope (one-way / nested / crossed) -> EMS engine.
#
# Each engine keeps its OWN singular guard (they do different jobs; not unified):
#   * EMS path  -> reintervals' Hessian -> bounded-Satterthwaite guard (the
#     PI-width guard pinned by test-singular-guard.R).
#   * slope path -> the guard carried by vfun_extract (its own
#     Satterthwaite -> containment df fallback + singular/boundary flags) plus
#     the singular -> GPQ engine switch in .choose_slope_engine().
#
# The Satterthwaite df FALLBACK is also NOT shared: the EMS fallback
# (.satterthwaite_df, vc-extract.R) and the slope fallback (.vfun_df_containment,
# vfun-extract.R) are distinct functions validated against distinct oracles.
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

## The single detector for "fixed-slope model with single-observation /
## unbalanced lots" -- the structure on which (a) ti_lmm re-routes the TOLERANCE
## interval to the Montes ANOVA closed form, and (b) the EMS CI/PI/new-group-mean
## drift anti-conservative (REML precision-weighting). TRUE iff: no random slope,
## exactly one fixed covariate (a fixed-slope model), and NOT clean balanced.
## Shared by ti_lmm (routing) and the coverage verbs (warning) so the two can
## never diverge.
.is_single_obs_m1 <- function(model) {
  if (.has_random_slope(model)) return(FALSE)
  covars <- setdiff(names(lme4::fixef(model)), "(Intercept)")
  if (length(covars) != 1L) return(FALSE)
  !.m1_is_clean_balanced(model, covars)
}

## Target-aware slope-engine selection (the GPQ routing policy).
## Returns "closed" or "gpq"; errors only on the one unsafe forced combination.
## Pure function of (target, singular, engine) so the (target x singular) routing
## matrix is testable WITHOUT a fit -- the singular-vs-closed routing rule is
## exercised directly.
.choose_slope_engine <- function(target, singular, engine) {
  if (engine == "closed" && target == "true_value") {
    stop("The closed form has no validated batch-mean (true_value) variant; it ",
         "is precisely the route the dispatcher avoids for batch-mean (verified ",
         "non-singular breakdown). Use engine = \"gpq\" or ",
         "\"auto\".", call. = FALSE)
  }
  switch(engine,
    gpq    = "gpq",
    closed = "closed",
    auto   = if (target == "true_value") "gpq"           # batch-mean: GPQ always
             else if (singular)          "gpq"           # future-obs: GPQ if singular
             else                        "closed"        # future-obs: closed otherwise
  )
}

## Normalize .gpq_core output to the slope-path row schema (matches ti_vfun's
## columns so a dispatch call returns one uniform frame; the cross-engine output contract is handled at the unify step). `boundary` is NA: the GPQ method has no fitted-
## boundary notion (frac_boundary is an internal MC diagnostic, not reported).
.gpq_rows <- function(g, design, sides, singular) {
  data.frame(
    time = g$time, type = "TI", estimate = g$estimate,
    lower = g$lower, upper = g$upper,
    content = g$content, conf = g$conf, df = g$df_resid, sides = sides,
    design = design, method = "gpq",
    singular = singular, boundary = NA,
    ## GPQ audit fields carried through (-> attr(.,"diagnostics") at unify)
    tau2_med = g$tau2_med, s2e_med = g$s2e_med, frac_boundary = g$frac_boundary,
    K = g$K, df_resid = g$df_resid, phi1 = g$phi1, c_resid = g$c_resid,
    stringsAsFactors = FALSE
  )
}

## slope-engine route (random slope present): pick the engine target-aware, then
## run it at each requested time point. Keeps the slope engine's own singular/df
## handling; never touches the EMS guard.
.ti_slope_path <- function(model, newdata, content, conf, target, sides,
                           engine, seed, M) {
  fe <- lme4::fixef(model)
  time_var <- setdiff(names(fe), "(Intercept)")
  if (length(time_var) != 1L) {
    stop("random-slope intervals require a single time covariate (~ 1 + time); ",
         "multiple fixed effects with a random slope are a planned v2 extension.",
         call. = FALSE)
  }
  t0s <- if (is.null(newdata)) stats::model.frame(model)[[time_var]][1] else
    newdata[[time_var]]
  if (is.null(t0s) || !length(t0s)) {
    stop("Could not determine the evaluation time point(s) for this interval; ",
         "supply `newdata` with a numeric '", time_var,
         "' column giving the time(s) to evaluate at.", call. = FALSE)
  }

  singular <- isTRUE(tryCatch(lme4::isSingular(model), error = function(e) NA))
  chosen   <- .choose_slope_engine(target, singular, engine)
  design   <- .slope_design(model)

  if (chosen == "gpq") {
    yts <- .model_yts(model, time_var)
    c_resid <- if (target == "observable") 1 else 0       # future-obs vs batch-mean
    g <- .gpq_core(yts$y, yts$time, yts$batch, t0 = t0s,
                    content = content, conf = conf, sides = sides,
                    c_resid = c_resid, M = M, seed = seed)
    return(.gpq_rows(g, design = design, sides = sides, singular = singular))
  }

  ## closed-form MLS (observable target). true_value never reaches here: auto
  ## routes batch-mean to GPQ, and engine = "closed" + true_value already errored.
  rows <- lapply(t0s, function(t0) {
    comp <- vfun_extract(model, t0 = t0, target = target)
    ti_vfun(comp, content = content, conf = conf, sides = sides)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

## Fixed-slope single-obs / staggered route -> the Montes ANOVA/MLS closed form.
.ti_m1_montes_path <- function(model, newdata, content, conf, sides, time_name) {
  yts <- .model_yts(model, time_name)
  t0s <- if (is.null(newdata)) yts$time[1] else newdata[[time_name]]
  if (is.null(t0s) || !length(t0s)) {
    stop("Could not determine the evaluation time point(s) for this interval; ",
         "supply `newdata` with a numeric '", time_name,
         "' column giving the time(s) to evaluate at.", call. = FALSE)
  }
  .montes_core(yts$y, yts$time, yts$batch, t0 = t0s,
               content = content, conf = conf, sides = sides)
}

#' Tolerance interval for a linear mixed model
#'
#' One of the four interval verbs, and the dispatched tolerance interval. It
#' classifies the fitted model and routes by random/fixed structure, the
#' tolerance target (`over`), and the singular flag to
#' the engine that is exact/validated for that case --- all transparently:
#'
#' * **Random slope.** `over = "observation"` uses the REML closed-form MLS
#'   interval when the fit is non-singular and the GPQ interval when singular;
#'   `over = "group_mean"` always uses GPQ.
#' * **Random intercept + fixed slope.** Clean balanced -> EMS engine; single-
#'   observation / staggered / unbalanced -> the Montes ANOVA/MLS closed form.
#' * **Random intercept, no fixed slope** -> EMS engine.
#'
#' For slope models the between-group variance `V_G(t0)` is **quadratic in the
#' time covariate**, so evaluating at a `t0` outside the observed time range is
#' an extrapolation of that quadratic and should be done knowingly (the interval
#' is computed without a guard or warning). The same applies to the other slope
#' verbs ([ci_lmm()], [pi_lmm()], [new_group_mean_lmm()]).
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
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`.
#' @param method `"auto"` (default) chooses the algorithm by structure; one of
#'   `"ems-mls"`, `"reml-mls"`, `"anova-mls"`, `"gpq"` forces it and errors if
#'   incompatible with the model.
#' @param control [lmm_interval_control()] tuning (GPQ seed / M).
#'
#' @inherit ci_lmm return
#' @seealso [ci_lmm()], [pi_lmm()], [new_group_mean_lmm()], [lmm_predict()],
#'   [reintervals-models].
#' @examples
#' # random intercepts, no slope -> EMS engine
#' fi <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
#' ti_lmm(fi, newdata = data.frame(Days = 0))
#' # random slope, future-observation -> REML closed form
#' fs <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy)
#' ti_lmm(fs, newdata = data.frame(Days = 5))
#' # random slope, batch-mean -> GPQ
#' ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean",
#'        control = lmm_interval_control(seed = 1))
#' @export
ti_lmm <- function(model, newdata = NULL, P = 0.95,
                   over = c("observation", "group_mean"),
                   alpha = NULL, level = NULL,
                   sides = c("two.sided", "lower", "upper"),
                   method = "auto", control = lmm_interval_control()) {
  .assert_lmm(model)
  over   <- match.arg(over)
  sides  <- match.arg(sides)
  method <- match.arg(method, c("auto", "ems-mls", "reml-mls", "anova-mls", "gpq"))
  ## boundary: public alpha/level -> internal `conf` (gamma); public `P` is the
  ## content proportion passed straight to the engines as `content`.
  conf   <- 1 - .resolve_alpha(alpha, level)
  if (!is.numeric(P) || length(P) != 1L || is.na(P) || P <= 0 || P >= 1) {
    stop("`P` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  seed <- control$seed
  M    <- control$M
  target <- if (over == "observation") "observable" else "true_value"

  ## ---- random slope present -> target-aware slope path -----------------------
  if (.has_random_slope(model)) {
    engine <- switch(method,
      auto       = "auto",
      `reml-mls` = "closed",
      gpq        = "gpq",
      stop(sprintf("method = \"%s\" is not available for a %s model; use ",
                   method, .model_class(model)),
           "\"auto\", \"reml-mls\", or \"gpq\".", call. = FALSE))
    out <- .ti_slope_path(model, newdata = newdata, content = P, conf = conf,
                          target = target, sides = sides, engine = engine,
                          seed = seed, M = M)
    return(.unify_intervals(out, model, kind = "ti", P_value = P))
  }

  ## ---- random intercept only -------------------------------------------------
  fe <- lme4::fixef(model)
  covars <- setdiff(names(fe), "(Intercept)")
  time_name <- if (length(covars) == 1L) covars else NA_character_
  is_m1 <- !is.na(time_name)

  ## batch-mean tolerance is GPQ-only (random slope); DNE for random intercept / fixed-slope
  if (over == "group_mean") {
    avail <- if (is_m1) {
      "ti_lmm(over = \"observation\") or new_group_mean_lmm() (the prediction interval)"
    } else {
      "ti_lmm(over = \"observation\")"
    }
    stop(sprintf(paste0("The batch-mean tolerance interval (over = ",
         "\"group_mean\") is not implemented for a %s model; available here: %s."),
         .model_class(model), avail), call. = FALSE)
  }

  if (method %in% c("reml-mls", "gpq")) {
    stop(sprintf("method = \"%s\" is not available for a %s model; use ",
                 method, .model_class(model)),
         "\"auto\", \"ems-mls\", or \"anova-mls\".", call. = FALSE)
  }

  ## Fixed-slope single-obs / staggered -> Montes ANOVA; else EMS. `method` can
  ## force. Uses the shared .is_single_obs_m1() detector (identical to the inline
  ## `is_m1 && !.m1_is_clean_balanced()` it replaces) so the routing here and the
  ## coverage-verb warning can never diverge.
  use_montes <- .is_single_obs_m1(model)
  if (method == "anova-mls") {
    if (!is_m1) {
      stop("method = \"anova-mls\" needs a fixed-slope (random intercept + fixed slope) model.", call. = FALSE)
    }
    use_montes <- TRUE
  }
  if (method == "ems-mls") use_montes <- FALSE

  if (use_montes) {
    out <- .ti_m1_montes_path(model, newdata = newdata, content = P, conf = conf,
                              sides = sides, time_name = time_name)
    return(.unify_intervals(out, model, kind = "ti", P_value = P))  # carries anova-mls
  }

  if (sides != "two.sided") {
    stop("The EMS engine ships two-sided tolerance intervals only; one-sided ",
         "applies to the slope path.", call. = FALSE)
  }
  out <- .ti_ems(model, newdata = newdata, level = P, conf = conf, target = target)
  .unify_intervals(out, model, kind = "ti", P_value = P, method_tag = "ems-mls")
}
