# =============================================================================
# Slope-engine CI / PI / new-group-mean interval constructors. The TI branch is
# ti_vfun(); these are the three coverage-based interval branches for the slope
# models, the closed-form analogues of the EMS ci_mean / pi_newobs:
#
#   CI              : mu +/- t(nu_F) * sqrt(V_F)
#   PI              : mu +/- t(nu_T) * sqrt(V_F + V_T)        (V_T = V_G + s2e)
#   new_group_mean  : mu +/- t(nu_G) * sqrt(V_F + V_G)        (PI minus residual;
#                     the CInew quantity -- a prediction interval for one new
#                     group's mean; coverage = `level`, NO content P).
#
# All inputs are already on the vfun_components container (mean, var_mean = V_F,
# V_G, V_T, dfs[c("F","T","G")]); vfun_extract computes those df with the same
# resolve_df (Satterthwaite -> containment) logic, so these constructors are
# pure arithmetic off the container -- method-pure, no lme4, exactly like ti_vfun.
#
# These constructors are validated against the external Montes/GPQ-source oracle
# (test-slope-oracle.R) on committed components.
#
# V_G handling: vfun_extract floors V_G at 0 (the signed-s01 covariance can
# legitimately make the quadratic form dip below 0 at an extreme t0); the floor
# is flagged via vg_floored. See ?reintervals-models.
# =============================================================================

## One Student-t interval off a vfun_components container, for whichever of the
## three "level-based" slope intervals is requested. 
##   kind = "CI"             -> se = sqrt(V_F),        df = nu_F
##   kind = "PI"             -> se = sqrt(V_F + V_T),  df = nu_T
##   kind = "new_group_mean" -> se = sqrt(V_F + V_G),  df = nu_G
.vfun_interval <- function(comp, kind = c("CI", "PI", "new_group_mean"),
                           level = 0.95,
                           sides = c("two.sided", "lower", "upper")) {
  if (!inherits(comp, "vfun_components")) {
    stop("`comp` must be a `vfun_components` object (see ?vfun_extract).",
         call. = FALSE)
  }
  kind  <- match.arg(kind)
  sides <- match.arg(sides)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }

  mu  <- comp$mean
  V_F <- comp$var_mean
  V_G <- comp$V_G
  V_T <- comp$V_T

  pick <- switch(kind,
    CI             = list(df = comp$dfs[["F"]], var = V_F),
    PI             = list(df = comp$dfs[["T"]], var = V_F + V_T),
    new_group_mean = list(df = comp$dfs[["G"]], var = V_F + V_G)
  )
  nu <- pick$df
  se <- sqrt(pick$var)

  ## a non-finite/non-positive df degrades to ~normal.
  alpha <- 1 - level
  nu_use <- if (!is.finite(nu) || nu <= 0) 1e6 else nu
  tcrit <- if (sides == "two.sided") stats::qt(1 - alpha / 2, nu_use) else
    stats::qt(level, nu_use)
  hw <- tcrit * se

  data.frame(
    time = comp$t0, type = kind, estimate = mu,
    lower = if (sides == "upper") -Inf else mu - hw,
    upper = if (sides == "lower")  Inf else mu + hw,
    se = se, level = level, df = nu, sides = sides,
    design = comp$type, method = "reml-mls",        # harmonized engine tag
    singular = isTRUE(comp$singular), boundary = isTRUE(comp$boundary),
    vg_floored = isTRUE(attr(comp, "vg_floored")),  # V_G<0 floored at this t0
    stringsAsFactors = FALSE
  )
}

## EMS-side new-group-mean (ENW for random intercept and fixed-slope): the Francq
## prediction interval for the unobservable level "true value" -- pi_newobs on a
## vc_extract(target = "true_value") container (between-level components only).
## Reuses the existing, EMS-validated pi_newobs true_value path.
.ngm_ems <- function(model, newdata = NULL, level = 0.95) {
  ## Mirror intervals_lmm's per-row loop, but on the between-level "true value"
  ## (target = "true_value") -> the EMS new-group-mean. Reuses the frozen EMS
  ## building blocks (.theta_cov / .fixed / re_components / pi_newobs); numbers
  ## are pi_newobs's. Returns a tidy native data.frame (NOT an re_interval),
  ## fixing the recon S5 list-vs-data.frame inconsistency.
  desc  <- classify_design(lme4::getME(model, "flist"),
                           length(lme4::getME(model, "y")))
  theta <- .theta_cov(model)
  dc    <- .design_components(desc, theta$components, "true_value",
                             .counts_if_unbalanced(model, desc))
  fx    <- .fixed(model, newdata)
  rows <- lapply(seq_len(nrow(fx$L)), function(i) {
    comp <- re_components(
      components = theta$components, dfs = c(ci = Inf, pi = theta$df_pi),
      mean = fx$fit[i], var_mean = fx$var_fix[i], coefs = dc$coefs,
      ems = dc$ems, target = "true_value", design = desc$type
    )
    it <- pi_newobs(comp, level = level)
    data.frame(estimate = it$estimate, lower = it$lower, upper = it$upper,
               level = level, df = it$df, sides = "two.sided",
               stringsAsFactors = FALSE)
  })
  out <- cbind(fx$newdata, do.call(rbind, rows)); rownames(out) <- NULL
  attr(out, "singular") <- isTRUE(theta$singular)
  out
}

## Internal dispatch for the three coverage intervals, wrapped by the ci_lmm /
## pi_lmm / new_group_mean_lmm verbs. Routes
## per the verified capability matrix: random-intercept / fixed-slope -> EMS (CI/PI via
## intervals_lmm, new_group_mean via .ngm_ems); random-slope -> the slope ports.
.interval_lmm <- function(model, newdata = NULL,
                          kind = c("CI", "PI", "new_group_mean"),
                          level = 0.95, sides = c("two.sided", "lower", "upper")) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted `lmerMod` (an lmer fit).", call. = FALSE)
  }
  kind  <- match.arg(kind)
  sides <- match.arg(sides)

  if (.has_random_slope(model)) {
    fe <- lme4::fixef(model)
    time_var <- setdiff(names(fe), "(Intercept)")
    if (length(time_var) != 1L) {
      stop("random-slope intervals require a single continuous covariate (e.g. ~ 1 + time); ",
         "multiple fixed effects with a random slope are a planned v2 extension.",
         call. = FALSE)
    }
    t0s <- if (is.null(newdata)) stats::model.frame(model)[[time_var]][1] else
      newdata[[time_var]]
    rows <- lapply(t0s, function(t0) {
      comp <- vfun_extract(model, t0 = t0, target = "observable")
      .vfun_interval(comp, kind = kind, level = level, sides = sides)
    })
    out <- do.call(rbind, rows); rownames(out) <- NULL
    return(out)
  }

  ## random intercept (incl. fixed-slope): EMS engine.
  if (kind == "new_group_mean") {
    if (sides != "two.sided") {
      stop("The EMS new-group-mean interval is two-sided only in this stage.",
           call. = FALSE)
    }
    return(.ngm_ems(model, newdata = newdata, level = level))
  }
  if (sides != "two.sided") {
    stop("The EMS CI/PI are two-sided only; one-sided applies to the slope path.",
         call. = FALSE)
  }
  intervals_lmm(model, newdata = newdata, level = level, conf = level,
                which = kind)
}
