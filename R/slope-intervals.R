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

## TRUE when the random-intercept new-group-mean interval uses the anova-mls
## core (2a / 2b): one grouping factor, and either an intercept-only fixed part
## or a single numeric (fixed-slope) covariate. Shared with the coverage verb's
## `method` check so routing and validation cannot diverge.
.ngm_uses_anova <- function(model) {
  if (.has_random_slope(model)) return(FALSE)
  length(lme4::getME(model, "cnms")) == 1L &&
    (.is_oneway_intercept_only(model) || !is.na(.m1_time_name(model)))
}

## Random-intercept new-group-mean interval. Target: the TRUE mean of a new
## innermost cell -- a new group (one-way), a new B-within-A cell (nested), a
## new (a, b) cell (crossed) -- or, for the fixed-slope design, of a new lot at
## t0.
## Constructions (validated by simulation in the paper repo's defect study,
## E5a / E5e / E6, except 2e):
##   2a fixed slope, one grouping factor, any balance -> Fa (anova-mls core):
##      W = S2_L (1 + 1/I) + S2_E (t0 - tbar*)^2 / S_tt, Satterthwaite df
##   2b intercept-only one-way -> the slope-free branch: W = S2_L (1 + 1/I),
##      df = I - 1 (== the classical Ca interval on balanced data)
##   2c nested -> Gm: W = V_F + sum_{j != e} k_j EMS_j, df = A - 1 (top factor)
##   2d crossed with interaction -> Gv: same W, V_F folded into the A and B
##      coefficients in proportion to k_A EMS_A : k_B EMS_B, Satterthwaite df on
##      all terms
##   2e one-way with other fixed effects -> Gm with the one-way k, df = A - 1
##   2f anything else (no closed-form EMS decomposition) -> error
## EMS_j are the REML-implied expected mean squares (.design_components; the
## synthesis k for unbalanced nested / crossed). Two-sided only.
.ngm_ri <- function(model, newdata = NULL, level = 0.95,
                   sides = c("two.sided", "lower", "upper")) {
  sides <- match.arg(sides)
  desc  <- .design_of(model)
  ## two-sided: t_{1 - alpha/2}; one-sided: t_{1 - alpha} on the same pivot
  tcrit <- function(df) if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df)
                        else stats::qt(level, df)
  value_rows <- function(est, W, df, rule, method) {
    if (any(!is.finite(W) | W <= 0)) {
      stop("The new-group-mean variance W is not positive for this fit (e.g. ",
           "all group means are identical, so the between-group mean square is ",
           "zero); the interval is not available.", call. = FALSE)
    }
    hw <- tcrit(df) * sqrt(W)
    data.frame(type = "new_group_mean", estimate = est,
               lower = if (sides == "upper") -Inf else est - hw,
               upper = if (sides == "lower")  Inf else est + hw,
               level = level, df = df, df_rule = rule, W = W,
               sides = sides, method = method, stringsAsFactors = FALSE)
  }

  if (.ngm_uses_anova(model)) {
    ## 2a / 2b: the unweighted lot-mean (anova-mls) core
    time_name <- .m1_time_name(model)
    has_time <- !is.na(time_name)
    lot <- lme4::getME(model, "flist")[[1]]
    time <- if (has_time) as.numeric(stats::model.frame(model)[[time_name]]) else NULL
    t0 <- if (!has_time) NULL else if (is.null(newdata)) time[1] else {
      .check_newdata(model, newdata); as.numeric(newdata[[time_name]])
    }
    core <- .montes_core(as.numeric(lme4::getME(model, "y")), time, lot, t0)
    out <- value_rows(core$estimate, attr(core, "ngm_W"), attr(core, "ngm_df"),
                      if (has_time) "anova-satterthwaite" else "anova-single",
                      "anova-mls")
    if (has_time) out <- cbind(time = t0, out)
    attr(out, "singular") <- .is_singular(model)
    return(out)
  }
  if (!desc$type %in% c("oneway", "nested", "crossed")) {
    stop(sprintf(paste0("The new-group-mean interval is not available for a '%s' ",
                        "design: it has no closed-form variance decomposition. ",
                        "Supported: one-way, two-fold nested, and two crossed ",
                        "factors with their interaction."), desc$type),
         call. = FALSE)
  }
  ## 2c / 2d / 2e: EMS route. Covariates beyond the intercept enter through V_F
  ## only. UNVALIDATED: the one-way design with other fixed effects (2e), and
  ## nested / crossed designs with covariates, were not part of the validation
  ## simulations (which used intercept-only nested / crossed fits); 2e is
  ## checked by simulation separately.
  dc <- .design_components_fit(model, desc, .var_components(model), "observable")
  if (is.null(dc$ems)) {
    stop("The new-group-mean interval is not available: no closed-form ",
         "variance decomposition for this design.", call. = FALSE)
  }
  ms <- dc$ems$ms; k <- dc$ems$k; dfj <- dc$ems$df
  j <- seq_len(length(ms) - 1L)            # every term but the residual (last)
  S <- sum(k[j] * ms[j])
  fx <- .fixed(model, newdata)
  VF <- fx$var_fix
  W  <- VF + S
  if (desc$type == "crossed") {
    ## Gv: V_F attached to the A and B terms (positions 1, 2) in proportion to
    ## k_A EMS_A : k_B EMS_B; W = sum_j c_j EMS_j is unchanged.
    pA <- k[1] * ms[1] / (k[1] * ms[1] + k[2] * ms[2])
    df <- vapply(VF, function(v) {
      co <- k[j]
      co[1] <- co[1] + pA * v / ms[1]
      co[2] <- co[2] + (1 - pA) * v / ms[2]
      sum(co * ms[j])^2 / sum((co * ms[j])^2 / dfj[j])
    }, numeric(1))
    rule <- "satterthwaite-vf"
  } else {
    df <- rep(dfj[[1]], length(VF))        # Gm: top-factor stratum df (A - 1,
                                           # less any between-group fixed columns)
    rule <- "top-factor"
  }
  out <- cbind(fx$newdata, value_rows(fx$fit, W, df, rule, "ems-mls"))
  rownames(out) <- NULL
  attr(out, "singular") <- .is_singular(model)
  out
}

## Internal dispatch for the three coverage intervals, wrapped by the ci_lmm /
## pi_lmm / new_group_mean_lmm verbs. Routes
## per the verified capability matrix: random-intercept / fixed-slope -> EMS (CI/PI via
## intervals_lmm, new_group_mean via .ngm_ems); random-slope -> the slope ports.
.interval_lmm <- function(model, newdata = NULL,
                          kind = c("CI", "PI", "new_group_mean"),
                          level = 0.95, sides = c("two.sided", "lower", "upper"),
                          ci_df = c("satterthwaite", "kenward-roger")) {
  if (!inherits(model, "lmerMod")) {
    stop("`model` must be a fitted `lmerMod` (an lmer fit).", call. = FALSE)
  }
  kind  <- match.arg(kind)
  sides <- match.arg(sides)
  ci_df <- match.arg(ci_df)

  if (.has_random_slope(model)) {
    if (kind == "CI" && ci_df == "kenward-roger") {
      stop("ci_df = \"kenward-roger\" applies to random-intercept models; the ",
           "random-slope confidence interval uses its own generalized ",
           "Satterthwaite df.", call. = FALSE)
    }
    fe <- lme4::fixef(model)
    time_var <- setdiff(names(fe), "(Intercept)")
    if (length(time_var) != 1L) {
      stop("random-slope intervals require a single continuous covariate (e.g. ~ 1 + time); ",
         "multiple fixed effects with a random slope are not currently supported.",
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
    return(.ngm_ri(model, newdata = newdata, level = level, sides = sides))
  }
  out <- intervals_lmm(model, newdata = newdata, level = level, conf = level,
                       which = kind, ci_method = ci_df)
  .one_sided_t(out, level = level, sides = sides)
}

## One-sided EMS CI / PI: t_{1 - alpha} on the existing symmetric t pivot (same
## centre, same df, same standard error, recovered from the two-sided
## half-width as hw / t_{1 - alpha/2}). Two-sided input passes through, tagged.
.one_sided_t <- function(out, level, sides) {
  out$sides <- sides
  if (sides == "two.sided") return(out)
  se <- (out$upper - out$lower) / 2 / stats::qt(1 - (1 - level) / 2, out$df)
  q1 <- stats::qt(level, out$df)
  if (sides == "lower") {
    out$lower <- out$estimate - q1 * se; out$upper <- Inf
  } else {
    out$upper <- out$estimate + q1 * se; out$lower <- -Inf
  }
  out
}
