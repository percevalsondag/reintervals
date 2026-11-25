# =============================================================================
# Slope-engine CI / PI / new-group-mean interval constructors (method "reml"):
# the three coverage-based interval branches for the slope models, the
# REML-based analogues of the EMS ci_mean / pi_newobs. (The slope tolerance
# interval is GPQ, gpq.R; the closed-form `ti_vfun` / "reml-mls" TI was removed
# in 1.0.0.)
#
#   CI              : mu +/- t(nu_F) * sqrt(V_F)
#   PI              : mu +/- t(nu_T) * sqrt(V_F + V_T)        (V_T = V_G + s2e)
#   new_group_mean  : mu +/- t(nu_G) * sqrt(V_F + V_G)        (PI minus residual;
#                     the CInew quantity -- a prediction interval for one new
#                     group's mean; coverage = `level`, NO content P).
#
# The PI has the validation arms slope_pi_method = "fv-anova" / "hybrid"
# (.pi_fv_anova below, on the per-lot OLS lines of the fv-anova NGM arm).
#
# All inputs are already on the vfun_components container (mean, var_mean = V_F,
# V_G, V_T, dfs[c("F","T","G")]); vfun_extract computes those df with the same
# resolve_df (Satterthwaite -> containment) logic, so these constructors are
# pure arithmetic off the container -- method-pure, no lme4.
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
    .abort("`comp` must be a `vfun_components` object (see ?vfun_extract).",
           class = "reintervals_error_input")
  }
  kind  <- match.arg(kind)
  sides <- match.arg(sides)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        level <= 0 || level >= 1) {
    .abort("`level` must be a single number strictly between 0 and 1.",
           class = "reintervals_error_input")
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

  ## df source, from vfun_extract's per-df fallback flags: TRUE = the
  ## Satterthwaite df was unavailable and containment was used; NA = containment
  ## was forced (df_method = "containment"); no attribute = a hand-built
  ## container (source unknown).
  key <- c(CI = "F", PI = "T", new_group_mean = "G")[[kind]]
  fb  <- attr(comp, "df_fallback")
  rule <- if (is.null(fb) || !key %in% names(fb)) NA_character_
          else if (isFALSE(fb[[key]])) "satterthwaite" else "containment"
  ## a non-default slope_ci_df arm replaced the CI df: report the arm
  if (kind == "CI" && !is.null(attr(comp, "ci_arm"))) rule <- attr(comp, "ci_arm")
  ## A non-finite or non-positive df (only reachable from a hand-built container;
  ## vfun_extract always returns a finite df) used to be replaced by 1e6, i.e. a
  ## normal quantile. It now takes the slope engine's own containment df, the
  ## same rule vfun_extract falls back to (.vfun_df_containment, vfun-extract.R):
  ## B - k_re (B - 2 with a random slope) for the between-group pieces (CI, NGM;
  ## the PI when V_G / V_T >= 0.5), the residual df otherwise. The design sizes
  ## it needs are carried by vfun_extract as attr "containment".
  if (!is.finite(nu) || nu <= 0) {
    parts <- attr(comp, "containment")
    if (is.null(parts)) {
      .abort("The ", kind, " df is not finite and the container carries no ",
             "design sizes to form the containment df; build it with ",
             "vfun_extract().", class = "reintervals_error_numerical")
    }
    w_between <- switch(kind, CI = 1, new_group_mean = 1,
                        PI = if (V_T > 0) V_G / V_T else 0)
    nu <- .vfun_df_containment(parts, w_between)
    rule <- "containment"
  }
  alpha <- 1 - level
  tcrit <- if (sides == "two.sided") stats::qt(1 - alpha / 2, nu) else
    stats::qt(level, nu)
  hw <- tcrit * se

  data.frame(
    time = comp$t0, type = kind, estimate = mu,
    lower = if (sides == "upper") -Inf else mu - hw,
    upper = if (sides == "lower")  Inf else mu + hw,
    se = se, level = level, df = nu, df_rule = rule, sides = sides,
    design = comp$type, method = "reml",            # harmonized engine tag
    singular = if (is.na(comp$singular)) NA else isTRUE(comp$singular), boundary = isTRUE(comp$boundary),
    vg_floored = isTRUE(attr(comp, "vg_floored")),  # V_G<0 floored at this t0
    vc_cov_route = if (is.null(attr(comp, "vc_cov_route"))) NA_character_
                   else attr(comp, "vc_cov_route"),  # Cov(phi) route (.vfun_vc_cov)
    ## why Cov(phi) is unavailable (NA when it was obtained, or unknown)
    vc_cov_reason = if (is.null(attr(comp, "vc_cov_reason"))) NA_character_
                    else attr(comp, "vc_cov_reason"),
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
  value_rows <- function(est, W, df, rule, method, reason = NA_character_) {
    if (any(!is.finite(W) | W <= 0)) {
      .abort("The new-group-mean variance W is not positive for this fit (e.g. ",
             "all group means are identical, so the between-group mean square is ",
             "zero); the interval is not available.", class = "reintervals_error_numerical")
    }
    hw <- tcrit(df) * sqrt(W)
    data.frame(type = "new_group_mean", estimate = est,
               lower = if (sides == "upper") -Inf else est - hw,
               upper = if (sides == "lower")  Inf else est + hw,
               level = level, df = df, df_rule = rule, W = W,
               df_fallback_reason = reason,   # NA where the df did not fall back
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
    .abort(sprintf(paste0("The new-group-mean interval is not available for a '%s' ",
                          "design: it has no closed-form variance decomposition. ",
                          "Supported: one-way, two-fold nested, and two crossed ",
                          "factors with their interaction."), desc$type),
           class = "reintervals_error_design")
  }
  ## 2c / 2d / 2e: EMS route. Covariates beyond the intercept enter through V_F
  ## only. UNVALIDATED: the one-way design with other fixed effects (2e), and
  ## nested / crossed designs with covariates, were not part of the validation
  ## simulations (which used intercept-only nested / crossed fits); 2e is
  ## checked by simulation separately.
  dc <- .design_components_fit(model, desc, .var_components(model), "observable")
  if (is.null(dc$ems)) {
    .abort("The new-group-mean interval is not available: no closed-form ",
           "variance decomposition for this design.", class = "reintervals_error_design")
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
      if (!is.finite(pA)) return(NA_real_)
      co <- k[j]
      co[1] <- co[1] + pA * v / ms[1]
      co[2] <- co[2] + (1 - pA) * v / ms[2]
      sum(co * ms[j])^2 / sum((co * ms[j])^2 / dfj[j])
    }, numeric(1))
    ## Undefined Satterthwaite term (k_A EMS_A + k_B EMS_B not positive, so pA
    ## is NaN/Inf, or a non-finite / non-positive df): never return a NaN
    ## limit. Fall back to the containment df min(A, B) - 1, reported per row.
    ## limit. A finite df below 1 falls back too (the default-rule alignment:
    ## non-finite or < 1 -> containment, as .ci_ddf).
    bad  <- !is.finite(df) | df < 1
    rule <- ifelse(bad, "containment", "satterthwaite-vf")
    row_reason <- ifelse(!bad, NA_character_,
                         if (!is.finite(pA)) "k_A EMS_A + k_B EMS_B is not positive"
                         else ifelse(is.finite(df), "df below 1", "non-finite df"))
    fb_reason <- if (!any(bad)) NULL else unique(row_reason[bad])
    if (any(bad)) df[bad] <- .ddf_containment(model)
  } else {
    df <- rep(dfj[[1]], length(VF))        # Gm: top-factor stratum df (A - 1,
                                           # less any between-group fixed columns)
    rule <- "top-factor"
    fb_reason <- NULL
    row_reason <- NA_character_
  }
  out <- cbind(fx$newdata, value_rows(fx$fit, W, df, rule, "ems-mls", row_reason))
  rownames(out) <- NULL
  attr(out, "singular") <- .is_singular(model)
  attr(out, "df_fallback_reason") <- fb_reason
  out
}

## Internal dispatch for the three coverage intervals, wrapped by the ci_lmm /
## pi_lmm / new_group_mean_lmm verbs. Routes
## per the verified capability matrix: random-intercept / fixed-slope -> EMS (CI/PI via
## intervals_lmm, new_group_mean via .ngm_ems); random-slope -> the slope ports.
.interval_lmm <- function(model, newdata = NULL,
                          kind = c("CI", "PI", "new_group_mean"),
                          level = 0.95, sides = "two.sided",
                          control = lmm_interval_control(),
                          applied = .applied_arms(model, kind, control)) {
  if (!inherits(model, "lmerMod")) {
    .abort("`model` must be a fitted `lmerMod` (an lmer fit).",
           class = "reintervals_error_input")
  }
  kind  <- match.arg(kind)
  ## vector level / sides (B3) reach this function only for the GPQ
  ## new-group-mean arm, which evaluates every combination on one draw set;
  ## every other route is called once per combination by .coverage_verb
  sides <- .resolve_sides(sides)
  ## df arms (lmm_interval_control): each option is read only where it applies
  ci_df       <- .arm_value(control, "ci_df")
  slope_ci_df <- .arm_value(control, "slope_ci_df")
  pi_df       <- .arm_value(control, "pi_df")

  if (.has_random_slope(model)) {
    fe <- lme4::fixef(model)
    time_var <- setdiff(names(fe), "(Intercept)")
    if (length(time_var) != 1L) {
      .abort("random-slope intervals require a single continuous covariate (e.g. ~ 1 + time); ",
             "multiple fixed effects with a random slope are not currently supported.",
             class = "reintervals_error_unsupported")
    }
    t0s <- if (is.null(newdata)) stats::model.frame(model)[[time_var]][1] else
      newdata[[time_var]]
    ## slope_ci_df arms other than the default replace the CI df (and, for
    ## Kenward-Roger, the variance) of the container; PI / NGM are unaffected.
    arm <- if (kind == "CI") slope_ci_df else "hessian-satterthwaite"
    kr  <- if (arm == "kenward-roger") .kr_parts(model, option = "slope_ci_df") else NULL
    ## slope_ngm_method arms other than the default "reml" (one call over all
    ## t0: the GPQ arm uses common random numbers across t0)
    slope_ngm <- if (kind == "new_group_mean" && "slope_ngm_method" %in% names(applied))
      applied[["slope_ngm_method"]] else "reml"
    if (slope_ngm != "gpq") .assert_scalar_combo(level, sides)
    if (slope_ngm != "reml") {
      return(.ngm_rs_arm(model, t0s, time_var, level = level, sides = sides,
                         arm = slope_ngm, control = control))
    }
    ## slope_pi_method (PI only): "reml" (default, the code below), "fv-anova",
    ## or "hybrid" = "fv-anova" on a singular fit (lme4::isSingular() TRUE or
    ## the check failed, read as singular as everywhere else), "reml" otherwise
    slope_pi <- if (kind == "PI" && "slope_pi_method" %in% names(applied))
      applied[["slope_pi_method"]] else "reml"
    ## hybrid rows carry slope_pi_route_reason: NA, or why the fv-anova PI was
    ## infeasible on this singular fit and "reml" was used instead
    route_reason <- NULL
    if (slope_pi == "hybrid") {
      route_reason <- NA_character_
      slope_pi <- if (isFALSE(.is_singular(model))) "reml" else "fv-anova"
      if (slope_pi == "fv-anova") {
        why <- .fv_pi_infeasible(.model_yts(model, time_var))
        if (is.null(why)) {
          out <- tryCatch(
            .pi_rs_fv_anova(model, t0s, time_var, level = level, sides = sides,
                            opt = "slope_pi_method = \"hybrid\" (fv-anova on this singular fit)"),
            reintervals_error_numerical = function(e) e)
          if (!inherits(out, "error")) {
            out$slope_pi_route_reason <- NA_character_
            return(out)
          }
          why <- "zero prediction variance"
        }
        ## never an error where "reml" computes: fall back, with a warning
        route_reason <- paste0("fv-anova infeasible: ", why)
        .warn("slope_pi_method = \"hybrid\": the fit is singular but the fv-anova ",
              "prediction interval is infeasible (", why, "); the \"reml\" ",
              "interval is used instead.",
              class = "reintervals_warning_hybrid_fallback",
              data = list(reason = route_reason))
        slope_pi <- "reml"
      }
    }
    if (slope_pi == "fv-anova") {
      return(.pi_rs_fv_anova(model, t0s, time_var, level = level, sides = sides,
                             opt = "slope_pi_method = \"fv-anova\""))
    }
    reasons <- character(0)
    rows <- lapply(t0s, function(t0) {
      comp <- vfun_extract(model, t0 = t0, guard_t0 = t0s)
      if (arm != "hessian-satterthwaite") comp <- .slope_ci_arm(comp, model, arm, kr)
      row <- .vfun_interval(comp, kind = kind, level = level, sides = sides)
      ## why the Hessian-Satterthwaite df fell back (read by .warn_df_fallback)
      row_reason <- NA_character_
      if (identical(row$df_rule, "containment")) {
        key <- c(CI = "F", PI = "T", new_group_mean = "G")[[kind]]
        why <- attr(comp, "df_fallback_why")[[key]]
        ## no covariance: its reason (vc_cov_reason; "singular fit" or
        ## "singularity check failed; treated as singular" on a singular fit)
        row_reason <-
          if (is.null(comp$cov)) {
            r <- attr(comp, "vc_cov_reason")
            if (!is.null(r) && !is.na(r)) r
            else "non-PD or unavailable covariance of the variance parameters"
          }
          else if (!is.null(why) && !is.na(why)) why else "non-finite df"
        reasons <<- c(reasons, row_reason)
      }
      ## the PI rows carry the reason per row (NA where the df did not fall
      ## back), as the random-intercept rows do
      if (kind == "PI") row$df_fallback_reason <- row_reason
      row
    })
    out <- do.call(rbind, rows); rownames(out) <- NULL
    if (length(reasons)) attr(out, "df_fallback_reason") <- unique(reasons)
    if (!is.null(route_reason)) out$slope_pi_route_reason <- route_reason
    return(out)
  }

  ## random intercept (incl. fixed-slope): EMS engine.
  .assert_scalar_combo(level, sides)
  if (kind == "new_group_mean") {
    ngm <- if ("ngm_method" %in% names(applied)) applied[["ngm_method"]] else "package"
    if (ngm != "package") {
      return(.ngm_ri_arm(model, newdata = newdata, level = level, sides = sides,
                         arm = ngm))
    }
    return(.ngm_ri(model, newdata = newdata, level = level, sides = sides))
  }
  out <- intervals_lmm(model, newdata = newdata, level = level, conf = level,
                       which = kind, ci_method = ci_df, pi_method = pi_df)
  .one_sided_t(out, level = level, sides = sides)
}

## slope_pi_method = "fv-anova" (and "hybrid" on a singular fit): a prediction
## interval for one new observation at t0 from a new lot, on the per-lot OLS
## lines of the fv-anova new-group-mean arm (.lotline, ngm-arms.R). With K lots,
## S2L the variance of the K fitted values at t0, S2E the pooled within-lot
## residual variance on df_E = n - 2K df:
##   center  Ybar(t0), the mean of the K fitted values (as the NGM arm)
##   W       S2L (1 + 1/K) + S2E   (the NGM arm's W plus the residual variance)
##   df      W^2 / [ (S2L (1 + 1/K))^2 / df_L + S2E^2 / df_E ]
##           (Satterthwaite on the two independent mean squares), with df_L the
##           effective df of S2L (.fv_df_lots): K - 1 in a balanced design;
##           smaller when the lots' leverages at t0 differ
## In a balanced design E[W] - Var(Y_new - Ybar) = sigma^2 h(t0) >= 0 (h the
## common leverage at t0): conservative, like the NGM arm, and W > 0 whenever
## S2E > 0. Uses the raw data only, not the fitted random-effects structure.
## No df fallback. One-sided: t_{1 - alpha}. df_E = 0 (every lot has exactly
## two observations) is an error: there is no residual variance to add.
## Effective (Satterthwaite) df of S2L = Y'AY, A = (I - J/K) / (K - 1), when
## the K per-lot fitted values at t0 are independent with variances
## v_k = max(VG_hat, 0) + a_k S2E (a_k the leverage of lot k at t0, VG_hat =
## S2L - abar S2E): df_L = tr(AV)^2 / tr((AV)^2), V = diag(v). K - 1 exactly
## when the leverages are equal (max |a_k - abar| <= 1e-12 abar; the balanced
## value, analytically) or every v_k is 0.
.fv_df_lots <- function(L) {
  K <- L$K; a <- L$a
  if (max(abs(a - L$abar)) <= 1e-12 * L$abar) return(K - 1)
  v <- max(L$S2L - L$abar * L$S2E, 0) + a * L$S2E
  if (all(v == 0)) return(K - 1)
  AV <- (diag(K) - 1 / K) %*% diag(v, K) / (K - 1)
  sum(diag(AV))^2 / sum(AV * t(AV))
}

.pi_fv_anova <- function(y, time, batch, t0s, level = 0.95,
                         sides = c("two.sided", "lower", "upper"),
                         opt = "slope_pi_method = \"fv-anova\"") {
  sides <- match.arg(sides)
  rows <- lapply(t0s, function(t0) {
    L <- .lotline(y, time, batch, t0, what = "prediction", opt = opt)
    if (L$df_E < 1L) {
      .abort(sprintf(paste0("%s needs a within-lot residual variance: every lot has ",
                            "exactly two observations (n = 2K), so it has no ",
                            "residual df."), opt),
             class = "reintervals_error_time_points")
    }
    W_L <- L$S2L * (1 + 1 / L$K)           # the fv-anova NGM variance
    W <- W_L + L$S2E
    df_L <- .fv_df_lots(L)
    df <- W^2 / (W_L^2 / df_L + L$S2E^2 / L$df_E)
    if (!is.finite(W) || W <= 0 || !is.finite(df) || df <= 0) {
      .abort("The fv-anova prediction variance at t0 = ", t0, " is zero (equal per-lot ",
             "fitted values and no residual variation); the interval is not available.",
             class = "reintervals_error_numerical")
    }
    tq <- if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df) else
      stats::qt(level, df)
    hw <- tq * sqrt(W)
    data.frame(time = t0, estimate = L$Ybar,
               lower = if (sides == "upper") -Inf else L$Ybar - hw,
               upper = if (sides == "lower")  Inf else L$Ybar + hw,
               df = df, W = W, se = sqrt(W), K = L$K, df_resid = L$df_E,
               fv_df_lots = df_L, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

## Why the fv-anova prediction interval cannot be formed on these data (the
## refusals of .lotline / .pi_fv_anova), or NULL when it can. Read by the
## "hybrid" route, which then falls back to "reml" instead of erroring.
.fv_pi_infeasible <- function(yts) {
  ok <- !is.na(yts$y) & !is.na(yts$time) & !is.na(yts$batch)
  b <- droplevels(factor(yts$batch)[ok]); tt <- yts$time[ok]
  K <- nlevels(b)
  if (K < 2L) return("fewer than 2 lots")
  n_t <- tapply(tt, b, function(v) length(unique(v)))
  if (any(n_t < 2L))
    return(sprintf("%d lot(s) with fewer than two distinct time points", sum(n_t < 2L)))
  if (sum(ok) - 2L * K < 1L) return("no residual df (n = 2K)")
  NULL
}

## The fv-anova prediction interval in the slope-path row schema. `method` is
## "fv-anova", the construction used (the rows of the default "reml" PI carry
## "reml"); df_fallback_reason is NA (no fallback).
.pi_rs_fv_anova <- function(model, t0s, time_var, level, sides, opt) {
  yts <- .model_yts(model, time_var)
  f <- .pi_fv_anova(yts$y, yts$time, yts$batch, t0s, level = level, sides = sides,
                    opt = opt)
  data.frame(time = f$time, type = "PI", estimate = f$estimate,
             lower = f$lower, upper = f$upper, se = f$se, level = level,
             df = f$df, df_rule = "fv-anova-satterthwaite", W = f$W, sides = sides,
             design = .slope_design(model), method = "fv-anova",
             singular = .is_singular(model), boundary = NA, vg_floored = NA,
             K = f$K, df_resid = f$df_resid, fv_df_lots = f$fv_df_lots,
             df_fallback_reason = NA_character_,
             stringsAsFactors = FALSE)
}

## Internal guard (B3): a scalar engine route received a vector level / sides.
.assert_scalar_combo <- function(level, sides) {
  if (length(level) != 1L || length(sides) != 1L) {
    .abort("Internal: this interval route takes one level and one side per call.",
           class = "reintervals_error_internal")
  }
  invisible(TRUE)
}

## Random-slope CI df arm (slope_ci_df) applied to a vfun_components container
## at its t0: replaces dfs["F"] and, for Kenward-Roger, var_mean, and marks the
## F df source with the arm (attr "ci_arm", read by .vfun_interval as df_rule).
##   "k-minus-1"     df = K - 1 (K = number of lots), the unadjusted V_F: the
##                   exact df for balanced random-coefficient models (GLS =
##                   mean of the per-lot OLS lines; Chow & Shao 1991).
##   "kenward-roger" pbkrtest::Lb_ddf with the KR-adjusted covariance
##                   pbkrtest::vcovAdj(); V_F = l' Vadj l (A3 pairing).
.slope_ci_arm <- function(comp, model, arm, kr = NULL) {
  if (arm == "k-minus-1") {
    comp$dfs[["F"]] <- max(1, attr(comp, "containment")$B - 1)
  } else if (arm == "kenward-roger") {
    fe <- lme4::fixef(model)                 # (Intercept), time: l = (1, t0)
    l <- matrix(ifelse(names(fe) == "(Intercept)", 1, comp$t0), nrow = 1L)
    ## A validation arm never falls back to another df (step 13 would credit the
    ## fallback's coverage to Kenward-Roger): a failed, non-finite or
    ## non-positive KR df, or a non-positive KR variance, is an error.
    ddf <- tryCatch(as.numeric(pbkrtest::Lb_ddf(l, kr$V0, kr$Vadj)),
                    error = function(e) NA_real_)
    var_kr <- as.numeric(l %*% kr$Vadj_mat %*% t(l))
    if (!isTRUE(length(ddf) == 1L && is.finite(ddf) && ddf > 0) ||
          !isTRUE(is.finite(var_kr) && var_kr > 0)) {
      .abort(sprintf(paste0("slope_ci_df = \"kenward-roger\": the Kenward-Roger df ",
                            "or variance could not be computed at t0 = %g for this ",
                            "model (no fallback is applied to a validation arm); ",
                            "use the default slope_ci_df = \"hessian-satterthwaite\"."),
                     comp$t0),
             class = "reintervals_error_numerical")
    }
    comp$dfs[["F"]] <- ddf
    comp$var_mean   <- var_kr
  } else {
    .abort("Internal: unknown slope_ci_df arm '", arm, "'.",
           class = "reintervals_error_internal")
  }
  attr(comp, "ci_arm") <- arm
  comp
}

## One-sided EMS CI / PI: t_{1 - alpha} on the existing symmetric t pivot (same
## center, same df, same standard error, recovered from the two-sided
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
