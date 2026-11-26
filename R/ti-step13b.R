# =============================================================================
# Two-sided MLS tolerance interval on random-intercept models (ems-mls and
# anova-mls routes): the step-13b predicate quantities, the extrapolation
# warning and the GPQ reroute (lmm_interval_control(reroute = TRUE)). Under
# validation (step 13b); the thresholds are provisional.
#
# Per row (two-sided MLS rows; NA on one-sided rows):
#   ti_R    the variance of the predicted mean at the evaluation point over the
#           estimated total variance. ems-mls: var_mean / sigma2_t, var_mean =
#           x0' V x0 (fixed-effect variance of the predicted mean) and sigma2_t
#           = sum_j k_j EMS_j, the two inputs of ti_francq(); anova-mls:
#           1 / n_E = VhatYh / Vhat_Y (the Montes / Hoffman-Kringle core,
#           VhatYh = S2_L / I + S2_E (t0 - tbar*)^2 / S_tt).
#   ti_df_T the Satterthwaite df of sigma2_t from the EMS decomposition,
#           (sum_j k_j EMS_j)^2 / sum_j (k_j EMS_j)^2 / df_j, on either route.
# Both are scale-free: a change of units of the response multiplies var_mean,
# every EMS_j and (S2_L, S2_E) by the same factor.
# =============================================================================

## R (ems-mls form, one per evaluation point) and df_T_hat of a
## random-intercept fit, from its EMS decomposition, exactly as .ti_ems()
## builds the ti_francq() inputs. NULL when the decomposition is unavailable
## (no closed form, or any error): the quantities are then NA and the warning
## and the reroute cannot fire.
.ti_mls_predicates <- function(model, newdata) {
  tryCatch(suppressWarnings({
    desc <- .design_of(model)
    comps <- .var_components(model)
    dc <- .design_components_fit(model, desc, comps, "observable")
    if (is.null(dc$ems)) return(NULL)
    fx <- .fixed(model, newdata)
    per <- vapply(seq_len(nrow(fx$L)), function(i) {
      comp <- re_components(components = comps, dfs = c(pi = Inf), mean = fx$fit[i],
                            var_mean = fx$var_fix[i], coefs = dc$coefs, ems = dc$ems,
                            target = "observable", design = desc$type)
      a <- comp$ems$k * comp$ems$ms
      c(R = comp$var_mean / sum(a),
        dfT = if (anyNA(a)) NA_real_ else sum(a)^2 / sum(a^2 / comp$ems$df))
    }, numeric(2))
    list(R_ems = unname(per["R", ]), df_T = unname(per["dfT", 1L]), fx = fx)
  }), error = function(e) NULL)
}

## The step-13b TI columns for the rows of a route they do not apply to (GPQ
## routes: random slope, method = "gpq").
.ti_step13b_na <- function(out) {
  out$ti_R <- NA_real_; out$ti_df_T <- NA_real_; out$ti_warn_fired <- NA
  out$rerouted <- NA; out$reroute_reason <- "not_applicable"
  out
}

## Add ti_R, ti_df_T, ti_warn_fired, rerouted and reroute_reason to the
## native rows `out` of an MLS tolerance-interval call (ems-mls or anova-mls,
## one block of n_eval rows per (sides, P, level) combination, evaluation point
## fastest); signal the extrapolation warning; with reroute = TRUE, replace the
## qualifying two-sided rows by the random-intercept GPQ interval with
## gpq_two_sided = "content" (same seed / M contract as an explicit
## method = "gpq" call).
.ti_step13b <- function(out, model, newdata, P_row, n_combos, control, target) {
  n <- nrow(out)
  two <- out$sides == "two.sided"
  meth <- if ("method" %in% names(out)) out$method else rep("ems-mls", n)
  n_eval <- n %/% n_combos
  ev <- rep_len(seq_len(n_eval), n)
  pr <- .ti_mls_predicates(model, newdata)
  nE <- attr(out, "n_E")
  R_anova <- if (is.null(nE)) rep(NA_real_, n) else 1 / rep_len(nE, n)
  R_ems <- if (is.null(pr)) rep(NA_real_, n) else pr$R_ems[ev]
  R <- ifelse(meth == "anova-mls", R_anova, R_ems)
  dfT <- rep(if (is.null(pr)) NA_real_ else pr$df_T, n)
  R[!two] <- NA_real_; dfT[!two] <- NA_real_

  ## extrapolation warning (provisional thresholds; on by default)
  wR <- .ctl_value(control, "ti_warn_R"); wD <- .ctl_value(control, "ti_warn_dfT")
  fired <- R > wR & dfT >= wD
  fired[!two] <- NA
  if (any(fired %in% TRUE)) .warn_ti_mls_extrapolation(model, newdata, R, dfT, ev,
                                                       fired %in% TRUE, wR, wD)

  rerouted <- ifelse(two, FALSE, NA)
  reason <- ifelse(two, NA_character_, "not_applicable")
  if (isTRUE(.ctl_value(control, "reroute")) && any(two)) {
    rR <- .ctl_value(control, "ti_reroute_R"); rD <- .ctl_value(control, "ti_reroute_dfT")
    .warn("reroute = TRUE is a validation arm (under validation, step 13b): ",
          "two-sided ems-mls / anova-mls tolerance-interval rows with R > ", rR,
          " and df_T_hat >= ", rD, " are computed with method = \"gpq\" and ",
          "gpq_two_sided = \"content\" on fits eligible for the random-intercept GPQ. ",
          "The default is reroute = FALSE.",
          class = "reintervals_warning_unvalidated",
          data = list(option = "reroute", value = TRUE))
    cand <- two & (R > rR & dfT >= rD) %in% TRUE
    if (any(cand)) {
      if (!is.null(.gpq_ri_ineligible(model))) {
        reason[cand] <- "not GPQ-eligible"
      } else {
        out <- .ti_reroute_rows(out, model, which(cand), ev, P_row, control, target,
                                meth, fx = if (is.null(pr)) .fixed(model, newdata) else pr$fx)
        rerouted[cand] <- TRUE
        reason[cand] <- "R > ti_reroute_R and df_T_hat >= ti_reroute_dfT"
      }
    }
  }
  out$ti_R <- R; out$ti_df_T <- dfT; out$ti_warn_fired <- fired
  out$rerouted <- rerouted; out$reroute_reason <- reason
  out
}

## The rerouted rows `idx` of `out`, recomputed by the random-intercept GPQ
## (gpq_two_sided = "content"): one GPQ call over the evaluation points, P and
## levels of those rows (common random numbers, so each row equals the explicit
## method = "gpq" call at that point, P and level with the same seed and M).
## The GPQ route's own unvalidated warnings (method, gpq_two_sided) are muffled:
## the reroute warning covers them.
.ti_reroute_rows <- function(out, model, idx, ev, P_row, control, target, meth, fx) {
  pts <- sort(unique(ev[idx]))
  P_u <- unique(P_row[idx]); conf_u <- unique(out$conf[idx])
  g <- withCallingHandlers(
    .ti_ri_gpq_path(model, newdata = fx$newdata[pts, , drop = FALSE], content = P_u,
                    conf = conf_u, target = target, sides = "two.sided",
                    seed = control$seed, M = control$M, two_sided = "content"),
    reintervals_warning_unvalidated = function(w) {
      if (isTRUE(w$option %in% c("method", "gpq_two_sided"))) invokeRestart("muffleWarning")
    })
  key <- function(p, P, conf) paste(p, sprintf("%.17g", P), sprintf("%.17g", conf))
  gk <- key(rep_len(seq_along(pts), nrow(g)), g$content, g$conf)
  m <- match(key(match(ev[idx], pts), P_row[idx], out$conf[idx]), gk)
  if (anyNA(m)) {
    .abort("Internal: a rerouted row has no GPQ counterpart.",
           class = "reintervals_error_internal")
  }
  out$method <- meth
  out$method[idx] <- "gpq"
  for (cn in c("estimate", "lower", "upper")) out[[cn]][idx] <- g[[cn]][m]
  if ("singular" %in% names(out)) out$singular[idx] <- g$singular[m]
  if ("mls_bound" %in% names(out)) out$mls_bound[idx] <- NA_character_
  if (!"df" %in% names(out)) out$df <- NA_real_
  out$df[idx] <- g$df[m]
  if (!"gpq_two_sided" %in% names(out)) out$gpq_two_sided <- "not_applicable"
  out$gpq_two_sided[idx] <- g$gpq_two_sided[m]
  for (cn in c("tau2_med", "s2e_med", "frac_boundary", "n_na_draws", "mc_se_lower",
               "mc_se_upper", "K", "df_resid", "phi1", "c_resid", "limitation_flag",
               "route_reason")) {
    v <- g[[cn]]
    out[[cn]] <- v[NA_integer_][rep(1L, nrow(out))]   # NA of the column's type
    out[[cn]][idx] <- v[m]
  }
  out$route_reason[idx] <- "reroute"
  for (a in c("gpq_seed", "rng_kind", "mc_sections", "M")) attr(out, a) <- attr(g, a)
  out
}

## One reintervals_warning_limitation per call when the two-sided MLS
## tolerance interval is in the extrapolation region of the step-13 training
## analysis (R > ti_warn_R and df_T_hat >= ti_warn_dfT). The coverage ranges
## are per-cell coverage on firing replicates at the default thresholds
## (cells with >= 200 firing replicates; ems-mls and anova-mls pooled).
.warn_ti_mls_extrapolation <- function(model, newdata, R, dfT, ev, fired, wR, wD) {
  pts <- sort(unique(ev[fired]))
  tn <- .m1_time_name(model)
  fmt <- function(v) paste(format(signif(v, 3), trim = TRUE), collapse = ", ")
  where <- if (!is.na(tn)) {
    nd <- if (is.null(newdata)) stats::model.frame(model)[1, , drop = FALSE] else newdata
    paste0("t0 = ", fmt(nd[[tn]][pts]))
  } else paste0("evaluation point(s) ", paste(pts, collapse = ", "), " of newdata")
  Rf <- tapply(R[fired], ev[fired], max)
  .warn("Two-sided tolerance interval at ", where, ": the predicted mean is ",
        "uncertain relative to the estimated total variance (R = ", fmt(Rf), "), ",
        "and the total variance is well-determined (df about ",
        fmt(unique(dfT[fired])), "). In simulation the interval covered 0.76-0.88 ",
        "(gamma = 0.90) and 0.84-0.95 (gamma = 0.95) in such settings, below the ",
        "nominal confidence. [provisional until step 13b: step-13 training data, ",
        "thresholds R > 0.25 and df_T_hat >= 40; this call R > ", wR,
        " and df_T_hat >= ", wD, "]",
        class = "reintervals_warning_limitation",
        data = list(limitation = "ti_mls_extrapolation", R = unname(Rf),
                    df_T = unique(dfT[fired]), provisional = TRUE))
}
