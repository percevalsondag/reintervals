# =============================================================================
# New-group-mean validation arms (B2; prereg H3 and H4). Selected by
# lmm_interval_control(ngm_method = , slope_ngm_method = ); the defaults
# ("package", "reml") are the 0.2.0 constructions in slope-intervals.R and are
# not routed through this file.
#
# Random intercept (ngm_method):
#   "broad-inference"  W = V_F + sum of the between-group variance components
#                      (REML, as lme4 returns them), Satterthwaite df
#                      2 W^2 / (g' A g) on the variance scale (SAS PROC MIXED
#                      broad-inference prediction of a new group, DDFM=SATTERTH;
#                      McLean, Sanders & Stroup 1991). No df fallback: a
#                      non-finite or < 1 df gives NA limits, ngm_df_failed.
#   "hts"              mu +/- t_{I-2} sqrt(sigma2_L + V_F) (Higgins, Thompson &
#                      Spiegelhalter 2009, with the REML between-group variance
#                      in place of the DerSimonian-Laird tau^2); intercept-only
#                      one-way, I >= 3.
#   "classical"        Ybar +/- t_{I-1} S_L sqrt(1 + 1/I) on the group means
#                      (Hahn & Meeker 1991); intercept-only one-way, balanced.
# Random slope (slope_ngm_method):
#   "broad-inference"  the default "reml" variance W = V_F + V_G(t0), with the
#                      Satterthwaite df of the whole W, 2 W^2 / (g' A g), g the
#                      gradient of V_F + V_G in the variance-covariance
#                      parameters, A their observed-information covariance
#                      (SAS DDFM=SATTERTH broad inference of a new lot). Same
#                      no-fallback rule as the random-intercept arm.
#   "gpq"              generalized prediction interval on the GPQ pivots
#                      (.gpq_ngm_core, gpq.R).
#   "fv-anova"         Ybar_fit +/- t_{K-1} S_L sqrt(1 + 1/K) on the per-lot OLS
#                      fitted values at t0 (package-original; the lot-line
#                      reference is lotline() in the simstudy pilot).
# =============================================================================

## Does the random-intercept arm `arm` apply to `model`? An arm that does not
## apply is a no-op (the default runs, ngm_arm = "not_applicable").
.ngm_arm_applies <- function(model, arm) {
  if (.has_random_slope(model)) return(FALSE)
  switch(arm,
    package = TRUE,
    ## every random-intercept design the package's NGM accepts
    "broad-inference" = .ngm_uses_anova(model) ||
      isTRUE(tryCatch(.design_of(model)$type %in% c("oneway", "nested", "crossed"),
                      error = function(e) FALSE)),
    hts = .is_oneway_intercept_only(model) &&
      nlevels(lme4::getME(model, "flist")[[1]]) >= 3L,
    classical = .is_oneway_intercept_only(model) && {
      n_i <- table(lme4::getME(model, "flist")[[1]])
      length(n_i) >= 2L && length(unique(as.integer(n_i))) == 1L
    },
    FALSE)
}

## Random-intercept NGM by a non-default arm. Rows and eval columns follow the
## default route (.ngm_ri): the time column for the fixed-slope design, none for
## the intercept-only one-way model, the newdata columns otherwise.
.ngm_ri_arm <- function(model, newdata = NULL, level = 0.95,
                        sides = c("two.sided", "lower", "upper"), arm) {
  sides <- match.arg(sides)
  tcrit <- function(df) if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df)
                        else stats::qt(level, df)
  uses_anova <- .ngm_uses_anova(model)
  time_name <- if (uses_anova) .m1_time_name(model) else NA_character_
  oneway_io <- .is_oneway_intercept_only(model)
  fx <- .fixed(model, if (oneway_io) NULL else newdata)
  method <- if (uses_anova) "anova-mls" else "ems-mls"

  if (arm == "broad-inference") {
    theta <- .theta_cov(model)
    k <- theta$parts$k
    between <- sum(theta$components[-k])     # REML, bounded at 0 by lme4
    a01 <- c(rep(1, k - 1L), 0)              # dW / d(variance components) beyond V_F
    W <- fx$var_fix + between
    df <- vapply(seq_along(W), function(i) {
      if (is.null(theta$cov_theta)) return(NA_real_)
      cg <- .contrast_var_grad(fx$L[i, ], theta$parts)
      if (is.null(cg)) return(NA_real_)
      g <- cg$grad + a01
      v <- as.numeric(t(g) %*% theta$cov_theta %*% g)
      if (!is.finite(v) || v <= 0) return(NA_real_)
      2 * W[i]^2 / v
    }, numeric(1))
    est <- fx$fit
    failed <- !is.finite(df) | df < 1
    rule <- "satterthwaite-broad"
  } else if (arm == "hts") {
    I <- nlevels(lme4::getME(model, "flist")[[1]])
    s2L <- .var_components(model)[[1]]
    W <- s2L + fx$var_fix
    df <- I - 2
    est <- fx$fit
    failed <- FALSE
    rule <- "hts-i-minus-2"
  } else if (arm == "classical") {
    g <- lme4::getME(model, "flist")[[1]]
    gm <- as.numeric(tapply(as.numeric(lme4::getME(model, "y")), g, mean))
    I <- length(gm)
    W <- stats::var(gm) * (1 + 1 / I)
    df <- I - 1
    est <- mean(gm)
    failed <- FALSE
    rule <- "classical-i-minus-1"
  } else {
    .abort("Internal: unknown ngm_method arm '", arm, "'.",
           class = "reintervals_error_internal")
  }
  if (any(!is.finite(W) | W <= 0)) {
    .abort("The new-group-mean variance W is not positive for this fit; the ",
           "ngm_method = \"", arm, "\" interval is not available.",
           class = "reintervals_error_numerical")
  }
  hw <- ifelse(failed, NA_real_, tcrit(pmax(df, 1)) * sqrt(W))
  rows <- data.frame(type = "new_group_mean", estimate = est,
                     lower = if (sides == "upper") -Inf else est - hw,
                     upper = if (sides == "lower")  Inf else est + hw,
                     level = level, df = df, df_rule = rule, W = W, se = sqrt(W),
                     sides = sides, method = method, stringsAsFactors = FALSE)
  if (arm == "broad-inference") {
    rows$ngm_df_failed <- failed
    ## NA limits on both sides when the df failed (no fallback)
    rows$lower[failed] <- NA_real_; rows$upper[failed] <- NA_real_
  }
  out <- if (oneway_io) rows
         else if (!is.na(time_name)) cbind(time = fx$newdata[[time_name]], rows)
         else cbind(fx$newdata, rows)
  rownames(out) <- NULL
  attr(out, "singular") <- .is_singular(model)
  out
}

## Per-lot OLS lines at t0 (base R; the package port of the simstudy pilot's
## lotline(), commit f4ce68e): the fitted values Y_k(t0), their mean Ybar and
## variance S2L, the mean leverage abar and the pooled residual variance S2E.
## Every lot needs two distinct time points.
.lotline <- function(y, time, batch, t0) {
  batch <- factor(batch)
  ok <- !is.na(y) & !is.na(time) & !is.na(batch)
  y <- y[ok]; time <- time[ok]; batch <- droplevels(batch[ok])
  K <- nlevels(batch); n <- length(y)
  if (K < 2L) {
    .abort(sprintf("The fv-anova new-group-mean interval needs >= 2 lots; got %d.", K),
           class = "reintervals_error_groups")
  }
  n_t <- tapply(time, batch, function(v) length(unique(v)))
  few <- names(n_t)[n_t < 2L]
  if (length(few)) {
    .abort(sprintf(paste0("slope_ngm_method = \"fv-anova\" fits a line per lot and ",
                          "needs at least two distinct time points per lot; %d lot(s) ",
                          "have fewer: %s."), length(few),
                   paste(utils::head(few, 5L), collapse = ", ")),
           class = "reintervals_error_time_points", data = list(lots = few))
  }
  Yk <- a <- numeric(K); sse <- 0
  for (k in seq_len(K)) {
    i <- which(as.integer(batch) == k); Z <- cbind(1, time[i])
    ZtZi <- solve(crossprod(Z)); b <- ZtZi %*% crossprod(Z, y[i])
    z0 <- c(1, t0); Yk[k] <- sum(z0 * b); a[k] <- drop(t(z0) %*% ZtZi %*% z0)
    sse <- sse + sum((y[i] - Z %*% b)^2)
  }
  list(K = K, Ybar = mean(Yk), S2L = stats::var(Yk),
       S2E = if (n > 2L * K) sse / (n - 2L * K) else NA_real_, abar = mean(a))
}

## slope_ngm_method = "fv-anova": Ybar_fit +/- t_{K-1} S_L sqrt(1 + 1/K), one
## row per t0. One-sided: t_{1 - alpha}.
.ngm_fv_anova <- function(y, time, batch, t0s, level = 0.95,
                          sides = c("two.sided", "lower", "upper")) {
  sides <- match.arg(sides)
  rows <- lapply(t0s, function(t0) {
    L <- .lotline(y, time, batch, t0)
    W <- L$S2L * (1 + 1 / L$K)
    df <- L$K - 1
    if (!is.finite(W) || W <= 0) {
      .abort("The per-lot fitted values at t0 = ", t0, " are all equal, so the ",
             "fv-anova new-group-mean variance is zero; the interval is not available.",
             class = "reintervals_error_numerical")
    }
    tq <- if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df) else
      stats::qt(level, df)
    hw <- tq * sqrt(W)
    data.frame(time = t0, estimate = L$Ybar,
               lower = if (sides == "upper") -Inf else L$Ybar - hw,
               upper = if (sides == "lower")  Inf else L$Ybar + hw,
               df = df, W = W, se = sqrt(W), K = L$K, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

## Satterthwaite df of the whole random-slope prediction variance W = V_F + V_G
## at t0, for parameters `phi` with covariance `cov` (the observed-information
## inverse at phi). V_F(phi) is parsed$vF_fun; its gradient is numerical, that
## of V_G analytic. Returns c(W, df); df is NA when it cannot be formed.
## Shared by slope_ngm_method = "broad-inference" and its SAS-estimate check.
.ngm_rs_broad_df <- function(parsed, phi, cov, t0, W = NULL) {
  x0 <- stats::setNames(c(1, t0), c("(Intercept)", parsed$fe_time))
  if (is.null(W)) W <- parsed$vF_fun(phi, x0) + .vfun_value_VG(parsed$type, phi, t0)
  gF <- numDeriv::grad(function(p) parsed$vF_fun(p, x0), phi)
  gG <- .vfun_grad_VG(parsed$type, t0)[names(phi)]
  c(W = W, df = .vfun_df_satterthwaite(W, gF + gG, cov))
}

## slope_ngm_method = "broad-inference": the default "reml" rows (same center
## and W, so the same SE) with the df replaced by the Satterthwaite df of W.
## No fallback: a non-finite or < 1 df gives NA limits and ngm_df_failed = TRUE.
.ngm_rs_broad <- function(model, t0s, level, sides) {
  parsed <- .vfun_parse(model)
  rows <- lapply(t0s, function(t0) {
    comp <- vfun_extract(model, t0 = t0)
    row <- .vfun_interval(comp, kind = "new_group_mean", level = level, sides = sides)
    df <- .ngm_rs_broad_df(parsed, comp$phi, comp$cov, t0,
                           W = comp$var_mean + comp$V_G)[["df"]]
    failed <- !is.finite(df) || df < 1
    row$df <- df; row$df_rule <- "satterthwaite-broad"; row$ngm_df_failed <- failed
    if (failed) {
      row$lower <- NA_real_; row$upper <- NA_real_
    } else {
      tq <- if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df) else
        stats::qt(level, df)
      hw <- tq * row$se
      row$lower <- if (sides == "upper") -Inf else row$estimate - hw
      row$upper <- if (sides == "lower")  Inf else row$estimate + hw
    }
    row
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

## Random-slope NGM by a non-default arm, in the slope-path row schema.
.ngm_rs_arm <- function(model, t0s, time_var, level, sides, arm, control) {
  if (arm == "broad-inference") return(.ngm_rs_broad(model, t0s, level, sides))
  yts <- .model_yts(model, time_var)
  singular <- .is_singular(model)   # NA when the check failed
  if (arm == "gpq") {
    g <- .gpq_ngm_core(yts$y, yts$time, yts$batch, t0 = t0s, level = level,
                       sides = sides, M = control$M, seed = control$seed)
    out <- data.frame(time = g$time, type = "new_group_mean", estimate = g$estimate,
                      lower = g$lower, upper = g$upper, level = level,
                      df = NA_real_, df_rule = "gpq", sides = sides,
                      design = .slope_design(model), method = "reml",
                      singular = singular, boundary = NA, vg_floored = NA,
                      tau2_med = g$tau2_med, s2e_med = g$s2e_med,
                      frac_boundary = g$frac_boundary, n_na_draws = g$n_na_draws,
                      mc_se_lower = g$mc_se_lower, mc_se_upper = g$mc_se_upper,
                      K = g$K, df_resid = g$df_resid, phi1 = g$phi1,
                      stringsAsFactors = FALSE)
    for (a in c("gpq_seed", "rng_kind", "mc_sections")) attr(out, a) <- attr(g, a)
    return(out)
  }
  if (arm == "fv-anova") {
    f <- .ngm_fv_anova(yts$y, yts$time, yts$batch, t0s, level = level, sides = sides)
    return(data.frame(time = f$time, type = "new_group_mean", estimate = f$estimate,
                      lower = f$lower, upper = f$upper, se = f$se, level = level,
                      df = f$df, df_rule = "k-minus-1", W = f$W, sides = sides,
                      design = .slope_design(model), method = "reml",
                      singular = singular, boundary = NA, vg_floored = NA,
                      K = f$K, stringsAsFactors = FALSE))
  }
  .abort("Internal: unknown slope_ngm_method arm '", arm, "'.",
         class = "reintervals_error_internal")
}
