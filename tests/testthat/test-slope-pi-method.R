# lmm_interval_control(slope_pi_method = ): the random-slope prediction
# interval. "reml" (default) is the 0.2.0 / pre-option interval, unchanged;
# "fv-anova" is the fv-anova new-group-mean construction plus the pooled
# within-lot residual variance, with a Satterthwaite df on the two mean squares;
# "hybrid" is "fv-anova" on a singular fit (or a failed singularity check) and
# "reml" otherwise. Rows carry the construction used in `method`, the option in
# diagnostics$slope_pi_arm, and df_fallback_reason per row.

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))
dg <- function(x) attr(x, "diagnostics")
ctl <- function(m) lmm_interval_control(slope_pi_method = m)
## muffle the classes these calls signal by design (the unvalidated-arm warning,
## the singular fit's df fallback, extrapolation); any other warning surfaces
quiet_pi <- function(expr) {
  withCallingHandlers(expr,
    reintervals_warning_unvalidated   = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"))
}

## the singular || fixture of test-vc-cov-reason.R / test-vfun-extract.R
sing_rs <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  lmer_q(y ~ time + (1 | batch) + (0 + time | batch), dn, control = lc)
}
sing_nd <- data.frame(time = c(0, 6, 12, 24))
ss_fit <- function() lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
ss_nd <- data.frame(Days = c(0, 4.5, 9))

## tr(AV)^2 / tr((AV)^2) written out from the per-lot OLS lines
df_lots_ref <- function(y, time, lot, t0) {
  lot <- factor(lot); K <- nlevels(lot); n <- length(y)
  fit <- a <- numeric(K); sse <- 0
  for (k in seq_len(K)) {
    i <- lot == levels(lot)[k]; Z <- cbind(1, time[i]); m <- stats::lm.fit(Z, y[i])
    fit[k] <- sum(c(1, t0) * m$coefficients); sse <- sse + sum(m$residuals^2)
    a[k] <- drop(t(c(1, t0)) %*% solve(crossprod(Z)) %*% c(1, t0))
  }
  S2E <- sse / (n - 2 * K); S2L <- stats::var(fit)
  v <- max(S2L - mean(a) * S2E, 0) + a * S2E
  A <- (diag(K) - matrix(1 / K, K, K)) / (K - 1); AV <- A %*% diag(v)
  sum(diag(AV))^2 / sum(diag(AV %*% AV))
}

## the fv-anova PI recomputed from its formula, off the raw data
fv_pi_ref <- function(y, time, batch, t0, level = 0.95, sides = "two.sided") {
  batch <- factor(batch); K <- nlevels(batch); n <- length(y)
  fit <- vapply(levels(batch), function(b) {
    i <- batch == b; sum(c(1, t0) * stats::coef(stats::lm(y[i] ~ time[i])))
  }, 0)
  sse <- sum(vapply(levels(batch), function(b) {
    i <- batch == b; sum(stats::resid(stats::lm(y[i] ~ time[i]))^2)
  }, 0))
  S2E <- sse / (n - 2 * K); WL <- stats::var(fit) * (1 + 1 / K); W <- WL + S2E
  df <- W^2 / (WL^2 / df_lots_ref(y, time, batch, t0) + S2E^2 / (n - 2 * K))
  q <- if (sides == "two.sided") stats::qt(1 - (1 - level) / 2, df) else stats::qt(level, df)
  list(estimate = mean(fit), lower = mean(fit) - q * sqrt(W),
       upper = mean(fit) + q * sqrt(W), W = W, df = df)
}

test_that("slope_pi_method: match.arg values, default 'reml', error class as the other options", {
  expect_identical(lmm_interval_control()$slope_pi_method, "reml")
  expect_identical(lmm_interval_control(slope_pi_method = "hyb")$slope_pi_method, "hybrid")
  e1 <- tryCatch(lmm_interval_control(slope_pi_method = "nope"), error = identity)
  e2 <- tryCatch(lmm_interval_control(slope_ngm_method = "nope"), error = identity)
  expect_identical(class(e1), class(e2))
  expect_match(conditionMessage(e1), "should be one of", fixed = TRUE)
  expect_match(conditionMessage(e1), "“reml”|\"reml\"")
})

test_that("default and explicit slope_pi_method = 'reml' are identical() (singular and non-singular fits)", {
  skip_if_not_installed("lme4")
  fs <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fs)), "fixture fit is not singular here")
  fits <- list(singular = list(fs, sing_nd), interior = list(ss_fit(), ss_nd))
  for (nm in names(fits)) {
    fm <- fits[[nm]][[1]]; nd <- fits[[nm]][[2]]
    for (sd in list("two.sided", c("lower", "upper"))) {
      a <- quiet_pi(pi_lmm(fm, nd, sides = sd))
      b <- expect_no_warning(quiet_pi(pi_lmm(fm, nd, sides = sd, control = ctl("reml"))))
      expect_identical(a, b, info = nm)
      expect_identical(a$method, rep("reml", nrow(a)), info = nm)
      expect_identical(dg(a)$slope_pi_arm, rep("reml", nrow(a)), info = nm)
    }
  }
})

test_that("random-slope PI rows carry df_fallback_reason (singular: 'singular fit'; interior: NA)", {
  skip_if_not_installed("lme4")
  fs <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fs)), "fixture fit is not singular here")
  p <- quiet_pi(pi_lmm(fs, sing_nd))
  expect_identical(dg(p)$df_rule, rep("containment", 4))
  expect_identical(dg(p)$df_fallback_reason, rep("singular fit", 4))
  ## the reason is also the one the df-fallback warning reports
  w <- tryCatch(pi_lmm(fs, sing_nd), reintervals_warning_df_fallback = function(w) w)
  expect_identical(w$reason, "singular fit")
  q <- quiet_pi(pi_lmm(ss_fit(), ss_nd))
  expect_identical(dg(q)$df_fallback_reason, rep(NA_character_, 3))
  ## vectorised: row-aligned
  v <- quiet_pi(pi_lmm(fs, sing_nd, level = c(0.9, 0.95), sides = c("two.sided", "lower")))
  expect_identical(dg(v)$df_fallback_reason, rep("singular fit", 16))
  ## not added to the random-slope CI or new-group-mean rows
  expect_null(dg(quiet_pi(ci_lmm(fs, sing_nd)))$df_fallback_reason)
  ## a non-singular fit without a covariance: the covariance reason
  local_mocked_bindings(.is_pd_cov = function(M, tol = 1e-8) FALSE)
  r <- quiet_pi(pi_lmm(ss_fit(), ss_nd))
  expect_identical(dg(r)$df_fallback_reason, rep("Hessian not PD or not invertible", 3))
})

test_that("singular fit: 'hybrid' equals 'fv-anova', tagged 'fv-anova', finite, wider than the fv-anova NGM", {
  skip_if_not_installed("lme4")
  fs <- sing_rs()
  skip_if_not(isTRUE(lme4::isSingular(fs)), "fixture fit is not singular here")
  for (sd in c("two.sided", "lower", "upper")) {
    fv <- quiet_pi(pi_lmm(fs, sing_nd, sides = sd, control = ctl("fv-anova")))
    hy <- quiet_pi(pi_lmm(fs, sing_nd, sides = sd, control = ctl("hybrid")))
    ## identical but for the option recorded in slope_pi_arm
    expect_identical(dg(fv)$slope_pi_arm, rep("fv-anova", 4))
    expect_identical(dg(hy)$slope_pi_arm, rep("hybrid", 4))
    a <- fv; b <- hy
    attr(a, "diagnostics")$slope_pi_arm <- NULL; attr(b, "diagnostics")$slope_pi_arm <- NULL
    expect_identical(dg(b)$slope_pi_route_reason, rep(NA_character_, 4))
    attr(b, "diagnostics")$slope_pi_route_reason <- NULL
    expect_identical(a, b, info = sd)
    expect_identical(hy$method, rep("fv-anova", 4))
    expect_identical(dg(hy)$df_rule, rep("fv-anova-satterthwaite", 4))
    expect_identical(dg(hy)$df_fallback_reason, rep(NA_character_, 4))
    expect_true(all(hy$singular))
    lim <- c(if (sd != "upper") hy$lower, if (sd != "lower") hy$upper)
    expect_true(all(is.finite(lim)), info = sd)
    ngm <- quiet_pi(new_group_mean_lmm(fs, sing_nd, sides = sd,
                    control = lmm_interval_control(slope_ngm_method = "fv-anova")))
    expect_identical(hy$estimate, ngm$estimate)            # the same center
    if (sd != "upper") expect_true(all(hy$lower < ngm$lower), info = sd)
    if (sd != "lower") expect_true(all(hy$upper > ngm$upper), info = sd)
    ## W is the NGM arm's W plus the residual variance
    expect_true(all(dg(hy)$W > dg(ngm)$W))
  }
})

test_that("non-singular fit: 'hybrid' equals 'reml' (also without a covariance: containment kept)", {
  skip_if_not_installed("lme4")
  fm <- ss_fit()
  expect_false(lme4::isSingular(fm))
  strip <- function(x) {
    attr(x, "diagnostics")$slope_pi_arm <- NULL
    attr(x, "diagnostics")$slope_pi_route_reason <- NULL; x
  }
  re <- quiet_pi(pi_lmm(fm, ss_nd))
  hy <- quiet_pi(pi_lmm(fm, ss_nd, control = ctl("hybrid")))
  expect_identical(strip(hy), strip(re))
  expect_identical(hy$method, rep("reml", 3))
  expect_identical(dg(hy)$slope_pi_arm, rep("hybrid", 3))
  ## non-singular, covariance unavailable: stays "reml" with the containment df
  local_mocked_bindings(.is_pd_cov = function(M, tol = 1e-8) FALSE)
  re2 <- quiet_pi(pi_lmm(fm, ss_nd))
  hy2 <- quiet_pi(pi_lmm(fm, ss_nd, control = ctl("hybrid")))
  expect_identical(strip(hy2), strip(re2))
  expect_identical(dg(hy2)$df_rule, rep("containment", 3))
  expect_warning(withCallingHandlers(pi_lmm(fm, ss_nd, control = ctl("hybrid")),
                   reintervals_warning_unvalidated = function(w) invokeRestart("muffleWarning")),
                 class = "reintervals_warning_df_fallback")
})

test_that("a failed singularity check is read as singular: 'hybrid' takes 'fv-anova'", {
  skip_if_not_installed("lme4")
  fm <- ss_fit()
  ref <- quiet_pi(pi_lmm(fm, ss_nd, control = ctl("fv-anova")))
  local_mocked_bindings(isSingular = function(x, ...) stop("boom"), .package = "lme4")
  hy <- quiet_pi(pi_lmm(fm, ss_nd, control = ctl("hybrid")))
  expect_identical(hy$method, rep("fv-anova", 3))
  expect_identical(hy$lower, ref$lower)
  expect_identical(hy$upper, ref$upper)
  expect_true(all(is.na(hy$singular)))
})

test_that("'fv-anova' recomputes from its formula (sides, levels, balanced and unbalanced)", {
  skip_if_not_installed("lme4")
  ss <- lme4::sleepstudy
  ub <- ss[!(ss$Subject %in% c("308", "309", "310") & ss$Days > 5), ]
  for (d in list(ss, ub)) {
    fm <- lmer_q(Reaction ~ Days + (Days | Subject), d)
    for (lv in c(0.9, 0.95)) for (sd in c("two.sided", "lower", "upper")) {
      p <- quiet_pi(pi_lmm(fm, ss_nd, level = lv, sides = sd, control = ctl("fv-anova")))
      for (i in seq_len(nrow(ss_nd))) {
        r <- fv_pi_ref(d$Reaction, d$Days, d$Subject, ss_nd$Days[i], level = lv, sides = sd)
        info <- paste(nrow(d), lv, sd, i)
        expect_equal(p$estimate[i], r$estimate, tolerance = 1e-10, info = info)
        expect_equal(dg(p)$W[i], r$W, tolerance = 1e-10, info = info)
        expect_equal(dg(p)$df[i], r$df, tolerance = 1e-10, info = info)
        if (sd != "upper") expect_equal(p$lower[i], r$lower, tolerance = 1e-10, info = info)
        else expect_identical(p$lower[i], -Inf)
        if (sd != "lower") expect_equal(p$upper[i], r$upper, tolerance = 1e-10, info = info)
        else expect_identical(p$upper[i], Inf)
      }
    }
  }
  ## vectorised call = the scalar calls, row by row
  fm <- lmer_q(Reaction ~ Days + (Days | Subject), ss)
  v <- quiet_pi(pi_lmm(fm, ss_nd, level = c(0.9, 0.95), sides = c("lower", "two.sided"),
                       control = ctl("fv-anova")))
  s <- quiet_pi(pi_lmm(fm, ss_nd, level = 0.95, sides = "lower", control = ctl("fv-anova")))
  expect_identical(v$lower[7:9], s$lower)
  expect_identical(v$method, rep("fv-anova", 12))
})

test_that("non-default slope_pi_method warns as a validation arm; no-op elsewhere", {
  skip_if_not_installed("lme4")
  fm <- ss_fit()
  for (m in c("fv-anova", "hybrid"))
    expect_warning(pi_lmm(fm, ss_nd, control = ctl(m)), m, fixed = TRUE,
                   class = "reintervals_warning_unvalidated")
  ## other random-slope verbs: not applied, recorded "not_applicable"
  a <- quiet_pi(ci_lmm(fm, ss_nd, control = ctl("fv-anova")))
  expect_identical(a, quiet_pi(ci_lmm(fm, ss_nd)))
  expect_identical(dg(a)$slope_pi_arm, rep("not_applicable", 3))
  ## random-intercept PI: not applied (pi_df governs it)
  fr <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  b <- expect_no_warning(pi_lmm(fr, ss_nd, control = ctl("fv-anova")))
  expect_identical(b, pi_lmm(fr, ss_nd))
  expect_identical(dg(b)$slope_pi_arm, rep("not_applicable", 3))
  ## lmm_predict passes it through
  lp <- quiet_pi(lmm_predict(fm, "prediction", ss_nd, control = ctl("fv-anova")))
  expect_identical(lp$method, rep("fv-anova", 3))
})

test_that("'fv-anova' refusals: lots with one time point; no residual df (n = 2K)", {
  skip_if_not_installed("lme4")
  ss <- lme4::sleepstudy
  one <- ss[!(ss$Subject == "308" & ss$Days > 0), ]
  fm <- lmer_q(Reaction ~ Days + (Days | Subject), one)
  expect_error(quiet_pi(pi_lmm(fm, ss_nd, control = ctl("fv-anova"))),
               "slope_pi_method = \"fv-anova\" fits a line per lot",
               class = "reintervals_error_time_points")
  two <- ss[ss$Days %in% c(0, 9), ]
  f2 <- lmer_q(Reaction ~ Days + (Days || Subject), two)
  expect_error(quiet_pi(pi_lmm(f2, ss_nd, control = ctl("fv-anova"))),
               "no residual df", class = "reintervals_error_time_points")
  ## the default "reml" still computes on the first fit
  expect_true(all(is.finite(quiet_pi(pi_lmm(fm, ss_nd))$lower)))
})

## ---- effective df of S2_L (fv_df_lots) and the hybrid fallback ---------------

## the unbalanced K = 6 fixture: three lots on 0..24, three on 0..12
unbal_rs <- function(seed = 3) {
  set.seed(seed); full <- c(0, 3, 6, 9, 12, 18, 24); short <- c(0, 3, 6, 9, 12)
  do.call(rbind, lapply(1:6, function(k) {
    tt <- if (k <= 3) full else short
    data.frame(lot = factor(k), time = tt,
               y = 100 + rnorm(1) + (-0.3 + rnorm(1, 0, 0.05)) * tt + rnorm(length(tt)))
  }))
}
test_that("balanced designs: fv_df_lots is K - 1 and the df is the b8ee4e6 expression, bit for bit", {
  skip_if_not_installed("lme4")
  fits <- list(list(ss_fit(), ss_nd, lme4::sleepstudy$Reaction, lme4::sleepstudy$Days,
                    lme4::sleepstudy$Subject))
  fs <- sing_rs()
  if (isTRUE(lme4::isSingular(fs))) {
    mf <- stats::model.frame(fs)
    fits[[2]] <- list(fs, sing_nd, mf$y, mf$time, mf$batch)
  }
  for (f in fits) {
    for (m in c("fv-anova", "hybrid")) {
      p <- quiet_pi(pi_lmm(f[[1]], f[[2]], control = ctl(m)))
      if (p$method[1] != "fv-anova") next          # hybrid on the interior fit is "reml"
      K <- nlevels(factor(f[[5]]))
      expect_identical(dg(p)$fv_df_lots, rep(K - 1, nrow(p)))
      for (i in seq_len(nrow(p))) {
        L <- .lotline(f[[3]], f[[4]], f[[5]], f[[2]][[1]][i])
        WL <- L$S2L * (1 + 1 / L$K); W <- WL + L$S2E
        expect_identical(dg(p)$df[i], W^2 / (WL^2 / (L$K - 1) + L$S2E^2 / L$df_E))
      }
    }
  }
})

test_that("unbalanced design: fv_df_lots < K - 1 and equals tr(AV)^2 / tr((AV)^2)", {
  skip_if_not_installed("lme4")
  d <- unbal_rs()
  fm <- lmer_q(y ~ time + (time | lot), d)
  nd <- data.frame(time = c(24, 12, 6))
  p <- quiet_pi(pi_lmm(fm, nd, control = ctl("fv-anova")))
  for (i in seq_len(nrow(nd))) {
    ref <- df_lots_ref(d$y, d$time, d$lot, nd$time[i])
    expect_equal(dg(p)$fv_df_lots[i], ref, tolerance = 1e-10)
    expect_lt(dg(p)$fv_df_lots[i], 6 - 1)
    ## the PI df uses it; W and the center are those of b8ee4e6
    L <- .lotline(d$y, d$time, d$lot, nd$time[i])
    WL <- L$S2L * (1 + 1 / 6); W <- WL + L$S2E
    expect_identical(dg(p)$W[i], W)
    expect_identical(p$estimate[i], L$Ybar)
    expect_equal(dg(p)$df[i], W^2 / (WL^2 / ref + L$S2E^2 / L$df_E), tolerance = 1e-10)
  }
  ## the edge cases of .fv_df_lots
  L <- list(K = 4, a = c(0.2, 0.5, 0.9, 1.4), abar = 0.75, S2L = 0, S2E = 0)
  expect_identical(.fv_df_lots(L), 3)                    # every v_k = 0 -> K - 1
  L$S2E <- 1                                             # VG_hat < 0: v_k = a_k S2E
  v <- L$a; ref <- sum(v - mean(v))^2 / sum(((diag(4) - 1 / 4) %*% diag(v))^2 * 1)
  M <- (diag(4) - 1 / 4) %*% diag(v)
  expect_equal(.fv_df_lots(L), sum(diag(M))^2 / sum(diag(M %*% M)), tolerance = 1e-12)
  L$a <- rep(0.5, 4); L$abar <- 0.5
  expect_identical(.fv_df_lots(L), 3)                    # equal leverages
})

test_that("Monte Carlo: two-sided 95% coverage of the fv-anova PI, K - 1 vs effective df (3 x 0..24 + 3 x {0,3,6}, t0 = 24)", {
  skip_on_cran()
  full <- c(0, 3, 6, 9, 12, 18, 24); short <- c(0, 3, 6); t0 <- 24; K <- 6
  lot <- factor(rep(seq_len(K), c(7, 7, 7, 3, 3, 3))); tt <- c(rep(full, 3), rep(short, 3))
  cover <- function(s2_0, s2_1, s2e, R = 2000) {
    set.seed(20261125)
    sdn <- sqrt(s2_0 + t0^2 * s2_1 + s2e); mu0 <- 100 - 0.5 * t0
    out <- vapply(seq_len(R), function(r) {
      b0 <- stats::rnorm(K, 0, sqrt(s2_0)); b1 <- stats::rnorm(K, 0, sqrt(s2_1))
      y <- 100 + b0[lot] + (-0.5 + b1[lot]) * tt + stats::rnorm(length(tt), 0, sqrt(s2e))
      p <- .pi_fv_anova(y, tt, lot, t0)
      L <- .lotline(y, tt, lot, t0); WL <- L$S2L * (1 + 1 / K); W <- WL + L$S2E
      hk <- stats::qt(0.975, W^2 / (WL^2 / (K - 1) + L$S2E^2 / L$df_E)) * sqrt(W)
      c(stats::pnorm(L$Ybar + hk, mu0, sdn) - stats::pnorm(L$Ybar - hk, mu0, sdn),
        stats::pnorm(p$upper, mu0, sdn) - stats::pnorm(p$lower, mu0, sdn))
    }, numeric(2))
    rowMeans(out)
  }
  for (s in list(c(1, 0.01, 1), c(1, 0.05, 0.1))) {
    cv <- cover(s[1], s[2], s[3])
    info <- sprintf("s2_0 = %g, s2_1 = %g, s2e = %g: K - 1 %.4f, effective %.4f",
                    s[1], s[2], s[3], cv[1], cv[2])
    expect_gte(cv[2], 0.94)
    expect_gte(cv[2], cv[1])        # df_L <= K - 1: never narrower
    message(info)
  }
})

test_that("hybrid: an infeasible fv-anova PI on a singular fit falls back to 'reml' with a classed warning", {
  skip_if_not_installed("lme4")
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  dn <- dn[!(dn$batch == levels(factor(dn$batch))[1] & dn$time > 0), ]   # one release-only lot
  fm <- lmer_q(y ~ time + (1 | batch) + (0 + time | batch), dn, control = lc)
  skip_if_not(isTRUE(lme4::isSingular(fm)), "fixture fit is not singular here")
  re <- quiet_pi(pi_lmm(fm, sing_nd))
  w <- NULL
  hy <- withCallingHandlers(quiet_pi(pi_lmm(fm, sing_nd, control = ctl("hybrid"))),
    reintervals_warning_hybrid_fallback = function(cnd) { w <<- cnd; invokeRestart("muffleWarning") })
  expect_s3_class(w, c("reintervals_warning_hybrid_fallback", "reintervals_warning_fallback",
                       "reintervals_warning"))
  expect_match(conditionMessage(w), "1 lot(s) with fewer than two distinct time points", fixed = TRUE)
  reason <- "fv-anova infeasible: 1 lot(s) with fewer than two distinct time points"
  expect_identical(w$reason, reason)
  expect_identical(hy$method, rep("reml", 4))
  expect_true(all(is.finite(c(hy$lower, hy$upper))))
  expect_identical(dg(hy)$slope_pi_route_reason, rep(reason, 4))
  expect_identical(hy$lower, re$lower); expect_identical(hy$upper, re$upper)
  ## the explicit arm keeps its error
  expect_error(quiet_pi(pi_lmm(fm, sing_nd, control = ctl("fv-anova"))),
               class = "reintervals_error_time_points")
  ## n = 2K: no residual df -> the same fallback
  d2 <- dn[dn$time %in% c(0, 24), ]
  d2 <- d2[d2$batch %in% names(which(table(d2$batch) == 2)), ]
  f2 <- lmer_q(y ~ time + (1 | batch) + (0 + time | batch), d2, control = lc)
  if (isTRUE(lme4::isSingular(f2))) {
    expect_warning(quiet_pi(pi_lmm(f2, sing_nd, control = ctl("hybrid"))),
                   "no residual df", class = "reintervals_warning_hybrid_fallback")
    h2 <- suppressWarnings(pi_lmm(f2, sing_nd, control = ctl("hybrid")))
    expect_identical(h2$method, rep("reml", 4))
    expect_identical(dg(h2)$slope_pi_route_reason,
                     rep("fv-anova infeasible: no residual df (n = 2K)", 4))
  }
  ## a feasible singular hybrid call: NA route reason
  ok <- quiet_pi(pi_lmm(sing_rs(), sing_nd, control = ctl("hybrid")))
  expect_identical(dg(ok)$slope_pi_route_reason, rep(NA_character_, 4))
})
