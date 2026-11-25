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
  df <- W^2 / (WL^2 / (K - 1) + S2E^2 / (n - 2 * K))
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
  strip <- function(x) { attr(x, "diagnostics")$slope_pi_arm <- NULL; x }
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
