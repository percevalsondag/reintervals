# ti_lmm(method = "gpq") on random-intercept models (step 13b; under
# validation): the random-intercept sub-class of Oliva-Aviles C, Hauser P
# (2025), Technometrics 67(2):193-202, doi:10.1080/00401706.2024.2407324,
# sections 4.1 and 4.3. Behaviour, eligibility, seeds and a Monte-Carlo sanity
# study. The published checks are in test-oracle-lli2005.R.

gri_q <- function(e) suppressWarnings(suppressMessages(e))

## value of `expr` and the classes of its warnings (each muffled)
gri_warn <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, class(cnd)[1]); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = unique(w))
}

## random intercept + fixed slope, unbalanced, with single-observation lots
gri_data <- function(seed = 3) {
  set.seed(seed)
  nk <- c(5, 1, 4, 7, 1, 3, 6, 2)
  tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- data.frame(lot = factor(rep(seq_along(nk), nk)))
  d$time <- unlist(lapply(nk, function(k) tt[seq_len(k)]))
  d$y <- 100 - 0.3 * d$time + rnorm(length(nk), 0, 1.2)[d$lot] + rnorm(nrow(d), 0, 0.8)
  d
}
gri_ctl <- lmm_interval_control(seed = 42L, M = 2000L)

test_that("gpq on an eligible random-intercept fit: all sides, both targets, GPQ diagnostics, one unvalidated warning", {
  d <- gri_data()
  fm <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d))
  nd <- data.frame(time = c(0, 12, 36))
  for (ov in c("observation", "group_mean")) for (s in c("two.sided", "lower", "upper")) {
    info <- paste(ov, s)
    rw <- gri_warn(ti_lmm(fm, nd, over = ov, sides = s, method = "gpq", control = gri_ctl))
    r <- rw$value
    ## the unvalidated warning (plus the extrapolation warning for t0 = 36)
    expect_setequal(rw$warnings, c("reintervals_warning_unvalidated",
                                   "reintervals_warning_extrapolation"))
    expect_identical(r$method, rep("gpq", 3), info = info)
    expect_identical(r$design, rep("random_intercept_fixed_slope", 3), info = info)
    expect_identical(r$time, nd$time, info = info)
    expect_identical(r$sides, rep(s, 3), info = info)
    dg <- attr(r, "diagnostics")
    expect_true(all(c("gpq_seed", "rng_kind", "mc_se_lower", "mc_se_upper", "M",
                      "n_na_draws", "mc_sections", "phi1", "K", "df_resid",
                      "c_resid", "route_reason") %in% names(dg)), info = info)
    expect_identical(dg$gpq_seed, 42L)
    expect_identical(dg$M, 2000L)
    expect_identical(dg$route_reason, rep("method", 3))
    expect_identical(dg$c_resid, rep(if (ov == "observation") 1 else 0, 3))
    expect_identical(dg$df_resid, rep(nrow(d) - 1L - 8L, 3))
    expect_true(all(r$lower < r$estimate & r$estimate < r$upper), info = info)
    if (s == "lower") expect_true(all(r$upper == Inf & is.na(dg$mc_se_upper)), info = info)
    if (s == "upper") expect_true(all(r$lower == -Inf & is.na(dg$mc_se_lower)), info = info)
  }
  ## the batch-mean target (c_resid = 0) is narrower than the observation target
  a <- gri_q(ti_lmm(fm, nd, method = "gpq", control = gri_ctl))
  b <- gri_q(ti_lmm(fm, nd, over = "group_mean", method = "gpq", control = gri_ctl))
  expect_true(all(b$upper - b$lower < a$upper - a$lower))
  ## the one-way model (q = 0) and a model with a factor covariate (q = 2)
  hk <- gri_q(ti_lmm(lme4::lmer(conc ~ 1 + (1 | run), data = hk_assay), method = "gpq",
                     control = gri_ctl))
  expect_identical(hk$method, "gpq")
  d$site <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  fx <- gri_q(lme4::lmer(y ~ time + site + (1 | lot), data = d))
  rx <- gri_q(ti_lmm(fx, data.frame(time = c(6, 6), site = factor(c("a", "c"), levels = c("a", "b", "c"))), method = "gpq",
                     control = gri_ctl))
  expect_identical(attr(rx, "diagnostics")$df_resid, rep(nrow(d) - 3L - 8L, 2))
  expect_true(all(is.finite(c(rx$lower, rx$upper))))
})

test_that("the default route and the MLS routes are unchanged; group_mean still errors off the gpq route", {
  d <- gri_data()
  fm <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d))
  nd <- data.frame(time = 12)
  au <- gri_q(ti_lmm(fm, nd))
  expect_identical(au$method, "anova-mls")            # not clean balanced
  expect_identical(au, gri_q(ti_lmm(fm, nd, method = "anova-mls")))
  for (m in c("auto", "ems-mls", "anova-mls"))
    expect_error(ti_lmm(fm, nd, over = "group_mean", method = m),
                 class = "reintervals_error_unsupported")
  fb <- gri_q(lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy))
  expect_identical(gri_q(ti_lmm(fb, data.frame(Days = 3)))$method, "ems-mls")
  expect_identical(gri_q(ti_lmm(fit_52()))$method, "anova-mls")
})

test_that("eligibility: one random factor, W of full column rank, positive residual df", {
  set.seed(1)
  ## nested and crossed random factors -> reintervals_error_method
  ne <- data.frame(a = factor(rep(1:4, each = 6)), b = factor(rep(1:12, each = 2)))
  ne$y <- rnorm(4)[ne$a] + rnorm(12)[ne$b] + rnorm(24)
  fn <- gri_q(lme4::lmer(y ~ 1 + (1 | a/b), data = ne))
  expect_error(ti_lmm(fn, method = "gpq"), class = "reintervals_error_method")
  expect_error(ti_lmm(fn, method = "gpq"), "single random factor")
  cr <- expand.grid(r = 1:2, a = factor(1:3), b = factor(1:4))
  cr$y <- rnorm(3)[cr$a] + rnorm(4)[cr$b] + rnorm(nrow(cr))
  fc <- gri_q(lme4::lmer(y ~ 1 + (1 | a) + (1 | b), data = cr))
  expect_error(ti_lmm(fc, method = "gpq"), "single random factor")
  ## a between-lot covariate lies in the column space of Z: condition (b) fails
  d <- gri_data()
  d$grp <- c(0, 1)[as.integer(d$lot) %% 2 + 1]
  fbw <- gri_q(lme4::lmer(y ~ time + grp + (1 | lot), data = d))
  expect_error(ti_lmm(fbw, data.frame(time = 3, grp = 1), method = "gpq"),
               class = "reintervals_error_method")
  expect_error(ti_lmm(fbw, data.frame(time = 3, grp = 1), method = "gpq"), "full column rank")
  ## prior weights
  fw <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d, weights = rep(c(1, 2), length.out = nrow(d))))
  expect_error(ti_lmm(fw, data.frame(time = 3), method = "gpq"), class = "reintervals_error_method")
  ## single-observation lots are allowed; n - q - K < 1 is refused cleanly
  d0 <- data.frame(lot = factor(c(1, 1, 2, 3, 4)), time = c(0, 6, 0, 3, 9),
                   y = c(10, 9.4, 10.3, 9.9, 9.1))
  f0 <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d0))
  expect_error(ti_lmm(f0, data.frame(time = 3), method = "gpq"),
               class = "reintervals_error_residual_df")
  d1 <- rbind(d0, data.frame(lot = factor(1), time = 12, y = 8.9))
  f1 <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d1))
  r1 <- gri_q(ti_lmm(f1, data.frame(time = 3), method = "gpq", control = gri_ctl))
  expect_identical(attr(r1, "diagnostics")$df_resid, 1L)
  expect_true(is.finite(r1$lower) && is.finite(r1$upper))
})

test_that("gpq random intercept: the same seed gives identical() results; a scalar call equals its row of a vectorised call", {
  d <- gri_data()
  fm <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d))
  nd <- data.frame(time = c(0, 18, 30))
  a <- gri_q(ti_lmm(fm, nd, method = "gpq", control = gri_ctl))
  b <- gri_q(ti_lmm(fm, nd, method = "gpq", control = gri_ctl))
  expect_identical(a, b)
  ## seed = NULL: one seed drawn from the caller's stream; set.seed reproduces
  ctl0 <- lmm_interval_control(M = 2000L)
  set.seed(9); n1 <- gri_q(ti_lmm(fm, nd, method = "gpq", control = ctl0))
  set.seed(9); n2 <- gri_q(ti_lmm(fm, nd, method = "gpq", control = ctl0))
  expect_identical(n1, n2)
  ## vectorised: rows t0 fastest, then sides, then P, then level
  P <- c(0.9, 0.99); lev <- c(0.9, 0.95); sd <- c("two.sided", "lower", "upper")
  for (ov in c("observation", "group_mean")) {
    v <- gri_q(ti_lmm(fm, nd, P = P, level = lev, sides = sd, over = ov, method = "gpq",
                      control = gri_ctl))
    vd <- attr(v, "diagnostics")
    expect_identical(nrow(v), 3L * 3L * 2L * 2L)
    i <- 0L
    for (l in lev) for (p in P) for (s in sd) for (k in seq_len(nrow(nd))) {
      i <- i + 1L
      sc <- gri_q(ti_lmm(fm, nd[k, , drop = FALSE], P = p, level = l, sides = s, over = ov,
                         method = "gpq", control = gri_ctl))
      sdg <- attr(sc, "diagnostics")
      info <- paste(ov, l, p, s, k)
      expect_identical(c(v$time[i], v$P[i], v$level[i]), c(nd$time[k], p, l), info = info)
      expect_identical(v$sides[i], s, info = info)
      expect_identical(c(v$estimate[i], v$lower[i], v$upper[i]),
                       c(sc$estimate, sc$lower, sc$upper), info = info)
      expect_identical(c(vd$mc_se_lower[i], vd$mc_se_upper[i]),
                       c(sdg$mc_se_lower, sdg$mc_se_upper), info = info)
    }
  }
})

## Monte-Carlo sanity study (report only; no coverage threshold is asserted).
## 1000 data sets per setting, M = 2000, fixed seeds; two-sided, observation
## target, (P 0.95, gamma 0.90) and (P 0.90, gamma 0.95); gpq versus ems-mls
## (the auto route on these clean balanced designs). Content of an interval:
## Phi((U - mu)/sd) - Phi((L - mu)/sd), sd^2 = sigma_a^2 + sigma_e^2.
gri_mc_cov <- function(K, rho, t0, reps, M = 2000L, seed = 20261125L) {
  tt <- c(0, 3, 6, 9, 12, 18, 24); s2a <- rho / (1 - rho); s2e <- 1
  mu <- function(t) 100 - 0.5 * t
  P <- c(0.95, 0.90); lev <- c(0.90, 0.95)
  set.seed(seed)
  res <- vector("list", reps)
  for (i in seq_len(reps)) {
    d <- expand.grid(time = tt, lot = factor(seq_len(K)))
    d$y <- mu(d$time) + stats::rnorm(K, 0, sqrt(s2a))[d$lot] + stats::rnorm(nrow(d), 0, sqrt(s2e))
    fm <- gri_q(lme4::lmer(y ~ time + (1 | lot), data = d))
    nd <- data.frame(time = t0)
    g <- gri_q(ti_lmm(fm, nd, P = P, level = lev, method = "gpq",
                      control = lmm_interval_control(seed = 1000L + i, M = M)))
    e <- gri_q(ti_lmm(fm, nd, P = P, level = lev, method = "ems-mls"))
    sd_t <- sqrt(s2a + s2e); m0 <- mu(t0)
    cont <- function(x) stats::pnorm((x$upper - m0) / sd_t) - stats::pnorm((x$lower - m0) / sd_t)
    keep <- (g$P == 0.95 & g$level == 0.90) | (g$P == 0.90 & g$level == 0.95)
    res[[i]] <- data.frame(P = g$P[keep], level = g$level[keep],
      cov_gpq = cont(g)[keep] >= g$P[keep], cov_ems = cont(e)[keep] >= e$P[keep],
      w_gpq = (g$upper - g$lower)[keep], w_ems = (e$upper - e$lower)[keep])
  }
  x <- do.call(rbind, res)
  a <- stats::aggregate(cbind(cov_gpq, cov_ems) ~ P + level, x, mean)
  b <- stats::aggregate(cbind(w_gpq, w_ems) ~ P + level, x, stats::median)
  out <- merge(a, b); out$K <- K; out$rho <- rho; out$t0 <- t0; out$reps <- reps
  out
}

test_that("Monte-Carlo sanity study: gpq vs ems-mls coverage and median width (report only)", {
  skip_on_cran()
  r <- rbind(gri_mc_cov(12, 0.05, 48, 1000), gri_mc_cov(3, 0.95, 24, 1000))
  cat("\ngpq (random intercept) vs ems-mls, two-sided, observation target, 1000 reps, M = 2000:\n")
  print(r, row.names = FALSE)
  expect_identical(nrow(r), 4L)
  expect_true(all(r$cov_gpq >= 0 & r$cov_gpq <= 1))
})
