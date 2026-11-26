# lmm_interval_control(gpq_two_sided): the two-sided GPQ tolerance-interval
# construction. "approx" (default) is Liao, Lin & Iyer (2005) eq. (2) =
# Oliva-Aviles & Hauser (2025) eq. (1), T +/- z_{(1+beta)/2} R_{zeta,gamma};
# "content" (package-original, under validation, step 13b) solves the exact
# content criterion of Liao, Lin & Iyer (2005) eq. (1) draw by draw. The
# default's bit-identity against the pre-option code was checked on 39 saved
# objects outside the test suite (step-13b report); here the default is pinned
# to the eq. (1) formula and to an explicit "approx".

with_warn2 <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd; invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w,
       options = vapply(w, function(x) if (is.null(x$option)) NA_character_ else x$option, ""))
}

## ---- an independent implementation (bisection on the content difference) ----
## h_k: Phi((T + h - R_k)/s_k) - Phi((T - h - R_k)/s_k) = beta; h_k = |T - R_k|
## when v_k = 0. Plain bisection on [0, |T - R_k| + 10 s_k], 300 halvings.
indep_h <- function(T, R, v, beta) {
  d <- abs(T - R); s <- sqrt(v); h <- d; pos <- s > 0
  lo <- rep(0, sum(pos)); hi <- d[pos] + 10 * s[pos]
  cont <- function(x) pnorm((T + x - R[pos]) / s[pos]) - pnorm((T - x - R[pos]) / s[pos])
  for (it in 1:300) {
    mid <- (lo + hi) / 2; up <- cont(mid) >= beta
    hi[up] <- mid[up]; lo[!up] <- mid[!up]
  }
  h[pos] <- (lo + hi) / 2
  h
}
## limits and MC SE of the "content" interval from the engine's draws (kept
## draws only; MC section of each kept draw by its index before the NA drop)
indep_content <- function(dr, beta, gamma, c_resid) {
  R <- dr$Rtheta; v <- pmax(0, dr$Rtau2 + c_resid * dr$Rs2e)
  lim <- function(i) {
    T <- median(R[i]); h <- quantile(indep_h(T, R[i], v[i], beta), gamma, names = FALSE, type = 7)
    c(T - h, T + h)
  }
  L <- lim(seq_along(R))
  secs <- vapply(seq_len(dr$S), function(s) lim(which(dr$sec == s)), numeric(2))
  list(lower = L[1], upper = L[2], mc_se = apply(secs, 1, sd) / sqrt(dr$S))
}

ss <- lme4::sleepstudy

test_that("default is 'approx': identical to an explicit 'approx' and to the eq. (1) formula", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject), ss))
  nd <- data.frame(Days = c(0, 4.5, 9))
  a <- ti_lmm(fm, nd, P = c(0.9, 0.99), sides = c("two.sided", "lower"),
              control = lmm_interval_control(seed = 4, M = 2000))
  b <- ti_lmm(fm, nd, P = c(0.9, 0.99), sides = c("two.sided", "lower"),
              control = lmm_interval_control(seed = 4, M = 2000, gpq_two_sided = "approx"))
  expect_identical(a, b)
  expect_identical(lmm_interval_control()$gpq_two_sided, "approx")
  dg <- attr(a, "diagnostics")
  expect_identical(dg$gpq_two_sided,
                   ifelse(a$sides == "two.sided", "approx", "not_applicable"))
  ## eq. (1) from the draws
  pre <- reintervals:::.gpq_prepare(ss$Reaction, ss$Days, ss$Subject)
  piv <- reintervals:::.gpq_pivots(pre$df_resid, pre$K, 2000L, seed = 4L)
  for (k in seq_len(nrow(a))[a$sides == "two.sided"]) {
    dr <- reintervals:::.gpq_theta_draws(pre, piv, a$Days[k])
    T <- median(dr$Rtheta)
    Rz <- quantile(sqrt(pmax(0, dr$Rtau2 + dr$Rs2e) + dr$Rs2_theta), 0.95, names = FALSE, type = 7)
    z <- qnorm((1 + a$P[k]) / 2)
    expect_equal(c(a$lower[k], a$upper[k]), c(T - z * Rz, T + z * Rz), tolerance = 1e-14)
  }
})

test_that("'content' agrees with an independent implementation on the same draws (1e-10), random slope and random intercept", {
  skip_if_not_installed("lme4")
  M <- 2000L; seed <- 21L
  ## random slope (sleepstudy), both targets
  pre <- reintervals:::.gpq_prepare(ss$Reaction, ss$Days, ss$Subject)
  piv <- reintervals:::.gpq_pivots(pre$df_resid, pre$K, M, seed = seed)
  for (cr in c(1, 0)) for (P in c(0.9, 0.99)) for (t0 in c(0, 4.5, 20)) {
    r <- with_warn2(ti_gpq_raw(ss$Reaction, ss$Days, ss$Subject, t0 = t0, P = P, level = 0.95,
                               c_resid = cr, control = lmm_interval_control(seed = seed, M = M,
                                                                            gpq_two_sided = "content")))$value
    ref <- indep_content(reintervals:::.gpq_theta_draws(pre, piv, t0), P, 0.95, cr)
    info <- sprintf("slope c=%g P=%g t0=%g", cr, P, t0)
    expect_equal(c(r$lower, r$upper), c(ref$lower, ref$upper), tolerance = 1e-10, info = info)
    dg <- attr(r, "diagnostics")
    expect_equal(c(dg$mc_se_lower, dg$mc_se_upper), ref$mc_se, tolerance = 1e-8, info = info)
    expect_identical(dg$gpq_two_sided, "content")
  }
  ## random intercept with a fixed slope (unbalanced, single-observation lots)
  set.seed(5)
  nk <- c(7, 1, 5, 7, 2, 7); tt <- c(0, 3, 6, 9, 12, 18, 24)
  d <- data.frame(lot = factor(rep(seq_along(nk), nk)))
  d$time <- unlist(lapply(nk, function(k) tt[seq_len(k)]))
  d$y <- 100 - 0.5 * d$time + rnorm(6, 0, 0.7)[d$lot] + rnorm(nrow(d))
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (1 | lot), d)))
  pre <- reintervals:::.gpq_prepare_ri(d$y, matrix(d$time, dimnames = list(NULL, "time")), d$lot)
  piv <- reintervals:::.gpq_pivots(pre$df_resid, pre$K, M, seed = seed)
  for (ov in c("observation", "group_mean")) {
    cr <- if (ov == "observation") 1 else 0
    r <- with_warn2(ti_lmm(fm, data.frame(time = c(12, 48)), P = c(0.9, 0.99), level = c(0.9, 0.95),
                           over = ov, method = "gpq",
                           control = lmm_interval_control(seed = seed, M = M,
                                                          gpq_two_sided = "content")))$value
    dg <- attr(r, "diagnostics")
    for (k in seq_len(nrow(r))) {
      ref <- indep_content(reintervals:::.gpq_theta_draws(pre, piv, r$time[k]), r$P[k], r$level[k], cr)
      info <- sprintf("ri %s row %d", ov, k)
      expect_equal(c(r$lower[k], r$upper[k]), c(ref$lower, ref$upper), tolerance = 1e-10, info = info)
      expect_equal(c(dg$mc_se_lower[k], dg$mc_se_upper[k]), ref$mc_se, tolerance = 1e-8, info = info)
    }
  }
})

test_that("the per-draw solver: correct content, monotone in beta, the bracket, and the boundaries", {
  h_of <- reintervals:::.gpq_content_h
  set.seed(8)
  R <- rnorm(500, 0, 2); v <- rexp(500) * 3; T <- 0.3
  d <- abs(T - R); s <- sqrt(v)
  prev <- rep(0, 500)
  for (b in c(0.5, 0.9, 0.95, 0.99, 0.999, 1 - 1e-6, 1 - 1e-10)) {
    h <- h_of(T, R, v, b)
    z <- qnorm((1 + b) / 2)
    ## miss mass (upper-tail form, accurate for beta near 1) equals 1 - beta
    miss <- pnorm(-(d + h) / s) + pnorm((d - h) / s)
    expect_true(all(abs(miss - (1 - b)) <= 1e-9 * (1 - b)), info = sprintf("beta %g", b))
    expect_true(all(h >= z * s * (1 - 1e-12) & h <= (d + z * s) * (1 + 1e-12)), info = sprintf("bracket %g", b))
    expect_true(all(h > prev), info = sprintf("monotone %g", b))
    prev <- h
  }
  ## v = 0: the point mass at R, h = |T - R|; mixed with positive v
  expect_identical(h_of(1, c(3, -2, 1), c(0, 0, 0), 0.95), c(2, 3, 0))
  hm <- h_of(1, c(3, -2), c(0, 4), 0.95)
  expect_identical(hm[1], 2)
  expect_equal(pnorm((1 + hm[2] + 2) / 2) - pnorm((1 - hm[2] + 2) / 2), 0.95, tolerance = 1e-12)
  ## centred draw: h = z_{(1+beta)/2} s; far draw: h = d + z_beta s (one tail)
  expect_equal(h_of(0, 0, 4, 0.95), qnorm(0.975) * 2, tolerance = 1e-14)
  expect_equal(h_of(0, 1e6, 1, 0.95), 1e6 + qnorm(0.95), tolerance = 1e-14)
  ## no NA / non-finite on extreme scales
  expect_true(all(is.finite(h_of(0, c(1e-200, 1e200), c(1e-300, 1e300), 0.99))))
})

test_that("'content' vectorised calls equal scalar calls; seed reproducibility", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject), ss))
  nd <- data.frame(Days = c(2, 12))
  ctl <- lmm_interval_control(seed = 99, M = 2000, gpq_two_sided = "content")
  q <- function(...) suppressWarnings(ti_lmm(fm, nd, ..., control = ctl))
  v <- q(P = c(0.9, 0.99), level = c(0.9, 0.95), sides = c("two.sided", "upper"))
  for (i in seq_len(nrow(v))) {
    s <- q(P = v$P[i], level = v$level[i], sides = v$sides[i])
    k <- match(v$Days[i], s$Days)
    expect_identical(c(v$lower[i], v$upper[i]), c(s$lower[k], s$upper[k]), info = paste("row", i))
    expect_identical(attr(v, "diagnostics")$mc_se_upper[i], attr(s, "diagnostics")$mc_se_upper[k])
  }
  ## the same seed gives identical results; seed = NULL is reproducible from gpq_seed
  expect_identical(q(P = 0.95), q(P = 0.95))
  set.seed(123)
  a <- suppressWarnings(ti_lmm(fm, nd, control = lmm_interval_control(M = 2000, gpq_two_sided = "content")))
  b <- suppressWarnings(ti_lmm(fm, nd, control = lmm_interval_control(
    seed = attr(a, "diagnostics")$gpq_seed, M = 2000, gpq_two_sided = "content")))
  expect_identical(c(a$lower, a$upper), c(b$lower, b$upper))
  ## one-sided rows are unchanged by the option
  ap <- suppressWarnings(ti_lmm(fm, nd, sides = c("lower", "upper"), control = lmm_interval_control(seed = 99, M = 2000)))
  co <- suppressWarnings(ti_lmm(fm, nd, sides = c("lower", "upper"), control = ctl))
  expect_identical(c(ap$lower, ap$upper), c(co$lower, co$upper))
})

test_that("warnings and the gpq_two_sided diagnostics field on every TI route", {
  skip_if_not_installed("lme4")
  fm <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject), ss))
  ctl <- lmm_interval_control(seed = 1, M = 2000, gpq_two_sided = "content")
  ## random slope, two-sided: one gpq_two_sided warning
  r <- with_warn2(ti_lmm(fm, data.frame(Days = 3), control = ctl))
  expect_identical(sum(r$options %in% "gpq_two_sided"), 1L)
  w <- r$warnings[[which(r$options == "gpq_two_sided")]]
  expect_s3_class(w, "reintervals_warning_unvalidated")
  expect_identical(w$value, "content")
  expect_identical(attr(r$value, "diagnostics")$gpq_two_sided, "content")
  ## one-sided only: not applied, no warning
  r <- with_warn2(ti_lmm(fm, data.frame(Days = 3), sides = "lower", control = ctl))
  expect_false("gpq_two_sided" %in% r$options)
  expect_identical(attr(r$value, "diagnostics")$gpq_two_sided, "not_applicable")
  ## the default never warns about it
  r <- with_warn2(ti_lmm(fm, data.frame(Days = 3), control = lmm_interval_control(seed = 1, M = 2000)))
  expect_false("gpq_two_sided" %in% r$options)
  ## ti_gpq_raw
  r <- with_warn2(ti_gpq_raw(ss$Reaction, ss$Days, ss$Subject, t0 = 3, control = ctl))
  expect_identical(sum(r$options %in% "gpq_two_sided"), 1L)
  expect_identical(attr(r$value, "diagnostics")$gpq_two_sided, "content")
  ## random intercept: method = "gpq" applies it (two warnings: method and option);
  ## the closed-form routes record not_applicable and do not warn about it
  fr <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 | Subject), ss))
  r <- with_warn2(ti_lmm(fr, data.frame(Days = 3), method = "gpq", control = ctl))
  expect_setequal(r$options, c("method", "gpq_two_sided"))
  expect_identical(attr(r$value, "diagnostics")$gpq_two_sided, "content")
  for (m in c("auto", "ems-mls")) {
    r <- with_warn2(ti_lmm(fr, data.frame(Days = 3), P = c(0.9, 0.95), method = m, control = ctl))
    expect_false("gpq_two_sided" %in% r$options, info = m)
    expect_identical(unique(attr(r$value, "diagnostics")$gpq_two_sided), "not_applicable", info = m)
  }
  expect_error(lmm_interval_control(gpq_two_sided = "exact"))
})

test_that("'content' on an extrapolated random-intercept data set: wider than eq. (1), same centre (regression pin)", {
  skip_if_not_installed("lme4")
  ## a regression pin on one data set from the step-13 failing G4 cell
  ## (K 12, rho 0.05, full schedule, t0 48): the content interval is wider than eq. (1)
  set.seed(20261126)
  d <- expand.grid(time = c(0, 3, 6, 9, 12, 18, 24), lot = factor(1:12))
  d$y <- 100 - 2 * d$time + rnorm(12, 0, sqrt(0.05 / 0.95))[d$lot] + rnorm(nrow(d))
  fm <- suppressWarnings(suppressMessages(lme4::lmer(y ~ time + (1 | lot), d)))
  f <- function(o) suppressWarnings(ti_lmm(fm, data.frame(time = 48), P = 0.9, method = "gpq",
         control = lmm_interval_control(seed = 3, M = 4000, gpq_two_sided = o)))
  a <- f("approx"); co <- f("content")
  expect_equal(a$estimate, co$estimate)
  expect_gt(co$upper - co$lower, a$upper - a$lower)
  expect_equal(c(a$lower, a$upper), c(0.8576379362, 5.8862383263), tolerance = 1e-9)
  expect_equal(c(co$lower, co$upper), c(0.6039554744, 6.1399207881), tolerance = 1e-9)
})
