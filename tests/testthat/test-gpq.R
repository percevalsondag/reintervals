# A1: GPQ numeric validation against the published Oliva-Aviles & Hauser (2025)
# Section 5 number, and A6: seed-before-call reproducibility.
#
# The published functional is the content variance tau^2(t0) = z' V z at t0 = 20,
# kappa = 1, which equals 8.439388 (the paper's 8.439). We assert it two ways:
#   (1) EXACT -- the closed-form value, reproduced exactly by the engine's V_G
#       machinery (deterministic, platform-robust, tight tolerance);
#   (2) ENGINE-IN-EXPECTATION -- the GPQ Monte-Carlo median recovers it on
#       average over seeded datasets simulated from the DGP. (A single dataset's
#       tau2_med is a noisy estimate, sd ~ 2.6 at K = 30; the MEAN over reps is
#       the engine check, hence the wider stated tolerance.)

# oracle: oliva2025-tau2 (fixtures/oracles.csv)
test_that("A1: reintervals reproduces the published closed-form tau^2 = 8.439 EXACTLY", {
  dgp <- oliva_hauser_dgp
  expect_equal(dgp$tau2_at_t0, 8.439388, tolerance = 1e-5)   # the committed value
  # the engine's own V_G(t0) = z' Sigma z must reproduce it from the DGP phi
  phi <- c(s2_0 = dgp$s2_0, s01 = dgp$s01, s2_1 = dgp$s2_1, s2e = dgp$s2e)
  expect_equal(.vfun_value_VG("M2c", phi, dgp$t0), 8.439388, tolerance = 1e-3)
})

test_that("A1: the GPQ engine recovers tau^2 = 8.439 in expectation (seeded reps)", {
  skip_if_not_installed("lme4")
  dgp <- oliva_hauser_dgp
  R <- 40L
  tau2_med <- vapply(seq_len(R), function(s) {
    d <- vfun_sim(B = 30, times = dgp$times, beta0 = dgp$beta0, beta1 = dgp$beta1,
                  s2_0 = dgp$s2_0, s2_1 = dgp$s2_1, s01 = dgp$s01, s2e = dgp$s2e,
                  seed = 1000 + s)
    g <- quiet_extrapolation(ti_gpq_raw(d$y, d$time, d$batch, t0 = dgp$t0, P = 0.90,
                    level = 0.95, c_resid = 0, control = lmm_interval_control(seed = s, M = 2000)))
    attr(g, "diagnostics")$tau2_med
  }, numeric(1))
  # mean over reps ~ the truth; tolerance 1.5 absorbs per-sample noise + any
  # cross-version RNG drift while still rejecting a broken engine.
  expect_equal(mean(tau2_med), 8.439388, tolerance = 1.5)
})

test_that("A6: GPQ is reproducible under an ambient set.seed() before the call (raw engine)", {
  d <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18), beta0 = 100, beta1 = -0.03,
                s2_0 = 1.5, s2_1 = 0.01, s01 = 0.07, s2e = 1, seed = 3)
  set.seed(123); a <- quiet_extrapolation(ti_gpq_raw(d$y, d$time, d$batch, t0 = 20, c_resid = 0))
  set.seed(123); b <- quiet_extrapolation(ti_gpq_raw(d$y, d$time, d$batch, t0 = 20, c_resid = 0))
  expect_identical(a, b)            # ambient-seed reproducibility (no control seed)
})

test_that("A6: GPQ is reproducible under an ambient set.seed() through ti_lmm(method='gpq')", {
  skip_if_not_installed("lme4")
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  fs <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
                                    data = lme4::sleepstudy, control = lc))
  ctl <- lmm_interval_control(M = 2000)   # reproducibility test, not an oracle: small M
  set.seed(123); a <- ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean", control = ctl)
  set.seed(123); b <- ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean", control = ctl)
  expect_identical(a, b)
})

## ---------------------------------------------------------------------------
## REGRESSION PIN, not an oracle: detects unintended changes to the GPQ draw
## scheme. Call: sleepstudy Reaction ~ Days + (1 + Days | Subject), newdata
## Days = c(0, 5, 9), P = 0.99, level = 0.95, seed = 20261107, M = 2000 (small M
## is fine: this pins the Monte-Carlo stream, it does not check coverage).
## First pinned on dev (A1 tree, 2026-11-07).
## re-pinned in A2: common random numbers (one pivot draw set per call, reused
## for every t0). The Days = 0 values are unchanged (the draw order U, Q, Z is
## the A1 order, and Days = 0 was the first t0, so it saw the same draws); the
## Days = 5 and Days = 9 values changed because under A1 they used the 2nd and
## 3rd per-t0 draw sets. The MC SE (MC sections) consumes no random numbers.
## ---------------------------------------------------------------------------
test_that("REGRESSION PIN (not an oracle): random-slope GPQ TI limits, sleepstudy", {
  skip_if_not_installed("lme4")
  fs <- suppressWarnings(suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy)))
  ctl <- lmm_interval_control(seed = 20261107, M = 2000)
  nd  <- data.frame(Days = c(0, 5, 9))
  pin <- list(
    observation = list(
      two.sided = list(lower = c(132.2476126096, 142.0696865491, 111.3936166104),
                       upper = c(370.7005974070, 465.6071263883, 580.1235511015)),
      lower     = list(lower = c(141.9646754089, 155.1699903682, 131.2061556454),
                       upper = c(Inf, Inf, Inf)),
      upper     = list(lower = c(-Inf, -Inf, -Inf),
                       upper = c(360.0519545146, 451.1151381177, 559.3199285048))),
    group_mean = list(
      two.sided = list(lower = c(152.7000817467, 156.3612792006, 121.5178799511),
                       upper = c(350.2481282699, 451.3155337367, 569.9992877608)),
      lower     = list(lower = c(161.0066973224, 168.1496405247, 139.5336678269),
                       upper = c(Inf, Inf, Inf)),
      upper     = list(lower = c(-Inf, -Inf, -Inf),
                       upper = c(340.6157438629, 438.1431594977, 549.7941022778))))
  for (ov in names(pin)) for (sd in names(pin[[ov]])) {
    r <- ti_lmm(fs, newdata = nd, P = 0.99, level = 0.95, over = ov, sides = sd,
                control = ctl)
    info <- paste(ov, sd)
    expect_identical(r$method, rep("gpq", 3), info = info)
    expect_equal(r$lower, pin[[ov]][[sd]]$lower, tolerance = 1e-8, info = info)
    expect_equal(r$upper, pin[[ov]][[sd]]$upper, tolerance = 1e-8, info = info)
  }
})

## ---- PENDING oracle: oliva2025-worked-example (fixtures/oracles.csv) ----------
## Oliva-Aviles & Hauser (2025) worked example values. To collect: transcribe the
## example data and the printed TI (and any printed pivot summaries) to
## fixtures/oliva2025_worked_example/ (data.csv; limits.csv: t0, P, level, side,
## lower, upper, dp; `#` header with the reference and section, transcriber, date).
## The engine is Monte Carlo, so compare the limits at the larger of the printed
## precision and 3 x the MC SE of each limit (diagnostics mc_se_lower /
## mc_se_upper, A2), at
## M = 20000 with a recorded seed; the paper's own Monte-Carlo error, if printed,
## widens the tolerance accordingly and is recorded in the registry notes.
## test_that("oliva2025-worked-example: GPQ TI reproduces the published example", {
##   fx <- oracle_fixture("oliva2025-worked-example")
##   fm <- lme4::lmer(y ~ t + (1 + t | lot), data = fx$data)
##   for (i in seq_len(nrow(fx$limits))) {
##     L <- fx$limits[i, ]
##     r <- ti_lmm(fm, newdata = data.frame(t = L$t0), P = L$P, level = L$level,
##                 sides = L$side, control = lmm_interval_control(seed = 20261109L))
##     d <- attr(r, "diagnostics")
##     expect_lt(abs(r$lower - L$lower), max(0.5 * 10^-L$dp, 3 * d$mc_se_lower))
##     expect_lt(abs(r$upper - L$upper), max(0.5 * 10^-L$dp, 3 * d$mc_se_upper))
##   }
## })
