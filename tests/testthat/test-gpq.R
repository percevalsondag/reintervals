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
    g <- ti_gpq_raw(d$y, d$time, d$batch, t0 = dgp$t0, P = 0.90, level = 0.95,
                    c_resid = 0, control = lmm_interval_control(seed = s, M = 2000))
    attr(g, "diagnostics")$tau2_med
  }, numeric(1))
  # mean over reps ~ the truth; tolerance 1.5 absorbs per-sample noise + any
  # cross-version RNG drift while still rejecting a broken engine.
  expect_equal(mean(tau2_med), 8.439388, tolerance = 1.5)
})

test_that("A6: GPQ is reproducible under an ambient set.seed() before the call (raw engine)", {
  d <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18), beta0 = 100, beta1 = -0.03,
                s2_0 = 1.5, s2_1 = 0.01, s01 = 0.07, s2e = 1, seed = 3)
  set.seed(123); a <- ti_gpq_raw(d$y, d$time, d$batch, t0 = 20, c_resid = 0)
  set.seed(123); b <- ti_gpq_raw(d$y, d$time, d$batch, t0 = 20, c_resid = 0)
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
## scheme; expected to change deliberately in A2 (common random numbers), and
## that change must be visible in the diff.
## Values computed on dev (A1 tree, 2026-11-07) with the call below: sleepstudy
## Reaction ~ Days + (1 + Days | Subject), newdata Days = c(0, 5, 9), P = 0.99,
## level = 0.95, seed = 20261107, M = 2000 (small M is fine: this pins the
## Monte-Carlo stream, it does not check coverage).
## ---------------------------------------------------------------------------
test_that("REGRESSION PIN (not an oracle): random-slope GPQ TI limits, sleepstudy", {
  skip_if_not_installed("lme4")
  fs <- suppressWarnings(suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy)))
  ctl <- lmm_interval_control(seed = 20261107, M = 2000)
  nd  <- data.frame(Days = c(0, 5, 9))
  pin <- list(
    observation = list(
      two.sided = list(lower = c(132.2476126096, 137.8980914938, 108.7649061215),
                       upper = c(370.7005974070, 470.2129931313, 582.5487016444)),
      lower     = list(lower = c(141.9646754089, 151.2000857935, 129.2630089864),
                       upper = c(Inf, Inf, Inf)),
      upper     = list(lower = c(-Inf, -Inf, -Inf),
                       upper = c(360.0519545146, 454.6152803259, 562.3157002639))),
    group_mean = list(
      two.sided = list(lower = c(152.7000817467, 151.6316182215, 118.4618669310),
                       upper = c(350.2481282699, 456.4794664036, 572.8517408348)),
      lower     = list(lower = c(161.0066973224, 164.7261376006, 138.2779902132),
                       upper = c(Inf, Inf, Inf)),
      upper     = list(lower = c(-Inf, -Inf, -Inf),
                       upper = c(340.6157438629, 441.0548808166, 553.0228005650))))
  for (ov in names(pin)) for (sd in names(pin[[ov]])) {
    r <- ti_lmm(fs, newdata = nd, P = 0.99, level = 0.95, over = ov, sides = sd,
                control = ctl)
    info <- paste(ov, sd)
    expect_identical(r$method, rep("gpq", 3), info = info)
    expect_equal(r$lower, pin[[ov]][[sd]]$lower, tolerance = 1e-8, info = info)
    expect_equal(r$upper, pin[[ov]][[sd]]$upper, tolerance = 1e-8, info = info)
  }
})
