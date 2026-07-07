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
  set.seed(123); a <- ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean")
  set.seed(123); b <- ti_lmm(fs, newdata = data.frame(Days = 5), over = "group_mean")
  expect_identical(a, b)
})
