# The slope-extraction path builds a populated vfun_components from a
# real lmer fit. It is a SEPARATE path from vc_extract (EMS) and is validated
# here against an external oracle (Montes / MBL Appendix A), not just itself.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)

test_that("builds a valid, populated vfun_components from an M2c fit", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18, 24),
                s2_0 = 2.0, s2_1 = 0.01, s01 = -0.05, s2e = 0.7, seed = 3)
  fm <- suppressWarnings(lme4::lmer(y ~ time + (1 + time | batch), data = d,
                                    REML = TRUE, control = lc))
  ve <- vfun_extract(fm, t0 = 24, target = "observable")

  expect_s3_class(ve, "vfun_components")
  expect_identical(ve$type, "M2c")
  expect_identical(names(ve$phi), c("s2_0", "s01", "s2_1", "s2e"))
  expect_true(all(is.finite(c(ve$mean, ve$var_mean, ve$V_G, ve$V_T))))
  expect_true(ve$n_E > 0)
  expect_true(all(is.finite(ve$dfs)) && all(ve$dfs > 0))
  expect_identical(names(ve$dfs), c("F", "T", "G", "e"))
  expect_identical(attr(ve, "target"), "observable")
  # healthy fit -> a 4x4 Cov(phi) reconstructed from the Hessian
  expect_true(is.matrix(ve$cov) && all(dim(ve$cov) == 4L))
})

test_that("V_G is the t0-quadratic z0' Sigma z0", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18, 24),
                s2_0 = 2.0, s2_1 = 0.01, s01 = -0.05, s2e = 0.7, seed = 3)
  fm <- suppressWarnings(lme4::lmer(y ~ time + (1 + time | batch), data = d,
                                    REML = TRUE, control = lc))
  for (t0 in c(0, 12, 24)) {
    ve <- vfun_extract(fm, t0 = t0)
    p <- ve$phi
    hand <- max(0, p[["s2_0"]] + 2 * t0 * p[["s01"]] + t0^2 * p[["s2_1"]])
    expect_equal(ve$V_G, unname(hand))
    expect_equal(ve$V_T, unname(hand + p[["s2e"]]))
  }
})

test_that("ORACLE: balanced M1 matches the Montes route exactly (REML == ANOVA)", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 8, times = c(0, 3, 6, 9, 12), beta0 = 100, beta1 = -0.25,
                s2_0 = 2.0, s2_1 = 0, s01 = 0, s2e = 0.7, seed = 11)
  fm <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch), data = d,
                                    REML = TRUE, control = lc))
  ve <- vfun_extract(fm, t0 = 12, target = "observable")
  mr <- montes_ref(d$y, d$time, d$batch, t0 = 12)
  # on balanced data the REML and ANOVA/MLS routes coincide
  expect_equal(ve$V_T,      mr$Vhat_Y,    tolerance = 1e-4)
  expect_equal(ve$var_mean, mr$Vhat_Yhat, tolerance = 1e-4)
  expect_equal(ve$n_E,      mr$n_E,       tolerance = 1e-4)
})

test_that("ORACLE: montes_ref reproduces MBL Appendix A published TI", {
  d <- mbl_appendix_a
  mr <- montes_ref(d$Y, d$Month, d$Lot, t0 = 12)
  expect_equal(mr$lower, 57.25, tolerance = 0.01)   # published
  expect_equal(mr$upper, 75.41, tolerance = 0.01)
  expect_equal(mr$n_E,   6.729, tolerance = 0.01)
  expect_equal(mr$U,    10.820, tolerance = 0.01)
})

test_that("MBL (unbalanced + single-obs lots): n_E matches Montes; variances REML-narrower", {
  skip_if_not_installed("lme4")
  d <- mbl_appendix_a
  fm <- suppressWarnings(lme4::lmer(Y ~ Month + (1 | Lot), data = d,
                                    REML = TRUE, control = lc))
  ve <- vfun_extract(fm, t0 = 12, target = "observable")
  mr <- montes_ref(d$Y, d$Month, d$Lot, t0 = 12)
  # the scale-invariant ratio agrees closely with the oracle ...
  expect_equal(ve$n_E, mr$n_E, tolerance = 0.01)            # ~0.2% in practice
  # ... while the variance scale is REML-narrower (REML-narrower, ~2-3%):
  expect_lt(ve$V_T, mr$Vhat_Y)                              # narrower, not equal
  expect_lt(abs(ve$V_T - mr$Vhat_Y) / mr$Vhat_Y, 0.05)
  expect_lt(ve$var_mean, mr$Vhat_Yhat)
})

test_that("boundary PREDICATE keys on variances only, never on the signed s01 (deterministic)", {
  # Pure function of (phi, tol) -- no fit, no lme4 optimizer, platform-robust.
  healthy <- c(s2_0 = 2.0, s01 = -0.05, s2_1 = 0.01, s2e = 0.7)
  expect_false(.vfun_is_boundary(healthy, 1e-4))       # all variances comfortably > 0
  # a strongly NEGATIVE s01 must NOT trip it (the inherited-safeguard bug)
  big_neg_cov <- c(s2_0 = 2.0, s01 = -0.13, s2_1 = 0.01, s2e = 0.7)
  expect_false(.vfun_is_boundary(big_neg_cov, 1e-4))
  # a small-but-positive slope variance DOES trip it
  near_bd <- c(s2_0 = 2.0, s01 = 0, s2_1 = 1e-6, s2e = 0.7)
  expect_true(.vfun_is_boundary(near_bd, 1e-3))
  # M1: tiny between-batch variance
  expect_true(.vfun_is_boundary(c(s2_0 = 1e-6, s2e = 5.0), 1e-3))
  expect_false(.vfun_is_boundary(c(s2_0 = 2.0, s2e = 0.7), 1e-3))
})

test_that("boundary flag: healthy fit with negative s01 is NOT flagged near-boundary", {
  skip_if_not_installed("lme4")
  # healthy M2c with NEGATIVE s01: variances large, so boundary must be FALSE
  # regardless of where lme4 lands -- this is robust (not a near-boundary fit).
  dh <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18, 24),
                 s2_0 = 2.0, s2_1 = 0.01, s01 = -0.05, s2e = 0.7, seed = 3)
  fh <- suppressWarnings(lme4::lmer(y ~ time + (1 + time | batch), data = dh,
                                    REML = TRUE, control = lc))
  vh <- vfun_extract(fh, t0 = 24, boundary_tol = 1e-4)
  expect_false(vh$boundary)
  expect_lt(vh$phi[["s01"]], 0)            # negative covariance present, but not flagged
  # the flag the fit produces is exactly what the pure predicate says
  expect_identical(vh$boundary, .vfun_is_boundary(vh$phi, 1e-4))
})

test_that("target selects which content variance n_E is formed from", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 12, times = c(0, 3, 6, 9, 12, 18, 24),
                s2_0 = 2.0, s2_1 = 0.01, s01 = -0.05, s2e = 0.7, seed = 3)
  fm <- suppressWarnings(lme4::lmer(y ~ time + (1 + time | batch), data = d,
                                    REML = TRUE, control = lc))
  vo <- vfun_extract(fm, t0 = 24, target = "observable")
  vt <- vfun_extract(fm, t0 = 24, target = "true_value")
  expect_equal(vo$n_E, vo$V_T / vo$var_mean)
  expect_equal(vt$n_E, vt$V_G / vt$var_mean)
  expect_identical(attr(vt, "target"), "true_value")
})

test_that("singular fit: cov is NULL and df falls back to containment (finite)", {
  skip_if_not_installed("lme4")
  # force a degenerate slope structure -> singular fit (probe-confirmed config)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  fn <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch) + (0 + time | batch),
                                    data = dn, REML = TRUE, control = lc))
  skip_if_not(lme4::isSingular(fn), "fit did not land singular in this environment")
  ve <- vfun_extract(fn, t0 = 12)
  expect_true(ve$singular)
  expect_null(ve$cov)                                  # Hessian route abandoned
  expect_true(all(is.finite(ve$dfs)) && all(ve$dfs > 0))  # containment fallback (safeguard 3)
})

test_that("input guards: one grouping factor, merMod, scalar finite t0", {
  skip_if_not_installed("lme4")
  d <- vfun_sim(B = 8, times = c(0, 3, 6, 9, 12), s2_0 = 2, s2_1 = 0, s2e = 0.7, seed = 5)
  fm <- suppressWarnings(lme4::lmer(y ~ time + (1 | batch), data = d,
                                    REML = TRUE, control = lc))
  expect_error(vfun_extract(fm, t0 = c(0, 12)), "single finite")
  expect_error(vfun_extract(fm, t0 = Inf), "single finite")
  expect_error(vfun_extract(lm(y ~ time, data = d), t0 = 12), "merMod")
})
