# Random-intercept confidence interval, default Satterthwaite df
# (ci_df = "satterthwaite"), against lmerTest (Kuznetsova, Brockhoff &
# Christensen 2017): oracle lmertest-satterthwaite-ci.
#
# lmerTest is independent software: its Satterthwaite df comes from its own
# numerical Hessian of the REML deviance in the (theta, sigma) parameterisation
# and the Jacobian of vcov(beta) there, whereas reintervals differentiates the
# REML log-likelihood in the variance components (sigma^2 scale). The two
# derivatives are taken independently, so the df agree only to numerical
# precision; the SE is sqrt(l vcov(beta) l') on both sides.
#
# Reference call: lmerTest::contest1D(as_lmerModLmerTest(fit), l,
# ddf = "Satterthwaite"), l the fixed-effect row of the mean at the newdata
# point (the intercept, or (1, t0) on the fixed-slope fits).
#
# Fits (all non-singular; singular fits are excluded because the default rule
# then falls back to the containment df by design): the five random-intercept
# fits of fixtures/sim_designs/, two balanced one-way fits (lme4::Dyestuff,
# hk_assay), two unbalanced one-way fits (the same data with rows dropped) and
# one fixed-slope fit (lme4::sleepstudy, Reaction ~ Days + (1 | Subject)).
#
# Measured on 2026-11-26 (lmerTest 3.2.1): max relative difference in df
# 3.6e-5 (hk_assay, whose variance components are of order 1e-3; every other
# fit <= 1e-6), in SE 9e-15. Tolerance: df 1e-4 relative (the agreement two
# independent numerical Hessians support), SE 1e-8 relative.

lmertest_fits <- function() {
  sd <- function(f) utils::read.csv(test_path("fixtures", "sim_designs", f),
                                    stringsAsFactors = TRUE)
  dye <- lme4::Dyestuff
  hk <- hk_assay
  mo <- c(0, 12, 24, 36, 48)
  list(
    OW_BAL    = list(y ~ 1 + (1 | grp), sd("oneway_bal.csv"), NULL),
    OW_UNBAL  = list(y ~ 1 + (1 | grp), sd("oneway_unbal.csv"), NULL),
    NEST_BAL  = list(y ~ 1 + (1 | a) + (1 | a:b), sd("nested_bal.csv"), NULL),
    FS_BAL    = list(y ~ month + (1 | lot), sd("stab_bal.csv"), data.frame(month = mo)),
    FS_UNBAL  = list(y ~ month + (1 | lot), sd("stab_unbal.csv"), data.frame(month = mo)),
    DYE       = list(Yield ~ 1 + (1 | Batch), dye, NULL),
    HK        = list(conc ~ 1 + (1 | run), hk, NULL),
    DYE_UNBAL = list(Yield ~ 1 + (1 | Batch), dye[-c(2, 7, 8, 14, 21, 22, 23), ], NULL),
    HK_UNBAL  = list(conc ~ 1 + (1 | run), hk[-c(3, 6, 9, 17), ], NULL),
    SLEEP_FS  = list(Reaction ~ Days + (1 | Subject), lme4::sleepstudy,
                     data.frame(Days = c(0, 4.5, 9)))
  )
}

# oracle: lmertest-satterthwaite-ci
test_that("Satterthwaite CI df and SE match lmerTest on non-singular random-intercept fits", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("lmerTest")
  fits <- lmertest_fits()
  for (tag in names(fits)) {
    f <- fits[[tag]]
    fm <- lme4::lmer(f[[1]], data = f[[2]], REML = TRUE)
    expect_false(lme4::isSingular(fm), info = tag)
    nd <- f[[3]]
    L <- if (is.null(nd)) matrix(1, 1, 1) else cbind(1, nd[[1]])
    r <- quiet_extrapolation(ci_lmm(fm, newdata = nd))
    dg <- attr(r, "diagnostics")
    se <- (r$upper - r$lower) / 2 / stats::qt(0.975, dg$df)
    lt <- lmerTest::as_lmerModLmerTest(fm)
    ref <- do.call(rbind, lapply(seq_len(nrow(L)), function(i)
      lmerTest::contest1D(lt, L[i, ], ddf = "Satterthwaite")))
    expect_identical(dg$df_method, rep("satterthwaite", nrow(L)), info = tag)
    expect_equal(r$estimate, ref$Estimate, tolerance = 1e-8, info = tag)
    expect_equal(se, ref[["Std. Error"]], tolerance = 1e-8, info = tag)
    expect_equal(dg$df, ref$df, tolerance = 1e-4, info = tag)
  }
})
