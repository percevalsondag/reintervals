# Oracle: random-intercept confidence interval (ci_lmm) against SAS PROC MIXED.
#
# Fixtures: tests/testthat/fixtures/sas_ngm/ (copied from reintervals-simstudy
# d3d8a13, oracles/sas_ngm/part1; see the README there). The oracle rows are the
# ##MEAN lines (OUTPM=, the marginal mean x0'b with its SE and DDFM df) of
# ngm_sas_check_satterth.results.txt.
#
# What is compared, and why:
#   * Satterthwaite (ci_df = "satterthwaite") vs SAS DDFM=SATTERTH on every
#     ##MEAN row: estimate (absolute 1e-6), SE (relative 1e-4), df (absolute 0.05).
#   * Kenward-Roger: reintervals' "kenward-roger" is Kenward & Roger (1997) as
#     implemented in pbkrtest (SE from the adjusted covariance vcovAdj(), df from
#     Lb_ddf()). SAS DDFM=KR is NOT the KR oracle. On these fixtures it reports the
#     Satterthwaite df unchanged, with an adjusted SE, and on the unbalanced
#     one-way fit its adjusted SE also differs from pbkrtest's:
#       OW_UNBAL  SE  0.5217893 (reintervals = pbkrtest)  vs  0.5210774 (SAS KR)
#                 df  4.783     (reintervals = pbkrtest)  vs  5.150     (SAS KR = SAS SATTERTH)
#     The cause of the SE difference has not been established. So:
#       - vs SAS DDFM=KR: the SE only, and only on the BALANCED fits (OW_BAL,
#         NEST_BAL, FS_BAL), where the KR and Satterthwaite df coincide;
#       - vs pbkrtest on all five fits: SE = sqrt(l vcovAdj l') and
#         df = pbkrtest::get_Lb_ddf(model, l), tolerance 1e-8. This guards the
#         pairing of the KR df with the KR-adjusted SE (the 0.2.0 bug paired it
#         with the unadjusted SE).

sas_dir <- function() test_path("fixtures", "sas_ngm")

## the ##MEAN rows (fit, ddfm, t0, pred, se, df, ...); parser in helper-oracles.R
sas_mean_rows <- function() oracle_fixture("sas-part1-ci-satterth")$mean

## PROC MIXED model -> lmer equivalent (see the fixture README)
sas_specs <- list(
  OW_BAL   = list(y ~ 1 + (1 | grp),            "oneway_bal.csv"),
  OW_UNBAL = list(y ~ 1 + (1 | grp),            "oneway_unbal.csv"),
  NEST_BAL = list(y ~ 1 + (1 | a) + (1 | a:b),  "nested_bal.csv"),
  FS_BAL   = list(y ~ month + (1 | lot),        "stab_bal.csv"),
  FS_UNBAL = list(y ~ month + (1 | lot),        "stab_unbal.csv")
)
sas_balanced <- c("OW_BAL", "NEST_BAL", "FS_BAL")

sas_fit <- function(tag) {
  d <- utils::read.csv(file.path(sas_dir(), sas_specs[[tag]][[2]]), stringsAsFactors = TRUE)
  lme4::lmer(sas_specs[[tag]][[1]], data = d, REML = TRUE)
}

## ci_lmm at the SAS rows of one fit; SE recovered from the two-sided 95% width
pkg_ci <- function(fm, t0, ci_df) {
  nd <- if (all(is.na(t0))) NULL else data.frame(month = t0)
  ## the fixtures are evaluated up to month 48, outside some fits' time range
  r  <- quiet_extrapolation(ci_lmm(fm, newdata = nd, control = lmm_interval_control(ci_df = ci_df)))
  df <- attr(r, "diagnostics")$df
  list(estimate = r$estimate, df = df,
       se = (r$upper - r$lower) / 2 / stats::qt(0.975, df),
       df_method = attr(r, "diagnostics")$df_method)
}

## the fixed-effect row l at each SAS row
l_rows <- function(fm, t0) {
  if (all(is.na(t0))) return(matrix(1, nrow = 1))
  cbind(1, t0)
}

# oracle: sas-part1-ci-satterth
test_that("Satterthwaite CI reproduces SAS DDFM=SATTERTH (estimate, SE, df) on every ##MEAN row", {
  skip_if_not_installed("lme4")
  sas <- sas_mean_rows()
  for (tag in names(sas_specs)) {
    s  <- sas[sas$fit == tag & sas$ddfm == "SATTERTH", ]
    expect_gt(nrow(s), 0)
    fm <- sas_fit(tag)
    p  <- pkg_ci(fm, s$t0, "satterthwaite")
    info <- paste(tag, "t0 =", s$t0)
    expect_identical(p$df_method, rep("satterthwaite", nrow(s)), info = tag)
    expect_true(all(abs(p$estimate - s$pred) <= 1e-6), info = paste(info, collapse = "; "))
    expect_true(all(abs(p$se / s$se - 1) <= 1e-4),     info = paste(info, collapse = "; "))
    expect_true(all(abs(p$df - s$df) <= 0.05),         info = paste(info, collapse = "; "))
  }
})

# oracle: sas-part1-ci-kr-balanced
test_that("Kenward-Roger CI SE reproduces SAS DDFM=KR on the balanced fits (SE only)", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  sas <- sas_mean_rows()
  for (tag in sas_balanced) {
    s  <- sas[sas$fit == tag & sas$ddfm == "KR", ]
    expect_gt(nrow(s), 0)
    p  <- pkg_ci(sas_fit(tag), s$t0, "kenward-roger")
    expect_true(all(abs(p$estimate - s$pred) <= 1e-6), info = tag)
    expect_true(all(abs(p$se / s$se - 1) <= 1e-4), info = tag)
  }
})

# oracle: pbkrtest-kr-sas-fits
test_that("Kenward-Roger CI pairs the pbkrtest KR df with the KR-adjusted SE on all five fits", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  sas <- sas_mean_rows()
  for (tag in names(sas_specs)) {
    t0 <- sas$t0[sas$fit == tag & sas$ddfm == "KR"]
    fm <- sas_fit(tag)
    p  <- pkg_ci(fm, t0, "kenward-roger")
    L  <- l_rows(fm, t0)
    vadj <- as.matrix(pbkrtest::vcovAdj(fm))
    se_ref <- sqrt(diag(L %*% vadj %*% t(L)))
    df_ref <- vapply(seq_len(nrow(L)), function(i)
      pbkrtest::get_Lb_ddf(fm, L[i, , drop = FALSE]), numeric(1))
    expect_identical(p$df_method, rep("kenward-roger", nrow(L)), info = tag)
    expect_equal(p$se, se_ref, tolerance = 1e-8, info = tag)
    expect_equal(p$df, df_ref, tolerance = 1e-8, info = tag)
  }
})

# oracle: sas-part1-ci-kr-balanced (the documented OW_UNBAL disagreement)
test_that("documented SAS difference: OW_UNBAL Kenward-Roger (pbkrtest) vs SAS DDFM=KR", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  ## pins the numbers quoted in the header and the fixture README
  s <- sas_mean_rows()
  s <- s[s$fit == "OW_UNBAL", ]
  expect_equal(s$df[s$ddfm == "KR"], s$df[s$ddfm == "SATTERTH"])   # SAS KR df == Satterthwaite df
  expect_equal(s$se[s$ddfm == "KR"], 0.5210774, tolerance = 1e-6)
  p <- pkg_ci(sas_fit("OW_UNBAL"), NA, "kenward-roger")
  expect_equal(p$se, 0.5217893, tolerance = 1e-6)
  expect_equal(p$df, 4.783, tolerance = 1e-3)
})
