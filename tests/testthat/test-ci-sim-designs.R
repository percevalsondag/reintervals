# Random-intercept confidence interval (ci_lmm) on the simulated designs of
# fixtures/sim_designs/ (data and models in the README there): one-way balanced
# and unbalanced, two-fold nested, fixed-slope balanced and staggered.
#
# What is checked, and against what:
#   * Satterthwaite (ci_df = "satterthwaite"): REGRESSION PINS on the package's
#     own estimate, SE and df on every evaluation point. These cases have no
#     external reference value in the repo; the pins guard against unintended
#     changes only.
#   * Kenward-Roger (ci_df = "kenward-roger"): reintervals' "kenward-roger" is
#     Kenward & Roger (1997) as implemented in pbkrtest (SE from the adjusted
#     covariance vcovAdj(), df from Lb_ddf()).
#       - vs pbkrtest on all five fits (oracle pbkrtest-kr-sim-fits):
#         SE = sqrt(l vcovAdj l') and df = pbkrtest::get_Lb_ddf(model, l),
#         tolerance 1e-8. This guards the pairing of the KR df with the
#         KR-adjusted SE (the 0.2.0 bug paired it with the unadjusted SE).
#       - REGRESSION PINS on the package's own estimate and SE on the balanced
#         fits, and on its SE and df on the unbalanced one-way fit (no external
#         reference value in the repo).

sim_dir <- function() test_path("fixtures", "sim_designs")

sim_specs <- list(
  OW_BAL   = list(y ~ 1 + (1 | grp),            "oneway_bal.csv"),
  OW_UNBAL = list(y ~ 1 + (1 | grp),            "oneway_unbal.csv"),
  NEST_BAL = list(y ~ 1 + (1 | a) + (1 | a:b),  "nested_bal.csv"),
  FS_BAL   = list(y ~ month + (1 | lot),        "stab_bal.csv"),
  FS_UNBAL = list(y ~ month + (1 | lot),        "stab_unbal.csv")
)
sim_balanced <- c("OW_BAL", "NEST_BAL", "FS_BAL")
## evaluation points: the intercept (one-way, nested) or months 0-48
sim_t0 <- function(tag) if (startsWith(tag, "FS")) c(0, 12, 24, 36, 48) else NA

sim_fit <- function(tag) {
  d <- utils::read.csv(file.path(sim_dir(), sim_specs[[tag]][[2]]), stringsAsFactors = TRUE)
  lme4::lmer(sim_specs[[tag]][[1]], data = d, REML = TRUE)
}

## ci_lmm at the evaluation points of one fit; SE recovered from the two-sided
## 95% width
pkg_ci <- function(fm, t0, ci_df) {
  nd <- if (all(is.na(t0))) NULL else data.frame(month = t0)
  ## the evaluation points go up to month 48, outside some fits' time range
  r  <- quiet_extrapolation(ci_lmm(fm, newdata = nd, control = lmm_interval_control(ci_df = ci_df)))
  df <- attr(r, "diagnostics")$df
  list(estimate = r$estimate, df = df,
       se = (r$upper - r$lower) / 2 / stats::qt(0.975, df),
       df_method = attr(r, "diagnostics")$df_method)
}

## the fixed-effect row l at each evaluation point
l_rows <- function(fm, t0) {
  if (all(is.na(t0))) return(matrix(1, nrow = 1))
  cbind(1, t0)
}

## Package values (regression pins, not a reference), computed 2026-11-25;
## tolerance 1e-6 relative, as the other pins of the suite.
satterth_pins <- list(
  OW_BAL = list(estimate = 100.5905167, se = 0.7292111667, df = 4.999999944),
  OW_UNBAL = list(estimate = 99.82253922, se = 0.514949035, df = 5.150227412),
  NEST_BAL = list(estimate = 50.7056875, se = 0.4336089642, df = 3.000006595),
  FS_BAL = list(estimate = c(99.70832395, 96.74157934, 93.77483473, 90.80809012, 87.84134551),
    se = c(0.2995745264, 0.2612423776, 0.3270210246, 0.4536853632, 0.6040782688),
    df = c(12.13211213, 7.12762997, 16.66327621, 39.97497621, 52.55420258)),
  FS_UNBAL = list(estimate = c(99.57729627, 96.88963332, 94.20197037, 91.51430742, 88.82664448),
    se = c(0.3051555672, 0.2758635954, 0.3919596431, 0.5707143206, 0.7696415693),
    df = c(11.6722114, 7.885605411, 25.2840532, 43.68616654, 45.91820503)))
kr_balanced_pins <- list(
  OW_BAL = list(estimate = 100.5905167, se = 0.7292111667),
  NEST_BAL = list(estimate = 50.7056875, se = 0.4336089642),
  FS_BAL = list(estimate = c(99.70832395, 96.74157934, 93.77483473, 90.80809012, 87.84134551),
    se = c(0.2995745264, 0.2612423776, 0.3270210246, 0.4536853632, 0.6040782688)))

test_that("Satterthwaite CI on the simulated designs: regression pins (estimate, SE, df)", {
  skip_if_not_installed("lme4")
  for (tag in names(sim_specs)) {
    pin <- satterth_pins[[tag]]
    p  <- pkg_ci(sim_fit(tag), sim_t0(tag), "satterthwaite")
    expect_identical(p$df_method, rep("satterthwaite", length(pin$se)), info = tag)
    expect_equal(p$estimate, pin$estimate, tolerance = 1e-6, info = tag)
    expect_equal(p$se, pin$se, tolerance = 1e-6, info = tag)
    expect_equal(p$df, pin$df, tolerance = 1e-6, info = tag)
  }
})

test_that("Kenward-Roger CI SE on the balanced simulated fits: regression pins", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  for (tag in sim_balanced) {
    pin <- kr_balanced_pins[[tag]]
    p  <- pkg_ci(sim_fit(tag), sim_t0(tag), "kenward-roger")
    expect_equal(p$estimate, pin$estimate, tolerance = 1e-6, info = tag)
    expect_equal(p$se, pin$se, tolerance = 1e-6, info = tag)
  }
})

# oracle: pbkrtest-kr-sim-fits
test_that("Kenward-Roger CI pairs the pbkrtest KR df with the KR-adjusted SE on all five fits", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  for (tag in names(sim_specs)) {
    t0 <- sim_t0(tag)
    fm <- sim_fit(tag)
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

test_that("OW_UNBAL Kenward-Roger CI: regression pin on the package's SE and df", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  ## package values (= pbkrtest; see the test above), pinned
  p <- pkg_ci(sim_fit("OW_UNBAL"), NA, "kenward-roger")
  expect_equal(p$se, 0.5217893, tolerance = 1e-6)
  expect_equal(p$df, 4.783, tolerance = 1e-3)
})
