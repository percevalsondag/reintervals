# A single non-intercept fixed term is the M1 (fixed-slope) time design ONLY if
# it maps to a NUMERIC model-frame column. A 2-level factor or a logical
# covariate also yields one non-intercept coefficient ("grpb", "flagTRUE") but
# is not a time design: it must follow the EMS path on all four verbs, with no
# fixed-slope single-observation warning. A numeric covariate is unchanged.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)

## 6 lots x 4 obs (grp alternating a/b within lot); `single_obs` appends 4
## release-only lots with one observation each.
covar_data <- function(single_obs = FALSE, seed = 7) {
  set.seed(seed)
  A <- if (single_obs) 10L else 6L
  lot <- c(rep(seq_len(6), each = 4), if (single_obs) 7:10)
  d <- data.frame(lot = factor(lot),
                  grp = factor(rep_len(c("a", "b"), length(lot)), levels = c("a", "b")),
                  time = c(rep(c(0, 3, 6, 9), 6), if (single_obs) rep(0, 4)))
  d$flag <- d$grp == "b"
  d$y <- 100 + 0.8 * (d$grp == "b") - 0.2 * d$time +
    rnorm(A, 0, 1.5)[d$lot] + rnorm(nrow(d), 0, 0.7)
  d
}

## all four verbs run, route to ems-mls, and emit no warning at all
expect_ems_all_verbs <- function(fm, nd) {
  expect_false(lme4::isSingular(fm))                 # no singular-fit warning either
  expect_false(.m1_not_clean_balanced(fm))
  expect_false(.m1_has_single_obs_lot(fm))
  expect_true(is.na(.m1_time_name(fm)))
  for (verb in list(ci_lmm, pi_lmm, new_group_mean_lmm)) {
    expect_warning(r <- verb(fm, newdata = nd), NA)
    expect_identical(r$method, "ems-mls")
    expect_true(all(is.finite(c(r$lower, r$upper))))
  }
  expect_warning(r <- ti_lmm(fm, newdata = nd), NA)
  expect_identical(r$method, "ems-mls")
  expect_true(all(is.finite(c(r$lower, r$upper))))
}

test_that("2-level factor covariate, all lots >= 2 obs: not M1, EMS on all four verbs", {
  skip_if_not_installed("lme4")
  d <- covar_data(single_obs = FALSE)
  fm <- lme4::lmer(y ~ grp + (1 | lot), data = d, control = lc)
  expect_ems_all_verbs(fm, nd = data.frame(grp = factor("b", levels = c("a", "b"))))
})

test_that("2-level factor covariate with single-observation lots: not M1, EMS on all four verbs", {
  skip_if_not_installed("lme4")
  d <- covar_data(single_obs = TRUE)
  expect_true(any(table(d$lot) < 2L))
  fm <- lme4::lmer(y ~ grp + (1 | lot), data = d, control = lc)
  expect_ems_all_verbs(fm, nd = data.frame(grp = factor("b", levels = c("a", "b"))))
})

test_that("logical covariate: not M1, EMS on all four verbs", {
  skip_if_not_installed("lme4")
  d <- covar_data(single_obs = TRUE)
  fm <- lme4::lmer(y ~ flag + (1 | lot), data = d, control = lc)
  expect_ems_all_verbs(fm, nd = data.frame(flag = TRUE))
})

test_that("anova-mls cannot be forced on a factor-covariate model", {
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(y ~ grp + (1 | lot), data = covar_data(single_obs = TRUE),
                   control = lc)
  expect_error(ti_lmm(fm, method = "anova-mls"), "needs a fixed-slope")
})

test_that("numeric covariate with single-observation lots still routes to anova-mls (unchanged)", {
  skip_if_not_installed("lme4")
  d <- covar_data(single_obs = TRUE)
  fm <- lme4::lmer(y ~ time + (1 | lot), data = d, control = lc)
  expect_identical(.m1_time_name(fm), "time")
  expect_true(.m1_not_clean_balanced(fm))
  expect_true(.m1_has_single_obs_lot(fm))
  nd <- data.frame(time = 6)
  r <- ti_lmm(fm, newdata = nd)
  expect_identical(r$method, "anova-mls")
  expect_warning(ci_lmm(fm, newdata = nd), "anti-conservative")
})

test_that("numeric covariate CONSTANT within groups (between-group) is not M1: EMS on all four verbs, no warning", {
  skip_if_not_installed("lme4")
  # the 0.2.0 validation's V5(b) design: z_i ~ U(0, 10) per group, y ~ z + (1 | g);
  # before the fix it was misclassified as a fixed slope and the new-group mean
  # and TI errored ("no within-lot time variation")
  set.seed(5e8 + 1e7 * 2 + 1e5 * 2 + 1e4 * 2 + 1)
  A <- 5; g <- factor(rep(seq_len(A), each = 4)); a <- rnorm(A); z <- runif(A, 0, 10)
  d <- data.frame(g = g, z = z[g], y = 100 + 0.5 * z[g] + a[g] + rnorm(length(g)))
  fm <- lme4::lmer(y ~ z + (1 | g), data = d, control = lc)
  expect_ems_all_verbs(fm, nd = data.frame(z = 5))
  expect_identical(attr(new_group_mean_lmm(fm, newdata = data.frame(z = 5)), "diagnostics")$df_rule,
                   "top-factor")                     # Gm, the one-way-with-covariates rule
})

test_that("staggered fixed-slope design, every lot >= 2 obs: no CI/PI warning, TI still anova-mls", {
  skip_if_not_installed("lme4")
  set.seed(21)
  times <- list(c(0, 3, 6, 9), c(0, 3, 6), c(1, 4, 7, 10), c(0, 6, 12), c(2, 5, 8), c(0, 3, 6, 9, 12))
  lot <- factor(rep(seq_along(times), lengths(times)))
  d <- data.frame(lot = lot, time = unlist(times))
  d$y <- 100 - 0.5 * d$time + rnorm(nlevels(lot), 0, 1.5)[lot] + rnorm(nrow(d), 0, 0.7)
  fm <- lme4::lmer(y ~ time + (1 | lot), data = d, control = lc)
  expect_identical(.m1_time_name(fm), "time")
  expect_true(.m1_not_clean_balanced(fm))            # routing predicate
  expect_false(.m1_has_single_obs_lot(fm))           # warning predicate
  nd <- data.frame(time = 6)
  expect_warning(ci_lmm(fm, newdata = nd), NA)
  expect_warning(pi_lmm(fm, newdata = nd), NA)
  expect_identical(ti_lmm(fm, newdata = nd)$method, "anova-mls")
  expect_identical(new_group_mean_lmm(fm, newdata = nd)$method, "anova-mls")
})

test_that("MBL (single-observation lots) still warns on CI and PI and routes the TI to anova-mls", {
  skip_if_not_installed("lme4")
  fm <- lme4::lmer(Y ~ Month + (1 | Lot), data = mbl_appendix_a, control = lc)
  nd <- data.frame(Month = 12)
  expect_true(.m1_has_single_obs_lot(fm))
  expect_warning(ci_lmm(fm, newdata = nd), "single-observation lots")
  expect_warning(pi_lmm(fm, newdata = nd), "single-observation lots")
  expect_identical(ti_lmm(fm, newdata = nd)$method, "anova-mls")
})
