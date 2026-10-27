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
  expect_false(.is_single_obs_m1(fm))
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
  expect_true(.is_single_obs_m1(fm))
  nd <- data.frame(time = 6)
  r <- ti_lmm(fm, newdata = nd)
  expect_identical(r$method, "anova-mls")
  expect_warning(ci_lmm(fm, newdata = nd), "anti-conservative")
})
