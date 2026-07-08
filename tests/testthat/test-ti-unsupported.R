# ti_lmm on a design with no closed-form EMS decomposition (crossed factors with
# no interaction term; >= 4 grouping factors) must signal a clean error -- NOT
# return NA bounds. The coverage verbs (CI / PI / new_group_mean), which do not
# need the EMS decomposition, must stay available with finite bounds.

lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)

# Two crossed factors, NO interaction term -> classify_design() = "crossed-no-
# interaction" -> dc$ems is NULL (from the scenario harness construction).
crossed_noint_fit <- function() {
  set.seed(1)
  g <- expand.grid(a = factor(1:3), b = factor(1:4), rep = 1:3)
  g$y <- rnorm(3, 0, 2)[g$a] + rnorm(4, 0, 2)[g$b] + rnorm(36)
  suppressWarnings(lme4::lmer(y ~ 1 + (1 | a) + (1 | b), g, control = lc))
}

test_that("ti_lmm errors (no NA bounds) on a crossed-no-interaction design", {
  fit <- crossed_noint_fit()
  expect_error(ti_lmm(fit), "not available")
})

test_that("CI / PI / new_group_mean stay available with finite bounds there", {
  fit <- crossed_noint_fit()
  for (verb in list(ci_lmm, pi_lmm, new_group_mean_lmm)) {
    out <- verb(fit)
    expect_true(all(is.finite(out$lower)))
    expect_true(all(is.finite(out$upper)))
  }
})
