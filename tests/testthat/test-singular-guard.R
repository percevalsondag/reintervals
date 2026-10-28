# Without the containment-df guard, a singular fit sends df_PI -> 0 and the PI
# width -> ~10^244. The guard must keep the width finite and flag the event.

test_that("singular fit: PI width stays finite and the fallback attributes are set", {
  skip_if_not_installed("lme4")
  # Force a singular fit: identical run means -> zero between-run variance.
  d <- sqrt(0.001253)
  y <- as.vector(sapply(rep(0.981, 6), function(m) m + c(-d, 0, d)))
  dat <- data.frame(y = y, run = factor(rep(1:6, each = 3)))
  m <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ 1 + (1 | run), data = dat, REML = TRUE)
  ))
  expect_true(lme4::isSingular(m))

  comp <- vc_extract(m)
  expect_true(attr(comp, "singular"))
  expect_true(attr(comp, "pi_df_fallback"))          # df came from the guard
  expect_true(is.finite(comp$dfs[["pi"]]))

  pi <- pi_newobs(comp, level = 0.95)
  width <- pi$upper - pi$lower
  expect_true(is.finite(width))
  expect_lt(width, 1)                                # not the ~10^244 explosion
})

test_that("singular fit: all four verbs report singular = TRUE (CI included)", {
  skip_if_not_installed("lme4")
  d <- sqrt(0.001253)
  y <- as.vector(sapply(rep(0.981, 6), function(m) m + c(-d, 0, d)))
  dat <- data.frame(y = y, run = factor(rep(1:6, each = 3)))
  m <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ 1 + (1 | run), data = dat, REML = TRUE)
  ))
  expect_true(lme4::isSingular(m))

  expect_true(all(ci_lmm(m)$singular))
  expect_true(all(pi_lmm(m)$singular))
  # new_group_mean: every run mean here is exactly 0.981, so the unweighted
  # between-group mean square is 0 and the anova-mls W is 0 -> a clear error,
  # not a zero-width interval. On an ordinary singular fit (MSA < MSE) it
  # returns an interval flagged singular.
  expect_error(new_group_mean_lmm(m), "not positive")
  set.seed(4); g <- factor(rep(1:5, each = 3))
  m2 <- suppressMessages(lme4::lmer(y ~ 1 + (1 | g), data = data.frame(g = g, y = 10 + rnorm(15))))
  expect_true(lme4::isSingular(m2))
  expect_true(all(new_group_mean_lmm(m2)$singular))
  expect_true(all(suppressWarnings(ti_lmm(m))$singular))   # TI also warns on singular
})
