# Without the bounded-Satterthwaite guard, a singular fit sends df_PI -> 0 and the PI
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
