# Uses the shared `fixture_52()` (Francq Section 5.2) from helper-fixtures.R.

# Independent re-derivation of Eq. 26 from the EMS slot, used to check that the
# constructor implements the documented formula (technical record Part XI).
ti_by_hand <- function(comp, content, conf) {
  e <- comp$ems
  a <- 1 - conf
  h <- e$df / qchisq(a, e$df) - 1
  inner <- sum(h^2 * e$k^2 * e$ms^2)
  s2t <- sum(e$k * e$ms)
  hw <- qnorm((1 + content) / 2) * sqrt(comp$var_mean + s2t) *
    sqrt(1 + sqrt(inner) / s2t)
  c(lower = comp$mean - hw, upper = comp$mean + hw)
}

test_that("ti_francq implements the documented Eq. 26 (matches re-derivation)", {
  comp <- fixture_52()
  ti <- ti_francq(comp, content = 0.90, conf = 0.95)
  hand <- ti_by_hand(comp, 0.90, 0.95)
  expect_s3_class(ti, "re_interval")
  expect_equal(ti$lower, unname(hand[["lower"]]))
  expect_equal(ti$upper, unname(hand[["upper"]]))
})

test_that("ti_francq reproduces Francq Section 5.2: 95/90 TI ~ [0.845, 1.117]", {
  # "95% TI_90" denotes 95% content (beta) at 90% confidence
  # (gamma) -- the 95/90 pharma convention (maintainer-confirmed). With the
  # source-endorsed formula (a = 1 - conf, z_{(1+content)/2}, lower-tail qchisq)
  # this gives [0.8455, 1.1165], matching the paper's printed [0.845, 1.117] to
  # ~5e-4. The residual is component rounding: the available split nails CI/PI
  # but the TI is sensitive to it (maintainer accepted ~1e-3 tolerance).
  ti <- ti_francq(fixture_52(), content = 0.95, conf = 0.90)
  expect_equal(ti$lower, 0.845, tolerance = 1e-3)
  expect_equal(ti$upper, 1.117, tolerance = 1e-3)
  expect_identical(ti$type, "TI")
  expect_identical(ti$content, 0.95)
  expect_identical(ti$conf, 0.90)
  expect_identical(ti$sides, "two")
  expect_true(is.na(ti$df))
  expect_identical(ti$method, "francq")
})

test_that("ti_francq is symmetric about the estimate", {
  ti <- ti_francq(fixture_52(), content = 0.95, conf = 0.95)
  expect_equal(ti$estimate - ti$lower, ti$upper - ti$estimate)
  expect_equal(ti$estimate, 0.981)
})

test_that("TI is wider than the PI at the same nominal level", {
  comp <- fixture_52()
  ti <- ti_francq(comp, content = 0.95, conf = 0.95)
  pi <- pi_newobs(comp, level = 0.95)
  expect_gt(ti$upper - ti$lower, pi$upper - pi$lower)
})

test_that("ti_francq widens with content and with confidence", {
  comp <- fixture_52()
  w <- function(ct, cf) {
    ti <- ti_francq(comp, content = ct, conf = cf)
    ti$upper - ti$lower
  }
  expect_gt(w(0.99, 0.95), w(0.90, 0.95))   # more content -> wider
  expect_gt(w(0.95, 0.99), w(0.95, 0.90))   # more confidence -> wider
})

test_that("smaller EMS values give a narrower TI", {
  comp <- fixture_52()
  small <- comp
  small$ems$ms <- comp$ems$ms / 4
  ti_full  <- ti_francq(comp,  content = 0.95, conf = 0.95)
  ti_small <- ti_francq(small, content = 0.95, conf = 0.95)
  expect_lt(ti_small$upper - ti_small$lower, ti_full$upper - ti_full$lower)
})

test_that("ti_francq rejects target = 'true_value' for the fixed-slope design (does not silently use observable)", {
  tv <- re_components(
    components = c(run = 0.000681, Residual = 0.001253),
    dfs        = c(ci = 5, pi = 12.484),
    mean = 0.981, var_mean = 0.01353^2,
    coefs  = c(run = 1, Residual = 0),
    target = "true_value",
    ems = list(
      ms = c(run = 0.003296, Residual = 0.001253),
      k  = c(run = 1 / 3,    Residual = 2 / 3),
      df = c(run = 5,        Residual = 12)
    ),
    design = "oneway"
  )
  expect_error(ti_francq(tv), "true_value")
  expect_error(ti_francq(tv), "not supported")
})

test_that("ti_francq errors when no EMS decomposition is available", {
  comp <- re_components(
    components = c(run = 0.000681, Residual = 0.001253),
    dfs        = c(ci = 5, pi = 12.484),
    mean = 0.981, var_mean = 0.01353^2     # ems = NULL
  )
  expect_error(ti_francq(comp), "no expected-mean-square decomposition")
})

test_that("ti_francq rejects one-sided requests and bad levels", {
  comp <- fixture_52()
  expect_error(ti_francq(comp, sides = "lower"), "two-sided")
  expect_error(ti_francq(comp, sides = "upper"), "two-sided")
  expect_error(ti_francq(comp, content = 1.0), "between 0 and 1")
  expect_error(ti_francq(comp, conf = 0), "between 0 and 1")
  expect_error(ti_francq(list(mean = 1)), "re_components")
})
