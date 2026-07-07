test_that("ems_oneway reproduces the Section 5.2 balanced EMS decomposition", {
  desc <- list(type = "oneway", A = 6, N = 18, n0 = 3, balanced = TRUE, grp = "run")
  out <- ems_oneway(desc, c(run = 0.000681, Residual = 0.001253), "observable")
  expect_equal(out$ems$ms[["run"]], 0.001253 + 3 * 0.000681)   # EMS_A = 0.003296
  expect_equal(out$ems$ms[["Residual"]], 0.001253)
  expect_equal(unname(out$ems$k), c(1 / 3, 2 / 3))
  expect_equal(unname(out$ems$df), c(5, 12))
  expect_equal(sum(out$ems$k * out$ems$ms), 0.001934)          # = total variance
  expect_equal(unname(out$coefs), c(1, 1))
})

test_that("true_value target zeroes only the residual coef", {
  desc <- list(type = "oneway", A = 6, N = 18, n0 = 3, grp = "run")
  out <- ems_oneway(desc, c(run = 0.000681, Residual = 0.001253), "true_value")
  expect_equal(out$coefs[["Residual"]], 0)
  expect_equal(out$coefs[["run"]], 1)
})

test_that("ems_oneway feeds ti_francq to reproduce the Section 5.2 TI", {
  desc <- list(type = "oneway", A = 6, N = 18, n0 = 3, grp = "run")
  out <- ems_oneway(desc, c(run = 0.000681, Residual = 0.001253), "observable")
  comp <- re_components(
    components = c(run = 0.000681, Residual = 0.001253),
    dfs = c(ci = 5, pi = 12.484), mean = 0.981, var_mean = 0.01353^2,
    ems = out$ems, design = "oneway"
  )
  ti <- ti_francq(comp, content = 0.95, conf = 0.90)
  expect_equal(ti$lower, 0.845, tolerance = 1e-3)
  expect_equal(ti$upper, 1.117, tolerance = 1e-3)
})

test_that("classify_design detects a balanced one-way", {
  d <- classify_design(list(run = factor(rep(1:6, each = 3))), 18)
  expect_identical(d$type, "oneway")
  expect_equal(d$A, 6)
  expect_equal(d$n0, 3)
  expect_true(d$balanced)
})

test_that("classify_design computes n0 for the Section 5.3 unbalanced one-way", {
  # runs 1,2,6 -> 2 reps; runs 3,4,5 -> 3 reps; N = 15
  fl <- list(run = factor(c(rep(1, 2), rep(2, 2), rep(3, 3),
                            rep(4, 3), rep(5, 3), rep(6, 2))))
  d <- classify_design(fl, 15)
  expect_identical(d$type, "oneway")
  expect_false(d$balanced)
  expect_equal(round(d$n0, 2), 2.48)            # paper Section 5.3 n0
})
