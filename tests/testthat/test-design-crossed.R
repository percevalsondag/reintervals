test_that("ems_crossed builds the balanced crossed-with-interaction EMS decomposition", {
  desc <- list(type = "crossed", A = 3, B = 4, n = 2, N = 24, balanced = TRUE,
               grp_A = "a", grp_B = "b", grp_AB = "a:b")
  cv <- c(a = 1.5, b = 1.0, "a:b" = 0.5, Residual = 0.25)
  out <- ems_crossed(desc, cv, "observable")
  n <- 2
  aa <- 3
  bb <- 4
  expect_equal(out$ems$ms[["a:b"]], 0.25 + n * 0.5)
  expect_equal(out$ems$ms[["a"]], 0.25 + n * 0.5 + n * bb * 1.5)
  expect_equal(out$ems$ms[["b"]], 0.25 + n * 0.5 + n * aa * 1.0)
  expect_equal(out$ems$ms[["Residual"]], 0.25)
  expect_equal(out$ems$df[["a"]], 2)                            # A-1
  expect_equal(out$ems$df[["b"]], 3)                            # B-1
  expect_equal(out$ems$df[["a:b"]], 6)                          # (A-1)(B-1)
  expect_equal(out$ems$df[["Residual"]], 12)                    # AB(n-1) = 3*4*1
  expect_equal(sum(out$ems$k * out$ems$ms), sum(cv), tolerance = 1e-10)
})

test_that("classify_design detects a balanced two-crossed-with-interaction design", {
  a <- factor(rep(rep(1:3, each = 4), 2))
  b <- factor(rep(rep(1:4, times = 3), 2))
  ab <- factor(paste(a, b, sep = ":"))
  d <- classify_design(list(a = a, b = b, ab = ab), 24)
  expect_identical(d$type, "crossed")
  expect_equal(d$A, 3)
  expect_equal(d$B, 4)
  expect_equal(d$n, 2)
  expect_true(d$balanced)
})

test_that(".design_components synthesizes ems for unbalanced crossed (M2) when counts given", {
  desc <- list(type = "crossed", balanced = FALSE,
               grp_A = "a", grp_B = "b", grp_AB = "a:b")
  comps <- c(a = 2, b = 1, "a:b" = 0.5, Residual = 1)
  dc <- .design_components(desc, comps, "observable", counts = build_44_counts())
  expect_false(is.null(dc$ems))
  expect_equal(sum(dc$ems$k), 1)                          # synthesis invariant
  # without the count table (edge) it still declines
  expect_null(.design_components(desc, comps, "observable")$ems)
})
