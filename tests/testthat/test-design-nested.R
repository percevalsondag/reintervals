test_that("ems_nested builds the balanced nested EMS decomposition", {
  desc <- list(type = "nested", A = 4, B = 2, n = 3, N = 24, balanced = TRUE,
               grp_alpha = "a", grp_beta = "b")
  cv <- c(a = 2, b = 1, Residual = 0.5)
  out <- ems_nested(desc, cv, "observable")
  n <- 3
  bb <- 2
  expect_equal(out$ems$ms[["b"]], 0.5 + n * 1)                  # EMS_B
  expect_equal(out$ems$ms[["a"]], 0.5 + n * 1 + n * bb * 2)     # EMS_A
  expect_equal(out$ems$ms[["Residual"]], 0.5)
  expect_equal(out$ems$df[["a"]], 3)                            # A-1
  expect_equal(out$ems$df[["b"]], 4)                            # A(B-1)
  expect_equal(out$ems$df[["Residual"]], 16)                    # AB(n-1)
  expect_equal(sum(out$ems$k * out$ems$ms), sum(cv), tolerance = 1e-10)
})

test_that("classify_design detects a balanced two-fold nested design", {
  # 4 coarse 'a' (6 obs each); 8 fine 'ab' (2 per a, 3 reps each)
  fl <- list(a = factor(rep(1:4, each = 6)), ab = factor(rep(1:8, each = 3)))
  d <- classify_design(fl, 24)
  expect_identical(d$type, "nested")
  expect_equal(d$A, 4)
  expect_equal(d$B, 2)
  expect_equal(d$n, 3)
  expect_true(d$balanced)
})

test_that(".design_components synthesizes ems for unbalanced nested (M2) when counts given", {
  desc <- list(type = "nested", balanced = FALSE, grp_alpha = "a", grp_beta = "b")
  comps <- c(a = 2, b = 1, Residual = 0.5)
  dc <- .design_components(desc, comps, "observable", counts = list(c(2, 1), c(1, 2, 1)))
  expect_false(is.null(dc$ems))
  expect_equal(sum(dc$ems$k), 1)
  expect_equal(unname(dc$coefs), c(1, 1, 1))
  expect_null(.design_components(desc, comps, "observable")$ems)   # no counts -> declines
})
