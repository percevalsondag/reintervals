# The vfun_components container: these tests pin
# its construction/validation contract -- above all the sign rule that
# distinguishes it from the frozen re_components: variances are non-negative,
# but the intercept-slope covariance s01 MAY be negative (never truncated).

mk_m2c <- function(s01 = -0.05, ...) {
  s2_0 <- 2.0; s2_1 <- 0.01; s2e <- 0.7; t0 <- 24
  vg <- s2_0 + 2 * t0 * s01 + t0^2 * s2_1
  args <- list(
    type = "M2c",
    phi  = c(s2_0 = s2_0, s01 = s01, s2_1 = s2_1, s2e = s2e),
    mean = 98.5, var_mean = 0.12,
    V_G = vg, V_T = vg + s2e, n_E = 30,
    dfs = c(F = 8, T = 6.4, G = 5.1, e = 40), t0 = t0
  )
  modifyList(args, list(...))
}

test_that("constructs a valid M2c payload (baseline)", {
  vc <- do.call(vfun_components, mk_m2c())
  expect_s3_class(vc, "vfun_components")
  expect_identical(vc$type, "M2c")
  expect_identical(names(vc$phi), c("s2_0", "s01", "s2_1", "s2e"))
  expect_false(vc$singular)
  expect_null(vc$cov)
})

test_that("NEGATIVE (PSD-valid) s01 covariance is ACCEPTED (the point of replace)", {
  vc <- do.call(vfun_components, mk_m2c(s01 = -0.05))
  expect_identical(unname(vc$phi[["s01"]]), -0.05)
  # a strongly negative covariance that still keeps Sigma PSD
  # (corr = -0.13/sqrt(2*0.01) = -0.92): valid, and V_G stays >= 0.
  vc2 <- do.call(vfun_components, mk_m2c(s01 = -0.13))
  expect_identical(unname(vc2$phi[["s01"]]), -0.13)
  expect_gte(vc2$V_G, 0)
})

test_that("non-PD / boundary Sigma: a negative V_G is CARRIED, not rejected", {
  # s01 beyond the PSD bound (|corr|>1) makes V_G(t0)=z0'Sigma z0 < 0 at some t0.

  vc <- do.call(vfun_components, mk_m2c(s01 = -0.30))   # V_G = -6.64, V_T = -5.94
  expect_lt(vc$V_G, 0)
  expect_lt(vc$V_T, 0)
  expect_identical(unname(vc$phi[["s01"]]), -0.30)
  # NA / Inf are still rejected -- finiteness is the retained invariant.
  expect_error(do.call(vfun_components, mk_m2c(V_G = NA_real_)), "non-NA")
  expect_error(do.call(vfun_components, mk_m2c(V_T = Inf)), "finite")
})

test_that("negative VARIANCE components are REJECTED (not truncated)", {
  expect_error(
    do.call(vfun_components, mk_m2c(phi = c(s2_0 = -1, s01 = -0.05,
                                            s2_1 = 0.01, s2e = 0.7))),
    "non-negative"
  )
  expect_error(
    do.call(vfun_components, mk_m2c(phi = c(s2_0 = 2, s01 = -0.05,
                                            s2_1 = 0.01, s2e = -0.7))),
    "non-negative"
  )
})

test_that("M1 and M2i canonical phi layouts validate", {
  m1 <- vfun_components(
    type = "M1", phi = c(s2_0 = 0.5, s2e = 0.5),
    mean = 100, var_mean = 0.05, V_G = 0.5, V_T = 1.0, n_E = 20,
    dfs = c(F = 5, T = 7), t0 = 12
  )
  expect_identical(m1$type, "M1")

  m2i <- vfun_components(
    type = "M2i", phi = c(s2_0 = 0.5, s2_1 = 0.01, s2e = 0.5),
    mean = 100, var_mean = 0.05, V_G = 0.5, V_T = 1.0, n_E = 20,
    dfs = c(F = 5, T = 7, G = 6), t0 = 12
  )
  expect_identical(m2i$type, "M2i")
})

test_that("phi names/order must match the type", {
  # wrong order
  expect_error(
    do.call(vfun_components, mk_m2c(phi = c(s2_0 = 2, s2_1 = 0.01,
                                            s01 = -0.05, s2e = 0.7))),
    "names must be exactly"
  )
  # M1 phi handed to M2c
  expect_error(
    do.call(vfun_components, mk_m2c(phi = c(s2_0 = 2, s2e = 0.7))),
    "names must be exactly"
  )
  expect_error(
    vfun_components(type = "BOGUS", phi = c(s2_0 = 1, s2e = 1),
                    mean = 0, var_mean = 0, V_G = 0, V_T = 0, n_E = 1,
                    dfs = c(F = 1), t0 = 0),
    "single recognized model code"
  )
})

test_that("cov: NULL allowed, dimension and symmetry enforced", {
  good <- diag(4)
  expect_s3_class(do.call(vfun_components, mk_m2c(cov = good)), "vfun_components")
  # wrong dimension (3x3 against length-4 phi)
  expect_error(do.call(vfun_components, mk_m2c(cov = diag(3))), "to match `phi`")
  # non-symmetric
  ns <- diag(4); ns[1, 2] <- 5
  expect_error(do.call(vfun_components, mk_m2c(cov = ns)), "symmetric")
  # carries singular flag through
  vc <- do.call(vfun_components, mk_m2c(cov = NULL, singular = TRUE))
  expect_true(vc$singular)
})

test_that("scalar payload rules: var_mean>=0 enforced; V_G/V_T finite (signed); n_E>0 with Inf ok", {
  expect_error(do.call(vfun_components, mk_m2c(var_mean = -0.01)), "non-negative")  # V_F structural
  # V_G/V_T are signed carriers now: a negative value is ACCEPTED (not an error)
  expect_s3_class(do.call(vfun_components, mk_m2c(V_G = -1, V_T = -0.3)), "vfun_components")
  expect_error(do.call(vfun_components, mk_m2c(n_E = 0)), "positive")
  # n_E = Inf (V_F = 0 limit) is allowed
  expect_s3_class(do.call(vfun_components, mk_m2c(n_E = Inf)), "vfun_components")
})

test_that("dfs contract: named, positive, no NA, Inf allowed", {
  expect_error(do.call(vfun_components, mk_m2c(dfs = c(8, 6))), "named")
  expect_error(do.call(vfun_components, mk_m2c(dfs = c(F = NA, T = 6))), "NA")
  expect_error(do.call(vfun_components, mk_m2c(dfs = c(F = 0, T = 6))), "positive")
  expect_s3_class(do.call(vfun_components, mk_m2c(dfs = c(F = Inf, T = 6))),
                  "vfun_components")
})

test_that("boundary flag: defaults FALSE, carried through, validated logical", {
  expect_false(do.call(vfun_components, mk_m2c())$boundary)        # default
  expect_true(do.call(vfun_components, mk_m2c(boundary = TRUE))$boundary)
  expect_error(do.call(vfun_components, mk_m2c(boundary = NA)),
               "single TRUE/FALSE")
  expect_error(do.call(vfun_components, mk_m2c(boundary = "yes")),
               "single TRUE/FALSE")
  # singular and boundary are independent flags
  vc <- do.call(vfun_components, mk_m2c(singular = TRUE, boundary = FALSE))
  expect_true(vc$singular); expect_false(vc$boundary)
})

test_that("print returns its input invisibly and shows the NULL-cov note", {
  vc <- do.call(vfun_components, mk_m2c(cov = NULL, singular = TRUE))
  expect_output(print(vc), "vfun_components")
  expect_output(print(vc), "NULL")
  expect_invisible(print(vc))
})
