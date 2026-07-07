test_that("re_components builds a valid container with expected fields", {
  comp <- re_components(
    components = c(run = 0.00046, Residual = 0.00147),
    dfs        = c(run = 5, Residual = 12, ci = 5, pi = 12.484),
    mean       = 0.981,
    var_mean   = 0.01353^2,
    estimator  = "reml",
    design     = "oneway"
  )
  expect_s3_class(comp, "re_components")
  expect_identical(comp$mean, 0.981)
  expect_identical(comp$var_mean, 0.01353^2)
  expect_identical(comp$estimator, "reml")
  expect_identical(comp$design, "oneway")
  expect_named(comp$components, c("run", "Residual"))
})

test_that("coefs defaults to all-ones (observable target) and aligns to components", {
  comp <- re_components(
    components = c(a = 1, Residual = 2),
    dfs        = c(a = 3, Residual = 4, ci = 3, pi = 5),
    mean = 0, var_mean = 1
  )
  expect_identical(comp$coefs, c(a = 1, Residual = 1))

  # supplied out of order -> reordered to components order
  comp2 <- re_components(
    components = c(a = 1, Residual = 2),
    dfs        = c(a = 3, Residual = 4, ci = 3, pi = 5),
    coefs      = c(Residual = 0, a = 1),     # "true_value" target
    mean = 0, var_mean = 1
  )
  expect_identical(comp2$coefs, c(a = 1, Residual = 0))
})

test_that("negative ANOVA component estimates are truncated at 0 with a warning", {
  expect_warning(
    comp <- re_components(
      components = c(a = -0.5, Residual = 2),
      dfs        = c(a = 3, Residual = 4, ci = 3, pi = 5),
      mean = 0, var_mean = 1
    ),
    "truncated at 0"
  )
  expect_identical(unname(comp$components[["a"]]), 0)
})

test_that("re_components rejects malformed input", {
  ok_dfs <- c(a = 3, Residual = 4, ci = 3, pi = 5)
  expect_error(
    re_components(components = c(1, 2), dfs = ok_dfs, mean = 0, var_mean = 1),
    "named numeric"
  )
  expect_error(
    re_components(components = c(a = Inf, Residual = 1), dfs = ok_dfs,
                  mean = 0, var_mean = 1),
    "finite"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2),
                  dfs = c(a = 3, Residual = -1, ci = 3, pi = 5),
                  mean = 0, var_mean = 1),
    "positive"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2),
                  dfs = c(a = 3, Residual = NA, ci = 3, pi = 5),
                  mean = 0, var_mean = 1),
    "NA is not permitted"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2), dfs = ok_dfs,
                  coefs = c(b = 1, Residual = 0), mean = 0, var_mean = 1),
    "names must match"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2), dfs = ok_dfs,
                  mean = c(0, 1), var_mean = 1),
    "single finite numeric"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2), dfs = ok_dfs,
                  mean = 0, var_mean = -1),
    "non-negative"
  )
  expect_error(
    re_components(components = c(a = 1, Residual = 2), dfs = ok_dfs,
                  mean = 0, var_mean = 1, estimator = "mle"),
    "should be one of"
  )
})

test_that("Inf df is permitted (residual-only limit)", {
  comp <- re_components(
    components = c(Residual = 2),
    dfs        = c(Residual = Inf, ci = Inf, pi = Inf),
    mean = 0, var_mean = 1
  )
  expect_true(is.infinite(comp$dfs[["pi"]]))
})

test_that("ems defaults to NULL and target defaults to observable", {
  comp <- re_components(
    components = c(a = 1, Residual = 2),
    dfs        = c(ci = 3, pi = 5),
    mean = 0, var_mean = 1
  )
  expect_null(comp$ems)
  expect_identical(comp$target, "observable")
})

test_that("ems is accepted, validated, and aligned to ems$ms order", {
  comp <- re_components(
    components = c(a = 1, Residual = 2),
    dfs        = c(ci = 3, pi = 5),
    mean = 0, var_mean = 1,
    ems = list(
      ms = c(A = 4, Residual = 2),
      k  = c(Residual = 0.5, A = 0.25),   # supplied out of order
      df = c(A = 5, Residual = 12)
    )
  )
  expect_named(comp$ems, c("ms", "k", "df"))
  expect_named(comp$ems$k, c("A", "Residual"))   # realigned to ms order
  expect_identical(unname(comp$ems$k[["A"]]), 0.25)
})

test_that("re_components rejects malformed ems and bad target", {
  base <- list(components = c(a = 1, Residual = 2),
               dfs = c(ci = 3, pi = 5), mean = 0, var_mean = 1)
  call_ems <- function(ems) {
    do.call(re_components, c(base, list(ems = ems)))
  }
  expect_error(call_ems(list(ms = c(A = 1))), "\"ms\", \"k\", \"df\"")
  expect_error(
    call_ems(list(ms = c(A = 1), k = c(A = 1), df = c(B = 5))),
    "same names"
  )
  expect_error(
    call_ems(list(ms = c(A = 1), k = c(A = 1), df = c(A = -2))),
    "must be positive"
  )
  expect_error(
    call_ems(list(ms = c(A = -1), k = c(A = 1), df = c(A = 5))),
    "non-negative"
  )
  expect_error(
    do.call(re_components, c(base, list(target = "future"))),
    "should be one of"
  )
})
