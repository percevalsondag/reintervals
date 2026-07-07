# REML-slope interval arithmetic: regression check against committed reference values.
#
# Read this before changing the tolerance:
#   * This checks the INTERVAL ARITHMETIC (the method-pure constructors) against
#     committed fixed reference values in slope_oracle.rds, at 1e-9. No model is
#     fitted at test time, so the check is identical on every platform.
#   * The fit -> component EXTRACTION (vfun_extract) is a separate concern: any
#     lme4 wrapper inherits cross-platform REML/Hessian jitter (~1e-6, mostly in
#     the Hessian-derived df), so extraction is validated separately against the
#     Montes oracle (n_E etc.) in test-vfun-extract.R, not here.
#
# Why fixed components rather than a re-fit-and-compare: re-fitting on each
# platform lands on slightly different variance components, so end-to-end interval
# values drift ~1e-6 across platforms. That drift is lme4 fit nondeterminism, not
# a difference in this package's arithmetic. Running the constructors on committed
# extracted components isolates the arithmetic and keeps the tight 1e-9 gate.

test_that("A2 (CI gate): slope constructors reproduce oracle to 1e-9 (committed components, no re-fit)", {
  fx     <- readRDS(test_path("fixtures", "slope_oracle.rds"))
  comps  <- fx$comps
  oracle <- fx$oracle

  # constructor per interval kind (method-pure; consumes the committed component)
  build <- list(
    CI    = function(comp) .vfun_interval(comp, "CI"),
    PI    = function(comp) .vfun_interval(comp, "PI"),
    CInew = function(comp) .vfun_interval(comp, "new_group_mean"),
    TI    = function(comp) ti_vfun(comp, content = 0.99, conf = 0.95)
  )

  for (case in names(comps)) {
    comp <- comps[[case]]
    for (ty in names(build)) {
      ref <- oracle[oracle$case == case & oracle$type == ty, ]
      got <- build[[ty]](comp)
      expect_equal(got$lower, ref$lower, tolerance = 1e-9, info = paste(case, ty, "lower"))
      expect_equal(got$upper, ref$upper, tolerance = 1e-9, info = paste(case, ty, "upper"))
      expect_equal(got$df,    ref$df,    tolerance = 1e-9, info = paste(case, ty, "df"))
    }
  }
})
