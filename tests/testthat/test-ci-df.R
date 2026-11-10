# CI denominator df on the random-intercept (EMS) path: Satterthwaite with theta
# computed, SAS-style containment fallback when it is unavailable, Kenward-Roger
# only on request -- and never dependent on which packages are installed.

ci_nested <- function(seed, sa, sb) {           # A = 5, B = 3 within A, n = 2
  set.seed(seed); A <- 5; B <- 3; n <- 2
  a <- factor(rep(1:A, each = B * n)); b <- factor(rep(1:(A * B), each = n))
  data.frame(a = a, b = b, y = 10 + rnorm(A, 0, sa)[a] + rnorm(A * B, 0, sb)[b] + rnorm(A * B * n))
}
ci_crossed <- function(seed, sa, sb, sab) {     # A = 3 x B = 4, n = 2
  set.seed(seed); A <- 3; B <- 4; n <- 2
  g <- expand.grid(r = 1:n, a = factor(1:A), b = factor(1:B)); g$ab <- interaction(g$a, g$b)
  g$y <- 10 + rnorm(A, 0, sa)[g$a] + rnorm(B, 0, sb)[g$b] + rnorm(A * B, 0, sab)[g$ab] + rnorm(nrow(g))
  g
}
ci_oneway_singular <- function() {
  set.seed(4); g <- factor(rep(1:5, each = 3))
  data.frame(g = g, y = 10 + rnorm(15))
}
fit <- function(f, d) suppressMessages(lme4::lmer(f, data = d))
diag_of <- function(x) attr(x, "diagnostics")

## the designs: (formula, data, expected singular, expected containment df)
ci_cases <- function() list(
  oneway_ok   = list(conc ~ 1 + (1 | run), hk_assay, FALSE, 5),
  oneway_sing = list(y ~ 1 + (1 | g), ci_oneway_singular(), TRUE, 4),
  nested_ok   = list(y ~ 1 + (1 | a/b), ci_nested(1, 2, 1), FALSE, 4),
  nested_sing = list(y ~ 1 + (1 | a/b), ci_nested(1, 0, 0), TRUE, 4),
  crossed_ok  = list(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), ci_crossed(1, 2, 2, 1), FALSE, 2),
  crossed_sing = list(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), ci_crossed(1, 0, 0, 0), TRUE, 2)
)

test_that("containment df is the SAS rank rule: A - 1 (one-way, nested), min(A, B) - 1 (crossed)", {
  skip_if_not_installed("lme4")
  for (nm in names(ci_cases())) {
    cs <- ci_cases()[[nm]]
    expect_equal(.ddf_containment(fit(cs[[1]], cs[[2]])), cs[[4]], info = nm)
  }
  # a between-group fixed column costs one df: one-way with a run-level covariate
  d <- hk_assay; d$z <- as.numeric(d$run) %% 2
  expect_equal(.ddf_containment(fit(conc ~ z + (1 | run), d)), 4)
})

test_that("non-singular fits use Satterthwaite; singular fits fall back to containment", {
  skip_if_not_installed("lme4")
  for (nm in names(ci_cases())) {
    cs <- ci_cases()[[nm]]
    fm <- fit(cs[[1]], cs[[2]])
    expect_identical(lme4::isSingular(fm), cs[[3]], info = nm)
    x <- quiet_df_fallback(ci_lmm(fm))
    if (cs[[3]]) {
      expect_identical(diag_of(x)$df_method, "containment", info = nm)
      expect_equal(diag_of(x)$df, cs[[4]], info = nm)
    } else {
      expect_identical(diag_of(x)$df_method, "satterthwaite", info = nm)
      expect_equal(diag_of(x)$df,
                   .ci_ddf_satterthwaite(fm, 1, .theta_cov(fm)), info = nm)
    }
  }
})

# oracle: pbkrtest-kr-ci-cases (fixtures/oracles.csv)
test_that("Kenward-Roger is used only on request and reported as such", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  ## df AND SE: the KR df must be paired with the KR-adjusted covariance
  ## (vcovAdj), not the unadjusted vcov(). The SAS-fixture designs are covered in
  ## test-oracle-sas-ci.R; this adds the crossed design and the hk_assay one-way.
  for (nm in c("oneway_ok", "nested_ok", "crossed_ok")) {
    cs <- ci_cases()[[nm]]
    fm <- fit(cs[[1]], cs[[2]])
    k <- ci_lmm(fm, control = lmm_interval_control(ci_df = "kenward-roger"))
    expect_identical(diag_of(k)$df_method, "kenward-roger", info = nm)
    expect_equal(diag_of(k)$df, pbkrtest::get_Lb_ddf(fm, 1), info = nm)
    se_k <- (k$upper - k$lower) / 2 / stats::qt(0.975, diag_of(k)$df)
    expect_equal(se_k, sqrt(as.matrix(pbkrtest::vcovAdj(fm))[1, 1]),
                 tolerance = 1e-8, info = nm)
    expect_identical(diag_of(ci_lmm(fm))$df_method, "satterthwaite", info = nm)  # default
  }
})

test_that("results are identical with pbkrtest masked; KR then errors instead of falling back", {
  skip_if_not_installed("lme4")
  fits <- lapply(ci_cases(), function(cs) fit(cs[[1]], cs[[2]]))
  before <- lapply(fits, function(f) quiet_df_fallback(ci_lmm(f)))
  local_mocked_bindings(requireNamespace = function(package, ...) {
    if (identical(package, "pbkrtest")) FALSE else base::requireNamespace(package, ...)
  })
  after <- lapply(fits, function(f) quiet_df_fallback(ci_lmm(f)))
  for (nm in names(fits)) {
    expect_identical(after[[nm]], before[[nm]], info = nm)
  }
  expect_error(ci_lmm(fits$oneway_ok, control = lmm_interval_control(ci_df = "kenward-roger")),
               "requires the 'pbkrtest' package")
})

## B1 (Percy, 2026-11-08): an option on a model it does not apply to is a
## no-op, never an error (one control object across a simulation grid). Until
## B1, ci_df = "kenward-roger" on a random-slope model was an error; the
## random-slope KR arm is slope_ci_df = "kenward-roger" (test-df-arms.R).
test_that("ci_df = 'kenward-roger' is a no-op on random-slope models", {
  skip_if_not_installed("lme4")
  fs <- fit(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = 5)
  expect_no_warning(a <- ci_lmm(fs, newdata = nd,
                                control = lmm_interval_control(ci_df = "kenward-roger")))
  b <- ci_lmm(fs, newdata = nd)
  expect_identical(a, b)
  expect_identical(attr(a, "diagnostics")$ci_df_rule, "not_applicable")
})
