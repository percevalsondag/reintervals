# B1: the step-13 degrees-of-freedom arms as lmm_interval_control() options
# (H1 ci_df, H2 slope_ci_df, H6 pi_df). Not oracles, except where a published
# identity is named (the balanced random-coefficient K - 1 df). The default
# pins below were computed on the clean A2 commit 716e880 (CI / PI code
# identical to 2e09497) BEFORE the B1 edits, so they prove that no default
# output changed.

df_arm_fits <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  fit <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, control = lc)))
  set.seed(20261108)
  ow <- data.frame(g = factor(rep(1:6, each = 4)))
  ow$y <- 10 + rnorm(6, sd = 2)[ow$g] + rnorm(24)
  ne <- expand.grid(r = 1:2, b = 1:3, a = factor(1:5))
  ne$b <- factor(paste(ne$a, ne$b))
  ne$y <- 5 + rnorm(5, sd = 1.5)[ne$a] + rnorm(15)[ne$b] + rnorm(nrow(ne))
  cr <- expand.grid(r = 1:2, a = factor(1:4), b = factor(1:6))
  cr$y <- 3 + rnorm(4)[cr$a] + rnorm(6, sd = 1.2)[cr$b] + rnorm(24, sd = 0.5)[interaction(cr$a, cr$b)] + rnorm(nrow(cr))
  ss <- lme4::sleepstudy
  list(
    oneway  = list(fit = fit(y ~ 1 + (1 | g), ow), nd = data.frame(g = "1")),
    nested  = list(fit = fit(y ~ 1 + (1 | a/b), ne), nd = NULL),
    crossed = list(fit = fit(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), cr), nd = NULL),
    fixed_slope = list(fit = fit(Reaction ~ Days + (1 | Subject), ss), nd = data.frame(Days = c(0, 5, 9))),
    rs_un  = list(fit = fit(Reaction ~ Days + (1 + Days | Subject), ss), nd = data.frame(Days = c(0, 5, 9))),
    rs_ind = list(fit = fit(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), ss), nd = data.frame(Days = c(0, 5, 9))))
}

## muffle only the documented limitation warnings (few-groups PI, once per session)
quiet <- function(expr) withCallingHandlers(expr,
  reintervals_warning_limitation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))
## muffle the B1 unvalidated-arm warning (tested separately in (d))
quiet_arm <- function(expr) withCallingHandlers(expr,
  reintervals_warning_unvalidated = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_limitation  = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))

default_pins <- list(
  oneway_ci = list(lower = c(7.1170263152),
    upper = c(11.3684837973),
    df = c(4.9999993047)),
  oneway_pi = list(lower = c(3.6690447465),
    upper = c(14.8164653660),
    df = c(6.8511276696)),
  nested_ci = list(lower = c(1.7760443339),
    upper = c(7.1663806779),
    df = c(4.0000201885)),
  nested_pi = list(lower = c(-1.8939677150),
    upper = c(10.8363927268),
    df = c(5.6250680888)),
  crossed_ci = list(lower = c(1.2580973101),
    upper = c(3.5303280233),
    df = c(5.6170153722)),
  crossed_pi = list(lower = c(-0.8300071332),
    upper = c(5.6184324667),
    df = c(19.4080533733)),
  fixed_slope_ci = list(lower = c(231.2331986014, 284.6347295415, 325.4387722381),
    upper = c(271.5770110956, 322.8483397515, 365.7825847316),
    df = c(22.8101988892, 17.0671832940, 22.8101988967)),
  fixed_slope_pi = list(lower = c(151.7980294093, 204.3995316231, 246.0036030456),
    upper = c(351.0121802877, 403.0835376699, 445.2177539241),
    df = c(41.2819537637, 41.2819537637, 41.2819537637)),
  rs_un_ci = list(lower = c(237.0064386041, 283.5267835284, 314.7465565099),
    upper = c(265.8037710929, 323.9562857645, 376.4748004598),
    df = c(16.9995997496, 16.9997480698, 16.9998621234)),
  rs_un_pi = list(lower = c(178.0120114946, 205.2296389976, 206.1728748171),
    upper = c(324.7981982023, 402.2534302953, 485.0484821526),
    df = c(37.6672179695, 30.7765516191, 20.9635367973)),
  rs_ind_ci = list(lower = c(236.9483715713, 284.3790613100, 315.4443529741),
    upper = c(265.8618381257, 323.1040079829, 375.7770039956),
    df = c(18.1561946500, 25.0661100348, 19.7304903083)),
  rs_ind_pi = list(lower = c(177.6968544604, 207.9547161822, 208.9390259510),
    upper = c(325.1133552366, 399.5283531108, 482.2823310187),
    df = c(39.2948395272, 46.5864140224, 24.4929085003))
)

F_b1 <- df_arm_fits()
ri_designs <- c("oneway", "nested", "crossed", "fixed_slope")
rs_designs <- c("rs_un", "rs_ind")

test_that("B1 (a): defaults unchanged -- default == explicit default == pre-B1 pin", {
  explicit <- lmm_interval_control(ci_df = "satterthwaite",
                                   slope_ci_df = "hessian-satterthwaite",
                                   pi_df = "francq-satterthwaite")
  for (nm in names(F_b1)) for (vb in c("ci", "pi")) {
    f <- if (vb == "ci") ci_lmm else pi_lmm
    fx <- F_b1[[nm]]; info <- paste(nm, vb)
    expect_no_warning(a <- quiet(f(fx$fit, newdata = fx$nd)))      # no unvalidated warning
    b <- quiet(f(fx$fit, newdata = fx$nd, control = explicit))
    expect_identical(a, b, info = info)
    pin <- default_pins[[paste(nm, vb, sep = "_")]]
    expect_equal(a$lower, pin$lower, tolerance = 1e-10, info = info)
    expect_equal(a$upper, pin$upper, tolerance = 1e-10, info = info)
    expect_equal(attr(a, "diagnostics")$df, pin$df, tolerance = 1e-10, info = info)
    d <- attr(a, "diagnostics")
    rs <- nm %in% rs_designs
    expect_identical(unique(d$ci_df_rule),
                     if (vb == "ci" && !rs) "satterthwaite" else "not_applicable", info = info)
    expect_identical(unique(d$slope_ci_df_rule),
                     if (vb == "ci" && rs) "hessian-satterthwaite" else "not_applicable", info = info)
    expect_identical(unique(d$pi_df_rule),
                     if (vb == "pi" && !rs) "francq-satterthwaite" else "not_applicable", info = info)
  }
})

test_that("B1 (b): containment arms use the rank-rule df / K - 1 and recompute from qt()", {
  ## the 0(a) rule on these fixtures: A - 1 = 5, A - 1 = 4, min(A, B) - 1 = 3, K - 1 = 17
  expected <- c(oneway = 5, nested = 4, crossed = 3, fixed_slope = 17)
  q <- stats::qt(0.975, 1:200)
  for (nm in ri_designs) {
    fx <- F_b1[[nm]]; fit <- fx$fit
    expect_identical(.ddf_containment(fit), expected[[nm]], info = nm)
    fxd <- .fixed(fit, fx$nd)
    s2T <- sum(as.data.frame(lme4::VarCorr(fit))$vcov)
    ci <- quiet_arm(ci_lmm(fit, newdata = fx$nd,
                           control = lmm_interval_control(ci_df = "containment")))
    pi <- quiet_arm(pi_lmm(fit, newdata = fx$nd,
                           control = lmm_interval_control(pi_df = "containment")))
    dc <- attr(ci, "diagnostics"); dp <- attr(pi, "diagnostics")
    expect_identical(dc$df, rep(expected[[nm]], nrow(ci)), info = nm)
    expect_identical(dp$df, rep(expected[[nm]], nrow(pi)), info = nm)
    expect_identical(unique(dc$ci_df_rule), "containment", info = nm)
    expect_identical(unique(dp$pi_df_rule), "containment", info = nm)
    tq <- stats::qt(0.975, expected[[nm]])
    se_ci <- unname(sqrt(fxd$var_fix))              # unadjusted l C11 l'
    se_pi <- unname(sqrt(fxd$var_fix + s2T))        # same PI variance V_F + sigma2_T
    expect_equal(ci$lower, ci$estimate - tq * se_ci, tolerance = 1e-10, info = nm)
    expect_equal(ci$upper, ci$estimate + tq * se_ci, tolerance = 1e-10, info = nm)
    expect_equal(pi$lower, pi$estimate - tq * se_pi, tolerance = 1e-10, info = nm)
    expect_equal(pi$upper, pi$estimate + tq * se_pi, tolerance = 1e-10, info = nm)
  }
  for (nm in rs_designs) {
    fx <- F_b1[[nm]]; fit <- fx$fit
    K <- nlevels(lme4::getME(fit, "flist")[[1]])
    ci <- quiet_arm(ci_lmm(fit, newdata = fx$nd,
                           control = lmm_interval_control(slope_ci_df = "k-minus-1")))
    d <- attr(ci, "diagnostics")
    expect_identical(d$df, rep(K - 1, nrow(ci)), info = nm)
    expect_identical(d$df_rule, rep("k-minus-1", nrow(ci)), info = nm)
    expect_identical(unique(d$slope_ci_df_rule), "k-minus-1", info = nm)
    L <- cbind(1, fx$nd$Days)
    se <- sqrt(rowSums((L %*% as.matrix(stats::vcov(fit))) * L))   # unadjusted
    expect_equal(d$se, se, tolerance = 1e-10, info = nm)
    tq <- stats::qt(0.975, K - 1)
    expect_equal(ci$lower, ci$estimate - tq * se, tolerance = 1e-10, info = nm)
    expect_equal(ci$upper, ci$estimate + tq * se, tolerance = 1e-10, info = nm)
  }
})

test_that("B1 (c): slope_ci_df = 'kenward-roger' on random-slope fits; phase-0 equalities (UN only)", {
  skip_if_not_installed("pbkrtest")
  t0 <- c(0, 5, 9, 20); nd <- data.frame(Days = t0)
  kr  <- lmm_interval_control(slope_ci_df = "kenward-roger")
  km1 <- lmm_interval_control(slope_ci_df = "k-minus-1")
  ## UN fit (balanced sleepstudy): KR df == K - 1 == 17 at every t0, SE ratio 1
  un <- F_b1$rs_un$fit
  a <- quiet_arm(ci_lmm(un, newdata = nd, control = kr))
  b <- quiet_arm(ci_lmm(un, newdata = nd, control = km1))
  da <- attr(a, "diagnostics"); db <- attr(b, "diagnostics")
  expect_equal(da$df, rep(17, 4), tolerance = 1e-10)
  expect_equal(da$se / db$se, rep(1, 4), tolerance = 1e-10)   # KR-adjusted / unadjusted
  expect_equal(c(a$lower, a$upper), c(b$lower, b$upper), tolerance = 1e-8)
  expect_identical(da$df_rule, rep("kenward-roger", 4))
  ## || fit: KR runs and is finite; regression pin (values reported in B1 step 0)
  ind <- quiet_arm(ci_lmm(F_b1$rs_ind$fit, newdata = nd, control = kr))
  di <- attr(ind, "diagnostics")$df
  expect_true(all(is.finite(di) & di >= 1))
  expect_true(all(abs(di - c(18.19, 24.65, 19.60, 17.04)) < 1e-2),
              info = paste(signif(di, 6), collapse = ", "))
  ## unbalanced subset: drop the last four days of six subjects
  ss <- lme4::sleepstudy
  ss <- ss[!(ss$Subject %in% c("308", "309", "310", "330", "331", "332") & ss$Days > 5), ]
  ub <- suppressWarnings(suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = ss)))
  u <- quiet_arm(ci_lmm(ub, newdata = nd, control = kr))
  du <- attr(u, "diagnostics")$df
  expect_true(all(is.finite(du) & du >= 1))
  expect_true(all(is.finite(c(u$lower, u$upper))))
})

test_that("B1 (c'): a failed Kenward-Roger df on a random slope aborts; it never falls back", {
  ## A validation arm must not fall back to another df: step 13 would credit the
  ## fallback's coverage to Kenward-Roger. Force each failure mode of
  ## pbkrtest::Lb_ddf (error, NA, Inf, non-positive) and require the classed
  ## numerical error naming the arm.
  skip_if_not_installed("pbkrtest")
  kr <- lmm_interval_control(slope_ci_df = "kenward-roger")
  nd <- data.frame(Days = 5)
  for (bad in list(function(...) stop("forced"), function(...) NA_real_,
                   function(...) Inf, function(...) 0, function(...) -3)) {
    local({
      testthat::local_mocked_bindings(Lb_ddf = bad, .package = "pbkrtest")
      expect_error(quiet_arm(ci_lmm(F_b1$rs_ind$fit, newdata = nd, control = kr)),
                   class = "reintervals_error_numerical",
                   regexp = "slope_ci_df = \"kenward-roger\"", fixed = FALSE)
    })
  }
  ## unmocked again: the arm runs (the mock did not leak)
  ok <- quiet_arm(ci_lmm(F_b1$rs_ind$fit, newdata = nd, control = kr))
  expect_true(is.finite(attr(ok, "diagnostics")$df))
})

test_that("B1 (d): the unvalidated warning fires once per call for each applied non-default arm", {
  fx <- F_b1$oneway; fs <- F_b1$rs_un
  warns_once <- function(expr) {
    n <- 0L
    withCallingHandlers(expr,
      reintervals_warning_unvalidated = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") },
      reintervals_warning_limitation  = function(w) invokeRestart("muffleWarning"))
    n
  }
  expect_identical(warns_once(ci_lmm(fx$fit, newdata = fx$nd,
                   control = lmm_interval_control(ci_df = "containment"))), 1L)
  expect_identical(warns_once(pi_lmm(fx$fit, newdata = fx$nd,
                   control = lmm_interval_control(pi_df = "containment"))), 1L)
  expect_identical(warns_once(ci_lmm(fs$fit, newdata = fs$nd,
                   control = lmm_interval_control(slope_ci_df = "k-minus-1"))), 1L)
  expect_warning(ci_lmm(fs$fit, newdata = fs$nd,
                        control = lmm_interval_control(slope_ci_df = "k-minus-1")),
                 class = "reintervals_warning_unvalidated")
  ## every call warns (not once per session)
  expect_identical(warns_once(ci_lmm(fs$fit, newdata = fs$nd,
                   control = lmm_interval_control(slope_ci_df = "k-minus-1"))), 1L)
  ## defaults never warn
  expect_identical(warns_once(ci_lmm(fx$fit, newdata = fx$nd)), 0L)
  expect_identical(warns_once(pi_lmm(fx$fit, newdata = fx$nd)), 0L)
  expect_identical(warns_once(ci_lmm(fs$fit, newdata = fs$nd)), 0L)
  skip_if_not_installed("pbkrtest")
  ## ci_df = "kenward-roger" on a random intercept: a 0.2.0 option, no warning ...
  expect_identical(warns_once(ci_lmm(fx$fit, newdata = fx$nd,
                   control = lmm_interval_control(ci_df = "kenward-roger"))), 0L)
  ## ... slope_ci_df = "kenward-roger" warns
  expect_identical(warns_once(ci_lmm(fs$fit, newdata = fs$nd,
                   control = lmm_interval_control(slope_ci_df = "kenward-roger"))), 1L)
})

test_that("B1 (e): an option where it does not apply is a no-op, recorded as not_applicable", {
  fs <- F_b1$rs_un; fx <- F_b1$oneway
  arms <- lmm_interval_control(ci_df = "containment", slope_ci_df = "k-minus-1",
                               pi_df = "containment")
  chk <- function(a, b, info) {
    expect_identical(c(a$lower, a$upper), c(b$lower, b$upper), info = info)
    expect_identical(attr(a, "diagnostics")$df, attr(b, "diagnostics")$df, info = info)
  }
  ## non-default pi_df on a random-slope PI: same numbers, not_applicable, no warning
  expect_no_warning(a <- pi_lmm(fs$fit, newdata = fs$nd,
                                control = lmm_interval_control(pi_df = "containment")))
  chk(a, pi_lmm(fs$fit, newdata = fs$nd), "rs PI")
  expect_identical(unique(attr(a, "diagnostics")$pi_df_rule), "not_applicable")
  ## slope_ci_df on a random-intercept CI; ci_df on a random-slope CI
  expect_no_warning(a <- quiet(ci_lmm(fx$fit, newdata = fx$nd,
                    control = lmm_interval_control(slope_ci_df = "k-minus-1"))))
  chk(a, quiet(ci_lmm(fx$fit, newdata = fx$nd)), "ri CI")
  expect_identical(attr(a, "diagnostics")$slope_ci_df_rule, "not_applicable")
  expect_no_warning(a <- ci_lmm(fs$fit, newdata = fs$nd,
                                control = lmm_interval_control(ci_df = "containment")))
  chk(a, ci_lmm(fs$fit, newdata = fs$nd), "rs CI")
  ## every arm at once on the new-group mean (no df arm in B1) and on the PI of
  ## a random-slope model
  for (m in list(fx, fs)) {
    expect_no_warning(a <- quiet(new_group_mean_lmm(m$fit, newdata = m$nd, control = arms)))
    chk(a, quiet(new_group_mean_lmm(m$fit, newdata = m$nd)), "ngm")
    d <- attr(a, "diagnostics")
    expect_identical(unlist(lapply(d[c("ci_df_rule", "slope_ci_df_rule", "pi_df_rule")], unique),
                            use.names = FALSE), rep("not_applicable", 3))
  }
})
