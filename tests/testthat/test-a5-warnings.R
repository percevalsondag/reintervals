# A5 leftovers: the extrapolation warning, the default-df fallback warning and
# the slope-model newdata check. None of them changes a number: the warnings are
# raised alongside the unchanged interval, and the newdata check only rejects
# input that used to fail obscurely or be silently coerced.

lmer_a5 <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))

## every warning signalled by `expr`, by class (muffled); the value is attached
catch_warnings <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}
n_class <- function(res, class) sum(vapply(res$warnings, inherits, logical(1), class))
msg_class <- function(res, class) {
  hit <- Filter(function(w) inherits(w, class), res$warnings)
  if (length(hit)) conditionMessage(hit[[1]]) else NA_character_
}

ctl_ti <- lmm_interval_control(seed = 1, M = 2000)
ss <- lme4::sleepstudy
fit_rs <- function() lmer_a5(Reaction ~ Days + (1 + Days | Subject), ss)
fit_fs <- function() lmer_a5(Reaction ~ Days + (1 | Subject), ss)

## the five verbs, each called once
verbs <- list(
  ci  = function(m, nd) ci_lmm(m, newdata = nd),
  pi  = function(m, nd) pi_lmm(m, newdata = nd),
  ngm = function(m, nd) new_group_mean_lmm(m, newdata = nd),
  ti  = function(m, nd) ti_lmm(m, newdata = nd, control = ctl_ti),
  predict = function(m, nd) lmm_predict(m, "confidence", newdata = nd)
)

## ---- 1. extrapolation ----------------------------------------------------------
test_that("extrapolation warns at t0 = 48 on random- and fixed-slope fits, every verb", {
  skip_if_not_installed("lme4")
  for (nm in c("rs", "fs")) {
    m <- if (nm == "rs") fit_rs() else fit_fs()
    for (v in names(verbs)) {
      info <- paste(nm, v)
      out <- catch_warnings(verbs[[v]](m, data.frame(Days = 48)))
      expect_identical(n_class(out, "reintervals_warning_extrapolation"), 1L, info = info)
      msg <- msg_class(out, "reintervals_warning_extrapolation")
      expect_match(msg, "Days = 48", fixed = TRUE, info = info)
      expect_match(msg, "[0, 9]", fixed = TRUE, info = info)
      expect_match(msg, if (nm == "rs") "quadratic in t0" else "linear (fixed-slope)",
                   fixed = TRUE, info = info)
      ## inside the range: no extrapolation warning
      ins <- catch_warnings(verbs[[v]](m, data.frame(Days = 5)))
      expect_identical(n_class(ins, "reintervals_warning_extrapolation"), 0L, info = info)
    }
  }
})

test_that("extrapolation warns ONCE per call with several out-of-range t0, naming them all", {
  skip_if_not_installed("lme4")
  for (m in list(fit_rs(), fit_fs())) {
    for (v in names(verbs)) {
      out <- catch_warnings(verbs[[v]](m, data.frame(Days = c(-2, 5, 48, 60))))
      expect_identical(n_class(out, "reintervals_warning_extrapolation"), 1L, info = v)
      w <- Filter(function(w) inherits(w, "reintervals_warning_extrapolation"), out$warnings)[[1]]
      expect_identical(w$t0, c(-2, 48, 60), info = v)
      expect_identical(w$range, c(0, 9), info = v)
      expect_s3_class(w, "reintervals_warning")
      expect_match(conditionMessage(w), "Days = -2, 48, 60 are outside", fixed = TRUE)
    }
  }
})

test_that("the warning does not change the interval", {
  skip_if_not_installed("lme4")
  m <- fit_rs()
  a <- quiet_extrapolation(ci_lmm(m, newdata = data.frame(Days = c(5, 48))))
  b <- ci_lmm(m, newdata = data.frame(Days = 5))
  expect_identical(a$lower[1], b$lower); expect_identical(a$upper[1], b$upper)
})

test_that("ti_gpq_raw() and ti_anova_raw() warn on the range of their time argument", {
  out <- catch_warnings(ti_gpq_raw(ss$Reaction, ss$Days, ss$Subject, t0 = c(20, 30),
                                   control = ctl_ti))
  expect_identical(n_class(out, "reintervals_warning_extrapolation"), 1L)
  expect_match(msg_class(out, "reintervals_warning_extrapolation"), "t0 = 20, 30", fixed = TRUE)
  expect_identical(n_class(catch_warnings(ti_gpq_raw(ss$Reaction, ss$Days, ss$Subject,
                   t0 = 5, control = ctl_ti)), "reintervals_warning_extrapolation"), 0L)
  out <- catch_warnings(ti_anova_raw(ss$Reaction, ss$Days, ss$Subject, t0 = c(-1, 12)))
  expect_identical(n_class(out, "reintervals_warning_extrapolation"), 1L)
  expect_match(msg_class(out, "reintervals_warning_extrapolation"), "linear (fixed-slope)",
               fixed = TRUE)
  expect_identical(n_class(catch_warnings(ti_anova_raw(ss$Reaction, ss$Days, ss$Subject,
                   t0 = 9)), "reintervals_warning_extrapolation"), 0L)
  ## no time covariate: not applicable
  expect_identical(n_class(catch_warnings(ti_anova_raw(hk_assay$conc, time = NULL,
                   lot = hk_assay$run)), "reintervals_warning_extrapolation"), 0L)
})

test_that("random-intercept models without a time covariate never warn on extrapolation", {
  skip_if_not_installed("lme4")
  fh <- lmer_a5(conc ~ 1 + (1 | run), hk_assay)
  for (v in c("ci", "ngm", "ti")) {
    out <- catch_warnings(verbs[[v]](fh, NULL))
    expect_identical(n_class(out, "reintervals_warning_extrapolation"), 0L, info = v)
  }
})

## ---- 2. default-df fallback ----------------------------------------------------
## a singular random-intercept fit: identical run means -> zero between-run
## variance (test-singular-guard.R); 6 runs, so the few-groups PI limitation
## warning (<= 5 levels) does not fire
singular_ri <- function() {
  d <- sqrt(0.001253)
  y <- as.vector(sapply(rep(0.981, 6), function(m) m + c(-d, 0, d)))
  lmer_a5(y ~ 1 + (1 | run), data.frame(y = y, run = factor(rep(1:6, each = 3))))
}

test_that("df_fallback fires once on a singular random-intercept CI and PI", {
  skip_if_not_installed("lme4")
  m <- singular_ri()
  expect_true(lme4::isSingular(m))
  ci <- catch_warnings(ci_lmm(m))
  expect_identical(attr(ci$value, "diagnostics")$df_method, "containment")
  expect_identical(n_class(ci, "reintervals_warning_df_fallback"), 1L)
  w <- Filter(function(w) inherits(w, "reintervals_warning_df_fallback"), ci$warnings)[[1]]
  expect_identical(w$verb, "ci_lmm()"); expect_identical(w$rule, "containment")
  expect_identical(w$reason, "singular fit"); expect_equal(w$df, 5)
  expect_match(conditionMessage(w), "ci_lmm(): the default Satterthwaite df", fixed = TRUE)
  expect_match(conditionMessage(w), "df = 5", fixed = TRUE)
  ## the existing singular-fit warning is a TI-only warning: not on the CI / PI
  expect_identical(n_class(ci, "reintervals_warning_singular"), 0L)

  pi <- catch_warnings(pi_lmm(m))
  expect_true(attr(pi$value, "diagnostics")$pi_df_fallback)
  expect_identical(n_class(pi, "reintervals_warning_df_fallback"), 1L)
  expect_match(msg_class(pi, "reintervals_warning_df_fallback"),
               "pi_lmm(): the default Francq generalized-Satterthwaite df", fixed = TRUE)
})

test_that("df_fallback does not fire on a healthy fit or under a containment arm", {
  skip_if_not_installed("lme4")
  for (v in c("ci", "pi", "ngm")) {
    out <- catch_warnings(verbs[[v]](fit_fs(), data.frame(Days = 5)))
    expect_identical(n_class(out, "reintervals_warning_df_fallback"), 0L, info = v)
    out <- catch_warnings(verbs[[v]](fit_rs(), data.frame(Days = 5)))
    expect_identical(n_class(out, "reintervals_warning_df_fallback"), 0L, info = v)
  }
  m <- singular_ri()
  ci <- catch_warnings(ci_lmm(m, control = lmm_interval_control(ci_df = "containment")))
  expect_identical(n_class(ci, "reintervals_warning_df_fallback"), 0L)
  expect_identical(n_class(ci, "reintervals_warning_unvalidated"), 1L)
  pi <- catch_warnings(pi_lmm(m, control = lmm_interval_control(pi_df = "containment")))
  expect_identical(n_class(pi, "reintervals_warning_df_fallback"), 0L)
})

test_that("df_fallback fires on a random-slope fit whose Hessian-Satterthwaite df falls back", {
  skip_if_not_installed("lme4")
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  dn <- vfun_sim(B = 10, times = c(0, 3, 6, 9, 12, 18, 24), s2_0 = 2.0,
                 s2_1 = 1e-6, s01 = 0, s2e = 0.7, seed = 1)
  fn <- lmer_a5(y ~ time + (1 | batch) + (0 + time | batch), dn, control = lc)
  skip_if_not(is.null(vfun_extract(fn, t0 = 12)$cov), "the covariance did not fail here")
  nd <- data.frame(time = c(6, 12))
  for (v in c("ci", "pi", "ngm")) {
    out <- catch_warnings(verbs[[v]](fn, nd))
    expect_identical(attr(out$value, "diagnostics")$df_rule, rep("containment", 2), info = v)
    expect_identical(n_class(out, "reintervals_warning_df_fallback"), 1L, info = v)
    expect_match(msg_class(out, "reintervals_warning_df_fallback"),
                 "default Hessian-Satterthwaite df", fixed = TRUE, info = v)
    expect_match(msg_class(out, "reintervals_warning_df_fallback"),
                 "covariance of the variance parameters", fixed = TRUE, info = v)
  }
  ## a validation arm never falls back and does not warn
  k1 <- catch_warnings(ci_lmm(fn, nd, control = lmm_interval_control(slope_ci_df = "k-minus-1")))
  expect_identical(n_class(k1, "reintervals_warning_df_fallback"), 0L)
  ## mocked: the Satterthwaite df alone fails (the covariance is fine)
  fm <- fit_rs()
  local_mocked_bindings(.vfun_df_satterthwaite = function(value, grad_vec, vc_cov) NA_real_)
  out <- catch_warnings(pi_lmm(fm, data.frame(Days = 5)))
  expect_identical(n_class(out, "reintervals_warning_df_fallback"), 1L)
  expect_match(msg_class(out, "reintervals_warning_df_fallback"), "(non-finite df)", fixed = TRUE)
})

test_that("df_fallback fires on the crossed new-group-mean containment fallback", {
  skip_if_not_installed("lme4")
  set.seed(1); A <- 3; B <- 4; n <- 2
  g <- expand.grid(r = 1:n, a = factor(1:A), b = factor(1:B)); g$ab <- interaction(g$a, g$b)
  g$y <- 10 + rnorm(A, 0, 2)[g$a] + rnorm(B, 0, 2)[g$b] + rnorm(A * B, 0, 1)[g$ab] + rnorm(nrow(g))
  fm <- lmer_a5(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), g)
  expect_identical(n_class(catch_warnings(new_group_mean_lmm(fm)),
                           "reintervals_warning_df_fallback"), 0L)
  real <- .design_components_fit
  local_mocked_bindings(.design_components_fit = function(...) {
    dc <- real(...); dc$ems$k[1:2] <- 0; dc
  })
  out <- catch_warnings(new_group_mean_lmm(fm))
  expect_identical(attr(out$value, "diagnostics")$df_rule, "containment")
  expect_identical(n_class(out, "reintervals_warning_df_fallback"), 1L)
  expect_match(msg_class(out, "reintervals_warning_df_fallback"),
               "new_group_mean_lmm(): the default Satterthwaite (crossed Gv) df", fixed = TRUE)
  expect_match(msg_class(out, "reintervals_warning_df_fallback"),
               "k_A EMS_A + k_B EMS_B is not positive", fixed = TRUE)
  ## a non-default ngm_method arm: no fallback warning
  arm <- catch_warnings(new_group_mean_lmm(fm, control =
                          lmm_interval_control(ngm_method = "broad-inference")))
  expect_identical(n_class(arm, "reintervals_warning_df_fallback"), 0L)
})

## ---- 3. slope-model newdata check ---------------------------------------------
test_that("slope-model newdata: missing, character and NA time column error, classed, naming Days", {
  skip_if_not_installed("lme4")
  bad <- list(missing = data.frame(x = 5),
              character = data.frame(Days = "5", stringsAsFactors = FALSE),
              na = data.frame(Days = c(5, NA)))
  for (m in list(fit_rs(), fit_fs())) {
    for (cs in names(bad)) for (v in names(verbs)) {
      expect_error(verbs[[v]](m, bad[[cs]]), "Days", class = "reintervals_error_input",
                   info = paste(cs, v))
    }
  }
  expect_error(ci_lmm(fit_rs(), data.frame(Days = "5")), "is of class 'character', not numeric",
               fixed = TRUE)
  expect_error(ci_lmm(fit_fs(), data.frame(Days = c(5, NA))), "contains missing values (NA)",
               fixed = TRUE)
  expect_error(ci_lmm(fit_fs(), data.frame(Days = factor(5))), "is of class 'factor'",
               fixed = TRUE)
  expect_error(ci_lmm(fit_rs(), data.frame(Days = Inf)), "non-finite", fixed = TRUE)
})

test_that("slope-model newdata: extra columns are ignored", {
  skip_if_not_installed("lme4")
  for (m in list(fit_rs(), fit_fs())) {
    a <- ci_lmm(m, newdata = data.frame(Days = 5, other = "x"))
    b <- ci_lmm(m, newdata = data.frame(Days = 5))
    expect_identical(a$lower, b$lower); expect_identical(a$upper, b$upper)
  }
})

## ---- class tree -------------------------------------------------------------
test_that("df_fallback is a reintervals_warning_fallback; extrapolation is not", {
  skip_if_not_installed("lme4")
  expect_identical(.condition_chain("reintervals_warning_df_fallback"),
                   c("reintervals_warning_df_fallback", "reintervals_warning_fallback",
                     "reintervals_warning"))
  expect_identical(.condition_chain("reintervals_warning_extrapolation"),
                   c("reintervals_warning_extrapolation", "reintervals_warning"))
  ## a handler on the parent class catches the df-fallback warning
  caught <- NULL
  withCallingHandlers(ci_lmm(singular_ri()), reintervals_warning_fallback = function(w) {
    caught <<- class(w); invokeRestart("muffleWarning")
  })
  expect_true("reintervals_warning_df_fallback" %in% caught)
  expect_warning(ci_lmm(singular_ri()), class = "reintervals_warning_fallback")
})

## ---- #13: a failed lme4::isSingular() is not "non-singular" ---------------------
test_that("a failing isSingular() gives singular = NA, singular_check, and the singular branch", {
  skip_if_not_installed("lme4")
  fh <- lmer_a5(conc ~ 1 + (1 | run), hk_assay)   # healthy: Francq df is Satterthwaite
  fr <- fit_rs()
  expect_false(attr(pi_lmm(fh), "diagnostics")$pi_df_fallback)
  local_mocked_bindings(isSingular = function(x, ...) stop("boom"), .package = "lme4")
  expect_identical(.singular_status(fh), list(singular = NA, check = "failed: boom"))
  ## PI: the conservative (singular) branch -> containment df, and it says why
  pi <- catch_warnings(pi_lmm(fh))
  d <- attr(pi$value, "diagnostics")
  expect_true(is.na(pi$value$singular))
  expect_identical(d$singular_check, "failed: boom")
  expect_true(d$pi_df_fallback)
  expect_equal(d$df, .ddf_containment(fh))
  expect_match(msg_class(pi, "reintervals_warning_df_fallback"),
               "singularity check failed; treated as singular", fixed = TRUE)
  ## CI: Satterthwaite needs the Hessian, which the singular branch skips
  ci <- catch_warnings(ci_lmm(fh))
  expect_identical(attr(ci$value, "diagnostics")$df_method, "containment")
  expect_true(is.na(ci$value$singular))
  ## EMS TI: the singular warning fires (NA treated as singular)
  ti <- catch_warnings(ti_lmm(fh, method = "ems-mls"))
  expect_identical(n_class(ti, "reintervals_warning_singular"), 1L)
  expect_true(is.na(ti$value$singular))
  ## random slope: CI / NGM / GPQ TI report NA too
  for (r in list(ci_lmm(fr, data.frame(Days = 5)),
                 new_group_mean_lmm(fr, data.frame(Days = 5)),
                 ti_lmm(fr, data.frame(Days = 5), control = ctl_ti))) {
    expect_true(is.na(r$singular))
    expect_identical(attr(r, "diagnostics")$singular_check, "failed: boom")
  }
})

test_that("a working isSingular() adds no singular_check field", {
  skip_if_not_installed("lme4")
  expect_null(attr(ci_lmm(fit_fs()), "diagnostics")$singular_check)
  expect_identical(.singular_status(fit_fs()), list(singular = FALSE, check = NA_character_))
})

## ---- NaN guard in .unify_intervals -------------------------------------------------
test_that("a NaN estimate or limit aborts with reintervals_error_numerical, naming verb and row", {
  native <- data.frame(type = "CI", estimate = c(1, 2), lower = c(0, NaN),
                       upper = c(2, 3), level = 0.95)
  expect_error(.unify_intervals(native, NULL, "CI"), "ci_lmm(): the engine returned NaN for `lower` in row(s) 2",
               fixed = TRUE, class = "reintervals_error_numerical")
  native$lower[2] <- 1; native$estimate[1] <- NaN
  e <- tryCatch(.unify_intervals(native, NULL, "ti", verb = "ti_gpq_raw()"), error = identity)
  expect_s3_class(e, "reintervals_error_numerical")
  expect_identical(e$verb, "ti_gpq_raw()"); expect_identical(e$rows, 1L)
  expect_identical(e$column, "estimate")
  ## NA (not NaN) and +/-Inf pass
  native$estimate[1] <- 1; native$lower <- c(NA, -Inf)
  expect_s3_class(.unify_intervals(native, NULL, "CI", method_tag = "x"), "lmm_interval")
  ## end to end: a mocked engine row with a NaN limit
  skip_if_not_installed("lme4")
  local_mocked_bindings(.one_sided_t = function(out, level, sides) { out$upper[1] <- NaN; out })
  expect_error(pi_lmm(fit_fs(), data.frame(Days = c(1, 2))), "pi_lmm().*row\\(s\\) 1",
               class = "reintervals_error_numerical")
})

## ---- ems_synthesis input check ----------------------------------------------------
test_that("ems_synthesis() rejects a degenerate crossed count matrix before model.matrix()", {
  expect_error(ems_synthesis(matrix(1, 1, 1), "crossed"), "at least 2 levels of each factor",
               class = "reintervals_error_design")
  expect_error(ems_synthesis(matrix(2, 1, 3), "crossed"), class = "reintervals_error_design")
  expect_error(ems_synthesis(matrix(1, 2, 2), "crossed"), "replicated cell",
               class = "reintervals_error_design")
  expect_type(ems_synthesis(matrix(c(3, 4, 1, 1), 2, 2), "crossed")$k, "double")
})

## ---- #14: grouping_reason -------------------------------------------------------
test_that("a failed design classification records grouping_reason next to grouping = NA", {
  skip_if_not_installed("lme4")
  fh <- lmer_a5(conc ~ 1 + (1 | run), hk_assay)
  expect_null(attr(ci_lmm(fh), "diagnostics")$grouping_reason)
  local_mocked_bindings(.design_of = function(model) stop("cannot classify"))
  r <- ci_lmm(fh)
  expect_true(is.na(attr(r, "grouping")))
  expect_null(attributes(attr(r, "grouping")))
  expect_identical(attr(r, "diagnostics")$grouping_reason,
                   "classification failed: cannot classify")
})
