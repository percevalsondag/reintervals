# A5: classed conditions (R/conditions.R) and the formerly silent paths.

lmer_q <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))

## ---- the hierarchy itself ---------------------------------------------------
test_that("every class carries its ancestors up to the root", {
  expect_identical(.condition_chain("reintervals_error_groups"),
                   c("reintervals_error_groups", "reintervals_error_data", "reintervals_error"))
  expect_identical(.condition_chain("reintervals_warning_truncation"),
                   c("reintervals_warning_truncation", "reintervals_warning_fallback",
                     "reintervals_warning"))
  expect_error(.condition_chain("reintervals_error_bogus"), class = "reintervals_error_internal")
  e <- tryCatch(.abort("x", class = "reintervals_error_input"), error = identity)
  expect_s3_class(e, c("reintervals_error_input", "reintervals_error", "error", "condition"))
  expect_null(conditionCall(e))                                   # == stop(call. = FALSE)
  w <- tryCatch(.warn("y", class = "reintervals_warning_mc"), warning = identity)
  expect_s3_class(w, c("reintervals_warning_mc", "reintervals_warning", "warning"))
})

test_that("reserved classes (defunct, deprecated, mc) are defined and signalable", {
  expect_error(.abort("gone", class = "reintervals_error_defunct"),
               class = "reintervals_error_defunct")
  expect_warning(.warn("old", class = "reintervals_warning_deprecated"),
                 class = "reintervals_warning_deprecated")
  expect_warning(.warn("noisy", class = "reintervals_warning_mc"),
                 class = "reintervals_warning_mc")
})

## ---- one raising site per error class ----------------------------------------
test_that("reintervals_error_input", {
  expect_error(ci_lmm(1), class = "reintervals_error_input")
  expect_error(lmm_interval_control(M = 0), class = "reintervals_error_input")
  expect_error(ci_lmm(fit_52(), alpha = 0.05, level = 0.95), class = "reintervals_error_input")
})

test_that("reintervals_error_unsupported", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (1 + Days | Subject), lme4::sleepstudy)
  expect_error(vc_extract(fs), class = "reintervals_error_unsupported")
  fi <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  expect_error(ti_lmm(fi, over = "group_mean"), class = "reintervals_error_unsupported")
})

test_that("reintervals_error_design", {
  skip_if_not_installed("lme4")
  set.seed(2); d <- expand.grid(r = 1:2, a = factor(1:3), b = factor(1:4))
  d$y <- rnorm(nrow(d))
  fm <- lmer_q(y ~ 1 + (1 | a) + (1 | b), d)                      # crossed, no interaction
  expect_error(ti_lmm(fm), class = "reintervals_error_design")
  expect_error(new_group_mean_lmm(fm), class = "reintervals_error_design")
})

test_that("reintervals_error_method", {
  expect_error(ti_lmm(fit_52(), method = "gpq"), class = "reintervals_error_method")
  expect_error(ci_lmm(fit_52(), method = "gpq"), class = "reintervals_error_method")
})

test_that("reintervals_error_groups / time_points / residual_df (all reintervals_error_data)", {
  expect_error(ti_anova_raw(c(1, 2, 3), time = NULL, lot = c(1, 1, 1)),
               class = "reintervals_error_groups")
  y <- c(1, 2, 3, 4, 5, 6, 7, 8); b <- rep(1:2, each = 4)
  expect_error(ti_gpq_raw(y, time = rep(c(0, 0, 1, 2), 2) * c(rep(0, 4), rep(1, 4)),
                          batch = b, t0 = 1), class = "reintervals_error_time_points")
  expect_error(ti_gpq_raw(y[c(1, 2, 5, 6)], time = c(0, 1, 0, 1), batch = c(1, 1, 2, 2),
                          t0 = 1), class = "reintervals_error_residual_df")
  expect_error(ti_anova_raw(c(1, 2, 3), time = NULL, lot = c(1, 1, 1)),
               class = "reintervals_error_data")
})

test_that("reintervals_error_numerical", {
  comp <- re_components(components = c(g = 0, Residual = 0), dfs = c(pi = Inf),
                        mean = 0, var_mean = 0,
                        ems = list(ms = c(g = 0, Residual = 0), k = c(g = 0.5, Residual = 0.5),
                                   df = c(g = 3, Residual = 8)))
  expect_error(ti_francq(comp), class = "reintervals_error_numerical")
})

test_that("reintervals_error_dependency (pbkrtest masked)", {
  local_mocked_bindings(requireNamespace = function(package, ...) {
    if (identical(package, "pbkrtest")) FALSE else base::requireNamespace(package, ...)
  })
  expect_error(ci_lmm(fit_52(), control = lmm_interval_control(ci_df = "kenward-roger")),
               class = "reintervals_error_dependency")
})

test_that("reintervals_error_internal", {
  expect_error(.abort("boom", class = "reintervals_error_internal"),
               class = "reintervals_error_internal")
})

## ---- one raising site per warning class ---------------------------------------
test_that("reintervals_warning_truncation (a reintervals_warning_fallback)", {
  expect_warning(re_components(components = c(g = -1, Residual = 1), dfs = c(pi = 5),
                               mean = 0, var_mean = 0.1),
                 class = "reintervals_warning_truncation")
  expect_warning(.floor_vg(-1, t0 = 3), class = "reintervals_warning_fallback")
})

test_that("reintervals_warning_singular", {
  skip_if_not_installed("lme4")
  set.seed(4); g <- factor(rep(1:6, each = 3))
  f1 <- lmer_q(y ~ 1 + (1 | g), data.frame(g = g, y = 10 + rnorm(18)))
  skip_if_not(lme4::isSingular(f1))
  expect_warning(ti_lmm(f1, method = "ems-mls"), class = "reintervals_warning_singular")
})

test_that("reintervals_warning_limitation: the two documented 0.2.0 limitations", {
  skip_if_not_installed("lme4")
  ## few top-level groups (<= 5), PI only, random-intercept models only
  set.seed(1); g <- factor(rep(1:5, each = 4))
  d <- data.frame(g = g, y = 10 + rnorm(5, 0, 2)[g] + rnorm(20))
  f5 <- lmer_q(y ~ 1 + (1 | g), d)
  .reset_once_warnings()
  w <- tryCatch(pi_lmm(f5), warning = identity)
  expect_s3_class(w, "reintervals_warning_limitation")
  expect_match(conditionMessage(w), "has 5 levels")
  expect_identical(w$n_top_levels, 5)
  ## at most once per session: a second call is silent until the flag is reset
  .reset_once_warnings()
  expect_warning(pi_lmm(f5), class = "reintervals_warning_limitation")
  expect_warning(pi_lmm(f5), NA)
  expect_warning(pi_lmm(f5, level = 0.9), NA)
  .reset_once_warnings()
  expect_warning(pi_lmm(f5), class = "reintervals_warning_limitation")
  expect_warning(ci_lmm(f5), NA)                                     # CI does not warn
  expect_warning(pi_lmm(fit_52()), NA)                               # 6 runs: no warning
  ## fixed-slope single-observation lots (MBL)
  fm <- lmer_q(Y ~ Month + (1 | Lot), mbl_appendix_a)
  expect_warning(ci_lmm(fm, newdata = data.frame(Month = 12)),
                 class = "reintervals_warning_limitation")
})

## ---- surfaced silent paths ------------------------------------------------------
test_that("pi_df_fallback reaches attr(., 'diagnostics')", {
  skip_if_not_installed("lme4")
  set.seed(4); g <- factor(rep(1:6, each = 3))
  f1 <- lmer_q(y ~ 1 + (1 | g), data.frame(g = g, y = 10 + rnorm(18)))
  skip_if_not(lme4::isSingular(f1))
  expect_true(attr(quiet_df_fallback(pi_lmm(f1)), "diagnostics")$pi_df_fallback)
  expect_false(attr(pi_lmm(fit_52()), "diagnostics")$pi_df_fallback)
})

test_that("random-slope CI / PI / NGM report their df source as df_rule", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (1 + Days | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = 5)
  for (v in list(ci_lmm, pi_lmm, new_group_mean_lmm)) {
    expect_identical(attr(v(fs, newdata = nd), "diagnostics")$df_rule, "satterthwaite")
  }
  set.seed(11); d <- expand.grid(time = c(0, 3, 6, 9, 12), batch = factor(1:6))
  d$y <- 100 - 2 * d$time + rnorm(nrow(d))
  fz <- lmer_q(y ~ time + (time | batch), d)
  skip_if_not(lme4::isSingular(fz))
  expect_identical(attr(quiet_df_fallback(ci_lmm(fz, newdata = data.frame(time = 6))), "diagnostics")$df_rule,
                   "containment")
})

test_that("a non-finite slope df takes the slope engine's containment df (B - 2), not 1e6", {
  skip_if_not_installed("lme4")
  fs <- lmer_q(Reaction ~ Days + (1 + Days | Subject), lme4::sleepstudy)
  comp <- vfun_extract(fs, t0 = 5)
  parts <- attr(comp, "containment")
  expect_identical(c(parts$B, parts$k_re), c(18L, 2L))
  for (kind in c("CI", "new_group_mean")) {                         # between-group pieces
    cp <- comp; cp$dfs[[c(CI = "F", new_group_mean = "G")[[kind]]]] <- Inf
    r <- .vfun_interval(cp, kind = kind)
    expect_identical(r$df, 16, info = kind)                          # B - k_re = 18 - 2
    expect_identical(r$df, .vfun_df_between(parts), info = kind)     # the engine's helper
    expect_identical(r$df_rule, "containment", info = kind)
  }
  cp <- comp; cp$dfs[["F"]] <- Inf
  r <- .vfun_interval(cp, kind = "CI")
  expect_equal(r$upper - r$estimate, stats::qt(0.975, 16) * sqrt(comp$var_mean))
  cp <- comp; cp$dfs[["T"]] <- Inf                                   # PI: engine's w rule
  w <- comp$V_G / comp$V_T
  expect_identical(.vfun_interval(cp, kind = "PI")$df, .vfun_df_containment(parts, w))
  comp$dfs[["F"]] <- Inf
  attr(comp, "containment") <- NULL
  expect_error(.vfun_interval(comp, kind = "CI"), class = "reintervals_error_numerical")
})

test_that(".solve_tau2 returns NA (never an error or a number) for a draw without a solution", {
  ## lambda below phi1 by more than eps -> a + lambda s2e <= 0 at the lower bound
  r <- .solve_tau2(Qi = c(1, 2), lam = c(0.1, 0.5), s2e = c(1, 1), Qstar = c(3, 3),
                   phi1 = 0.5)
  expect_true(all(is.na(r)))
  ok <- .solve_tau2(Qi = c(1, 2), lam = c(0.5, 0.9), s2e = c(1, 1), Qstar = c(3, 3),
                    phi1 = 0.5)
  expect_true(all(is.finite(ok)))
})

test_that("GPQ drops NA draws, reports the count, and errors above 1%", {
  d <- lme4::sleepstudy
  real <- .solve_tau2
  na_frac <- 0.005
  local_mocked_bindings(.solve_tau2 = function(...) {
    out <- real(...); out[seq_len(floor(na_frac * length(out)))] <- NA; out
  })
  g <- ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5,
                  control = lmm_interval_control(seed = 1, M = 2000))
  expect_identical(attr(g, "diagnostics")$n_na_draws, 10L)
  expect_true(is.finite(g$lower) && is.finite(g$upper))
  na_frac <- 0.02
  expect_error(ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5,
                          control = lmm_interval_control(seed = 1, M = 2000)),
               class = "reintervals_error_numerical")
})

## ---- every user-facing error / warning in R/ is classed ---------------------------
test_that("R/ has no bare stop()/warning(); every .abort/.warn names a reintervals class", {
  rdir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir), "package sources not available (installed tests)")
  files <- setdiff(list.files(rdir, "\\.R$", full.names = TRUE),
                   file.path(rdir, "conditions.R"))
  src <- unlist(lapply(files, function(f) {
    l <- readLines(f); l <- sub("#.*$", "", l)     # crude comment strip (no '#' in our strings matters)
    paste(basename(f), seq_along(l), l, sep = ":")
  }))
  bare <- grep("(^|[^.[:alnum:]_])(stop|warning)\\(", sub("^[^:]*:[0-9]*:", "", src))
  expect_identical(src[bare], character(0))
  ## every helper call carries a class from the registered hierarchy
  txt <- paste(vapply(files, function(f) paste(readLines(f), collapse = "\n"), ""), collapse = "\n")
  calls <- regmatches(txt, gregexpr("\\.(abort|warn)\\(", txt))[[1]]
  classes <- regmatches(txt, gregexpr('class = "reintervals_(error|warning)[a-z_]*"', txt))[[1]]
  expect_gte(length(classes), length(calls))
  used <- unique(sub('class = "(.*)"', "\\1", classes))
  expect_true(all(used %in% names(.reintervals_condition_parents)),
              info = paste(setdiff(used, names(.reintervals_condition_parents)), collapse = ", "))
})
