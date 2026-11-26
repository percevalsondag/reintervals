# Step-13b validation arms (non-default; under validation):
#   pi_df = "capped" / "ratio-guard" (+ pi_df_ratio), ci_df = "capped-kr",
#   the two-sided MLS TI quantities ti_R / ti_df_T, the extrapolation warning
#   (ti_warn_R, ti_warn_dfT; on by default) and reroute = TRUE (ti_reroute_R,
#   ti_reroute_dfT).
# Every expected value below is computed in this file from the default
# interval's own centre and SE, the closed-form containment df, the REML
# variance components or the raw data, not from the arm's code path.

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))
## value and every warning (class, message, condition), each muffled
run_w <- function(expr) {
  .reset_once_warnings()
  w <- list()
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w[[length(w) + 1L]] <<- cnd; invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w, classes = vapply(w, function(x) class(x)[1], ""))
}
is_ti_warn <- function(w) inherits(w, "reintervals_warning_limitation") &&
  identical(w$limitation, "ti_mls_extrapolation")
n_ti_warn <- function(r) sum(vapply(r$warnings, is_ti_warn, logical(1)))
dg <- function(x) attr(x, "diagnostics")
qside <- function(s, level = 0.95) ifelse(s == "two.sided", 1 - (1 - level) / 2, level)
## half-width of an interval on its finite side(s)
hw_of <- function(x) ifelse(x$sides == "upper", x$upper - x$estimate, x$estimate - x$lower)
## limits from a centre, SE and df
lims <- function(est, se, df, sides, level = 0.95) {
  h <- se * stats::qt(qside(sides, level), df)
  list(lower = ifelse(sides == "upper", -Inf, est - h), upper = ifelse(sides == "lower", Inf, est + h))
}

## balanced fixed-slope stability design (EMS route): K = 6 lots, 5 times
fs_data <- function(seed = 5, K = 6, tt = c(0, 3, 6, 9, 12)) {
  set.seed(seed)
  d <- expand.grid(time = tt, lot = factor(seq_len(K)))
  d$y <- 100 - 0.4 * d$time + stats::rnorm(K, 0, 1.5)[d$lot] + stats::rnorm(nrow(d), 0, 0.8)
  d
}
## nested two-factor design (not GPQ-eligible)
nested_data <- function(seed = 1, A = 5, B = 3, n = 2) {
  set.seed(seed)
  a <- factor(rep(seq_len(A), each = B * n)); b <- factor(rep(seq_len(A * B), each = n))
  data.frame(a = a, b = b, y = 10 + stats::rnorm(A, 0, 2)[a] + stats::rnorm(A * B, 0, 1)[b] +
               stats::rnorm(A * B * n))
}
SIDES <- c("two.sided", "lower", "upper")

## ---- lmm_interval_control --------------------------------------------------
test_that("the step-13b controls: defaults and validation", {
  ctl <- lmm_interval_control()
  expect_identical(ctl$pi_df_ratio, 1.5)
  expect_identical(ctl$ti_warn_R, 0.25); expect_identical(ctl$ti_warn_dfT, 40)
  expect_false(ctl$reroute)
  expect_identical(ctl$ti_reroute_R, 0.25); expect_identical(ctl$ti_reroute_dfT, 20)
  expect_identical(ctl$pi_df, "francq-satterthwaite")
  expect_identical(ctl$ci_df, "satterthwaite")
  expect_identical(lmm_interval_control(pi_df = "cap")$pi_df, "capped")
  expect_identical(lmm_interval_control(pi_df = "ratio")$pi_df, "ratio-guard")
  expect_identical(lmm_interval_control(ci_df = "capped")$ci_df, "capped-kr")
  for (bad in list(1, 0.5, NA_real_, Inf, "2", c(2, 3)))
    expect_error(lmm_interval_control(pi_df_ratio = bad), class = "reintervals_error_input")
  expect_identical(lmm_interval_control(pi_df_ratio = 1 + 1e-9)$pi_df_ratio, 1 + 1e-9)
  for (nm in c("ti_warn_R", "ti_warn_dfT", "ti_reroute_R", "ti_reroute_dfT")) {
    for (bad in list(-1, NA_real_, "1", c(1, 2)))
      expect_error(do.call(lmm_interval_control, stats::setNames(list(bad), nm)),
                   class = "reintervals_error_input")
    expect_identical(do.call(lmm_interval_control, stats::setNames(list(Inf), nm))[[nm]], Inf)
  }
  for (bad in list(NA, "TRUE", 1, c(TRUE, FALSE)))
    expect_error(lmm_interval_control(reroute = bad), class = "reintervals_error_input")
  ## a hand-built control without the new fields gets the defaults
  expect_identical(.ctl_value(list(), "pi_df_ratio"), 1.5)
  expect_false(.ctl_value(list(), "reroute"))
})

## ---- PI: pi_df = "capped" / "ratio-guard" ----------------------------------
test_that("pi_df = 'capped': Francq centre and SE with df = min(df_F, df_C)", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 4.5, 9))
  base <- quiet_limitation(pi_lmm(fm, nd, sides = SIDES))
  df_F <- dg(base)$df
  df_C <- nlevels(lme4::sleepstudy$Subject) - 1        # one-way containment: K - 1 = 17
  se <- hw_of(base) / stats::qt(qside(base$sides), df_F)
  r <- run_w(pi_lmm(fm, nd, sides = SIDES, control = lmm_interval_control(pi_df = "capped")))
  x <- r$value
  expect_true("reintervals_warning_unvalidated" %in% r$classes)
  ref <- lims(base$estimate, se, pmin(df_F, df_C), base$sides)
  expect_identical(x$estimate, base$estimate)
  expect_equal(x$lower, ref$lower, tolerance = 1e-10)
  expect_equal(x$upper, ref$upper, tolerance = 1e-10)
  expect_identical(dg(x)$df, pmin(df_F, df_C))
  expect_identical(dg(x)$df_F, df_F); expect_identical(dg(x)$df_C, rep(df_C, 9))
  expect_identical(dg(x)$pi_guard_fired, rep(TRUE, 9))
  expect_identical(dg(x)$pi_df_rule, rep("capped", 9))
  expect_true(all(df_F > df_C))                          # the cap binds here
  ## default and containment rows report df_F / df_C, no guard
  expect_identical(dg(base)$df_C, rep(df_C, 9)); expect_identical(dg(base)$df_F, df_F)
  expect_identical(dg(base)$pi_guard_fired, rep(NA, 9))
})

test_that("pi_df = 'capped' / 'ratio-guard' keep the Francq df when it is below the containment df", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  base <- quiet_limitation(pi_lmm(fm, nd))
  local_mocked_bindings(.ddf_containment = function(model) 1000)   # df_C > df_F
  for (arm in c("capped", "ratio-guard")) {
    x <- suppressWarnings(pi_lmm(fm, nd, control = lmm_interval_control(pi_df = arm)))
    expect_identical(x$lower, base$lower, info = arm)
    expect_identical(x$upper, base$upper, info = arm)
    expect_identical(dg(x)$df, dg(base)$df, info = arm)
    expect_identical(dg(x)$df_C, c(1000, 1000), info = arm)
  }
  x <- suppressWarnings(pi_lmm(fm, nd, control = lmm_interval_control(pi_df = "capped")))
  expect_identical(dg(x)$pi_guard_fired, c(TRUE, TRUE))    # fires; min() leaves df_F
  x <- suppressWarnings(pi_lmm(fm, nd, control = lmm_interval_control(pi_df = "ratio-guard")))
  expect_identical(dg(x)$pi_guard_fired, c(FALSE, FALSE))
})

test_that("pi_df = 'ratio-guard': fires at df_F / df_C >= pi_df_ratio, equality included", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  base <- quiet_limitation(pi_lmm(fm, nd, sides = SIDES))
  df_F <- dg(base)$df[1]; df_C <- 17
  ratio <- df_F / df_C
  expect_gt(ratio, 1)
  se <- hw_of(base) / stats::qt(qside(base$sides), dg(base)$df)
  run_ratio <- function(r) suppressWarnings(pi_lmm(fm, nd, sides = SIDES,
    control = lmm_interval_control(pi_df = "ratio-guard", pi_df_ratio = r)))
  ## at equality: fires, df = df_C
  x <- run_ratio(ratio)
  expect_identical(dg(x)$pi_guard_fired, rep(TRUE, 6))
  ref <- lims(base$estimate, se, df_C, base$sides)
  expect_equal(x$lower, ref$lower, tolerance = 1e-10)
  expect_equal(x$upper, ref$upper, tolerance = 1e-10)
  expect_identical(dg(x)$df, rep(df_C, 6))
  ## just above: does not fire, the Francq interval unchanged
  x <- run_ratio(ratio * (1 + 1e-12))
  expect_identical(dg(x)$pi_guard_fired, rep(FALSE, 6))
  expect_identical(x$lower, base$lower); expect_identical(x$upper, base$upper)
  expect_identical(dg(x)$df, dg(base)$df)
  ## the default threshold 1.5 and the Francq fallback (singular fit: df_F = df_C)
  expect_identical(dg(run_ratio(1.5))$pi_guard_fired, rep(ratio >= 1.5, 6))
  set.seed(4); d <- data.frame(g = factor(rep(1:5, each = 4)), y = stats::rnorm(20))
  fs <- lmer_q(y ~ 1 + (1 | g), d)
  skip_if_not(isTRUE(lme4::isSingular(fs)))
  x <- suppressWarnings(pi_lmm(fs, control = lmm_interval_control(pi_df = "ratio-guard")))
  expect_identical(dg(x)$df_F, dg(x)$df_C)
  expect_false(dg(x)$pi_guard_fired)
})

## ---- CI: ci_df = "capped-kr" -----------------------------------------------
test_that("ci_df = 'capped-kr': KR SE with df = min(df_KR, df_C); cap_binds iff df_C < df_KR", {
  skip_if_not_installed("pbkrtest")
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 4.5, 30))
  kr <- quiet_limitation(ci_lmm(fm, nd, sides = SIDES,
                                control = lmm_interval_control(ci_df = "kenward-roger")))
  df_KR <- dg(kr)$df; df_C <- 17
  se <- hw_of(kr) / stats::qt(qside(kr$sides), df_KR)
  r <- run_w(ci_lmm(fm, nd, sides = SIDES, control = lmm_interval_control(ci_df = "capped-kr")))
  x <- r$value
  expect_true("reintervals_warning_unvalidated" %in% r$classes)
  ref <- lims(kr$estimate, se, pmin(df_KR, df_C), kr$sides)
  expect_identical(x$estimate, kr$estimate)
  expect_equal(x$lower, ref$lower, tolerance = 1e-10)
  expect_equal(x$upper, ref$upper, tolerance = 1e-10)
  expect_identical(dg(x)$df, pmin(df_KR, df_C))
  expect_identical(dg(x)$df_KR, df_KR); expect_identical(dg(x)$df_C, rep(df_C, 9))
  expect_identical(dg(x)$ci_cap_binds, df_C < df_KR)
  expect_identical(dg(x)$df_method, ifelse(df_C < df_KR, "containment", "kenward-roger"))
  expect_identical(dg(x)$ci_df_rule, rep("capped-kr", 9))
  ## other arms: df_KR is the KR df on the KR arm, NA elsewhere; no cap
  expect_identical(dg(kr)$df_KR, df_KR); expect_identical(dg(kr)$ci_cap_binds, rep(NA, 9))
  d0 <- dg(ci_lmm(fm, nd[1:2, , drop = FALSE]))
  expect_identical(d0$df_KR, c(NA_real_, NA_real_)); expect_identical(d0$df_C, c(17, 17))
})

test_that("ci_df = 'capped-kr' equals the KR interval when the cap does not bind (equality included)", {
  skip_if_not_installed("pbkrtest")
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = 0)
  kr <- quiet_limitation(ci_lmm(fm, nd, sides = SIDES,
                                control = lmm_interval_control(ci_df = "kenward-roger")))
  df_KR <- dg(kr)$df[1]
  for (dc in c(df_KR, 1000)) {          # df_C = df_KR exactly, and df_C > df_KR
    local_mocked_bindings(.ddf_containment = function(model) dc)
    x <- suppressWarnings(ci_lmm(fm, nd, sides = SIDES,
                                 control = lmm_interval_control(ci_df = "capped-kr")))
    expect_identical(x$lower, kr$lower, info = dc); expect_identical(x$upper, kr$upper, info = dc)
    expect_identical(dg(x)$ci_cap_binds, rep(FALSE, 3), info = dc)
    expect_identical(dg(x)$df_method, rep("kenward-roger", 3), info = dc)
  }
})

## ---- TI: ti_R and ti_df_T --------------------------------------------------
test_that("ti_R and ti_df_T on the EMS route match the balanced closed form", {
  d <- fs_data(); fm <- lmer_q(y ~ time + (1 | lot), d)
  nd <- data.frame(time = c(0, 12, 30))
  x <- quiet_limitation(ti_lmm(fm, nd, P = c(0.9, 0.99), level = c(0.9, 0.95), sides = SIDES))
  expect_identical(unique(x$method), "ems-mls")
  vc <- as.data.frame(lme4::VarCorr(fm)); s2L <- vc$vcov[1]; s2e <- vc$vcov[2]
  K <- 6; J <- 5; N <- 30
  X0 <- cbind(1, nd$time); var_mean <- diag(X0 %*% as.matrix(stats::vcov(fm)) %*% t(X0))
  sig2 <- s2L + s2e                         # sum k_j EMS_j = (s2e + J s2L) / J + (1 - 1/J) s2e
  R_ref <- var_mean / sig2
  dfT_ref <- sig2^2 / (((s2e + J * s2L) / J)^2 / (K - 1) + ((1 - 1 / J) * s2e)^2 / (N - K - 1))
  two <- x$sides == "two.sided"
  ev <- rep_len(1:3, nrow(x))
  expect_equal(dg(x)$ti_R[two], R_ref[ev[two]], tolerance = 1e-10)
  expect_equal(dg(x)$ti_df_T[two], rep(dfT_ref, sum(two)), tolerance = 1e-10)
  expect_true(all(is.na(dg(x)$ti_R[!two])) && all(is.na(dg(x)$ti_df_T[!two])))
  expect_identical(dg(x)$ti_warn_fired[!two], rep(NA, sum(!two)))
  expect_identical(dg(x)$rerouted, ifelse(two, FALSE, NA))
  expect_identical(dg(x)$reroute_reason, ifelse(two, NA_character_, "not_applicable"))
})

test_that("ti_R = 1 / n_E and ti_df_T on the anova-mls routes (Hoffman-Kringle, Montes)", {
  ## intercept-only one-way, balanced 6 x 3 (Hoffman-Kringle)
  fm <- lmer_q(conc ~ 1 + (1 | run), hk_assay)
  x <- ti_lmm(fm, P = 0.95, level = 0.9)
  expect_identical(x$method, "anova-mls")
  I <- 6; J <- 3
  ybar <- tapply(hk_assay$conc, hk_assay$run, mean)
  S2L <- stats::var(ybar)
  S2E <- sum((hk_assay$conc - ybar[hk_assay$run])^2) / (I * (J - 1))
  R_ref <- (S2L / I) / (S2L + S2E * (1 - 1 / J))
  vc <- as.data.frame(lme4::VarCorr(fm)); s2L <- vc$vcov[1]; s2e <- vc$vcov[2]
  a <- c((s2e + J * s2L) / J, (1 - 1 / J) * s2e)
  dfT_ref <- sum(a)^2 / sum(a^2 / c(I - 1, I * (J - 1)))
  expect_equal(dg(x)$ti_R, R_ref, tolerance = 1e-10)
  expect_equal(dg(x)$ti_df_T, dfT_ref, tolerance = 1e-10)
  expect_identical(dg(x)$ti_R, 1 / dg(x)$n_E)
  ## forced ems-mls on the same fit: R = var_mean / sigma2_t
  e <- ti_lmm(fm, P = 0.95, level = 0.9, method = "ems-mls")
  expect_equal(dg(e)$ti_R, stats::vcov(fm)[1, 1] / (s2L + s2e), tolerance = 1e-10)
  expect_equal(dg(e)$ti_df_T, dfT_ref, tolerance = 1e-10)
  ## staggered fixed slope (Montes route): R = VhatYh / Vhat_Y from the raw data
  d <- mbl_appendix_a; f2 <- lmer_q(Y ~ Month + (1 | Lot), d)
  nd <- data.frame(Month = c(0, 12, 36))
  y <- quiet_limitation(ti_lmm(f2, nd, P = 0.99, level = c(0.9, 0.95)))
  expect_identical(unique(y$method), "anova-mls")
  ref <- montes_ref(d$Y, d$Month, d$Lot, nd$Month)
  expect_equal(dg(y)$ti_R, rep(ref$Vhat_Yhat / ref$Vhat_Y, 2), tolerance = 1e-10)
})

test_that("ti_R and ti_df_T reproduce three step-13 training replicates (to 1e-10)", {
  skip_on_cran()
  skip_if_not_installed("digest")
  root <- path.expand(Sys.getenv("REINTERVALS_SIMSTUDY", "~/dev/reintervals-simstudy"))
  H <- file.path(root, "step13", "harness")
  D <- file.path(root, "step13", "analysis", "predicates", "data")
  skip_if_not(file.exists(file.path(D, "MANIFEST.csv")) && dir.exists(H),
              "step-13 training data not available")
  man <- utils::read.csv(file.path(D, "MANIFEST.csv"), stringsAsFactors = FALSE)
  cells <- utils::read.csv(file.path(H, "cells.csv"), stringsAsFactors = FALSE)
  env <- new.env()
  for (f in c("generators.R", "seeds.R")) sys.source(file.path(H, "R", f), envir = env)
  ## the harness switches the RNG kind: restore the caller's state afterwards
  old_kind <- RNGkind(); old_seed <- get0(".Random.seed", globalenv(), inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, globalenv())
  }, add = TRUE)
  pn <- function(s) as.numeric(strsplit(s, ";", fixed = TRUE)[[1]])
  ## one replicate (rep 1, every t0) of an ems-mls cell and two anova-mls cells
  for (id in c("G4_s12_rho=0.05_full_t12-24-48", "G4_s12_rho=0.05_staggered_t12-24-48",
               "G4_s12_rho=0.05_single_obs_t12-24-48")) {
    m <- man[man$family == "main" & man$cell_id == id, ]
    expect_identical(nrow(m), 1L)
    p <- file.path(D, m$file)
    ## the file is the one the MANIFEST lists (digest is not a declared
    ## dependency: reached through getExportedValue after the skip above)
    sha <- getExportedValue("digest", "digest")(file = p, algo = "sha256")
    expect_identical(sha, m$sha256)
    g <- readRDS(p)
    st <- env$cell_streams(id, 1, "main")
    expect_identical(env$stream_id(st$data[[1]]), g$fits$data_stream[g$fits$rep == 1][1])
    env$use_stream(st$data[[1]])
    gen <- env$generate(as.list(cells[cells$cell_id == id, ]))
    fm <- lmer_q(gen$f, gen$d, REML = TRUE)
    ## rep 1 of these cells is a singular fit (its warning is expected)
    x <- suppressWarnings(ti_lmm(fm, gen$nd, P = 0.95, level = 0.95))
    rows <- g$rows[g$rows$rep == 1 & g$rows$verb == "TI" & g$rows$side == "two.sided" &
                     g$rows$P == 0.95 & g$rows$level == 0.95, ]
    rows <- rows[match(gen$nd$Month, rows$t0), ]
    fits <- g$fits[g$fits$rep == 1, ]; fits <- fits[match(gen$nd$Month, fits$t0), ]
    ## the same dataset: the stored limits are reproduced exactly
    expect_identical(x$lower, rows$lower, info = id); expect_identical(x$upper, rows$upper, info = id)
    ## prep.R: R = R_ems (ems-mls) or 1 / n_E (anova-mls); df_T_hat from the
    ## stored EMS terms (10 significant digits)
    R_ref <- if (x$method[1] == "anova-mls") 1 / rows$dg_n_E else fits$R_ems
    dfT_ref <- mapply(function(ms, k, d) { a <- pn(k) * pn(ms); sum(a)^2 / sum(a^2 / pn(d)) },
                      fits$ems_ms, fits$ems_k, fits$ems_df)
    expect_equal(dg(x)$ti_R, R_ref, tolerance = 1e-10, info = id)
    expect_equal(dg(x)$ti_df_T, unname(dfT_ref), tolerance = 1e-10, info = id)
  }
})

## ---- TI: extrapolation warning ---------------------------------------------
test_that("the TI extrapolation warning: fires once per call when R > ti_warn_R and df_T_hat >= ti_warn_dfT", {
  d <- fs_data(); fm <- lmer_q(y ~ time + (1 | lot), d)
  nd <- data.frame(time = c(0, 12, 30))
  base <- run_w(ti_lmm(fm, nd, P = c(0.9, 0.99), sides = SIDES))
  R <- dg(base$value)$ti_R; dfT <- dg(base$value)$ti_df_T[1]
  two <- base$value$sides == "two.sided"
  ## thresholds straddling this fit's quantities, equality on df
  r_cut <- sort(unique(R[two]))[2]           # fires at the largest R only
  r <- run_w(ti_lmm(fm, nd, P = c(0.9, 0.99), sides = SIDES,
                    control = lmm_interval_control(ti_warn_R = r_cut, ti_warn_dfT = dfT)))
  expect_identical(n_ti_warn(r), 1L)
  fired <- dg(r$value)$ti_warn_fired
  expect_identical(fired, ifelse(two, R > r_cut, NA))
  w <- Filter(is_ti_warn, r$warnings)[[1]]
  expect_match(conditionMessage(w), "provisional until step 13b", fixed = TRUE)
  expect_match(conditionMessage(w), "t0 = 30", fixed = TRUE)
  expect_true(w$provisional)
  ## numbers are unchanged by the warning
  expect_identical(r$value$lower, base$value$lower); expect_identical(r$value$upper, base$value$upper)
  ## just above df_T_hat: no warning; Inf disables it
  r <- run_w(ti_lmm(fm, nd, control = lmm_interval_control(ti_warn_R = 0, ti_warn_dfT = dfT * (1 + 1e-12))))
  expect_identical(n_ti_warn(r), 0L)
  r <- run_w(ti_lmm(fm, nd, control = lmm_interval_control(ti_warn_R = Inf, ti_warn_dfT = 0)))
  expect_identical(n_ti_warn(r), 0L)
  ## one-sided calls never warn and carry NA
  r <- run_w(ti_lmm(fm, nd, sides = "upper", control = lmm_interval_control(ti_warn_R = 0, ti_warn_dfT = 0)))
  expect_identical(n_ti_warn(r), 0L)
  expect_identical(dg(r$value)$ti_warn_fired, rep(NA, 3))
  ## the defaults do not fire on this design (R at t0 = 30 below 0.25)
  expect_identical(n_ti_warn(base), 0L)
  expect_true(all(R[two] <= 0.25) || dfT < 40)
})

test_that("the GPQ routes carry NA / not_applicable step-13b TI fields and never warn", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  ctl <- lmm_interval_control(seed = 1, M = 2000, ti_warn_R = 0, ti_warn_dfT = 0, reroute = TRUE)
  r <- run_w(ti_lmm(fm, data.frame(Days = 30), method = "gpq", control = ctl))
  expect_identical(n_ti_warn(r), 0L)
  fs <- lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  r2 <- run_w(ti_lmm(fs, data.frame(Days = 5), control = ctl))
  for (x in list(r$value, r2$value)) {
    expect_identical(dg(x)$ti_R, NA_real_); expect_identical(dg(x)$ti_df_T, NA_real_)
    expect_identical(dg(x)$ti_warn_fired, NA); expect_identical(dg(x)$rerouted, NA)
    expect_identical(dg(x)$reroute_reason, "not_applicable")
  }
  ## reroute does not apply: no reroute warning
  rr <- vapply(r2$warnings, function(w) identical(w$option, "reroute"), logical(1))
  expect_false(any(rr))
})

## ---- TI: reroute -----------------------------------------------------------
test_that("reroute = TRUE: qualifying two-sided rows equal the explicit method = 'gpq', gpq_two_sided = 'content' call", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9, 30))
  P <- c(0.9, 0.99); lv <- c(0.9, 0.95)
  base <- quiet_limitation(ti_lmm(fm, nd, P = P, level = lv, sides = SIDES))
  R <- dg(base)$ti_R; dfT <- dg(base)$ti_df_T[1]
  two <- base$sides == "two.sided"
  r_cut <- sort(unique(R[two]))[2]           # reroutes the t0 = 30 rows only
  ctl <- lmm_interval_control(seed = 11, M = 2000, reroute = TRUE, ti_reroute_R = r_cut,
                              ti_reroute_dfT = dfT)
  r <- run_w(ti_lmm(fm, nd, P = P, level = lv, sides = SIDES, control = ctl))
  x <- r$value
  ## one reroute warning, the GPQ route's own unvalidated warnings muffled
  opts <- unlist(lapply(r$warnings, function(w)
    if (inherits(w, "reintervals_warning_unvalidated")) w$option))
  expect_identical(opts, "reroute")
  cand <- two & R > r_cut
  expect_identical(which(cand), which(two & rep_len(1:3, nrow(x)) == 3))
  expect_identical(x$method, ifelse(cand, "gpq", "ems-mls"))
  expect_identical(dg(x)$rerouted, ifelse(two, cand, NA))
  expect_identical(dg(x)$reroute_reason,
                   ifelse(cand, "R > ti_reroute_R and df_T_hat >= ti_reroute_dfT",
                          ifelse(two, NA_character_, "not_applicable")))
  expect_identical(dg(x)$gpq_seed, 11L)
  expect_identical(dg(x)$route_reason[cand], rep("reroute", sum(cand)))
  expect_identical(dg(x)$gpq_two_sided, ifelse(cand, "content", "not_applicable"))
  ## the other rows are the default rows, bit for bit
  expect_identical(x$lower[!cand], base$lower[!cand]); expect_identical(x$upper[!cand], base$upper[!cand])
  ## each rerouted row: the explicit scalar GPQ call with the same seed and M
  for (i in which(cand)) {
    e <- suppressWarnings(ti_lmm(fm, nd[3, , drop = FALSE], P = x$P[i], level = x$level[i],
                                 method = "gpq", control = lmm_interval_control(
                                   seed = 11, M = 2000, gpq_two_sided = "content")))
    expect_identical(c(x$estimate[i], x$lower[i], x$upper[i]), c(e$estimate, e$lower, e$upper))
    expect_identical(dg(x)$mc_se_lower[i], dg(e)$mc_se_lower)
  }
  ## R and df_T_hat are those of the MLS route, on rerouted rows too
  expect_identical(dg(x)$ti_R, R); expect_identical(dg(x)$ti_df_T, dg(base)$ti_df_T)
})

test_that("reroute = TRUE with seed = NULL: one seed from the caller's stream, reproducible from gpq_seed", {
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 30))
  ctl <- lmm_interval_control(M = 2000, reroute = TRUE, ti_reroute_R = 0, ti_reroute_dfT = 0)
  set.seed(42)
  x <- suppressWarnings(ti_lmm(fm, nd, control = ctl))
  expect_identical(dg(x)$rerouted, c(TRUE, TRUE))
  e <- suppressWarnings(ti_lmm(fm, nd, method = "gpq", control = lmm_interval_control(
    seed = dg(x)$gpq_seed, M = 2000, gpq_two_sided = "content")))
  expect_identical(x$lower, e$lower); expect_identical(x$upper, e$upper)
  set.seed(42)
  e2 <- suppressWarnings(ti_lmm(fm, nd, method = "gpq",
                                control = lmm_interval_control(M = 2000, gpq_two_sided = "content")))
  expect_identical(x$lower, e2$lower)
})

test_that("reroute = TRUE does not fire on an ineligible (nested) fit or below the thresholds; the warning still fires", {
  fn <- lmer_q(y ~ 1 + (1 | a/b), nested_data())
  ctl <- lmm_interval_control(reroute = TRUE, ti_reroute_R = 0, ti_reroute_dfT = 0,
                              ti_warn_R = 0, ti_warn_dfT = 0)
  base <- quiet_limitation(ti_lmm(fn))
  r <- run_w(ti_lmm(fn, control = ctl))
  expect_identical(r$value$method, "ems-mls")
  expect_identical(r$value$lower, base$lower); expect_identical(r$value$upper, base$upper)
  expect_false(dg(r$value)$rerouted)
  expect_identical(dg(r$value)$reroute_reason, "not GPQ-eligible")
  expect_identical(n_ti_warn(r), 1L)
  expect_true("reintervals_warning_unvalidated" %in% r$classes)
  ## eligible fit, below the reroute thresholds: nothing rerouted, warning kept
  fm <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(0, 9))
  b2 <- quiet_limitation(ti_lmm(fm, nd))
  r2 <- run_w(ti_lmm(fm, nd, control = lmm_interval_control(
    reroute = TRUE, ti_reroute_R = Inf, ti_warn_R = 0, ti_warn_dfT = 0)))
  expect_identical(r2$value$lower, b2$lower); expect_identical(r2$value$upper, b2$upper)
  expect_identical(dg(r2$value)$rerouted, c(FALSE, FALSE))
  expect_identical(dg(r2$value)$reroute_reason, c(NA_character_, NA_character_))
  expect_identical(n_ti_warn(r2), 1L)
  ## the default thresholds: R on sleepstudy at Days 0..9 is small
  r3 <- run_w(ti_lmm(fm, nd, control = lmm_interval_control(reroute = TRUE)))
  expect_identical(dg(r3$value)$rerouted, c(FALSE, FALSE))
  expect_true(all(dg(r3$value)$ti_R <= 0.25))
  ## one-sided rows are never rerouted
  r4 <- run_w(ti_lmm(fm, nd, sides = "lower", control = ctl))
  expect_identical(dg(r4$value)$rerouted, c(NA, NA))
  expect_identical(r4$value$method, c("ems-mls", "ems-mls"))
})

## ---- scale equivariance ----------------------------------------------------
## y -> c y and t -> c t at the transformed theta (no optimization; see
## test-scale-equivariance.R): limits scale by c (y) or are unchanged (t);
## every df, R, df_T_hat and flag is unchanged.
test_that("the step-13b arms are equivariant under y -> c y and invariant under t -> c t", {
  skip_on_cran()
  d <- fs_data()
  f1 <- lmer_q(y ~ time + (1 | lot), d)
  th <- lme4::getME(f1, "theta")
  refit <- function(dd) lmer_q(y ~ time + (1 | lot), dd, start = list(theta = th),
                               control = lme4::lmerControl(optimizer = NULL))
  nd <- data.frame(time = c(0, 12, 30))
  calls <- function(fm, nd) {
    tictl <- lmm_interval_control(seed = 3, M = 2000, reroute = TRUE, ti_reroute_R = 0.2,
                                  ti_reroute_dfT = 0, ti_warn_R = 0.2, ti_warn_dfT = 0)
    out <- list(
      pi_capped = run_w(pi_lmm(fm, nd, sides = SIDES, control = lmm_interval_control(pi_df = "capped"))),
      pi_ratio = run_w(pi_lmm(fm, nd, sides = SIDES,
                                control = lmm_interval_control(pi_df = "ratio-guard", pi_df_ratio = 1.1))),
      ti_default = run_w(ti_lmm(fm, nd, P = c(0.9, 0.99), sides = SIDES)),
      ti_reroute = run_w(ti_lmm(fm, nd, P = 0.99, sides = SIDES, control = tictl)))
    if (requireNamespace("pbkrtest", quietly = TRUE))
      out$ci_capped <- run_w(ci_lmm(fm, nd, sides = SIDES, control = lmm_interval_control(ci_df = "capped-kr")))
    out
  }
  base <- calls(f1, nd)
  expect_true(any(dg(base$ti_reroute$value)$rerouted %in% TRUE))
  expect_true(any(dg(base$ti_reroute$value)$rerouted %in% FALSE))
  expect_true(any(dg(base$pi_ratio$value)$pi_guard_fired %in% TRUE))
  rel <- function(a, b, tol) all((is.na(a) & is.na(b)) | (is.infinite(a) & a == b) |
                                   abs(a - b) <= tol * pmax(abs(b), .Machine$double.xmin))
  check <- function(s, b, k, info) {
    expect_identical(s$classes, b$classes, info = info)
    x <- s$value; y <- b$value
    gpq <- x$method == "gpq"
    expect_identical(x$method, y$method, info = info)
    for (cn in c("estimate", "lower", "upper"))
      expect_true(rel(x[[cn]], y[[cn]] * k, 1e-8), info = paste(info, cn))
    dx <- dg(x); dy <- dg(y)
    for (cn in c("df", "df_F", "df_C", "df_KR", "ti_R", "ti_df_T"))
      if (!is.null(dy[[cn]])) expect_true(rel(dx[[cn]], dy[[cn]], 1e-6), info = paste(info, cn))
    for (cn in c("pi_guard_fired", "ci_cap_binds", "ti_warn_fired", "rerouted", "reroute_reason",
                 "df_method", "gpq_seed"))
      expect_identical(dx[[cn]], dy[[cn]], info = paste(info, cn))
  }
  for (c in c(1e-4, 1e-2, 1e2, 1e4)) {
    dd <- d; dd$y <- dd$y * c
    sc <- calls(refit(dd), nd)
    for (v in names(base)) check(sc[[v]], base[[v]], c, sprintf("%s y x %g", v, c))
  }
  for (c in c(1e-2, 1 / 12, 12, 1e2)) {
    dd <- d; dd$time <- dd$time * c
    sc <- calls(refit(dd), transform(nd, time = time * c))
    for (v in names(base)) check(sc[[v]], base[[v]], 1, sprintf("%s t x %g", v, c))
  }
})
