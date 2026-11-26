# The documented limitation of one-sided GPQ tolerance limits (validation step
# 13, random-slope designs): a reintervals_warning_limitation, once per call,
# on one-sided GPQ tolerance-interval rows (a) at a confidence below 0.95,
# both targets, and (b) at 0.95 or above with over = "group_mean"; the per-row
# diagnostics field limitation_flag ("gpq_one_sided_lowconf",
# "gpq_one_sided_group_mean", NA). Two-sided rows never warn. The
# random-intercept GPQ route keeps its reintervals_warning_unvalidated as well.

lmer_q <- function(f, d, ...) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, ...)))
dg <- function(x) attr(x, "diagnostics")
ctl <- lmm_interval_control(seed = 1, M = 2000)
msg_low <- paste0("In validation (step 13, random-slope designs), one-sided GPQ limits at ",
                  "90% confidence covered 0.879 to 0.90 (median 0.894 to 0.903).")
msg_gm <- paste0("In validation (step 13), one-sided batch-mean GPQ limits at 95% ",
                 "confidence covered 0.936 to 0.95 (median 0.948).")
msg_two <- "two-sided GPQ intervals are not affected at these levels"

## every warning of a call (muffled), with its value
collect <- function(expr) {
  ws <- list()
  v <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w; invokeRestart("muffleWarning")
  })
  list(value = v, limitation = Filter(function(w) inherits(w, "reintervals_warning_limitation"), ws),
       unvalidated = Filter(function(w) inherits(w, "reintervals_warning_unvalidated"), ws))
}

rs_fit <- function() lmer_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
nd <- data.frame(Days = c(0, 5))

test_that("one-sided GPQ TI at a confidence below 0.95 warns once (both targets) with the validation figures", {
  skip_if_not_installed("lme4")
  fm <- rs_fit()
  for (ov in c("observation", "group_mean")) for (sd in c("lower", "upper")) {
    r <- collect(ti_lmm(fm, nd, over = ov, sides = sd, level = 0.90, control = ctl))
    info <- paste(ov, sd)
    expect_length(r$limitation, 1L)
    w <- r$limitation[[1]]
    expect_identical(w$limitation, "gpq_one_sided_lowconf", info = info)
    expect_match(conditionMessage(w), msg_low, fixed = TRUE, info = info)
    expect_match(conditionMessage(w), msg_two, ignore.case = TRUE, info = info)
    expect_identical(dg(r$value)$limitation_flag, rep("gpq_one_sided_lowconf", 2), info = info)
  }
})

test_that("one-sided batch-mean GPQ TI at 0.95 or above warns once with the validation figures", {
  skip_if_not_installed("lme4")
  fm <- rs_fit()
  for (lv in c(0.95, 0.99)) for (sd in c("lower", "upper")) {
    r <- collect(ti_lmm(fm, nd, over = "group_mean", sides = sd, level = lv, control = ctl))
    info <- paste(lv, sd)
    expect_length(r$limitation, 1L)
    w <- r$limitation[[1]]
    expect_identical(w$limitation, "gpq_one_sided_group_mean", info = info)
    expect_match(conditionMessage(w), msg_gm, fixed = TRUE, info = info)
    expect_match(conditionMessage(w), msg_two, ignore.case = TRUE, info = info)
    expect_identical(dg(r$value)$limitation_flag, rep("gpq_one_sided_group_mean", 2), info = info)
  }
  ## alpha = 0.05 is the same 0.95
  r <- collect(ti_lmm(fm, nd, over = "group_mean", sides = "lower", alpha = 0.05, control = ctl))
  expect_identical(dg(r$value)$limitation_flag, rep("gpq_one_sided_group_mean", 2))
})

test_that("no limitation warning: two-sided GPQ TI, and the observation target at 0.95", {
  skip_if_not_installed("lme4")
  fm <- rs_fit()
  for (ov in c("observation", "group_mean")) for (lv in c(0.90, 0.95)) {
    r <- collect(ti_lmm(fm, nd, over = ov, sides = "two.sided", level = lv, control = ctl))
    expect_length(r$limitation, 0L)
    expect_identical(dg(r$value)$limitation_flag, rep(NA_character_, 2), info = paste(ov, lv))
  }
  for (sd in c("lower", "upper")) for (lv in c(0.95, 0.99)) {
    r <- collect(ti_lmm(fm, nd, over = "observation", sides = sd, level = lv, control = ctl))
    expect_length(r$limitation, 0L)
    expect_identical(dg(r$value)$limitation_flag, rep(NA_character_, 2), info = paste(sd, lv))
  }
  ## the closed-form routes carry no flag and do not warn
  fi <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  r <- collect(ti_lmm(fi, nd, sides = "lower", level = 0.90))
  expect_length(r$limitation, 0L)
  expect_null(dg(r$value)$limitation_flag)
})

test_that("vectorised call: the flag is per row and each case warns once", {
  skip_if_not_installed("lme4")
  fm <- rs_fit()
  r <- collect(ti_lmm(fm, nd, over = "group_mean", sides = c("two.sided", "lower"),
                      level = c(0.90, 0.95), control = ctl))
  v <- r$value
  ## rows: t0 fastest, then sides, then level
  expect_identical(v$sides, rep(rep(c("two.sided", "lower"), each = 2), 2))
  expect_identical(dg(v)$limitation_flag,
                   c(NA, NA, "gpq_one_sided_lowconf", "gpq_one_sided_lowconf",
                     NA, NA, "gpq_one_sided_group_mean", "gpq_one_sided_group_mean"))
  expect_identical(vapply(r$limitation, `[[`, "", "limitation"),
                   c("gpq_one_sided_lowconf", "gpq_one_sided_group_mean"))
  ## the observation target: only the low-confidence rows
  o <- collect(ti_lmm(fm, nd, sides = c("upper", "two.sided"), level = c(0.90, 0.95), control = ctl))
  expect_identical(dg(o$value)$limitation_flag,
                   c("gpq_one_sided_lowconf", "gpq_one_sided_lowconf", NA, NA, NA, NA, NA, NA))
  expect_length(o$limitation, 1L)
})

test_that("random-intercept GPQ route: the limitation warning applies, and the unvalidated warning stays", {
  skip_if_not_installed("lme4")
  fi <- lmer_q(Reaction ~ Days + (1 | Subject), lme4::sleepstudy)
  r <- collect(ti_lmm(fi, data.frame(Days = 0), method = "gpq", sides = "lower",
                      level = 0.90, control = ctl))
  expect_length(r$limitation, 1L)
  expect_identical(r$limitation[[1]]$limitation, "gpq_one_sided_lowconf")
  expect_length(r$unvalidated, 1L)
  expect_identical(dg(r$value)$limitation_flag, "gpq_one_sided_lowconf")
  g <- collect(ti_lmm(fi, data.frame(Days = 0), method = "gpq", over = "group_mean",
                      sides = "upper", control = ctl))
  expect_identical(dg(g$value)$limitation_flag, "gpq_one_sided_group_mean")
  expect_length(g$limitation, 1L)
  t <- collect(ti_lmm(fi, data.frame(Days = 0), method = "gpq", control = ctl))
  expect_length(t$limitation, 0L)
  expect_length(t$unvalidated, 1L)
  expect_identical(dg(t$value)$limitation_flag, NA_character_)
})

test_that("ti_gpq_raw(): the same rule, c_resid = 0 being the batch-mean target", {
  d <- lme4::sleepstudy
  raw <- function(...) collect(ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5, control = ctl, ...))
  expect_identical(dg(raw(sides = "lower", level = 0.90)$value)$limitation_flag, "gpq_one_sided_lowconf")
  expect_identical(dg(raw(sides = "upper", c_resid = 0)$value)$limitation_flag, "gpq_one_sided_group_mean")
  r <- raw(sides = "lower")
  expect_identical(dg(r$value)$limitation_flag, NA_character_)
  expect_length(r$limitation, 0L)
  expect_length(raw(level = 0.90)$limitation, 0L)                    # two-sided
  expect_length(raw(sides = "upper", c_resid = 0)$limitation, 1L)
})

test_that(".gpq_limitation_flag(): the rule", {
  f <- .gpq_limitation_flag(conf = c(0.9, 0.95, 1 - 0.05, 0.99, 0.9, 0.95),
                            sides = c("lower", "lower", "upper", "upper", "two.sided", "two.sided"),
                            c_resid = c(1, 1, 0, 0, 0, 0))
  expect_identical(f, c("gpq_one_sided_lowconf", NA, "gpq_one_sided_group_mean",
                        "gpq_one_sided_group_mean", NA, NA))
})
