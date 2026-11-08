# A2: GPQ hardening (plumbing only). Seed contract, common random numbers (CRN),
# Monte-Carlo standard error by MC sections, the M >= 2000 check, and the
# per-lot time-point guard. None of these is an oracle: small M (2000-8000) is
# fine because they test the Monte-Carlo plumbing, not coverage.

lc_a2 <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
fs_a2 <- function() {
  suppressWarnings(suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy,
               control = lc_a2)))
}
ctl2k <- lmm_interval_control(M = 2000)        # not an oracle: small M

test_that("A2 (a): set.seed() before the call controls a seed = NULL GPQ result", {
  fs <- fs_a2(); nd <- data.frame(Days = c(0, 5, 9))
  set.seed(1); a <- ti_lmm(fs, newdata = nd, control = ctl2k)
  set.seed(1); b <- ti_lmm(fs, newdata = nd, control = ctl2k)
  expect_identical(a, b)
  # the recorded seed is the one integer drawn from the ambient stream ...
  set.seed(1); s <- sample.int(.Machine$integer.max, 1L)
  expect_identical(attr(a, "diagnostics")$gpq_seed, s)
  # ... and the call advances the caller's stream by exactly that one draw
  set.seed(1); invisible(ti_lmm(fs, newdata = nd, control = ctl2k)); after_call <- runif(3)
  set.seed(1); invisible(sample.int(.Machine$integer.max, 1L)); after_one <- runif(3)
  expect_identical(after_call, after_one)
  expect_identical(attr(a, "diagnostics")$rng_kind, RNGkind())
})

test_that("A2 (b): seed = gpq_seed from the diagnostics reproduces the limits exactly", {
  fs <- fs_a2(); nd <- data.frame(Days = c(0, 5, 9))
  set.seed(42); a <- ti_lmm(fs, newdata = nd, sides = "lower", control = ctl2k)
  d <- attr(a, "diagnostics")
  b <- ti_lmm(fs, newdata = nd, sides = "lower",
              control = lmm_interval_control(seed = d$gpq_seed, M = 2000))
  expect_identical(a$lower, b$lower); expect_identical(a$upper, b$upper)
  expect_identical(attr(b, "diagnostics")$gpq_seed, d$gpq_seed)
  # same contract under a non-default RNG kind: rng_kind records it
  old <- RNGkind("L'Ecuyer-CMRG"); withr::defer(RNGkind(old[1], old[2], old[3]))
  set.seed(7); x <- ti_lmm(fs, newdata = nd, control = ctl2k)
  dx <- attr(x, "diagnostics")
  expect_identical(dx$rng_kind[1], "L'Ecuyer-CMRG")
  y <- ti_lmm(fs, newdata = nd, control = lmm_interval_control(seed = dx$gpq_seed, M = 2000))
  expect_identical(x$lower, y$lower); expect_identical(x$upper, y$upper)
})

test_that("A2 (c): CRN -- the t0 = 5 limits do not depend on the other requested time points", {
  fs <- fs_a2()
  ctl <- lmm_interval_control(seed = 20261107, M = 2000)
  nds <- list(data.frame(Days = 5), data.frame(Days = c(0, 5, 9)),
              data.frame(Days = c(9, 5, 0)))
  for (ov in c("observation", "group_mean")) for (sd in c("two.sided", "lower", "upper")) {
    at5 <- lapply(nds, function(nd) {
      r <- ti_lmm(fs, newdata = nd, P = 0.99, over = ov, sides = sd, control = ctl)
      i <- which(r$Days == 5); dg <- attr(r, "diagnostics")
      c(r$lower[i], r$upper[i], r$estimate[i], dg$mc_se_lower[i], dg$mc_se_upper[i])
    })
    info <- paste(ov, sd)
    expect_identical(at5[[2]], at5[[1]], info = info)
    expect_identical(at5[[3]], at5[[1]], info = info)
  }
  # the raw engine too
  d <- lme4::sleepstudy
  r1 <- ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5, control = ctl)
  r3 <- ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = c(9, 5, 0), control = ctl)
  expect_identical(c(r3$lower[2], r3$upper[2]), c(r1$lower, r1$upper))
})

test_that("A2 (d): MC SE is finite/positive on finite sides, NA on infinite sides, ~1/sqrt(M)", {
  fs <- fs_a2(); nd <- data.frame(Days = c(0, 5, 9))
  ctl <- lmm_interval_control(seed = 1, M = 2000)
  for (sd in c("two.sided", "lower", "upper")) {
    dg <- attr(ti_lmm(fs, newdata = nd, sides = sd, control = ctl), "diagnostics")
    expect_identical(dg$mc_sections, 20L)
    if (sd == "upper") expect_true(all(is.na(dg$mc_se_lower)), info = sd)
    else expect_true(all(is.finite(dg$mc_se_lower) & dg$mc_se_lower > 0), info = sd)
    if (sd == "lower") expect_true(all(is.na(dg$mc_se_upper)), info = sd)
    else expect_true(all(is.finite(dg$mc_se_upper) & dg$mc_se_upper > 0), info = sd)
  }
  # 4x the draws -> about half the MC SE (two-sided, observation, Days = 5).
  # The MC SE is itself noisy (an sd over 20 MC sections, ~16% relative error),
  # so a single-seed ratio falls outside [1.3, 3.0] about 10-20% of the time
  # (checked over 40 seeds in A2). The bounds are therefore applied to the
  # geometric mean of the ratio over 5 fixed seeds, not to one seed.
  se_at <- function(M, seed) {
    dg <- attr(ti_lmm(fs, newdata = data.frame(Days = 5),
                      control = lmm_interval_control(seed = seed, M = M)), "diagnostics")
    c(dg$mc_se_lower, dg$mc_se_upper)
  }
  lr <- vapply(1:5, function(s) log(se_at(2000, s) / se_at(8000, s)), numeric(2))
  ratio <- exp(rowMeans(lr))
  expect_true(all(ratio >= 1.3 & ratio <= 3.0), info = paste(round(ratio, 3), collapse = ", "))
})

test_that("A2 (e): M below 2000 is an input error; M not divisible by 20 works", {
  expect_error(lmm_interval_control(M = 1999), class = "reintervals_error_input")
  expect_error(lmm_interval_control(M = 2000.5), class = "reintervals_error_input")
  expect_identical(lmm_interval_control(M = 2001)$M, 2001L)
  d <- lme4::sleepstudy
  g <- ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5,
                  control = lmm_interval_control(seed = 1, M = 2001))
  expect_true(is.finite(g$lower) && is.finite(g$upper))
  expect_true(all(is.finite(unlist(attr(g, "diagnostics")[c("mc_se_lower", "mc_se_upper")]))))
})

test_that("A2 (f): a lot with < 2 distinct time points aborts (raw engine and ti_lmm)", {
  d <- lme4::sleepstudy
  one_time <- c("308", "309", "310", "330", "331", "332", "333")   # 7 lots
  keep <- !(d$Subject %in% one_time) | d$Days == 0
  dd <- d[keep, ]
  err <- tryCatch(ti_gpq_raw(dd$Reaction, dd$Days, dd$Subject, t0 = 5, control = ctl2k),
                  error = identity)
  expect_s3_class(err, "reintervals_error_time_points")
  expect_s3_class(err, "reintervals_error_data")
  expect_match(conditionMessage(err), "at least two distinct time points per lot")
  expect_match(conditionMessage(err), "308, 309, 310, 330, 331 and 2 more", fixed = TRUE)
  expect_identical(err$lots, one_time)
  expect_match(conditionMessage(err), "y ~ time + (1 | lot)", fixed = TRUE)
  # repeated measurements at a single time point are still one time point
  d2 <- rbind(d, transform(d[d$Subject == "308" & d$Days == 0, ], Reaction = Reaction + 1))
  d2 <- d2[!(d2$Subject == "308") | d2$Days == 0, ]
  expect_error(ti_gpq_raw(d2$Reaction, d2$Days, d2$Subject, t0 = 5, control = ctl2k),
               class = "reintervals_error_time_points")
  # reachable through ti_lmm(): lme4 fits these data (it only checks the totals)
  fit <- suppressWarnings(suppressMessages(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = dd, control = lc_a2)))
  expect_error(ti_lmm(fit, newdata = data.frame(Days = 5), control = ctl2k),
               class = "reintervals_error_time_points")
})

test_that("A2 (g): ti_gpq_raw() returns the same diagnostics as the ti_lmm() GPQ route", {
  fs <- fs_a2(); d <- lme4::sleepstudy
  ctl <- lmm_interval_control(seed = 3, M = 2000)
  for (sd in c("two.sided", "lower")) {
    a <- ti_lmm(fs, newdata = data.frame(Days = c(2, 7)), method = "gpq", sides = sd,
                control = ctl)
    b <- ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = c(2, 7), sides = sd, control = ctl)
    expect_identical(names(attr(b, "diagnostics")), names(attr(a, "diagnostics")))
    expect_equal(attr(b, "diagnostics"), attr(a, "diagnostics"), info = sd)
    expect_identical(b$sides, rep(sd, 2))       # raw engine labels one-sided output
    expect_equal(c(b$lower, b$upper), c(a$lower, a$upper))
  }
})
