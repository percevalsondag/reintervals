# oracle: chow-shao1991-t23 (fixtures/oracles.csv)
#
# Chow S-C, Shao J (1991). Estimating drug shelf-life with random batches.
# Biometrics 47(3):1071-1079. Table 2 (p. 1078): potency (% of claim) of 24
# batches at 0, 12, 24, 36 months (fixtures/chow_shao1991/data.csv).
#
# Reference: on balanced data (every batch at the same times) with an interior
# REML fit of y ~ month + (month | batch), the GLS mean line at xi is the mean
# of the 24 per-batch OLS lines, its variance is their sample variance / 24, and
# the pivot is exactly t on K - 1 = 23 df. The exact limits are therefore
# closed form, computed in base R below. (The paper's own WLS bound uses an
# unstructured D and suggests t on k - 2 df; the gate uses the exact K - 1
# result for the REML model.)
#
# H2 gate (pre-registered 2026-11-22, before the package's results were
# computed): an arm misses if, at any xi and side, its limit is narrower than
# the exact t_23 limit (closer to the estimate) by more than 0.01 x the exact
# half-width; a wider limit is recorded and passes.

cs_data <- function() oracle_fixture("chow-shao1991-t23")
cs_xi <- c(0, 12, 24, 36, 48, 60)
cs_months <- c(0, 12, 24, 36)

## per-batch OLS lines -> exact t_23 reference at each xi (base R only)
cs_exact <- function(d, xi = cs_xi) {
  B <- t(sapply(split(d, d$batch), function(s) stats::coef(stats::lm(potency ~ month, s))))
  K <- nrow(B)
  do.call(rbind, lapply(xi, function(x) {
    fv <- B[, 1] + B[, 2] * x
    m <- mean(fv); se <- stats::sd(fv) / sqrt(K)
    q2 <- stats::qt(0.975, K - 1); q1 <- stats::qt(0.95, K - 1)
    data.frame(xi = x, sides = c("two.sided", "lower", "upper"), estimate = m, se = se,
               df = K - 1,
               lower = c(m - q2 * se, m - q1 * se, -Inf),
               upper = c(m + q2 * se, Inf, m + q1 * se),
               hw = c(q2, q1, q1) * se, stringsAsFactors = FALSE)
  }))
}

## The H2 gate: one row per finite limit (two-sided rows give a lower and an
## upper bound). narrower_by = how much closer to the estimate the package's
## limit is than the exact one, as a fraction of the exact half-width (negative
## = wider). miss when narrower_by > 0.01.
gate_chow_shao <- function(result, exact, tol = 0.01) {
  key <- function(x, s) paste(x, s)
  ex <- exact[match(key(result$xi, result$sides), key(exact$xi, exact$sides)), ]
  rows <- lapply(seq_len(nrow(result)), function(i) {
    b <- c(if (is.finite(ex$lower[i])) "lower", if (is.finite(ex$upper[i])) "upper")
    do.call(rbind, lapply(b, function(bd) {
      lim <- result[[bd]][i]; elim <- ex[[bd]][i]
      nb <- (if (bd == "lower") lim - elim else elim - lim) / ex$hw[i]
      data.frame(xi = result$xi[i], side = result$sides[i], bound = bd,
                 df = result$df[i], limit = lim, exact_limit = elim,
                 narrower_by = nb, verdict = if (nb > tol) "miss" else "pass",
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows)
}

## the package's ci_lmm at every xi x side for one slope_ci_df arm
cs_package <- function(fm, arm) {
  ctl <- lmm_interval_control(slope_ci_df = arm)
  do.call(rbind, lapply(c("two.sided", "lower", "upper"), function(sd) {
    r <- suppressWarnings(ci_lmm(fm, newdata = data.frame(month = cs_xi), level = 0.95,
                                 sides = sd, control = ctl))
    data.frame(xi = r$month, sides = sd, estimate = r$estimate,
               se = attr(r, "diagnostics")$se, df = attr(r, "diagnostics")$df,
               lower = r$lower, upper = r$upper, stringsAsFactors = FALSE)
  }))
}

## The identity is exact at the REML optimum. lme4's default nloptwrap
## tolerances stop short of it by enough to move the SE by up to 1e-5 relative
## (at xi = 60: 1.0e-5; residual variance 0.9145791 vs the pooled 0.9145833),
## above the 1e-6 tolerance; the same model with tight optimizer tolerances
## reaches it (SE within 5e-8). This converges the fit; it does not change the
## model or any package code.
cs_fit <- function(d) suppressWarnings(suppressMessages(
  lme4::lmer(potency ~ month + (month | batch), data = d, REML = TRUE,
             control = lme4::lmerControl(optimizer = "nloptwrap", optCtrl = list(
               xtol_abs = 1e-14, ftol_abs = 1e-15, xtol_rel = 1e-14,
               ftol_rel = 1e-15, maxeval = 1e5)))))

test_that("chow-shao1991-t23 (a): the transcription reproduces the paper's printed WLS and batch-21 OLS estimates", {
  d <- cs_data()
  expect_identical(dim(d), c(96L, 3L))
  expect_equal(sort(unique(d$batch)), 1:24)
  expect_equal(sort(unique(d$month)), cs_months)
  expect_true(all(table(d$batch) == 4L))
  d <- d[order(d$batch, d$month), ]
  Z <- do.call(rbind, split(d$potency, d$batch))      # 24 x 4 batch vectors
  k <- nrow(Z)
  X <- unname(cbind(1, cs_months))
  Di <- solve(stats::cov(Z))                          # divisor k - 1
  A <- solve(t(X) %*% Di %*% X)
  th <- unname(drop(A %*% t(X) %*% Di %*% colMeans(Z)))
  se <- unname(sqrt(diag(A / k)))
  expect_equal(round(th, 3), c(103.968, -0.199), tolerance = 0)
  expect_equal(round(se, 3), c(0.289, 0.009), tolerance = 0)
  f21 <- stats::lm(potency ~ month, d[d$batch == 21, ])
  expect_equal(round(unname(stats::coef(f21)), 3), c(102.500, -0.333), tolerance = 0)
  expect_equal(round(unname(sqrt(diag(stats::vcov(f21)))), 3), c(1.775, 0.079), tolerance = 0)
})

test_that("chow-shao1991-t23 (b): the exact t_23 identity's precondition holds (PD moment Sigma, non-singular REML fit)", {
  skip_if_not_installed("lme4")
  d <- cs_data()
  s <- split(d, d$batch)
  B <- t(sapply(s, function(x) stats::coef(stats::lm(potency ~ month, x))))
  s2 <- mean(sapply(s, function(x) summary(stats::lm(potency ~ month, x))$sigma^2))
  X <- cbind(1, cs_months)
  ev <- eigen(stats::cov(B) - s2 * solve(crossprod(X)), symmetric = TRUE)$values
  expect_true(all(ev > 0))
  expect_equal(ev, c(1.442, 0.00196), tolerance = 1e-3)
  expect_false(lme4::isSingular(cs_fit(d)))
  ex <- cs_exact(d)
  expect_true(all(ex$df == 23))
})

test_that("chow-shao1991-t23 (c, d): ci_lmm slope_ci_df arms vs the exact t_23 limits; H2 gate", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  d <- cs_data()
  fm <- cs_fit(d)
  ex <- cs_exact(d)
  for (arm in c("kenward-roger", "k-minus-1", "hessian-satterthwaite")) {
    r <- cs_package(fm, arm)
    e <- ex[match(paste(r$xi, r$sides), paste(ex$xi, ex$sides)), ]
    expect_equal(r$estimate, e$estimate, tolerance = 1e-6, info = arm)
    expect_equal(r$se, e$se, tolerance = 1e-6, info = arm)
    g <- gate_chow_shao(r, ex)
    if (arm != "hessian-satterthwaite") {
      expect_true(all(abs(r$df - 23) <= 1e-6), info = arm)
      expect_equal(r$lower, e$lower, tolerance = 1e-6, info = arm)
      expect_equal(r$upper, e$upper, tolerance = 1e-6, info = arm)
      expect_true(all(g$verdict == "pass"), info = arm)
    } else {
      ## the default arm: recorded, not asserted (H2 gate result reported)
      cat("\nchow-shao1991-t23 H2 gate, slope_ci_df = \"hessian-satterthwaite\":\n")
      tab <- transform(g, df = sprintf("%.7f", df), limit = sprintf("%.6f", limit),
                       exact_limit = sprintf("%.6f", exact_limit),
                       narrower_by = sprintf("%.2e", narrower_by))
      print(tab, row.names = FALSE)
      cat("verdict:", if (any(g$verdict == "miss")) "MISS" else "PASS",
          sprintf("(max narrower_by %.2e)\n", max(g$narrower_by)))
    }
  }
})

test_that("gate_chow_shao() flags a limit narrower than exact by more than 1% of the half-width", {
  ex <- data.frame(xi = 0, sides = c("two.sided", "lower"), estimate = 100, se = 1, df = 23,
                   lower = c(98, 98.5), upper = c(102, Inf), hw = c(2, 1.5))
  ok <- data.frame(xi = 0, sides = c("two.sided", "lower"), df = 23,
                   lower = c(97.5, 98.5), upper = c(102.01, Inf))           # wider / equal
  g <- gate_chow_shao(ok, ex)
  expect_identical(g$bound, c("lower", "upper", "lower"))
  expect_identical(g$verdict, rep("pass", 3))
  expect_equal(g$narrower_by, c(-0.25, -0.005, 0))
  bad <- ok; bad$upper[1] <- 101.97                                       # 1.5% narrower
  g <- gate_chow_shao(bad, ex)
  expect_identical(g$verdict, c("pass", "miss", "pass"))
  expect_equal(g$narrower_by[2], 0.015)
})
