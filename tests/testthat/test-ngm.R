# Random-intercept new-group-mean constructions. Each is checked against an
# independent reference written here from its defining formula (no package code),
# on datasets from the generators of the validation study (paper repo,
# scripts/defects/e5_fun.R, e5e_fun.R, e1.R):
#   2a Fa (fixed slope)   2b slope-free Fa / Ca (intercept-only one-way)
#   2c Gm (nested)        2d Gv (crossed with interaction)
hw_of <- function(x) (x$upper - x$lower) / 2
dg <- function(x) attr(x, "diagnostics")
lmer_q <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d)))

## ---- references -----------------------------------------------------------
ref_fa <- function(y, time, lot, t0, level = 0.95) {      # time = NULL: slope-free
  lot <- factor(lot); I <- nlevels(lot); N <- length(y); li <- as.integer(lot)
  yb <- tapply(y, lot, mean); ys <- mean(yb)
  if (is.null(time)) {
    S2L <- sum((yb - ys)^2) / (I - 1)
    W <- S2L * (1 + 1 / I); df <- I - 1; est <- ys
  } else {
    tb <- tapply(time, lot, mean); ts <- mean(tb)
    tc <- time - tb[li]; yc <- y - yb[li]
    Stt <- sum(tc^2); b <- sum(tc * yc) / Stt; r <- N - I - 1
    S2L <- sum((yb - ys + b * (ts - tb))^2) / (I - 1)
    S2E <- (sum(yc^2) - sum(tc * yc)^2 / Stt) / r
    a1 <- S2L * (1 + 1 / I); a2 <- S2E * (t0 - ts)^2 / Stt
    W <- a1 + a2; df <- W^2 / (a1^2 / (I - 1) + a2^2 / r); est <- ys + b * (t0 - ts)
  }
  list(est = est, hw = qt(1 - (1 - level) / 2, df) * sqrt(W), df = df, W = W)
}
ref_ca <- function(y, g, level = 0.95) {                  # E1's Ca, balanced one-way
  g <- factor(g); A <- nlevels(g); n <- length(y) / A
  MSA <- n * sum((tapply(y, g, mean) - mean(y))^2) / (A - 1)
  qt(1 - (1 - level) / 2, A - 1) * sqrt(MSA * (1 + 1 / A) / n)
}
blk <- function(f) { z <- model.matrix(~f - 1); z %*% (t(z) / colSums(z)) }
## E (component coefficients of each term's expected mean square), df, and k
## (E'k = 1) from the ANOVA projections: sequential for nested, sum-to-zero
## Type III for crossed.
ref_ems <- function(fa, fb, fab = NULL) {
  n <- length(fa)
  if (is.null(fab)) {
    Q <- list(A = blk(fa) - matrix(1 / n, n, n), B = blk(fb) - blk(fa))
    df <- c(A = nlevels(fa) - 1, B = nlevels(fb) - nlevels(fa))
    Z <- list(A = model.matrix(~fa - 1), B = model.matrix(~fb - 1))
  } else {
    x <- model.matrix(~fa * fb, contrasts.arg = list(fa = contr.sum, fb = contr.sum))
    asg <- attr(x, "assign"); wi <- solve(crossprod(x))
    Q <- lapply(list(A = which(asg == 1), B = which(asg == 2), AB = which(asg == 3)),
                function(f) { p <- wi[f, , drop = FALSE] %*% t(x); t(p) %*% solve(wi[f, f, drop = FALSE]) %*% p })
    df <- c(A = nlevels(fa) - 1, B = nlevels(fb) - 1, AB = (nlevels(fa) - 1) * (nlevels(fb) - 1))
    Z <- list(A = model.matrix(~fa - 1), B = model.matrix(~fb - 1), AB = model.matrix(~fab - 1))
  }
  Q$Residual <- diag(n) - (if (is.null(fab)) blk(fb) else x %*% wi %*% t(x))
  df <- c(df, Residual = n - sum(df) - 1); Z$Residual <- diag(n)
  E <- t(sapply(names(Q), function(j) sapply(names(Z), function(g) sum((Q[[j]] %*% Z[[g]]) * Z[[g]]) / df[[j]])))
  list(E = E, df = df, k = solve(t(E), rep(1, nrow(E))))
}
ref_g <- function(fm, em, s2, crossed) {
  ms <- as.vector(em$E %*% s2); j <- seq_len(length(ms) - 1L)
  VF <- vcov(fm)[1, 1]; k <- em$k; W <- VF + sum(k[j] * ms[j])
  if (!crossed) df <- em$df[[1]] else {
    co <- k[j]; pA <- k[1] * ms[1] / (k[1] * ms[1] + k[2] * ms[2])
    co[1] <- co[1] + pA * VF / ms[1]; co[2] <- co[2] + (1 - pA) * VF / ms[2]
    df <- sum(co * ms[j])^2 / sum((co * ms[j])^2 / em$df[j])
  }
  list(hw = qt(0.975, df) * sqrt(W), df = df, W = W)
}

## ---- generators (validation-study designs) -----------------------------------
gen_fixed_slope <- function(I, frac, ICC, seed) {           # e5_fun.R make_design
  set.seed(seed); ns <- round(frac * I); nm <- I - ns
  d <- rbind(data.frame(Lot = rep(sprintf("L%02d", seq_len(nm)), each = 5), Month = rep(c(0, 3, 6, 9, 12), nm)),
             data.frame(Lot = sprintf("S%02d", seq_len(ns)), Month = rep(0, ns)))
  d$Lot <- factor(d$Lot); s2a <- ICC / (1 - ICC)
  d$Y <- 100 - 2 * d$Month + rnorm(nlevels(d$Lot), 0, sqrt(s2a))[d$Lot] + rnorm(nrow(d))
  d
}
gen_oneway <- function(A, n, ICC, seed) {                     # e1.R (n = 0: unbalanced)
  set.seed(seed)
  ni <- if (n > 0) rep(n, A) else repeat { v <- sample(1:6, A, replace = TRUE); if (sum(v) > A && length(unique(v)) > 1) break }
  if (n == 0) ni <- v
  g <- factor(rep(seq_len(A), ni))
  data.frame(g = g, y = 100 + rnorm(A, 0, sqrt(ICC / (1 - ICC)))[g] + rnorm(length(g)))
}
gen_cells <- function(A, B, n, crossed) {                     # e5e_fun.R gen_cells_*
  if (n > 0) return(matrix(n, A, B))
  repeat { m <- matrix(sample(1:4, A * B, replace = TRUE), A, B)
           if (length(unique(as.vector(m))) > 1 && (crossed || sum(m) > A * B)) return(m) }
}
gen_two <- function(A, B, n, v, crossed, seed) {
  set.seed(seed); v <- v / sum(v); cells <- gen_cells(A, B, n, crossed)
  ci <- rep(seq_len(A * B), as.vector(cells))
  fa <- factor((ci - 1) %% A + 1)
  fb <- if (crossed) factor((ci - 1) %/% A + 1) else factor(ci)
  y <- 100 + rnorm(A, 0, sqrt(v[1]))[fa] +
    (if (crossed) rnorm(B, 0, sqrt(v[2]))[fb] + rnorm(A * B, 0, sqrt(v[3]))[ci] else rnorm(A * B, 0, sqrt(v[2]))[ci]) +
    rnorm(length(ci), 0, sqrt(v[length(v)]))
  data.frame(y = y, fa = fa, fb = fb, fab = factor(ci))
}

## ---- equality with the references ------------------------------------------
test_that("2a Fa: fixed slope (balanced and single-observation lots) equals the reference", {
  skip_if_not_installed("lme4")
  grid <- expand.grid(I = c(6, 11, 20), frac = c(0, 0.3, 0.6))[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 5), ]
  for (i in seq_len(nrow(grid))) {
    d <- gen_fixed_slope(grid$I[i], grid$frac[i], c(0.1, 0.5, 0.9)[i %% 3 + 1], seed = 700 + i)
    x <- new_group_mean_lmm(lmer_q(Y ~ Month + (1 | Lot), d), newdata = data.frame(Month = c(6, 12)))
    for (t in 1:2) {
      r <- ref_fa(d$Y, d$Month, d$Lot, c(6, 12)[t])
      expect_equal(hw_of(x)[t], r$hw, tolerance = 1e-8, info = paste(i, t))
      expect_equal(dg(x)$df[t], r$df, tolerance = 1e-8, info = paste(i, t))
      expect_equal(x$estimate[t], r$est, tolerance = 1e-10)
    }
    expect_identical(unique(x$method), "anova-mls")
    expect_identical(unique(dg(x)$df_rule), "anova-satterthwaite")
  }
})

test_that("2b: intercept-only one-way equals the slope-free reference, and Ca exactly when balanced", {
  skip_if_not_installed("lme4")
  cfg <- expand.grid(A = c(3, 5, 10), n = c(2, 5, 0))[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 3), ]
  for (i in seq_len(nrow(cfg))) {
    d <- gen_oneway(cfg$A[i], cfg$n[i], 0.5, seed = 900 + i)
    x <- new_group_mean_lmm(lmer_q(y ~ 1 + (1 | g), d))
    r <- ref_fa(d$y, NULL, d$g)
    expect_equal(hw_of(x), r$hw, tolerance = 1e-8, info = i)
    expect_equal(dg(x)$df, cfg$A[i] - 1)
    expect_identical(x$method, "anova-mls")
    expect_identical(dg(x)$df_rule, "anova-single")
    if (cfg$n[i] > 0) {
      expect_equal(hw_of(x), ref_ca(d$y, d$g), tolerance = 1e-12, info = i)   # == Ca (E1)
      expect_equal(x$estimate, mean(d$y), tolerance = 1e-12)
    }
  }
})

test_that("2c Gm: nested (balanced and unbalanced) equals the reference", {
  skip_if_not_installed("lme4")
  cfg <- expand.grid(A = c(3, 5, 10), B = c(2, 3), n = c(2, 0))[c(1, 2, 3, 4, 7, 8, 9, 10, 11, 12), ]
  for (i in seq_len(nrow(cfg))) {
    d <- gen_two(cfg$A[i], cfg$B[i], cfg$n[i], v = list(c(1, 1, 8), c(4, 2, 4), c(8, 1, 1))[[i %% 3 + 1]],
                 crossed = FALSE, seed = 1100 + i)
    fm <- lmer_q(y ~ 1 + (1 | fa) + (1 | fb), d)
    x <- new_group_mean_lmm(fm)
    vc <- as.data.frame(lme4::VarCorr(fm))
    s2 <- vc$vcov[match(c("fa", "fb", "Residual"), vc$grp)]
    r <- ref_g(fm, ref_ems(d$fa, d$fb), s2, crossed = FALSE)
    expect_equal(hw_of(x), r$hw, tolerance = 1e-8, info = i)
    expect_equal(dg(x)$df, r$df, tolerance = 1e-8, info = i)
    expect_equal(dg(x)$df, cfg$A[i] - 1)
    expect_identical(x$method, "ems-mls")
    expect_identical(dg(x)$df_rule, "top-factor")
  }
})

test_that("2d Gv: crossed with interaction (balanced and unbalanced) equals the reference", {
  skip_if_not_installed("lme4")
  cfg <- expand.grid(A = c(2, 3, 5), B = c(2, 3, 5), n = c(2, 0))[c(2, 4, 5, 6, 8, 11, 13, 14, 15, 18), ]
  for (i in seq_len(nrow(cfg))) {
    d <- gen_two(cfg$A[i], cfg$B[i], cfg$n[i], v = list(c(1, 1, 1, 7), c(3, 3, 2, 2), c(4, 1, 1, 4))[[i %% 3 + 1]],
                 crossed = TRUE, seed = 1300 + i)
    fm <- lmer_q(y ~ 1 + (1 | fa) + (1 | fb) + (1 | fab), d)
    x <- new_group_mean_lmm(fm)
    vc <- as.data.frame(lme4::VarCorr(fm))
    s2 <- vc$vcov[match(c("fa", "fb", "fab", "Residual"), vc$grp)]
    r <- ref_g(fm, ref_ems(d$fa, d$fb, d$fab), s2, crossed = TRUE)
    expect_equal(hw_of(x), r$hw, tolerance = 1e-8, info = i)
    expect_equal(dg(x)$df, r$df, tolerance = 1e-8, info = i)
    expect_identical(x$method, "ems-mls")
    expect_identical(dg(x)$df_rule, "satterthwaite-vf")
  }
})

## ---- routing ----------------------------------------------------------------
test_that("2e: one-way with other fixed effects uses Gm with the one-way k (df = A - 1)", {
  skip_if_not_installed("lme4")
  d <- gen_oneway(6, 4, 0.5, seed = 1); d$grp <- factor(rep(c("a", "b"), length.out = nrow(d)))
  fm <- lmer_q(y ~ grp + (1 | g), d)
  x <- new_group_mean_lmm(fm, newdata = data.frame(grp = factor("a", levels = c("a", "b"))))
  vc <- as.data.frame(lme4::VarCorr(fm)); s2a <- vc$vcov[vc$grp == "g"]; s2e <- vc$vcov[vc$grp == "Residual"]
  W <- vcov(fm)[1, 1] + (s2e + 4 * s2a) / 4                 # V_F + EMS_A / n0, n0 = 4
  expect_equal(hw_of(x), qt(0.975, 5) * sqrt(W), tolerance = 1e-8)
  expect_equal(dg(x)$df, 5)
  expect_identical(x$method, "ems-mls")
  expect_identical(dg(x)$df_rule, "top-factor")
})

test_that("2f: designs without a closed-form decomposition error", {
  skip_if_not_installed("lme4")
  set.seed(1); g <- expand.grid(a = factor(1:3), b = factor(1:4), rep = 1:3)
  g$y <- rnorm(3, 0, 2)[g$a] + rnorm(4, 0, 2)[g$b] + rnorm(36)
  expect_error(new_group_mean_lmm(lmer_q(y ~ 1 + (1 | a) + (1 | b), g)),
               "no closed-form variance decomposition")
})

test_that("method = is limited to the design's new-group-mean engine", {
  skip_if_not_installed("lme4")
  fm <- lmer_q(Y ~ Month + (1 | Lot), mbl_appendix_a)
  nd <- data.frame(Month = 12)
  expect_identical(new_group_mean_lmm(fm, nd, method = "anova-mls")$method, "anova-mls")
  expect_error(new_group_mean_lmm(fm, nd, method = "ems-mls"), "not available")
  d <- gen_two(3, 2, 2, c(1, 1, 8), crossed = FALSE, seed = 1)
  fn <- lmer_q(y ~ 1 + (1 | fa) + (1 | fb), d)
  expect_identical(new_group_mean_lmm(fn, method = "ems-mls")$method, "ems-mls")
  expect_error(new_group_mean_lmm(fn, method = "anova-mls"), "not available")
})

test_that("MBL (single-observation fixed slope) at Month 12: anova-mls, no anti-conservative warning", {
  skip_if_not_installed("lme4")
  fm <- lmer_q(Y ~ Month + (1 | Lot), mbl_appendix_a)
  expect_warning(x <- new_group_mean_lmm(fm, newdata = data.frame(Month = 12)), NA)
  expect_identical(x$method, "anova-mls")
  expect_identical(dg(x)$df_rule, "anova-satterthwaite")
  r <- ref_fa(mbl_appendix_a$Y, mbl_appendix_a$Month, mbl_appendix_a$Lot, 12)
  expect_equal(c(x$lower, x$upper), c(r$est - r$hw, r$est + r$hw), tolerance = 1e-10)
  expect_warning(ci_lmm(fm, newdata = data.frame(Month = 12)), "anti-conservative")   # CI still warns
})
