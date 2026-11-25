# GPQ tolerance interval on random-intercept models (ti_lmm(method = "gpq"),
# step 13b): the random-intercept sub-class of Oliva-Aviles C, Hauser P (2025),
# Technometrics 67(2):193-202, doi:10.1080/00401706.2024.2407324 (their
# sections 4.1 and 4.3), checked against Liao CT, Lin TY, Iyer HK (2005),
# Technometrics 47(3):323-335, doi:10.1198/004017005000000102.
#
# oracle: liao2005-sulfur      (fixtures/liao2005_sulfur; Example 2)
# oracle: liao2005-ex3-nist    (fixtures/nist_ex4; Example 3: canonical form
#                               is an oracle, the printed GPQ limits are a
#                               comparison only)
# oracle: liao2005-ex1-glucose (pending, out of scope: Example 1 is a general
#                               balanced model with several random factors,
#                               outside the Oliva-Aviles class; not collected)
#
# Equivalence with the random-slope GPQ path is NOT tested: the two paths fit
# different members of the class (s = 2, Z_k = [1 | t_k], no X0 versus s = 1,
# Z_k = 1, time in X0). The slope path needs >= 2 distinct time points per lot,
# and then its random-slope column is never degenerate, so no data set makes
# the two models (and their SSE df, n - 2K versus n - 1 - K) coincide. The
# engine is instead checked against an independent implementation of the
# paper's matrix algorithm (W, S, N, L, P and uniroot per draw) on the same
# pivot draws, for q = 0 (diagonal C) and q > 0 (full C).

lli_q <- function(e) suppressWarnings(suppressMessages(e))

## value of `expr` and the classes of its warnings (each muffled)
with_warn <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, class(cnd)[1]); invokeRestart("muffleWarning")
  })
  list(value = v, warnings = unique(w))
}

sulfur_data <- function() {
  d <- utils::read.csv(test_path("fixtures", "liao2005_sulfur", "data.csv"),
                       comment.char = "#", colClasses = c(bottle = "character"))
  d$bottle <- factor(d$bottle, levels = unique(d$bottle))
  d
}
sulfur_printed <- function() {
  p <- utils::read.csv(test_path("fixtures", "liao2005_sulfur", "printed.csv"), comment.char = "#")
  stats::setNames(p$value, p$quantity)
}

nist_ex4_data <- function() {
  l <- trimws(readLines(test_path("fixtures", "nist_ex4", "ex4.dat"), warn = FALSE))
  l <- l[nzchar(l) & !startsWith(l, "#")]
  d <- utils::read.table(text = l[4:75], col.names = c("level", "batch", "strength"))
  d$temp <- c(75, -67)[d$level]
  d$batch <- factor(d$batch)
  d
}

## ---- independent implementation of Oliva-Aviles & Hauser (2025) s3.4 --------
## Matrix form, straight from the paper: W = (X0, Z), SSE = Y'(I - P_W)Y,
## S = (1_K x0', I_K), NY = S (W'W)^{-1} W'Y, G-part C = S (W'W)^{-1} S',
## L, P, Q_i; per draw: uniroot on (-phi1 R_s2e, Inf), R_G = R_tau2 I + R_s2e C
## solved directly. Draws U, Q, Z under withr::with_seed(seed), in that order
## (the package's documented draw order). Shares no code with R/gpq.R.
paper_gpq <- function(y, X0, g, x0, c_resid, P, conf, sides, M, seed) {
  g <- factor(g); K <- nlevels(g); n <- length(y); q <- NCOL(X0)
  Z <- stats::model.matrix(~ g - 1)
  W <- cbind(X0, Z)
  WtWi <- solve(crossprod(W))
  sse <- sum((y - W %*% (WtWi %*% crossprod(W, y)))^2)
  S <- cbind(matrix(rep(x0, each = K), K, q), diag(K))
  NY <- as.numeric(S %*% WtWi %*% crossprod(W, y))
  C <- S %*% WtWi %*% t(S)
  L <- qr.Q(qr(diag(K) - 1 / K))[, seq_len(K - 1), drop = FALSE]
  e <- eigen(t(L) %*% C %*% L, symmetric = TRUE)
  phi1 <- min(eigen(C, symmetric = TRUE)$values)
  lam <- pmax(e$values, phi1)
  Qi <- as.numeric(t(e$vectors) %*% t(L) %*% NY)^2
  dr <- withr::with_seed(seed, list(U = stats::rchisq(M, n - q - K),
                                    Q = stats::rchisq(M, K - 1), Z = stats::rnorm(M)))
  Rth <- cv <- vth <- numeric(M)
  for (m in seq_len(M)) {
    s2 <- sse / dr$U[m]
    h <- function(a) sum(Qi / (a + lam * s2))
    lo <- -phi1 * s2
    h_lo <- h(lo + 1e-9 * max(s2, 1))
    t2 <- if (is.finite(h_lo) && dr$Q[m] >= h_lo) lo else
      stats::uniroot(function(a) h(a) - dr$Q[m], c(lo + 1e-12 * s2, lo + s2 + 1),
                     extendInt = "downX", tol = 1e-15 * max(s2, 1), maxiter = 5000)$root
    ## singular R_G at the clamp: the paper's right-hand limit, taken just
    ## inside the range (1e-10 x s2 phi1 above the bound)
    t2_eff <- if (t2 <= lo) lo + 1e-10 * s2 * phi1 else t2
    RG <- t2_eff * diag(K) + s2 * C
    Ri <- solve(RG, tol = 0)
    vth[m] <- 1 / sum(Ri)
    Rth[m] <- vth[m] * sum(Ri %*% NY) - dr$Z[m] * sqrt(vth[m])
    cv[m] <- max(0, t2 + c_resid * s2)
  }
  if (sides == "two.sided") {
    T <- stats::median(Rth); z <- stats::qnorm((1 + P) / 2)
    r <- stats::quantile(sqrt(cv + vth), conf, names = FALSE)
    c(T - z * r, T + z * r)
  } else if (sides == "lower") {
    c(stats::quantile(Rth - stats::qnorm(P) * sqrt(cv), 1 - conf, names = FALSE), Inf)
  } else {
    c(-Inf, stats::quantile(Rth + stats::qnorm(P) * sqrt(cv), conf, names = FALSE))
  }
}

test_that("liao2005-sulfur: the transcribed data reproduce the printed bottle means and sse", {
  d <- sulfur_data(); p <- sulfur_printed()
  expect_identical(dim(d), c(10L, 3L))
  expect_identical(nlevels(d$bottle), 6L)
  expect_identical(as.integer(table(d$bottle)), c(2L, 2L, 2L, 2L, 1L, 1L))
  ybar <- tapply(d$sulfur, d$bottle, mean)
  expect_equal(round(as.numeric(ybar), 3),
               unname(p[paste0("ybar_", levels(d$bottle))]), tolerance = 0)
  expect_equal(round(sum((d$sulfur - ybar[d$bottle])^2), 4), unname(p["sse"]), tolerance = 0)
})

test_that("liao2005-sulfur: two-sided (0.99, 0.95) GPQ TI for N(mu, sigma_alpha^2) reproduces [4.473, 4.870] within MC error", {
  d <- sulfur_data(); p <- sulfur_printed()
  fm <- lli_q(lme4::lmer(sulfur ~ 1 + (1 | bottle), data = d))
  M <- 200000L; seed <- 1L
  rw <- with_warn(ti_lmm(fm, P = 0.99, level = 0.95, over = "group_mean", method = "gpq",
                         control = lmm_interval_control(seed = seed, M = M)))
  expect_identical(rw$warnings, "reintervals_warning_unvalidated")
  r <- rw$value
  dg <- attr(r, "diagnostics")
  expect_identical(r$method, "gpq")
  expect_identical(dg$c_resid, 0)
  expect_identical(c(dg$K, dg$df_resid), c(6L, 4L))
  expect_equal(dg$phi1, 0.5)                         # 1 / max n_k
  ## our MC SD (20 MC sections) is below 0.002 at M = 200000
  expect_lt(max(dg$mc_se_lower, dg$mc_se_upper), 0.002)
  ## The printed limits carry the authors' own Monte-Carlo error (K = 10000
  ## runs). The MC SD scales as 1 / sqrt(M), so theirs is about
  ## mc_se x sqrt(M / 10000) (~0.0047; 200 runs of 10000 draws give 0.0049).
  ## Gate: 3 x the combined MC SD plus half a printed unit.
  sd_ours <- c(dg$mc_se_lower, dg$mc_se_upper)
  sd_theirs <- sd_ours * sqrt(M / 10000)
  got <- c(r$lower, r$upper); want <- unname(p[c("lower", "upper")])
  gate <- 3 * sqrt(sd_ours^2 + sd_theirs^2) + 0.0005
  expect_true(all(abs(got - want) <= gate),
              info = sprintf("got [%.4f, %.4f], printed [%.3f, %.3f], gate %.4f",
                             got[1], got[2], want[1], want[2], gate[1]))
  ## regression pin (seed 1, M = 200000)
  expect_equal(got, c(4.466148, 4.876889), tolerance = 1e-6)
  cat(sprintf(paste0("\nliao2005-sulfur: printed [%.3f, %.3f]; package [%.4f, %.4f] ",
                     "(M = %d, seed %d, MC SD %.5f / %.5f; deviation %.4f = %.1f x the ",
                     "authors' K = 10000 MC SD)\n"),
              want[1], want[2], got[1], got[2], M, seed, sd_ours[1], sd_ours[2],
              got[1] - want[1], (got[1] - want[1]) / sd_theirs[1]))
})

test_that("the package GPQ equals an independent implementation of the paper's algorithm on the same draws (q = 0 and q > 0)", {
  ## q = 0: the one-way sulfur data (C diagonal, the existing arithmetic)
  d <- sulfur_data()
  pre <- .gpq_prepare_ri(d$sulfur, matrix(0, nrow(d), 0), d$bottle)
  for (s in c("two.sided", "lower", "upper")) for (cr in c(0, 1)) {
    got <- .gpq_run(pre, list(numeric(0)), 0.9, 0.95, s, cr, M = 2000L, seed = 7L)
    want <- paper_gpq(d$sulfur, matrix(0, nrow(d), 0), d$bottle, numeric(0),
                      cr, 0.9, 0.95, s, 2000L, 7L)
    expect_equal(c(got$lower, got$upper), want, tolerance = 1e-7, info = paste(s, cr))
  }
  ## q = 1: NIST RECIPE ex4 (batch random intercept, fixed temperature), x0 off the data
  x <- nist_ex4_data()
  X0 <- matrix(x$temp, ncol = 1)
  pre <- .gpq_prepare_ri(x$strength, X0, x$batch)
  for (t0 in c(-67, 0, 120)) for (s in c("two.sided", "lower")) {
    got <- .gpq_run(pre, list(t0), 0.9, 0.95, s, 1, M = 2000L, seed = 3L)
    want <- paper_gpq(x$strength, X0, x$batch, t0, 1, 0.9, 0.95, s, 2000L, 3L)
    expect_equal(c(got$lower, got$upper), want, tolerance = 1e-7, info = paste(t0, s))
  }
  ## q = 2 with single-observation lots: time + a within-lot factor
  set.seed(4)
  nk <- c(4, 1, 3, 1, 5, 2)
  dd <- data.frame(lot = factor(rep(seq_along(nk), nk)))
  dd$time <- unlist(lapply(nk, function(k) sort(sample(c(0, 3, 6, 9, 12, 18), k))))
  dd$site <- factor(rep(c("a", "b"), length.out = nrow(dd)))
  dd$y <- 100 - 0.3 * dd$time + 0.5 * (dd$site == "b") +
    rnorm(length(nk), 0, 1)[dd$lot] + rnorm(nrow(dd), 0, 0.7)
  fm <- lli_q(lme4::lmer(y ~ time + site + (1 | lot), data = dd))
  X <- lme4::getME(fm, "X"); X0 <- X[, -1, drop = FALSE]
  pre <- .gpq_prepare_ri(dd$y, X0, dd$lot)
  expect_identical(pre$df_resid, nrow(dd) - 2L - 6L)
  for (x0 in list(c(24, 0), c(6, 1))) for (cr in c(0, 1)) {
    got <- .gpq_run(pre, list(x0), 0.95, 0.9, "two.sided", cr, M = 2000L, seed = 11L)
    want <- paper_gpq(dd$y, X0, dd$lot, x0, cr, 0.95, 0.9, "two.sided", 2000L, 11L)
    expect_equal(c(got$lower, got$upper), want, tolerance = 1e-7,
                 info = paste(paste(x0, collapse = ","), cr))
  }
})

## Liao, Lin & Iyer (2005), Example 3 (section 5.4): printed canonical-form inputs
## (LaMotte-McWhorter decomposition, X = (1, temp), Z = batch indicators) and
## lower (beta 0.90, gamma 0.95) limits for N(mu, sigma_alpha^2 + sigma_e^2)
lli_ex3 <- list(
  sse = 8663.776, sse_df = 63L,
  lambda = c(4.093, 5.407, 6.650, 11, 11.778), r = c(1L, 1L, 1L, 3L, 1L),
  q = c(1796.312, 1587.788, 11.492, 4538.396, 3231.697),
  t0 = c(-67, -50, -25, 0, 25, 50, 75),
  gpq = c(314.661, 312.216, 308.942, 306.056, 303.772, 302.392, 301.511),
  mls = c(286.895, 285.580, 283.558, 281.471, 279.336, 277.120, 274.784))

test_that("liao2005-ex3-nist: the printed canonical form is reproduced from the fixture, and it is the engine's Q pivot", {
  x <- nist_ex4_data(); n <- nrow(x)
  X <- cbind(1, x$temp); Z <- stats::model.matrix(~ batch - 1, x)
  H <- qr.Q(qr(X), complete = TRUE)[, 3:n]                 # H'H = I, HH' = I - P_X
  e <- eigen(t(H) %*% Z %*% t(Z) %*% H, symmetric = TRUE)
  ev <- e$values; Hy <- t(H) %*% x$strength
  ## eigenvalue 0 with multiplicity 63: its quadratic form is the SSE
  null <- which(ev < 1e-8)
  expect_identical(length(null), lli_ex3$sse_df)
  expect_equal(round(sum((t(e$vectors[, null]) %*% Hy)^2), 3), lli_ex3$sse, tolerance = 0)
  ## distinct non-zero eigenvalues, their multiplicities and q_i, at printed digits
  pos <- ev[-null]
  grp <- split(seq_along(pos), round(pos, 3))
  lam <- as.numeric(names(grp)); ord <- order(lam)
  expect_equal(lam[ord], lli_ex3$lambda, tolerance = 0)
  expect_identical(lengths(grp)[ord], lli_ex3$r, ignore_attr = TRUE)
  qi <- vapply(grp, function(i) sum((t(e$vectors[, i]) %*% Hy)^2), 0)[ord]
  expect_equal(round(unname(qi), 3), lli_ex3$q, tolerance = 0)
  ## the engine's pivot (Oliva-Aviles & Hauser eq. 5, B = L'CL) is the same
  ## chi-square: its eigenvalues are 1 / lambda_i and Q_i^OH = q_i / lambda_i,
  ## at any x0, so the two papers use the same model and the same Q pivot
  pre <- .gpq_prepare_ri(x$strength, matrix(x$temp, ncol = 1), x$batch)
  expect_equal(pre$sse, lli_ex3$sse, tolerance = 1e-7)
  expect_identical(pre$df_resid, lli_ex3$sse_df)
  L <- .ortho_complement_of_one(pre$K)
  for (t0 in c(-67, 30)) {
    dz <- .gpq_design_at(pre, t0)
    eb <- eigen(crossprod(L, dz$C %*% L), symmetric = TRUE)
    w <- as.numeric(crossprod(eb$vectors, crossprod(L, dz$NY)))^2
    expect_equal(sort(1 / eb$values), sort(pos), tolerance = 1e-8)
    oh_q <- tapply(w / eb$values, round(1 / eb$values, 3), sum)   # Q_i^OH x lambda_i
    expect_equal(as.numeric(oh_q), unname(qi), tolerance = 1e-8)
  }
})

test_that("liao2005-ex3-nist: package GPQ lower limits next to the printed GPQ and Vangel MLS (comparison; regression pin only)", {
  x <- nist_ex4_data()
  fm <- lli_q(lme4::lmer(strength ~ temp + (1 | batch), data = x))
  nd <- data.frame(temp = lli_ex3$t0)
  rw <- with_warn(ti_lmm(fm, nd, P = 0.90, level = 0.95, sides = "lower", method = "gpq",
                         control = lmm_interval_control(seed = 2005L, M = 20000L)))
  expect_true("reintervals_warning_unvalidated" %in% rw$warnings)
  r <- rw$value
  expect_identical(r$method, rep("gpq", 7))
  dg <- attr(r, "diagnostics")
  ## regression pin (seed 2005, M = 20000), not an oracle
  pin <- c(284.4342, 283.8185, 282.8228, 281.7851, 280.6848, 279.5615, 278.3568)
  expect_equal(r$lower, pin, tolerance = 1e-6)
  ## plug-in 10th percentile of N(mu, sigma_a^2 + sigma_e^2) at the REML fit:
  ## no valid 95% lower confidence limit for it can sit far above it
  vc <- as.data.frame(lme4::VarCorr(fm))$vcov
  plug <- unname(lme4::fixef(fm)[1] + lme4::fixef(fm)[2] * lli_ex3$t0) - stats::qnorm(0.9) * sqrt(sum(vc))
  tab <- data.frame(t0 = lli_ex3$t0, package_gpq = round(r$lower, 3),
                    mc_se = round(dg$mc_se_lower, 3), printed_gpq = lli_ex3$gpq,
                    printed_vangel_mls = lli_ex3$mls, plug_in_q10 = round(plug, 3))
  cat("\nliao2005-ex3-nist: lower (0.90, 0.95) limits\n")
  print(tab, row.names = FALSE)
  expect_true(all(lli_ex3$gpq > plug + 5))      # the printed GPQ limits exceed the plug-in quantile
  expect_true(all(r$lower < plug))
})

## ---- PENDING (not to be collected): liao2005-ex1-glucose ----------------------
## Liao, Lin & Iyer (2005) Example 1 (glucose meters): a general balanced mixed
## model with several random factors (meter, blood sample, strip lot), outside the
## Oliva-Aviles & Hauser class (one random factor). Recorded in oracles.csv so the
## decision is visible; there is no engine to compare with.
## test_that("liao2005-ex1-glucose: out of scope", {
##   ## would need a multi-random-factor GPQ engine
## })
