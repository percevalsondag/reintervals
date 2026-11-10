## ---------------------------------------------------------------------------
## gpq.R      Oliva-Aviles & Hauser (2025, Technometrics 67:2, 193-202)
##            generalized-pivotal-quantity (GPQ) (P-content, gamma-confidence)
##            tolerance interval for the RANDOM-INTERCEPT, RANDOM-SLOPE
##            single-temperature stability model (their Eq. 3):
##
##              y_kj = (eta1 + a_{k,1}) + (eta2 + a_{k,2}) t_kj + e_kj,
##              a_k = (a_{k,1}, a_{k,2})' ~ N(0, V)   [V correlated or diagonal],
##              e_kj ~ N(0, s2e),  batches independent.
##
## The paper's class (their Eq. 2) puts all fixed effects inside the column
## space of Z, so for the single-temperature stability model q = 0, s = 2,
## X0 is empty, A = I_2, and W = Z = block-diag of Z_k = [1 | t_k] (n x 2K).
## Then:
##   * W'W = blockdiag(Z_k'Z_k)            -> (W'W)^{-1} = blockdiag((Z_k'Z_k)^{-1})
##   * S   = blockdiag(z') with z = (1,t0) -> S(W'W)^{-1}S' = diag(a_k),
##           a_k = z'(Z_k'Z_k)^{-1} z   (the OLS leverage of batch k's line at t0)
##   * NY  = (z' bhat_k^OLS)_k            -> batch k's own OLS fitted value at t0
##   * SSE = sum_k within-batch residual SS,  U = SSE/s2e ~ chi^2(n - 2K)   (Eq. 4)
##
## So the three independent pivots of the paper specialize to:
##   U  ~ chi^2(n - 2K)                                    (-> GPQ for s2e)
##   Q  = sum_i Q_i / (tau2 + lambda_i s2e) ~ chi^2(K-1)   (-> GPQ for tau2)   (Eq. 5)
##   Z  ~ N(0,1)                                           (-> GPQ for theta)
## where lambda_i are the eigenvalues of B = L' diag(a_k) L (L = orthonormal
## complement of 1_K), Q_i = (P_i' L' NY)^2 with B = P diag(lambda) P', and
## phi1 = min_k a_k bounds how negative the tau2 realization may go.
##
## tau2 = z' V z is the ONLY functional of V the method needs -- individual
## variance components are never estimated. Hence its immunity to the
## correlated-slope identification pathology (correlation pinned to +/-1,
## singular REML fits) that degraded the REML closed-form interval removed in
## 1.0.0 (reml-mls): no V to identify.
##
## TARGET / RESIDUAL MULTIPLIER c:
##   c = 1  -> content distribution N(theta, tau2 + s2e): a FUTURE OBSERVATION.
##             This is ti_lmm(over = "observation"): content variance
##             V_T = V_G + s2e.
##   c = 0  -> content distribution N(theta, tau2): the BATCH-MEAN distribution
##             (the target of the paper's own Section 5 / shelf-life example;
##             ti_lmm(over = "group_mean")).
##
## Verified against the paper's Section 5 coverage study (see data-raw/gpq_coverage_study.R):
## with c = 0, t0 = 20, their DGP, the GPQ (P,gamma)-TI covers ~ gamma.
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## .gpq_core: the Oliva-Aviles & Hauser (2025) GPQ tolerance-interval Monte
## Carlo. Internal engine used by both the public raw-data wrapper ti_gpq_raw()
## and the slope dispatch (ti_lmm). Carries content/conf/M/seed internally; the
## public wrapper does the P/alpha/level/control harmonization. Provenance:
## method = "gpq".
##
## Structure (A2, plumbing only; the pivots, the content-variance construction
## and the quantile definitions are unchanged):
##   .gpq_prepare  data checks + per-lot OLS pieces (deterministic)
##   .gpq_pivots   the M pivot draws (U, Q, Z), drawn ONCE per call under one
##                 recorded integer seed
##   .gpq_eval_t0  everything t0-dependent, evaluated on those same draws
## Common random numbers: every t0 in the call, and both limits, use the same
## pivot draws, so a limit at a time point does not depend on which other time
## points were requested. The draw order (U, Q, Z) is the 0.2.0 / A1 order, so
## the FIRST t0 of a seeded call reproduces the A1 numbers exactly.
## ---------------------------------------------------------------------------
.gpq_core <- function(y, time, batch, t0, content = 0.99, conf = 0.95,
                      sides = c("two.sided", "lower", "upper"),
                      c_resid = 1, M = 20000L, seed = NULL) {
  sides <- match.arg(sides)
  pre <- .gpq_prepare(y, time, batch)          # errors here leave the RNG untouched
  piv <- .gpq_pivots(df_resid = pre$df_resid, K = pre$K, M = M, seed = seed)

  zP   <- stats::qnorm(content)             # one-sided content quantile
  zP2  <- stats::qnorm((1 + content) / 2)   # two-sided content quantile

  rows <- lapply(t0, function(tt)
    .gpq_eval_t0(pre, piv, tt, zP = zP, zP2 = zP2, conf = conf, sides = sides,
                 c_resid = c_resid, content = content))
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "gpq_seed")    <- piv$seed
  attr(out, "rng_kind")    <- piv$rng_kind
  attr(out, "mc_sections") <- .gpq_mc_sections
  out
}

## Number of MC sections for the Monte-Carlo standard error of each limit.
## ("Sections", not "batches": a batch is a stability lot in this package.)
.gpq_mc_sections <- 20L

## Advice appended to the reintervals_error_time_points message.
.time_points_advice <- paste0(
  "Refit without these lots, or, if the between-lot slope variance is negligible, ",
  "use a random-intercept model with a fixed slope (y ~ time + (1 | lot)), whose ",
  "tolerance interval (anova-mls) accepts single-observation lots.")

## Data checks + per-lot OLS pieces: (Z_k'Z_k)^{-1}, bhat_k, within-lot SSE.
## Deterministic; runs before any draw, so a data error never advances the
## caller's RNG stream.
.gpq_prepare <- function(y, time, batch) {
  batch <- factor(batch)
  ok <- !is.na(y) & !is.na(time) & !is.na(batch)
  y <- y[ok]; time <- time[ok]; batch <- droplevels(batch[ok])
  K <- nlevels(batch); n <- length(y)
  if (K < 2L) .abort(sprintf("ti_gpq_raw(): the GPQ tolerance interval needs >= 2 batches (grouping levels) to estimate the between-batch variance; got %d.", K),
                     class = "reintervals_error_groups")
  ## every lot needs >= 2 distinct time points for its own intercept + slope
  n_t <- tapply(time, batch, function(v) length(unique(v)))
  few <- names(n_t)[n_t < 2L]
  if (length(few)) {
    shown <- paste(few[seq_len(min(5L, length(few)))], collapse = ", ")
    if (length(few) > 5L) shown <- sprintf("%s and %d more", shown, length(few) - 5L)
    .abort(sprintf(paste0("A random-slope tolerance interval needs at least two ",
                          "distinct time points per lot; %d lot(s) have fewer: %s. %s"),
                   length(few), shown, .time_points_advice),
           class = "reintervals_error_time_points", data = list(lots = few))
  }
  s <- 2L                                   # random intercept + random slope
  df_resid <- n - s * K                     # n - 2K
  if (df_resid < 1L)
    .abort("ti_gpq_raw(): non-positive residual df (need more within-batch replication).",
           class = "reintervals_error_residual_df")

  bl  <- split(seq_len(n), batch)
  ZtZi <- vector("list", K)                 # (Z_k'Z_k)^{-1}, each 2x2
  bhat <- matrix(0, 2L, K)                  # OLS (intercept, slope) per batch
  sse  <- 0
  for (k in seq_len(K)) {
    idx <- bl[[k]]
    tk  <- time[idx]; yk <- y[idx]
    Zk  <- cbind(1, tk)                     # n_k x 2
    ZtZ <- crossprod(Zk)                    # 2x2
    ## numerical backstop (distinct but numerically coincident time points)
    if (qr(ZtZ)$rank < 2L)
      .abort(sprintf(paste0("A random-slope tolerance interval needs at least two ",
                            "distinct time points per lot; lot %s's time points are ",
                            "numerically indistinguishable. %s"), levels(batch)[k],
                     .time_points_advice),
             class = "reintervals_error_time_points", data = list(lots = levels(batch)[k]))
    ZtZinv <- solve(ZtZ)
    bk  <- ZtZinv %*% crossprod(Zk, yk)     # OLS coefficients
    ZtZi[[k]] <- ZtZinv
    bhat[, k] <- bk
    sse <- sse + sum((yk - Zk %*% bk)^2)    # within-batch residual SS
  }
  sse <- as.numeric(sse)
  ## every lot's data on its own straight line: the residual-variance pivot
  ## sse / U is then 0 in every draw. "Zero" is zero to round-off (relative to
  ## the total sum of squares), since the per-lot OLS fit leaves ~1e-28 residues
  ## even on exactly linear data.
  if (sse <= 1e-12 * max(sum((y - mean(y))^2), .Machine$double.xmin)) {
    .warn("The pooled within-lot residual sum of squares is zero (every lot's ",
          "data lie on a straight line), so the residual-variance pivot is ",
          "degenerate and the interval reflects only the between-lot variance.",
          class = "reintervals_warning_limitation")
  }
  list(K = K, n = n, df_resid = df_resid, ZtZi = ZtZi, bhat = bhat, sse = sse)
}

## The M pivot draws, once per call, in the fixed order U, Q, Z.
## seed = NULL: one integer seed is drawn from the caller's ambient stream
## (advancing it by exactly one draw), and the pivots are drawn under it.
## Either way the draws run inside withr::with_seed(), which restores the
## caller's RNG state afterwards, and the integer seed and the RNGkind() in
## effect during the draws are returned for the diagnostics.
## z_new = TRUE (the generalized prediction interval for a new lot's mean,
## new_group_mean_lmm(slope_ngm_method = "gpq")) draws one more standard-normal
## vector Z_new AFTER U, Q, Z, so U, Q, Z -- and every tolerance interval --
## are the same draws with or without it.
.gpq_pivots <- function(df_resid, K, M, seed = NULL, z_new = FALSE) {
  M <- as.integer(M)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.integer(seed)
  withr::with_seed(seed, {
    rng_kind <- RNGkind()
    U <- stats::rchisq(M, df_resid)
    Q <- stats::rchisq(M, K - 1L)
    Z <- stats::rnorm(M)
    Z_new <- if (z_new) stats::rnorm(M) else NULL
    list(U = U, Q = Q, Z = Z, Z_new = Z_new, M = M, seed = seed,
         rng_kind = rng_kind)
  })
}

## The t0-dependent pivots on the shared draws `piv` (A2 code, moved here
## unchanged in B2 so the new-group-mean GPQ uses the same realizations):
##   Rtheta  GPQ of the mean line mu(t0)                 (R_mu(t0))
##   Rtau2   GPQ of the between-lot variance z0' V z0    (R_VG(t0), c_resid = 0)
##   Rs2e    GPQ of the residual variance
## plus the draw bookkeeping: `keep` (draws with a valid tau2 solution), the MC
## section of each kept draw, the NA count and the boundary flags.
.gpq_theta_draws <- function(pre, piv, tt) {
  K <- pre$K; M <- piv$M
  z0 <- c(1, tt)
  a  <- vapply(seq_len(K), function(k) as.numeric(t(z0) %*% pre$ZtZi[[k]] %*% z0), 0) # leverages a_k
  NY <- as.numeric(crossprod(pre$bhat, z0)) # batch fitted values at t0 (K-vector)
  phi1 <- min(a)

  ## L = orthonormal complement of 1_K  (K x (K-1)); B = L' diag(a) L
  L  <- .ortho_complement_of_one(K)
  B  <- crossprod(L, a * L)                 # (K-1)x(K-1)
  eg <- eigen(B, symmetric = TRUE)
  lam <- pmax(eg$values, 0)                 # eigenvalues lambda_i (>=0)
  P  <- eg$vectors
  w  <- as.numeric(crossprod(P, crossprod(L, NY)))   # P' L' NY
  Qi <- w^2                                 # observed Q_i  (length K-1)

  ## ---- the M shared GPQ realizations ----
  Qstar <- piv$Q
  Zstar <- piv$Z
  Rs2e  <- pre$sse / piv$U                  # GPQ for s2e (M-vector)

  ## solve  sum_i Qi/(a_tau + lambda_i Rs2e) = Qstar  for a_tau in (-phi1 Rs2e, Inf)
  Rtau2 <- .solve_tau2(Qi, lam, Rs2e, Qstar, phi1)
  ## drop draws without a valid tau2 solution (NA); more than 1% is an error
  na <- is.na(Rtau2)
  n_na <- sum(na)
  if (n_na > 0.01 * M) {
    .abort(sprintf(paste0("ti_gpq_raw(): %d of %d GPQ draws (%.1f%%) have no valid ",
                          "between-batch variance solution at t0 = %g; more than ",
                          "1%% is treated as a numerical failure."),
                   n_na, M, 100 * n_na / M, tt),
           class = "reintervals_error_numerical")
  }
  ## MC section of each draw, by draw index BEFORE the NA drop; the remainder
  ## M %% S draws belong to no section (they still enter the limit itself)
  S <- .gpq_mc_sections
  m_sec <- M %/% S
  sec <- (seq_len(M) - 1L) %/% m_sec + 1L
  sec[sec > S] <- NA_integer_
  if (n_na > 0L) {
    Rtau2 <- Rtau2[!na]; Rs2e <- Rs2e[!na]; Zstar <- Zstar[!na]; sec <- sec[!na]
  }
  boundary <- Rtau2 <= (-phi1 * Rs2e + 1e-12 * pmax(Rs2e, 1))

  ## R_G = diag(Rtau2 + Rs2e a_k); weighted-mean center and its variance.
  ## g[m,k] = Rtau2[m] + Rs2e[m]*a[k]; clamp tiny denoms (right-hand limit at boundary).
  g <- outer(Rtau2, rep(1, K)) + outer(Rs2e, a)          # M x K
  g <- pmax(g, 1e-10 * pmax(outer(Rs2e, a), 1e-12))
  inv <- 1 / g
  den <- rowSums(inv)                                    # 1' RG^{-1} 1
  num <- as.numeric(inv %*% NY)                          # 1' RG^{-1} NY
  theta_tilde <- num / den
  Rs2_theta   <- 1 / den                                 # Var(theta-tilde) realization

  Rtheta <- theta_tilde - Zstar * sqrt(Rs2_theta)        # GPQ for theta
  list(Rtheta = Rtheta, Rtau2 = Rtau2, Rs2e = Rs2e, Rs2_theta = Rs2_theta,
       keep = !na, sec = sec,
       n_na = n_na, boundary = boundary, phi1 = phi1, K = K, M = M, S = S)
}

## Everything t0-dependent for the tolerance interval, on the shared pivot
## draws `piv`. Returns one row.
.gpq_eval_t0 <- function(pre, piv, tt, zP, zP2, conf, sides, c_resid, content) {
  dr <- .gpq_theta_draws(pre, piv, tt)
  Rtheta <- dr$Rtheta; Rtau2 <- dr$Rtau2; Rs2e <- dr$Rs2e; sec <- dr$sec
  S <- dr$S; K <- dr$K; phi1 <- dr$phi1; n_na <- dr$n_na; boundary <- dr$boundary
  Rs2_theta <- dr$Rs2_theta
  content_var <- pmax(0, Rtau2 + c_resid * Rs2e)         # tau2 + c s2e (>=0)

  if (sides == "two.sided") {
    Rzeta <- sqrt(content_var + Rs2_theta)               # ζ = sqrt(content var + center var)
  } else if (sides == "lower") {
    Rdelta <- Rtheta - zP * sqrt(content_var)            # GPQ for lower β-quantile
  } else {                                               # upper
    Rdelta <- Rtheta + zP * sqrt(content_var)            # GPQ for upper β-quantile
  }
  ## the limits from the draws indexed by i (all draws: the reported limit;
  ## one MC section: that section's replicate of it)
  limits <- function(i) {
    if (sides == "two.sided") {
      theta_hat <- stats::median(Rtheta[i])
      Rzeta_g   <- stats::quantile(Rzeta[i], conf, names = FALSE, type = 7)
      c(theta_hat - zP2 * Rzeta_g, theta_hat + zP2 * Rzeta_g)
    } else if (sides == "lower") {
      c(stats::quantile(Rdelta[i], 1 - conf, names = FALSE, type = 7), Inf)
    } else {
      c(-Inf, stats::quantile(Rdelta[i], conf, names = FALSE, type = 7))
    }
  }
  lim <- limits(seq_along(Rtheta))
  ## Monte-Carlo SE by MC sections: sd of the S section replicates / sqrt(S)
  sec_lim <- vapply(seq_len(S), function(s) limits(which(sec == s)), numeric(2))
  mc_se <- apply(sec_lim, 1L, stats::sd) / sqrt(S)
  mc_se[!is.finite(lim)] <- NA_real_

  data.frame(method = "gpq", time = tt,
             estimate = stats::median(Rtheta),
             lower = lim[1], upper = lim[2],
             mc_se_lower = mc_se[1], mc_se_upper = mc_se[2],
             tau2_med = stats::median(Rtau2),
             s2e_med = stats::median(Rs2e),
             frac_boundary = mean(boundary), n_na_draws = n_na,
             K = K, df_resid = pre$df_resid, phi1 = phi1,
             content = content, conf = conf, c_resid = c_resid,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## .gpq_ngm_core: generalized prediction interval for the mean of a NEW lot at
## t0, theta_new = mu(t0) + z0' b_new (new_group_mean_lmm, slope_ngm_method =
## "gpq"; a validation arm). On the tolerance interval's pivots (Oliva-Aviles &
## Hauser 2025) the generalized pivot of theta_new is (Krishnamoorthy & Mathew
## 2009, generalized prediction intervals)
##   R_theta_new = R_mu(t0) + Z_new * sqrt(R_VG(t0)),
## with R_mu = Rtheta, R_VG = max(0, Rtau2) (the c_resid = 0 content variance)
## and Z_new an independent N(0, 1) vector drawn after U, Q, Z. Limits are its
## quantiles: alpha/2 and 1 - alpha/2 two-sided, alpha (lower) or 1 - alpha
## (upper) one-sided. Same seed contract, common random numbers across t0, MC
## sections, NA-draw rule and per-lot time-point guard as the TI.
## ---------------------------------------------------------------------------
.gpq_ngm_core <- function(y, time, batch, t0, level = 0.95,
                          sides = c("two.sided", "lower", "upper"),
                          M = 20000L, seed = NULL) {
  sides <- match.arg(sides)
  pre <- .gpq_prepare(y, time, batch)          # errors here leave the RNG untouched
  piv <- .gpq_pivots(df_resid = pre$df_resid, K = pre$K, M = M, seed = seed,
                     z_new = TRUE)
  alpha <- 1 - level
  rows <- lapply(t0, function(tt) {
    dr <- .gpq_theta_draws(pre, piv, tt)
    Rnew <- dr$Rtheta + piv$Z_new[dr$keep] * sqrt(pmax(0, dr$Rtau2))
    limits <- function(i) {
      q <- function(p) stats::quantile(Rnew[i], p, names = FALSE, type = 7)
      switch(sides,
             two.sided = c(q(alpha / 2), q(1 - alpha / 2)),
             lower     = c(q(alpha), Inf),
             upper     = c(-Inf, q(1 - alpha)))
    }
    lim <- limits(seq_along(Rnew))
    sec_lim <- vapply(seq_len(dr$S), function(s) limits(which(dr$sec == s)), numeric(2))
    mc_se <- apply(sec_lim, 1L, stats::sd) / sqrt(dr$S)
    mc_se[!is.finite(lim)] <- NA_real_
    data.frame(time = tt, estimate = stats::median(dr$Rtheta),
               lower = lim[1], upper = lim[2],
               mc_se_lower = mc_se[1], mc_se_upper = mc_se[2],
               tau2_med = stats::median(dr$Rtau2), s2e_med = stats::median(dr$Rs2e),
               frac_boundary = mean(dr$boundary), n_na_draws = dr$n_na,
               K = dr$K, df_resid = pre$df_resid, phi1 = dr$phi1,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "gpq_seed")    <- piv$seed
  attr(out, "rng_kind")    <- piv$rng_kind
  attr(out, "mc_sections") <- .gpq_mc_sections
  out
}

## ---- helpers --------------------------------------------------------------

## orthonormal basis (K x (K-1)) for the orthogonal complement of 1_K.
.ortho_complement_of_one <- function(K) {
  ## Helmert-style: QR of the complement of the all-ones vector.
  A <- diag(K) - matrix(1 / K, K, K)        # centering projector (rank K-1)
  q <- qr(A)
  Q <- qr.Q(q)[, seq_len(K - 1L), drop = FALSE]
  Q
}

## Vectorized solve of  h(a) = sum_i Qi/(a + lambda_i * s2e) = Qstar  for each MC
## draw, with a in (-phi1*s2e, Inf). h is positive, strictly decreasing, convex
## on that range, so a monotone bracket + vectorized bisection is exact & fast.
## Boundary rule (paper Step 3): if Qstar >= h(-phi1*s2e), clamp a = -phi1*s2e.
.solve_tau2 <- function(Qi, lam, s2e, Qstar, phi1) {
  M <- length(s2e)
  lo <- -phi1 * s2e                          # per-draw lower bound on a
  ## h(a) for a vector of candidate a (length M), per-draw s2e:
  hfun <- function(a) {
    ## denom[m,i] = a[m] + lam[i]*s2e[m]
    denom <- outer(a, rep(1, length(lam))) + outer(s2e, lam)
    denom[denom <= 0] <- NA                  # guard; only happens just below lo
    rowSums(sweep(1 / denom, 2, Qi, `*`))
  }
  ## value of h at the lower bound (may be +Inf when lambda_1 == phi1)
  eps <- 1e-9 * pmax(s2e, 1)
  h_lo <- hfun(lo + eps)
  clamp <- is.finite(h_lo) & (Qstar >= h_lo) # interior root does NOT exist
  out <- lo
  interior <- !clamp
  ## A draw whose h is NA anywhere (a denominator a + lambda_i s2e <= 0, which
  ## can only arise from round-off pushing a lambda_i below phi1) has no valid
  ## solution: it is returned as NA, never as a number. The caller drops it and
  ## reports the count.
  bad <- is.na(h_lo)
  if (any(interior)) {
    blo <- (lo + eps)[interior]
    s2e_i <- s2e[interior]; Qs_i <- Qstar[interior]
    ## expand an upper bracket until h(bhi) < Qstar
    bhi <- pmax(blo + 1, blo + 10)
    h_at <- function(a, idx) {
      denom <- outer(a, rep(1, length(lam))) + outer(s2e_i[idx], lam)
      denom[denom <= 0] <- NA
      rowSums(sweep(1 / denom, 2, Qi, `*`))
    }
    bad_i <- bad[interior]
    need <- !bad_i; it <- 0L
    while (any(need) && it < 200L) {
      hv <- rep(Inf, length(blo))
      hv[need] <- h_at(bhi[need], which(need))
      done <- hv < Qs_i
      bad_i[is.na(done)] <- TRUE
      done[is.na(done)] <- TRUE               # stop expanding a bad draw
      need[done] <- FALSE
      bhi[need] <- blo[need] + (bhi[need] - blo[need]) * 2
      it <- it + 1L
    }
    ## vectorized bisection
    for (it in seq_len(80L)) {
      mid <- 0.5 * (blo + bhi)
      hv  <- h_at(mid, seq_along(mid))
      go_up <- hv > Qs_i                      # h decreasing: too high -> move lo up
      bad_i[is.na(go_up)] <- TRUE
      go_up[is.na(go_up)] <- FALSE
      blo[go_up] <- mid[go_up]
      bhi[!go_up] <- mid[!go_up]
    }
    out[interior] <- 0.5 * (blo + bhi)
    bad[interior] <- bad_i
  }
  out[bad] <- NA_real_
  out
}

#' GPQ tolerance interval from raw data (direct engine, advanced)
#'
#' The Oliva-Aviles & Hauser (2025) generalized-pivotal-
#' quantity tolerance interval for the random-intercept + random-slope stability
#' model, from raw `(y, time, batch)` vectors without a fitted `lmer` model. The
#' fit-based dispatcher [ti_lmm()] uses it for every random-slope tolerance
#' interval (both targets, all sides); use `ti_gpq_raw()` directly when you have
#' data but no model object.
#'
#' @param y Numeric response.
#' @param time Numeric time/age covariate (same length as `y`).
#' @param batch Batch factor (same length as `y`).
#' @param t0 Time point(s) at which to evaluate the interval.
#' @param P Content proportion (default 0.95; pharma stability typically 0.99).
#' @param alpha,level Confidence as a mutually-exclusive pair (give at most one);
#'   `level = 1 - alpha`, default `alpha = 0.05`.
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`.
#' @param c_resid Residual multiplier `c` in `N(theta, tau2 + c*s2e)`: `1`
#'   (default) a future observation, `0` the batch-mean.
#' @param control [lmm_interval_control()] tuning (the GPQ `seed` and `M`).
#'
#' @details
#' The seed handling, common random numbers across `t0`, and the Monte-Carlo
#' standard error are those of [ti_lmm()] (see its Details): one pivot draw set
#' per call, reused for every `t0` and both limits; the seed used is
#' `attr(., "diagnostics")$gpq_seed` (with `$rng_kind`), and
#' `lmm_interval_control(seed = <gpq_seed>)` with the same `M` under the same RNG
#' kind reproduces the result. `mc_se_lower` / `mc_se_upper` are the run-to-run
#' Monte-Carlo jitter of each limit (20 MC sections), not statistical
#' uncertainty. Every batch needs at least two distinct time points (error
#' class `reintervals_error_time_points`, naming the batches).
#'
#' @return An `lmm_interval` data frame, one row per `t0`: `time`, `type`,
#'   `estimate`, `lower`, `upper`, `level`, `P`, `sides`, `method = "gpq"`,
#'   `design = "random_slope_correlated"`, plus the flag columns. The
#'   diagnostics attribute carries the same fields as the [ti_lmm()] GPQ route
#'   (`df`, `tau2_med`, `s2e_med`, `frac_boundary`, `K`, `df_resid`, `phi1`,
#'   `c_resid`, `route_reason = "method"`, `n_na_draws`, `mc_se_lower`,
#'   `mc_se_upper`, `gpq_seed`, `rng_kind`, `mc_sections`).
#' @references Oliva-Aviles C, Hauser P (2025). *Technometrics* 67(2):193-202. \doi{10.1080/00401706.2024.2407324}
#' @seealso [ti_lmm()] (the dispatcher), [ti_anova_raw()] (the ANOVA raw engine),
#'   [reintervals-models].
#' @examples
#' fm <- lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = lme4::sleepstudy)
#' d <- lme4::sleepstudy
#' ti_gpq_raw(d$Reaction, d$Days, d$Subject, t0 = 5,
#'            control = lmm_interval_control(seed = 1, M = 2000))
#' @export
ti_gpq_raw <- function(y, time, batch, t0, P = 0.95, alpha = NULL, level = NULL,
                       sides = c("two.sided", "lower", "upper"),
                       c_resid = 1, control = lmm_interval_control()) {
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal `conf`; public `P` -> `content`.
  conf  <- 1 - .resolve_alpha(alpha, level)
  g <- .gpq_core(y, time, batch, t0, content = P, conf = conf, sides = sides,
                 c_resid = c_resid, M = control$M, seed = control$seed)
  .warn_extrapolation(t0, time, random = TRUE, var = "t0")
  ## the same row schema (and so the same diagnostics) as the ti_lmm() GPQ route;
  ## route_reason = "method": the GPQ engine was called explicitly
  out <- .gpq_rows(g, design = "random_slope_correlated", sides = sides,
                   singular = FALSE)
  out$route_reason <- "method"
  .unify_intervals(out, model = NULL, kind = "ti", P_value = P, verb = "ti_gpq_raw()")
}
