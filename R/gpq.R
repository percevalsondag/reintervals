## ---------------------------------------------------------------------------
## gpq.R  --  Oliva-Aviles & Hauser (2025, Technometrics 67:2, 193-202)
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
## singular REML fits) that degrades the closed form: no V to identify.
##
## TARGET / RESIDUAL MULTIPLIER c:
##   c = 1  -> content distribution N(theta, tau2 + s2e): a FUTURE OBSERVATION.
##             This matches our TI target (scored against V_T = V_G+s2e),
##             so c = 1 is the apples-to-apples head-to-head setting.
##   c = 0  -> content distribution N(theta, tau2): the BATCH-MEAN distribution
##             (the target of the paper's own Section 5 / shelf-life example).
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
## ---------------------------------------------------------------------------
.gpq_core <- function(y, time, batch, t0, content = 0.99, conf = 0.95,
                      sides = c("two.sided", "lower", "upper"),
                      c_resid = 1, M = 10000L, seed = NULL) {
  sides <- match.arg(sides)
  if (!is.null(seed)) { old <- .Random.seed_safe(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed) }

  batch <- factor(batch)
  ok <- !is.na(y) & !is.na(time) & !is.na(batch)
  y <- y[ok]; time <- time[ok]; batch <- droplevels(batch[ok])
  K <- nlevels(batch); n <- length(y)
  if (K < 2L) stop(sprintf("ti_gpq_raw(): the GPQ tolerance interval needs >= 2 batches (grouping levels) to estimate the between-batch variance; got %d.", K), call. = FALSE)
  s <- 2L                                   # random intercept + random slope
  df_resid <- n - s * K                     # n - 2K
  if (df_resid < 1L)
    stop("ti_gpq_raw(): non-positive residual df (need more within-batch replication).")

  ## ---- per-batch OLS pieces: (Z_k'Z_k)^{-1}, bhat_k, within-batch SSE ----
  bl  <- split(seq_len(n), batch)
  ZtZi <- vector("list", K)                 # (Z_k'Z_k)^{-1}, each 2x2
  bhat <- matrix(0, 2L, K)                  # OLS (intercept, slope) per batch
  sse  <- 0
  for (k in seq_len(K)) {
    idx <- bl[[k]]
    tk  <- time[idx]; yk <- y[idx]
    Zk  <- cbind(1, tk)                     # n_k x 2
    ZtZ <- crossprod(Zk)                    # 2x2
    if (qr(ZtZ)$rank < 2L)
      stop(sprintf("ti_gpq_raw(): batch %s has < 2 distinct time points.", levels(batch)[k]))
    ZtZinv <- solve(ZtZ)
    bk  <- ZtZinv %*% crossprod(Zk, yk)     # OLS coefficients
    ZtZi[[k]] <- ZtZinv
    bhat[, k] <- bk
    sse <- sse + sum((yk - Zk %*% bk)^2)    # within-batch residual SS
  }
  sse <- as.numeric(sse)

  zP   <- stats::qnorm(content)             # one-sided content quantile
  zP2  <- stats::qnorm((1 + content) / 2)   # two-sided content quantile

  ## ---- per-t0 GPQ Monte Carlo ----
  rows <- lapply(t0, function(tt) {
    z0 <- c(1, tt)
    a  <- vapply(seq_len(K), function(k) as.numeric(t(z0) %*% ZtZi[[k]] %*% z0), 0) # leverages a_k
    NY <- as.numeric(crossprod(bhat, z0))   # batch fitted values at t0 (K-vector)
    phi1 <- min(a)

    ## L = orthonormal complement of 1_K  (K x (K-1)); B = L' diag(a) L
    L  <- .ortho_complement_of_one(K)
    B  <- crossprod(L, a * L)               # (K-1)x(K-1)
    eg <- eigen(B, symmetric = TRUE)
    lam <- pmax(eg$values, 0)               # eigenvalues lambda_i (>=0)
    P  <- eg$vectors
    w  <- as.numeric(crossprod(P, crossprod(L, NY)))   # P' L' NY
    Qi <- w^2                               # observed Q_i  (length K-1)

    ## ---- M independent GPQ realizations (vectorized) ----
    Ustar <- stats::rchisq(M, df_resid)
    Qstar <- stats::rchisq(M, K - 1L)
    Zstar <- stats::rnorm(M)
    Rs2e  <- sse / Ustar                    # GPQ for s2e (M-vector)

    ## solve  sum_i Qi/(a_tau + lambda_i Rs2e) = Qstar  for a_tau in (-phi1 Rs2e, Inf)
    Rtau2 <- .solve_tau2(Qi, lam, Rs2e, Qstar, phi1)
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

    content_var <- pmax(0, Rtau2 + c_resid * Rs2e)         # tau2 + c s2e (>=0)
    Rtheta <- theta_tilde - Zstar * sqrt(Rs2_theta)        # GPQ for theta

    if (sides == "two.sided") {
      Rzeta <- sqrt(content_var + Rs2_theta)               # ζ = sqrt(content var + center var)
      theta_hat <- stats::median(Rtheta)
      Rzeta_g   <- stats::quantile(Rzeta, conf, names = FALSE, type = 7)
      lower <- theta_hat - zP2 * Rzeta_g
      upper <- theta_hat + zP2 * Rzeta_g
    } else if (sides == "lower") {
      Rdelta <- Rtheta - zP * sqrt(content_var)            # GPQ for lower β-quantile
      lower  <- stats::quantile(Rdelta, 1 - conf, names = FALSE, type = 7)
      upper  <- Inf
    } else {                                               # upper
      Rdelta <- Rtheta + zP * sqrt(content_var)            # GPQ for upper β-quantile
      lower  <- -Inf
      upper  <- stats::quantile(Rdelta, conf, names = FALSE, type = 7)
    }

    data.frame(method = "gpq", time = tt,
               estimate = stats::median(Rtheta),
               lower = lower, upper = upper,
               tau2_med = stats::median(Rtau2),
               s2e_med = stats::median(Rs2e),
               frac_boundary = mean(boundary),
               K = K, df_resid = df_resid, phi1 = phi1,
               content = content, conf = conf, c_resid = c_resid,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
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
    need <- rep(TRUE, length(blo)); it <- 0L
    while (any(need) && it < 200L) {
      hv <- rep(Inf, length(blo))
      hv[need] <- h_at(bhi[need], which(need))
      done <- hv < Qs_i
      need[done] <- FALSE
      bhi[need] <- blo[need] + (bhi[need] - blo[need]) * 2
      it <- it + 1L
    }
    ## vectorized bisection
    for (it in seq_len(80L)) {
      mid <- 0.5 * (blo + bhi)
      hv  <- h_at(mid, seq_along(mid))
      go_up <- hv > Qs_i                      # h decreasing: too high -> move lo up
      blo[go_up] <- mid[go_up]
      bhi[!go_up] <- mid[!go_up]
    }
    out[interior] <- 0.5 * (blo + bhi)
  }
  out
}

## RNG-state save/restore so a per-call seed doesn't disturb the caller's stream
.Random.seed_safe <- function() if (exists(".Random.seed", envir = .GlobalEnv))
  get(".Random.seed", envir = .GlobalEnv) else NULL
.restore_seed <- function(s) {
  if (is.null(s)) { if (exists(".Random.seed", envir = .GlobalEnv))
    rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", s, envir = .GlobalEnv)
}

#' GPQ tolerance interval from raw data (direct engine, advanced)
#'
#' The Oliva-Aviles & Hauser (2025) generalized-pivotal-
#' quantity tolerance interval for the random-intercept + random-slope stability
#' model, from raw `(y, time, batch)` vectors without a fitted `lmer` model. The
#' fit-based dispatcher [ti_lmm()] routes random-slope designs here automatically
#' (singular fits, and the batch-mean target `over = "group_mean"`); use
#' `ti_gpq_raw()` directly when you have data but no model object.
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
#' @return A data frame, one row per `t0`: `method = "gpq"`, `time`, `estimate`,
#'   `lower`, `upper`, `design = "random_slope_correlated"`, plus audit fields
#'   (`tau2_med`, `s2e_med`, `frac_boundary`, `K`, `df_resid`, `phi1`).
#' @references Oliva-Aviles C, Hauser S (2025). *Technometrics* 67(2):193-202. \doi{10.1080/00401706.2024.2407324}
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
  out <- .gpq_core(y, time, batch, t0, content = P, conf = conf, sides = sides,
                   c_resid = c_resid, M = control$M, seed = control$seed)
  out$design <- "random_slope_correlated"
  .unify_intervals(out, model = NULL, kind = "ti", P_value = P)
}
