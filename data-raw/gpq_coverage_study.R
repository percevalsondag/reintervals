# =============================================================================
# GPQ tolerance-interval coverage study
#
# PURPOSE. Measure the coverage of the GPQ (P, gamma) batch-mean tolerance
# interval before any coverage language is written in the paper/vignette. An
# earlier, lower-resolution look reported the GPQ TI coverage dipping to ~0.934
# in a lower-edge (small-tau^2) corner; this is a higher-resolution
# re-measurement to see whether that dip is real or a low-Monte-Carlo-draw
# (low-M) artifact.
#
# DESIGN CHOICES (stated explicitly):
#  * DGP: the committed, published Oliva-Aviles & Hauser (2025) Section-5 spec
#    `oliva_hauser_dgp` (random-slope stability model), evaluated at its t0 = 20.
#  * Target: the BATCH-MEAN tolerance interval (c_resid = 0) -- the case the GPQ
#    engine documents as covering ~gamma, and the claim the paper will make.
#    The population at t0 is the batch true-value distribution
#    theta(t0) ~ N(mu0, tau2), mu0 = beta0 + beta1*t0 = 99.4,
#    tau2 = z0' Sigma z0 with z0 = (1, t0).
#  * COVERAGE CRITERION (TI coverage, NOT interval-contains-a-point): for each
#    simulated dataset the interval [L, U] is computed, and the dataset COUNTS
#    as covered iff the TRUE content it captures,
#       content = Phi((U-mu0)/sigma) - Phi((L-mu0)/sigma),  sigma = sqrt(tau2),
#    is >= P. Coverage = fraction of datasets covered; nominal target = gamma.
#  * LOWER-EDGE SWEEP: the original grid is not committed, so we sweep the
#    variance ratio by scaling Sigma by `scale` in {0.1, 0.25, 0.5, 1, 2}
#    (residual s2e = 1 fixed), giving tau2(t0) in {0.84, 2.11, 4.22, 8.44, 16.88}.
#    Small scale = the "small tau^2" lower edge; scale >= 1 = interior.
#  * M ARTIFACT TEST: every condition is run at M = 4000 (the original low-M
#    resolution) AND M = 20000 (the paper resolution), so the table shows
#    whether any dip persists at high M or resolves.
#  * Fixed: K = 15 batches (a representative stability-study batch count; the
#    DGP does not fix K -- the K-scan during exploration showed coverage is
#    ~nominal across K in 8..50 at tau2 = 8.44, so K is not the dip driver),
#    n = 6 time points per batch (the DGP `times`), P = 0.95, gamma = 0.95,
#    two-sided.
#  * REPLICATION: N_sim = 3000 datasets per condition. Binomial SE of a
#    coverage p is sqrt(p(1-p)/N_sim); at p ~ 0.934, N_sim = 3000 gives
#    SE ~ 0.0045, enough to separate 0.934 from 0.95 (~3.5 SE).
#  * Seeded: condition c, replicate i -> set.seed(1000*c + i); GPQ control seed
#    = the same, so the whole study is reproducible.
#
# Slow (~minutes, parallelized). Run: Rscript data-raw/gpq_coverage_study.R
# Writes data-raw/gpq_coverage_results.csv.
# =============================================================================
suppressMessages(pkgload::load_all(quiet = TRUE))
library(parallel)

d   <- reintervals::oliva_hauser_dgp
Sig0 <- matrix(c(d$s2_0, d$s01, d$s01, d$s2_1), 2, 2)
z0  <- c(1, d$t0)
mu0 <- d$beta0 + d$beta1 * d$t0
ncores <- max(1L, detectCores() - 1L)

N_sim  <- 3000L
scales <- c(0.1, 0.25, 0.5, 1, 2)
Ms     <- c(4000L, 20000L)
K      <- 15L
P      <- 0.95
gamma  <- 0.95

cover_one <- function(seed, scale, M) {
  set.seed(seed)
  Sig <- scale * Sig0; Lc <- chol(Sig)
  tau2 <- as.numeric(t(z0) %*% Sig %*% z0); sig <- sqrt(tau2)
  x <- do.call(rbind, lapply(seq_len(K), function(k) {
    a <- as.numeric(crossprod(Lc, rnorm(2)))
    data.frame(batch = k, time = d$times,
               y = (d$beta0 + a[1]) + (d$beta1 + a[2]) * d$times +
                   rnorm(length(d$times), 0, sqrt(d$s2e)))
  }))
  g <- ti_gpq_raw(x$y, x$time, x$batch, t0 = d$t0, P = P, level = gamma,
                  c_resid = 0, control = lmm_interval_control(seed = seed, M = M))
  content <- pnorm(g$upper, mu0, sig) - pnorm(g$lower, mu0, sig)
  as.numeric(content >= P)
}

rows <- list(); ci <- 0L
for (M in Ms) for (s in scales) {
  ci <- ci + 1L
  base <- 1000L * ci
  cov <- unlist(mclapply(seq_len(N_sim), function(i) cover_one(base + i, s, M),
                         mc.cores = ncores))
  p <- mean(cov); se <- sqrt(p * (1 - p) / N_sim)
  tau2 <- as.numeric(t(z0) %*% (s * Sig0) %*% z0)
  rows[[length(rows) + 1L]] <- data.frame(
    M = M, scale = s, tau2 = round(tau2, 3), K = K, n_per = length(d$times),
    P = P, gamma = gamma, N_sim = N_sim,
    coverage = round(p, 4), se = round(se, 4))
  cat(sprintf("M=%5d scale=%.2f tau2=%6.3f  coverage=%.4f  se=%.4f\n",
              M, s, tau2, p, se))
}
res <- do.call(rbind, rows)
write.csv(res, "data-raw/gpq_coverage_results.csv", row.names = FALSE)
cat("\nWROTE data-raw/gpq_coverage_results.csv\n")
