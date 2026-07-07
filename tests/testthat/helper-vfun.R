# =============================================================================
# TEST-ONLY helpers for the slope-extraction tests. None of this is package API.
#   - vfun_sim()       : base-R stability-data generator (M1 / M2c / M2i).
#   - mbl_appendix_a() : the Montes-Burdick-Leblond (2019) Appendix A dataset.
#   - montes_ref()     : used here ONLY
#                        as the external ANOVA/MLS oracle to validate that the
#                        slope extraction reproduces a known reference (not just
#                        itself). It is NOT exported by the package.
# =============================================================================

## base-R random-slope/intercept generator 
## (lives only in the test scope).
vfun_sim <- function(B = 8, times = c(0, 3, 6, 9, 12, 18, 24),
                     beta0 = 100, beta1 = -0.25,
                     s2_0 = 0.5, s2_1 = 0.000625, s01 = 0, s2e = 0.5,
                     reps = 1L, seed = 1L) {
  set.seed(seed)
  time_list <- if (is.list(times)) times else rep(list(times), B)
  Sigma <- matrix(c(s2_0, s01, s01, s2_1), 2, 2)
  b <- if (s2_1 > 0 || s01 != 0) {
    L <- chol(Sigma + diag(1e-12, 2)); t(L) %*% matrix(stats::rnorm(2 * B), 2, B)
  } else {
    rbind(stats::rnorm(B, 0, sqrt(s2_0)), rep(0, B))
  }
  recs <- list(); ri <- 0L
  for (i in seq_len(B)) for (tt in time_list[[i]]) for (r in seq_len(reps)) {
    ri <- ri + 1L
    mu <- (beta0 + b[1, i]) + (beta1 + b[2, i]) * tt
    recs[[ri]] <- data.frame(batch = i, time = tt,
                             y = mu + stats::rnorm(1, 0, sqrt(s2e)))
  }
  d <- do.call(rbind, recs); d$batch <- factor(d$batch); d
}

## Montes-Burdick-Leblond (2019) Appendix A is now the committed package dataset
## `mbl_appendix_a` (see R/data.R); it is lazy-loaded (LazyData) and referenced by
## name in the tests. The previous test-helper that held the data literal is
## retired so the dataset has a single source of truth (data/mbl_appendix_a.rda).

montes_ref <- function(y, time, lot, t0, content = 0.99, conf = 0.95) {
  lot <- factor(lot)
  ok <- !is.na(y) & !is.na(time) & !is.na(lot)
  y <- y[ok]; time <- time[ok]; lot <- droplevels(lot[ok])
  I <- nlevels(lot); N <- length(y); Ji <- as.numeric(table(lot))
  li <- as.integer(lot)
  ybar_i <- tapply(y, lot, mean); tbar_i <- tapply(time, lot, mean)
  yc <- y - ybar_i[li]; tc <- time - tbar_i[li]
  S_tyw <- sum(tc * yc); S_ttw <- sum(tc * tc); S_yyw <- sum(yc * yc)
  beta_hat <- S_tyw / S_ttw
  ybar_star <- mean(ybar_i); tbar_star <- mean(tbar_i)
  JH <- I / sum(1 / Ji); s <- I - 1; r <- N - I - 1
  adj <- ybar_i - ybar_star + beta_hat * (tbar_star - tbar_i)
  S2_L <- sum(adj^2) / (I - 1)
  S2_E <- (S_yyw - S_tyw^2 / S_ttw) / r
  Vhat_Y <- S2_L + S2_E * (1 - 1 / JH)
  alpha <- 1 - conf
  C1 <- s / stats::qchisq(alpha, s) - 1
  C2 <- r / stats::qchisq(alpha, r) - 1
  U  <- Vhat_Y + sqrt(C1^2 * S2_L^2 + C2^2 * (1 - 1 / JH)^2 * S2_E^2)
  zP2 <- stats::qnorm((1 + content) / 2)
  Yhat   <- ybar_star + beta_hat * (t0 - tbar_star)
  VhatYh <- S2_L / I + S2_E * (t0 - tbar_star)^2 / S_ttw
  nE     <- Vhat_Y / VhatYh
  hw <- zP2 * sqrt(1 + 1 / nE) * sqrt(U)
  list(estimate = Yhat, lower = Yhat - hw, upper = Yhat + hw,
       beta_hat = beta_hat, Vhat_Y = Vhat_Y, Vhat_Yhat = VhatYh,
       n_E = nE, U = U)
}
