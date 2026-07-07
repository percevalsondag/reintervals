# Shared hand-entered fixtures from Francq, Lin & Hoyer (2019).
# The paper supplies the inputs (printed components, df) AND the outputs (printed
# CI/PI/TI bounds) independently, so feeding these to a constructor is a genuine,
# non-circular reproduction of the formula under test.

# --- Section 5.2: balanced one-way, 6 runs x 3 reps ------------------------
# EXACT printed values:
#   run component  = 0.000681,  residual = 0.001253  ->  sigma^2_T = 0.001934
#   df_PI          = 12.484
#   intercept      = 0.981 (SE 0.01353)  ->  var_mean = 0.01353^2,  CI df = 5
# Balanced one-way EMS (n0 = 3 reps, A = 6 runs, N = 18):
#   EMS_A = s2e + n0*s2a   (df r_A = A - 1 = 5,  k_A = 1/n0     = 1/3)
#   EMS_e = s2e            (df r_e = N - A = 12, k_e = 1 - 1/n0 = 2/3)
# so sum(k * EMS) = s2a + s2e = sigma^2_T, the standard identity.
# Deterministic balanced one-way data realizing the Section 5.2 components
# EXACTLY (run = 0.000681, residual = 0.001253; 6 runs x 3 reps). Because it is
# deterministic, the REML fit -> numerical Hessian reproduces 
# df_PI = 12.486 every time. Used by test-vc-extract.R.
# Not circular for df_PI: the construction realizes the (paper-reported)
# components; df_PI is an independent output of the Hessian pipeline.
build_52_data <- function() {
  s2run <- 0.000681
  s2res <- 0.001253
  n <- 3
  A <- 6
  d <- sqrt(s2res)                          # within {-d,0,d} has sample var = s2res
  varmean_runs <- (s2res + n * s2run) / n   # = MSB/3 = sample var of run means
  base <- c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)
  cc <- sqrt(varmean_runs / stats::var(base))
  runmean <- 0.981 + cc * base
  y <- as.vector(sapply(runmean, function(m) m + c(-d, 0, d)))
  data.frame(y = y,
             run = factor(rep(seq_len(A), each = n)),
             rep = rep(seq_len(n), A))
}

fit_52 <- function() {
  lme4::lmer(y ~ 1 + (1 | run), data = build_52_data(), REML = TRUE)
}

fixture_52 <- function() {
  s2a <- 0.000681
  s2e <- 0.001253
  n0  <- 3
  re_components(
    components = c(run = s2a, Residual = s2e),
    dfs        = c(ci = 5, pi = 12.484),
    mean       = 0.981,
    var_mean   = 0.01353^2,
    ems = list(
      ms = c(run = s2e + n0 * s2a, Residual = s2e),
      k  = c(run = 1 / n0,         Residual = 1 - 1 / n0),
      df = c(run = 5,              Residual = 12)
    ),
    target    = "observable",
    estimator = "reml",
    design    = "oneway"
  )
}
