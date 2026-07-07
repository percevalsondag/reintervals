## data-raw/make_fixtures.R -- regenerate the committed package datasets.
## Run with: source("data-raw/make_fixtures.R")  (writes to data/).

## --- Montes-Burdick-Leblond (2019) Appendix A: the M1 fixed-slope oracle ---
## 26 rows, 11 lots; mixes multi-timepoint stability lots (A-D) with release-only
## single-observation lots (E-K). Two-sided P=0.99/gamma=0.95 TI at t0=12
## reproduces the published (57.25, 75.41).
mbl_appendix_a <- data.frame(
  Lot   = c("A","A","A","A","A","A","B","B","B","B","B","C","C",
            "D","D","D","D","D","E","E","F","G","H","I","J","K"),
  Month = c(0,1,3,6,9,12, 0,1,3,6,9, 0,3, 0,1,3,6,9, 0,1, 0, 0, 0, 0, 0, 0),
  Y     = c(70,71,69,68,66,65, 71,69,68,66,66, 71,69,
            73,73,71,71,70, 72,73, 70, 69, 75, 75, 72, 72),
  stringsAsFactors = FALSE
)

## --- Oliva-Aviles & Hauser (2025) Section 5 DGP (the GPQ reference) ---
## The published simulation specification for the random-intercept + random-slope
## stability model y_kj = (eta1 + a_{k,1}) + (eta2 + a_{k,2}) t_kj + e_kj, with
## kappa = 1. The content variance the GPQ method targets is tau^2(t0) = z' V z,
## z = (1, t0); at t0 = 20 this equals 8.439388 (the paper's 8.439). Stored as a
## specification list (the paper reports a DGP/coverage study, not a single raw
## dataset), so the closed-form tau^2 is exactly reproducible and datasets can be
## simulated from it deterministically.
.s1sq <- 1.5; .s2sq <- 0.01; .rho <- 0.6
.s01  <- .rho * sqrt(.s1sq) * sqrt(.s2sq)
.z20  <- c(1, 20)
.Vk1  <- matrix(c(.s1sq, .s01, .s01, .s2sq), 2, 2)   # V at kappa = 1
oliva_hauser_dgp <- list(
  beta0 = 100, beta1 = -0.03,
  s2_0 = .s1sq, s01 = .s01, s2_1 = .s2sq, s2e = 1, kappa = 1,
  times = c(0, 3, 6, 9, 12, 18),
  t0 = 20,
  tau2_at_t0 = as.numeric(t(.z20) %*% .Vk1 %*% .z20),  # 8.439388 (paper: 8.439)
  citation = "Oliva-Aviles & Hauser (2025), Technometrics 67(2):193-202, Section 5"
)

save(mbl_appendix_a,   file = "data/mbl_appendix_a.rda",   version = 2)
save(oliva_hauser_dgp, file = "data/oliva_hauser_dgp.rda", version = 2)
cat(sprintf("wrote data/: mbl_appendix_a (%d rows), oliva_hauser_dgp (tau2=%.6f)\n",
            nrow(mbl_appendix_a), oliva_hauser_dgp$tau2_at_t0))
