# GPQ tolerance-interval coverage study — result

**Reproduce:** `Rscript data-raw/gpq_coverage_study.R` → `data-raw/gpq_coverage_results.csv`.
Seeded; `N_sim = 3000` datasets per condition; parallel over cores.

## What was measured

Coverage of the **GPQ batch-mean `(P, gamma)`-tolerance interval** (`c_resid = 0`)
for the committed, published **Oliva-Aviles & Hauser (2025) Section-5 DGP**
(`oliva_hauser_dgp`), evaluated at its `t0 = 20`.

* **Population at t0** (batch true values): `theta(t0) ~ N(mu0, tau2)`,
  `mu0 = beta0 + beta1*t0 = 99.4`, `tau2 = z0' Sigma z0`, `z0 = (1, t0)`.
* **Coverage criterion (TI coverage, NOT interval-contains-a-point):** a dataset
  counts as covered iff the true content its interval `[L, U]` captures,
  `content = Phi((U-mu0)/sigma) - Phi((L-mu0)/sigma)`, `sigma = sqrt(tau2)`,
  is `>= P`. Coverage = fraction of `N_sim` datasets covered. **Nominal = gamma.**
* **Lower-edge sweep:** `Sigma` scaled by `scale ∈ {0.1, 0.25, 0.5, 1, 2}`
  (residual variance `s2e = 1` fixed) → `tau2(t0) ∈ {0.84, 2.11, 4.22, 8.44, 16.88}`.
  Small `tau2` = the "small-tau^2" lower edge where a ~0.934 coverage dip had
  previously been reported.
* **M-artifact test:** every condition run at **M = 4000** (low/original resolution)
  and **M = 20000** (paper resolution).
* **Fixed:** `K = 15` batches, `n = 6` time points/batch, `P = 0.95`,
  `gamma = 0.95`, two-sided. (A K-scan in 8..50 at `tau2 = 8.44` found coverage
  ~nominal across K, so K is not the dip driver; `K = 15` is a representative
  stability-study size.)

## Coverage table (nominal = 0.95; binomial SE = sqrt(p(1-p)/3000))

| tau2 | M = 4000 | SE | M = 20000 | SE |
|---:|---:|---:|---:|---:|
| 0.84 | 0.9393 | 0.0044 | **0.9443** | 0.0042 |
| 2.11 | 0.9473 | 0.0041 | 0.9537 | 0.0038 |
| 4.22 | 0.9503 | 0.0040 | 0.9513 | 0.0039 |
| 8.44 | 0.9427 | 0.0042 | 0.9550 | 0.0038 |
| 16.88 | 0.9423 | 0.0043 | 0.9460 | 0.0041 |

## Finding (honest)

The original **~0.934 dip does not persist at the paper resolution — it was
largely a low-Monte-Carlo-draw (M ≈ 4000) artifact.** At **M = 20000**, coverage
is **at or near nominal across the whole tested variance-ratio range
(0.944–0.955 at gamma = 0.95)**: every condition is within **~1.3 binomial SE of
0.95**, i.e. not statistically distinguishable from nominal. The slightest
residual softness (**0.944**, ~1.3 SE below nominal) is at the **smallest variance
ratio (tau2 ≈ 0.84 with residual variance 1)**; the interior (tau2 ≈ 2–8) is
nominal or slightly above (0.951–0.955).

At low M (4000) the same conditions read ~0.94 (down to 0.939 at the small-tau2
edge), reproducing the order of the prior dip — confirming the dip is a
finite-M phenomenon that the paper's M resolves.

## One-line characterization for the paper / vignette (matched to the number)

> At the paper's Monte-Carlo resolution (M = 20000), the GPQ batch-mean
> (P, gamma)-tolerance interval attains its nominal confidence across the tested
> variance-ratio range — empirical coverage 0.944–0.955 at gamma = 0.95
> (N_sim = 3000), with no condition more than ~1.3 binomial standard errors below
> nominal. The sub-0.94 coverage seen at low Monte-Carlo resolution (M ≈ 4000) is
> an artifact of too few GPQ draws and resolves at high M; the slightest softness
> (~0.944) is at the smallest variance ratio (tau^2 ≈ 0.8).
