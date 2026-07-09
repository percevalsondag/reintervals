<!-- README.md is generated from README.Rmd. Please edit that file. -->



# reintervals

> Confidence, prediction, tolerance, and new-group-mean intervals for linear
> mixed models.

`reintervals` computes confidence (CI), prediction (PI), tolerance (TI), and
new-group-mean (a.k.a. CInew / new-batch-mean) intervals for linear mixed models
fitted with [lme4](https://cran.r-project.org/package=lme4), through one uniform
interface. It is a general mixed-model interval package, with particular support
for pharmaceutical stability / shelf-life and bioassay work.

## Why this package

Several tools already do part of the job well:

- **`predict.merMod()`** (lme4) and **`emmeans`** / **`ggeffects`** give point
  predictions and confidence intervals for the mean / marginal effects.
- **`merTools::predictInterval()`** gives simulation-based prediction intervals
  for `merMod` fits.
- **The `tolerance` package** gives tolerance intervals for many distributions
  and some regression settings (for simpler designs, without the mixed-model
  degrees-of-freedom treatment).

What none of them provides is **tolerance and new-group-mean (CInew) intervals
across both random-intercept and random-slope / fixed-slope (stability) designs**
— with the observed-information degrees of freedom of Francq, Lin & Hoyer
(2019), the unweighted-lot-mean ANOVA closed form of Montes, Burdick & Leblond
(2019), and the generalized-pivotal-quantity method of Oliva-Aviles & Hauser
(2025) — **behind one uniform interface.**

`reintervals` provides four interval types from four engines, selected
automatically for the fitted model (with a manual override), and returns them in
one tidy `lmm_interval` data frame. Every interval the package ships is validated
against its published source on continuous integration. **The methods are
published and cited; the contribution here is the validated, unified
implementation, not new methodology.**

## Installation

```r
# pre-CRAN: install from GitHub
# install.packages("remotes")
remotes::install_github("percevalsondag/reintervals")
```

(A CRAN release will be noted here once available.)

## Quick start


``` r
library(reintervals)

fm <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)

# 95% confidence interval for the mean at Days = 5
ci_lmm(fm, newdata = data.frame(Days = 5))
#> <lmm_interval> ci  (method = ems-mls, design = random_intercept_fixed_slope)
#>  Days type estimate    lower    upper level  P     sides singular boundary vg_floored
#>     5   ci 303.7415 284.6347 322.8483  0.95 NA two.sided    FALSE       NA         NA
#> diagnostics: df  (see attr(., "diagnostics"))

# 99%-content / 95%-confidence tolerance interval at Days = 5
ti_lmm(fm, newdata = data.frame(Days = 5), P = 0.99)
#> <lmm_interval> ti  (method = ems-mls, design = random_intercept_fixed_slope)
#>  Days type estimate    lower    upper level    P     sides singular boundary vg_floored
#>     5   ti 303.7415 142.9115 464.5716  0.95 0.99 two.sided    FALSE       NA         NA
```

Every verb returns the same `lmm_interval` data frame, so confidence,
prediction, tolerance, and new-group-mean results all share one column contract
(see `?ci_lmm` for the schema).

## The four verbs

| Verb | Computes |
|---|---|
| `ci_lmm()` | confidence interval for the mean |
| `pi_lmm()` | prediction interval for a future observation |
| `ti_lmm()` | tolerance interval (content `P`); `over = "group_mean"` for the batch-mean target |
| `new_group_mean_lmm()` | the CInew / new-batch-mean interval (prediction for one new group's mean) |

`lmm_predict(model, interval = ...)` is a `predict()`-style wrapper over the
four; `lmm_interval_control()` supplies the GPQ tuning (`seed`, `M`).

## Engine and scope map

The engine is chosen from the fitted model's structure (full mapping in
`?reintervals-models`):

| Model class | Engine (`method`) |
|---|---|
| random intercept (any fixed effects; one-way / nested / crossed) | expected-mean-square / MLS (`"ems-mls"`) |
| random intercept + fixed slope, clean balanced | `"ems-mls"` |
| random intercept + fixed slope, single-observation / unbalanced | unweighted lot-mean ANOVA (`"anova-mls"`) |
| random slope, future-observation | REML closed form (`"reml-mls"`); GPQ when singular |
| random slope, batch-mean | GPQ Monte Carlo (`"gpq"`) |

**v1 scope:** general fixed effects are supported for random-intercept models;
random-slope models are supported for a **single continuous covariate** (the stability
design). The tier-2 engines `ti_anova_raw()` and `ti_gpq_raw()` compute from raw
`(y, time, group)` vectors without a fitted model.

## Validation

Every shipped interval is checked against a published or external-source oracle
on CI: the Francq, Lin & Hoyer (2019) worked numbers; the Montes-Burdick-Leblond
(2019) Appendix A tolerance interval `(57.25, 75.41)`; the random-slope interval
arithmetic against committed reference values to ~`1e-9`; and the GPQ content
variance against the Oliva-Aviles & Hauser (2025) closed form (`tau^2 = 8.439`).

## Scope and roadmap

v1 covers general-fixed-effect random-intercept models and single-time-covariate
random-slope (stability) models. The planned v2 extension is multi-covariate
random-slope intervals; see `NEWS.md` for the change log.

## License and citation

MIT © Perceval Sondag. A JOSS DOI will be added here on acceptance; until then,
please cite the package (`citation("reintervals")`) and the underlying methods:

- Francq, B. G., Lin, D., & Hoyer, W. (2019). Confidence, prediction, and
  tolerance in linear mixed models. *Statistics in Medicine*, 38(30), 5603–5622.
  <doi:10.1002/sim.8386>
- Montes, R. O., Burdick, R. K., & Leblond, D. J. (2019). Simple approach to
  calculate random effects model tolerance intervals to set release and
  shelf-life specification limits of pharmaceutical products. *PDA Journal of
  Pharmaceutical Science and Technology*. <doi:10.5731/pdajpst.2018.008839>
- Oliva-Aviles, C., & Hauser, P. (2025). *Technometrics*, 67(2), 193–202.
  <doi:10.1080/00401706.2024.2407324>
