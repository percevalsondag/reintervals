## ---------------------------------------------------------------------------
## .montes_core: the Montes-Burdick-Leblond (2019, Appendix A) ANOVA/MLS
## tolerance-interval computation (unweighted lot-mean estimator). Internal
## engine used by both the public raw-data wrapper ti_anova_raw() and the
## dispatch (ti_lmm). Method-pure: computes from raw (y, time, lot), never lme4.
## Carries content/conf internally; the public wrapper does the P/alpha/level
## harmonization. Provenance: method = "anova-mls".
##
## Slope-free branch (time = NULL): the Hoffman & Kringle (2005, Eq. 3;
## restated in Hoffman & Kringle 2007, Eqs. 1-2) one-way random-effects
## two-sided tolerance interval. It is the same estimator with the fixed slope
## removed: beta_hat is dropped, the residual df is N - I (not N - I - 1), the
## centre is the unweighted mean of the group means, and VhatYh = S2_L / I.
## Everything else is unchanged, with the identities (HK notation on the right)
##   S2_L = MSBU / J_H,   Vhat_Y = sTOT2,   n_E = N_e,   C1 = H1,   C2 = H2,
## where MSBU = J_H * sum((ybar_i - ybar_star)^2) / (I - 1) is HK's
## unweighted-means between mean square and J_H = I / sum(1 / J_i).
## Provenance there: design = "random_intercept".
##
## New-group-mean quantities (attributes ngm_W, ngm_df; consumed by
## new_group_mean_lmm(), never by the TI). Target: the TRUE mean of a new lot at
## t0, theta_new = mu + beta t0 + a_new. With D = sum_i (tbar_i - tbar_star)^2,
##   Var(theta_new - Yhat) = s2a (1 + 1/I) + s2e / (I J_H)
##                           + s2e (t0 - tbar_star)^2 / S_tt,
##   E(S2_L) = s2a + s2e / J_H + s2e D / ((I - 1) S_tt),   E(S2_E) = s2e,
## so W = S2_L (1 + 1/I) + S2_E (t0 - tbar_star)^2 / S_tt has
## E(W) - Var(theta_new - Yhat) = s2e / J_H + (1 + 1/I) s2e D / ((I - 1) S_tt)
## >= 0.
## S2_L and S2_E are independent; df = Satterthwaite on the two positive terms
## (df I - 1 and r). Slope-free branch: W = S2_L (1 + 1/I) with df I - 1, which
## on balanced data is the classical t_{I-1} sqrt(MSA (1 + 1/I) / n) interval.
## Validated by simulation (the paper repo's defect study, E5a: no cell below
## 0.95 - 2 SE on balanced and single-observation designs).
## ---------------------------------------------------------------------------
.montes_core <- function(y, time, lot, t0, content = 0.99, conf = 0.95,
                         sides = c("two.sided", "lower", "upper")) {
  sides <- match.arg(sides)
  for (nm in c("content", "conf")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0 || v >= 1) {
      stop(sprintf("`%s` must be a single number strictly between 0 and 1.", nm),
           call. = FALSE)
    }
  }
  has_time <- !is.null(time)
  if (has_time) {
    if (is.null(t0) || !is.numeric(t0) || !length(t0) || any(!is.finite(t0))) {
      stop("`t0` must be a non-empty numeric vector of finite time point(s).",
           call. = FALSE)
    }
  } else {
    if (!is.null(t0)) {
      stop("`t0` applies only with a `time` covariate; with `time = NULL` (the ",
           "Hoffman-Kringle one-way interval) leave `t0` as NULL.", call. = FALSE)
    }
    if (sides != "two.sided") {
      stop("The Hoffman-Kringle (2005) one-way tolerance interval (time = NULL) ",
           "is two-sided only.", call. = FALSE)
    }
  }

  lot <- factor(lot)
  ok  <- !is.na(y) & !is.na(lot) & (if (has_time) !is.na(time) else TRUE)
  y <- y[ok]; lot <- droplevels(lot[ok]); if (has_time) time <- time[ok]
  I <- nlevels(lot); N <- length(y)
  if (I < 2L) stop(sprintf("ti_anova_raw(): the ANOVA tolerance interval needs >= 2 lots (grouping levels) to estimate the between-lot variance; got %d.", I), call. = FALSE)

  Ji <- as.numeric(table(lot)); li <- as.integer(lot)
  ybar_i <- tapply(y, lot, mean)
  yc <- y - ybar_i[li]; S_yyw <- sum(yc * yc)

  if (has_time) {
    ## within-lot (pooled) regression slope
    tbar_i <- tapply(time, lot, mean)
    tc <- time - tbar_i[li]
    S_tyw <- sum(tc * yc); S_ttw <- sum(tc * tc)
    if (!is.finite(S_ttw) || S_ttw <= 0) {
      stop("ti_anova_raw(): no within-lot time variation; the fixed-slope ANOVA ",
           "route needs >= 2 distinct time points within lots.", call. = FALSE)
    }
    beta_hat <- S_tyw / S_ttw
    tbar_star <- mean(tbar_i)
  } else {
    ## HK (2005) one-way: no slope; the lot-mean adjustment and the slope's
    ## residual df disappear.
    beta_hat <- 0; S_tyw <- 0; S_ttw <- NA_real_
    tbar_i <- rep(0, I); tbar_star <- 0
  }

  ybar_star <- mean(ybar_i)
  JH <- I / sum(1 / Ji)                     # harmonic-mean lot size (HK J_H)
  s  <- I - 1                               # between-lot df
  r  <- if (has_time) N - I - 1 else N - I  # residual df (slope-adjusted or not)
  if (r < 1L) {
    stop("ti_anova_raw(): non-positive residual df; need more within-lot ",
         "replication.", call. = FALSE)
  }

  ## between- and within-lot mean squares (both >= 0 by construction)
  adj  <- ybar_i - ybar_star + beta_hat * (tbar_star - tbar_i)
  S2_L <- sum(adj^2) / (I - 1)              # HK: MSBU / J_H
  S2_E <- (S_yyw - (if (has_time) S_tyw^2 / S_ttw else 0)) / r   # HK: MSE

  Vhat_Y <- S2_L + S2_E * (1 - 1 / JH)      # content variance (HK: sTOT2)
  alpha  <- 1 - conf
  C1 <- s / stats::qchisq(alpha, s) - 1     # HK: H1
  C2 <- r / stats::qchisq(alpha, r) - 1     # HK: H2
  U  <- Vhat_Y + sqrt(C1^2 * S2_L^2 + C2^2 * (1 - 1 / JH)^2 * S2_E^2)

  zc <- if (sides == "two.sided") stats::qnorm((1 + content) / 2) else
    stats::qnorm(content)
  df <- min(s, r)

  t_eval <- if (has_time) t0 else NA_real_
  nE_vec <- numeric(length(t_eval))
  ## new-group-mean variance W and its df (see header)
  c1 <- 1 + 1 / I
  c2 <- if (has_time) (t_eval - tbar_star)^2 / S_ttw else 0
  ngm_W  <- c1 * S2_L + c2 * S2_E
  ngm_df <- if (has_time) ngm_W^2 / ((c1 * S2_L)^2 / s + (c2 * S2_E)^2 / r) else
    rep(s, length(t_eval))
  rows <- lapply(seq_along(t_eval), function(i) {
    tt     <- t_eval[i]
    Yhat   <- if (has_time) ybar_star + beta_hat * (tt - tbar_star) else ybar_star
    VhatYh <- if (has_time) S2_L / I + S2_E * (tt - tbar_star)^2 / S_ttw else S2_L / I
    nE     <- Vhat_Y / VhatYh                # HK: N_e
    nE_vec[i] <<- nE
    hw <- zc * sqrt(1 + 1 / nE) * sqrt(U)
    row <- data.frame(
      type = "TI", estimate = Yhat,
      lower = if (sides == "upper") -Inf else Yhat - hw,
      upper = if (sides == "lower")  Inf else Yhat + hw,
      content = content, conf = conf, df = df, sides = sides,
      design = if (has_time) "random_intercept_fixed_slope" else "random_intercept",
      method = "anova-mls",
      singular = FALSE, boundary = FALSE, stringsAsFactors = FALSE
    )
    if (has_time) cbind(time = tt, row) else row
  })

  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (has_time) attr(out, "beta_hat") <- beta_hat
  attr(out, "n_E")      <- nE_vec
  attr(out, "U")        <- U
  attr(out, "Vhat_Y")   <- Vhat_Y
  attr(out, "C1")       <- C1
  attr(out, "C2")       <- C2
  attr(out, "ngm_W")    <- ngm_W
  attr(out, "ngm_df")   <- ngm_df
  out
}

#' ANOVA/MLS tolerance interval from raw data (direct engine, advanced)
#'
#' Tier-2 direct engine: the published closed-form (ANOVA mean-square / modified
#' large-sample) tolerance interval from raw vectors, without a fitted `lmer`
#' model. Two cases share one estimator (the *unweighted lot-mean* ANOVA):
#'
#' * **With a time covariate** (`time`, `t0` given): the Montes, Burdick &
#'   Leblond (2019, Appendix A) interval for the random-intercept + fixed-slope
#'   stability model. [ti_lmm()] routes single-observation / staggered
#'   fixed-slope designs here automatically.
#' * **Without one** (`time = NULL`, `t0 = NULL`): the Hoffman & Kringle (2005)
#'   two-sided interval for the one-way random-effects model, balanced or
#'   unbalanced (restated in Hoffman & Kringle 2007, Eqs. 1-2). It is the same
#'   computation with the slope removed (residual df `N - I`, centre = the
#'   unweighted mean of the group means). [ti_lmm()] uses it by default for an
#'   intercept-only model with one random-intercept factor. Two-sided only.
#'
#' Method-pure (never touches `lme4`). The between- and within-lot mean squares
#' are sums of squares over their degrees of freedom and so are non-negative by
#' construction; no ANOVA-estimate truncation.
#'
#' @param y Numeric response.
#' @param time Numeric time/age covariate (same length as `y`), or `NULL` for
#'   the one-way (Hoffman-Kringle) interval with no covariate.
#' @param lot Lot/batch (grouping) factor (same length as `y`).
#' @param t0 Time point(s) at which to evaluate the interval; must be `NULL`
#'   (the default) when `time = NULL`.
#' @param P Content proportion (default 0.95; pharma stability typically 0.99).
#' @param alpha,level Confidence as a mutually-exclusive pair (give at most one);
#'   `level = 1 - alpha`, default `alpha = 0.05`.
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`. One-sided
#'   intervals are available only with a time covariate.
#'
#' @return A data frame (one row per `t0`, or one row when `time = NULL`):
#'   `time` (only with a time covariate), `type = "TI"`, `estimate`, `lower`,
#'   `upper`, `content`, `conf`, `df`, `sides`, `design`
#'   (`"random_intercept_fixed_slope"`, or `"random_intercept"` when
#'   `time = NULL`), `method = "anova-mls"`, and `singular`/`boundary` flags.
#'   ANOVA intermediates are attached as attributes: `n_E`, `U`, `Vhat_Y`
#'   (Hoffman-Kringle's `N_e` and `sTOT2` when `time = NULL`), `C1`, `C2`
#'   (their `H1`, `H2`), and `beta_hat` (with a time covariate).
#' @references Montes RO, Burdick RK, Leblond DJ (2019). Simple approach to
#'   calculate random effects model tolerance intervals to set release and
#'   shelf-life specification limits of pharmaceutical products. *PDA Journal of Pharmaceutical Science and Technology*. \doi{10.5731/pdajpst.2018.008839}
#'
#'   Hoffman D, Kringle R (2005). Two-sided tolerance intervals for balanced and
#'   unbalanced random effects models. *Journal of Biopharmaceutical
#'   Statistics* 15(2):283-293. \doi{10.1081/BIP-200048826}
#'
#'   Hoffman D, Kringle R (2007). A total error approach for the validation of
#'   quantitative analytical methods. *Pharmaceutical Research*
#'   24(6):1157-1164. \doi{10.1007/s11095-007-9242-3}
#' @seealso [ti_lmm()] (the dispatcher), [ti_gpq_raw()] (the GPQ raw engine),
#'   [hk_assay], [reintervals-models].
#' @examples
#' lot <- rep(c("A", "B", "C", "D", "E", "F"), c(6, 5, 2, 5, 2, 1))
#' ti_anova_raw(
#'   y    = c(70, 71, 69, 68, 66, 65, 71, 69, 68, 66, 66, 71, 69,
#'            73, 73, 71, 71, 70, 72, 73, 70),
#'   time = c(0, 1, 3, 6, 9, 12, 0, 1, 3, 6, 9, 0, 3, 0, 1, 3, 6, 9, 0, 1, 0),
#'   lot  = lot, t0 = 12, P = 0.99
#' )
#' # one-way, no covariate: Hoffman & Kringle (2005), here (0.846, 1.114)
#' ti_anova_raw(hk_assay$conc, time = NULL, lot = hk_assay$run,
#'              P = 0.95, level = 0.90)
#' @export
ti_anova_raw <- function(y, time, lot, t0 = NULL, P = 0.95, alpha = NULL,
                         level = NULL, sides = c("two.sided", "lower", "upper")) {
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal `conf`; public `P` -> `content`.
  conf  <- 1 - .resolve_alpha(alpha, level)
  native <- .montes_core(y, time, lot, t0, content = P, conf = conf, sides = sides)
  .unify_intervals(native, model = NULL, kind = "ti", P_value = P)
}
