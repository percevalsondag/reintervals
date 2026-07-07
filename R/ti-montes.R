## ---------------------------------------------------------------------------
## .montes_core: the Montes-Burdick-Leblond (2019, Appendix A) ANOVA/MLS
## tolerance-interval computation (unweighted lot-mean estimator). Internal
## engine used by both the public raw-data wrapper ti_anova_raw() and the fixed-slope
## dispatch (ti_lmm). Method-pure: computes from raw (y, time, lot), never lme4.
## Carries content/conf internally; the public wrapper does the P/alpha/level
## harmonization. Provenance: design = "random_intercept_fixed_slope",
## method = "anova-mls".
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
  if (!is.numeric(t0) || !length(t0) || any(!is.finite(t0))) {
    stop("`t0` must be a non-empty numeric vector of finite time point(s).",
         call. = FALSE)
  }

  lot <- factor(lot)
  ok  <- !is.na(y) & !is.na(time) & !is.na(lot)
  y <- y[ok]; time <- time[ok]; lot <- droplevels(lot[ok])
  I <- nlevels(lot); N <- length(y)
  if (I < 2L) stop(sprintf("ti_anova_raw(): the ANOVA tolerance interval needs >= 2 lots (grouping levels) to estimate the between-lot variance; got %d.", I), call. = FALSE)

  Ji <- as.numeric(table(lot)); li <- as.integer(lot)
  ybar_i <- tapply(y, lot, mean); tbar_i <- tapply(time, lot, mean)

  ## within-lot (pooled) regression slope
  yc <- y - ybar_i[li]; tc <- time - tbar_i[li]
  S_tyw <- sum(tc * yc); S_ttw <- sum(tc * tc); S_yyw <- sum(yc * yc)
  if (!is.finite(S_ttw) || S_ttw <= 0) {
    stop("ti_anova_raw(): no within-lot time variation; the fixed-slope ANOVA ",
         "route needs >= 2 distinct time points within lots.", call. = FALSE)
  }
  beta_hat <- S_tyw / S_ttw

  ybar_star <- mean(ybar_i); tbar_star <- mean(tbar_i)
  JH <- I / sum(1 / Ji)                     # harmonic-mean lot size
  s  <- I - 1                               # between-lot df
  r  <- N - I - 1                           # residual df (slope-adjusted)
  if (r < 1L) {
    stop("ti_anova_raw(): non-positive residual df; need more within-lot ",
         "replication.", call. = FALSE)
  }

  ## between- and within-lot mean squares (both >= 0 by construction)
  adj  <- ybar_i - ybar_star + beta_hat * (tbar_star - tbar_i)
  S2_L <- sum(adj^2) / (I - 1)
  S2_E <- (S_yyw - S_tyw^2 / S_ttw) / r

  Vhat_Y <- S2_L + S2_E * (1 - 1 / JH)      # content variance
  alpha  <- 1 - conf
  C1 <- s / stats::qchisq(alpha, s) - 1
  C2 <- r / stats::qchisq(alpha, r) - 1
  U  <- Vhat_Y + sqrt(C1^2 * S2_L^2 + C2^2 * (1 - 1 / JH)^2 * S2_E^2)

  zc <- if (sides == "two.sided") stats::qnorm((1 + content) / 2) else
    stats::qnorm(content)
  df <- min(s, r)

  nE_vec <- numeric(length(t0))
  rows <- lapply(seq_along(t0), function(i) {
    tt     <- t0[i]
    Yhat   <- ybar_star + beta_hat * (tt - tbar_star)
    VhatYh <- S2_L / I + S2_E * (tt - tbar_star)^2 / S_ttw
    nE     <- Vhat_Y / VhatYh
    nE_vec[i] <<- nE
    hw <- zc * sqrt(1 + 1 / nE) * sqrt(U)
    data.frame(
      time = tt, type = "TI", estimate = Yhat,
      lower = if (sides == "upper") -Inf else Yhat - hw,
      upper = if (sides == "lower")  Inf else Yhat + hw,
      content = content, conf = conf, df = df, sides = sides,
      design = "random_intercept_fixed_slope", method = "anova-mls",
      singular = FALSE, boundary = FALSE, stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "beta_hat") <- beta_hat
  attr(out, "n_E")      <- nE_vec
  attr(out, "U")        <- U
  attr(out, "Vhat_Y")   <- Vhat_Y
  out
}

#' Montes ANOVA/MLS tolerance interval from raw data (direct engine, advanced)
#'
#' Tier-2 direct engine: the published closed-form (ANOVA mean-square / modified
#' large-sample) tolerance interval of Montes, Burdick & Leblond (2019, Appendix
#' A) for the random-intercept + fixed-slope stability model --- the
#' *unweighted lot-mean* estimator --- from raw `(y, time, lot)` vectors, without
#' a fitted `lmer` model. The fit-based dispatcher [ti_lmm()] routes single-
#' observation / staggered fixed-slope designs to this engine automatically; use
#' `ti_anova_raw()` directly when you have data but no model object.
#'
#' Method-pure (never touches `lme4`). The between- and within-lot mean squares
#' are sums of squares over their degrees of freedom and so are non-negative by
#' construction; no ANOVA-estimate truncation arises.
#'
#' @param y Numeric response.
#' @param time Numeric time/age covariate (same length as `y`).
#' @param lot Lot/batch factor (same length as `y`).
#' @param t0 Time point(s) at which to evaluate the interval.
#' @param P Content proportion (default 0.95; pharma stability typically 0.99).
#' @param alpha,level Confidence as a mutually-exclusive pair (give at most one);
#'   `level = 1 - alpha`, default `alpha = 0.05`.
#' @param sides `"two.sided"` (default), `"lower"`, or `"upper"`.
#'
#' @return A data frame, one row per `t0`: `time`, `type = "TI"`, `estimate`,
#'   `lower`, `upper`, `content`, `conf`, `df`, `sides`,
#'   `design = "random_intercept_fixed_slope"`, `method = "anova-mls"`, and
#'   `singular`/`boundary` flags. ANOVA intermediates `beta_hat`, `n_E`, `U`,
#'   `Vhat_Y` are attached as attributes.
#' @references Montes RO, Burdick RK, Leblond DJ (2019). Simple approach to
#'   calculate random effects model tolerance intervals to set release and
#'   shelf-life specification limits of pharmaceutical products. *PDA Journal of Pharmaceutical Science and Technology*. \doi{10.5731/pdajpst.2018.008839}
#' @seealso [ti_lmm()] (the dispatcher), [ti_gpq_raw()] (the GPQ raw engine),
#'   [reintervals-models].
#' @examples
#' lot <- rep(c("A", "B", "C", "D", "E", "F"), c(6, 5, 2, 5, 2, 1))
#' ti_anova_raw(
#'   y    = c(70, 71, 69, 68, 66, 65, 71, 69, 68, 66, 66, 71, 69,
#'            73, 73, 71, 71, 70, 72, 73, 70),
#'   time = c(0, 1, 3, 6, 9, 12, 0, 1, 3, 6, 9, 0, 3, 0, 1, 3, 6, 9, 0, 1, 0),
#'   lot  = lot, t0 = 12, P = 0.99
#' )
#' @export
ti_anova_raw <- function(y, time, lot, t0, P = 0.95, alpha = NULL, level = NULL,
                         sides = c("two.sided", "lower", "upper")) {
  sides <- match.arg(sides)
  ## boundary: public alpha/level -> internal `conf`; public `P` -> `content`.
  conf  <- 1 - .resolve_alpha(alpha, level)
  native <- .montes_core(y, time, lot, t0, content = P, conf = conf, sides = sides)
  .unify_intervals(native, model = NULL, kind = "ti", P_value = P)
}
