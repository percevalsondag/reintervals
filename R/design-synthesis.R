# =============================================================================
# Type III / sequential EMS synthesis for UNBALANCED designs (pure, base R).
#
# lme4 does not expose the expected-mean-square (EMS) linear-combination
# coefficients that the tolerance interval (Francq Eq. 26) needs. For balanced
# designs these are integers (design-*.R); for unbalanced designs they are
# non-integer and design-dependent. This module computes them by Hartley's
# method of synthesis on the ANOVA quadratic forms:
#
#   E[j, g] = tr(Q_j Z_g Z_g') / df_j
#
# where Q_j is the sum-of-squares quadratic form of ANOVA term j and Z_g is the
# incidence matrix of variance component g. The combination coefficients solve
#   E' k = c   (observable target c = 1)  ->  sigma^2_T = sum_j k_j EMS_j.
#
# Crossed-with-interaction uses Type III SS (sum-to-zero contrasts = SAS Type
# III); validated against Francq Section 4.4 (kA=0.128/kB=0.151/kAB=0.265/
# ke=0.456). Nested uses the sequential (hierarchical) block projections.
# Pure base R only -- no lme4/Matrix/numDeriv -- so it is unit-testable from a
# hand-entered cell-count table. The nested arm has no published worked oracle at all.
# =============================================================================

#' Synthesize EMS coefficients for an unbalanced design (pure)
#'
#' @param counts For `design = "crossed"`, an A-by-B matrix of cell counts
#'   `n_ij` (all cells must be filled). For `design = "nested"`, a list of length
#'   A whose i-th element is the vector of replicate counts for the fine levels
#'   nested within coarse level i.
#' @param design `"crossed"` or `"nested"`.
#' @return A list with `E` (EMS coefficient matrix, rows = ANOVA terms, cols =
#'   variance components, both ending in `"Residual"`), `df` (per-term degrees of
#'   freedom), `k` (observable combination coefficients; guaranteed `sum(k) == 1`),
#'   and `terms`.
#' @noRd
ems_synthesis <- function(counts, design = c("crossed", "nested")) {
  design <- match.arg(design)
  if (design == "crossed") .ems_synth_crossed(counts) else .ems_synth_nested(counts)
}

## tr(Q Z Z') without forming Z Z'
.tr_QZ <- function(q, z) sum((q %*% z) * z)

## Assemble E, df, k from the term quadratic forms Q, the component incidence
## matrices Zg, and the term df. Enforces the Sum(k) = 1 hard invariant.
.assemble_ems <- function(Q, Zg, dfF, terms) {
  E <- outer(terms, terms, Vectorize(function(f, g) .tr_QZ(Q[[f]], Zg[[g]]) / dfF[[f]]))
  dimnames(E) <- list(term = terms, component = terms)
  E[abs(E) < 1e-9] <- 0                                   # clean numerical zeros
  list(E = E, df = dfF, k = .solve_k_observable(E), terms = terms)
}

## Solve the observable combination coefficients k (E' k = 1) and enforce the
## hard Sum(k) = 1 invariant. Every EMS has residual-component coefficient 1, so
## sum(k) must be exactly 1; this caught the Francq Section 4.4 transcription typo
## (printed ke=0.464; forced 0.456) and will catch a malformed count table in
## production.
.solve_k_observable <- function(E) {
  k <- tryCatch(solve(t(E), rep(1, nrow(E))),
                error = function(e)
                  stop("EMS synthesis: the coefficient matrix is singular ",
                       "(degenerate design).", call. = FALSE))
  names(k) <- rownames(E)
  if (abs(sum(k) - 1) > 1e-6) {
    stop(sprintf(paste("EMS synthesis invariant violated: sum(k) = %.8f, must",
                       "be 1. The count table is malformed or the design is",
                       "unsupported."), sum(k)), call. = FALSE)
  }
  k
}

## ---- crossed-with-interaction: Type III SS (sum-to-zero contrasts) ----------
.ems_synth_crossed <- function(counts) {
  if (!is.matrix(counts) || any(counts < 1)) {
    stop("crossed synthesis needs an A-by-B matrix with every cell filled ",
         "(n_ij >= 1).", call. = FALSE)
  }
  rg <- expand.grid(a = seq_len(nrow(counts)), b = seq_len(ncol(counts)))
  rows <- rg[rep(seq_len(nrow(rg)), as.vector(counts)), ]
  af <- factor(rows$a)
  bf <- factor(rows$b)
  cellf <- factor(paste(rows$a, rows$b, sep = ":"))
  n <- nrow(rows)

  x <- stats::model.matrix(~ af * bf,
                           contrasts.arg = list(af = stats::contr.sum,
                                                bf = stats::contr.sum))
  asg <- attr(x, "assign")                               # 0 int, 1 af, 2 bf, 3 af:bf
  blocks <- list(A = which(asg == 1), B = which(asg == 2), AB = which(asg == 3))
  winv <- solve(crossprod(x))
  # Type III SS quadratic form Q_F = P' (W^-1_FF)^-1 P, P = W^-1_F. X'
  Q <- lapply(blocks, function(f) {
    p <- winv[f, , drop = FALSE] %*% t(x)
    t(p) %*% solve(winv[f, f, drop = FALSE]) %*% p
  })
  Q$Residual <- diag(n) - x %*% winv %*% t(x)
  dfF <- c(A = length(blocks$A), B = length(blocks$B), AB = length(blocks$AB),
           Residual = n - ncol(x))
  Zg <- list(A = stats::model.matrix(~ af - 1),
             B = stats::model.matrix(~ bf - 1),
             AB = stats::model.matrix(~ cellf - 1),
             Residual = diag(n))
  .assemble_ems(Q, Zg, dfF, c("A", "B", "AB", "Residual"))
}

## ---- nested: sequential (hierarchical) block projections --------------------
.ems_synth_nested <- function(counts) {
  if (!is.list(counts) || length(counts) < 2L) {
    stop("nested synthesis needs a list of >= 2 coarse levels, each a vector ",
         "of fine-level replicate counts.", call. = FALSE)
  }
  a_levels <- length(counts)
  k_fine <- sum(lengths(counts))
  n <- sum(unlist(counts))
  rows <- do.call(rbind, lapply(seq_len(a_levels), function(i) {
    do.call(rbind, lapply(seq_along(counts[[i]]), function(j) {
      data.frame(a = i, b = paste(i, j, sep = "."))[rep(1, counts[[i]][j]), ]
    }))
  }))
  af <- factor(rows$a)
  bf <- factor(rows$b)                                   # globally-unique fine levels

  h1 <- matrix(1 / n, n, n)
  ha <- .block_proj(af)                                  # onto coarse-group means
  hb <- .block_proj(bf)                                  # onto fine-group means
  Q <- list(A = ha - h1, B = hb - ha, Residual = diag(n) - hb)
  dfF <- c(A = a_levels - 1, B = k_fine - a_levels, Residual = n - k_fine)
  Zg <- list(A = stats::model.matrix(~ af - 1),
             B = stats::model.matrix(~ bf - 1),
             Residual = diag(n))
  .assemble_ems(Q, Zg, dfF, c("A", "B", "Residual"))
}

## orthogonal projection onto the group means of a factor
.block_proj <- function(f) {
  z <- stats::model.matrix(~ f - 1)
  z %*% (t(z) / colSums(z))
}

#' Build the `ems` list for an unbalanced design from a synthesis
#'
#' Turns an `ems_synthesis()` result plus the estimated variance components into
#' the `ems = list(ms, k, df)` slot an `re_components` / `ti_francq` expects:
#' `ms_j = sum_g E[j,g] sigma2[g]` (the estimated EMS values), `k` and `df` from
#' the synthesis.
#'
#' @param syn An `ems_synthesis()` result.
#' @param sigma2 Named numeric variance components, aligned to `syn$terms`
#'   (component order, residual last).
#' @return `list(ms, k, df)` with names `syn$terms`.
#' @noRd
ems_list_from_synthesis <- function(syn, sigma2) {
  comp <- sigma2[syn$terms]
  if (anyNA(comp)) {
    stop("`sigma2` is missing components for: ",
         paste(setdiff(syn$terms, names(sigma2)), collapse = ", "), call. = FALSE)
  }
  list(ms = stats::setNames(as.vector(syn$E %*% comp), syn$terms),
       k = syn$k,
       df = syn$df)
}

#' Unbalanced nested/crossed EMS adapter (pure)
#'
#' Bridges the design descriptor + cell-count table to the synthesis: maps the
#' grouping-named variance components to the synthesis canonical order, runs
#' `ems_synthesis()`, and returns the `ems`/`coefs` an `re_components` expects.
#' Pure --- the `counts` table is plain data extracted by the lme4 layer.
#'
#' @param desc Design descriptor (uses `grp_A`/`grp_B`/`grp_AB` for crossed,
#'   `grp_alpha`/`grp_beta` for nested).
#' @param components Named numeric variance components (grouping names + `Residual`).
#' @param target `"observable"` / `"true_value"` (selects `coefs`).
#' @param counts Cell-count table (matrix for crossed, list for nested).
#' @param design `"crossed"` or `"nested"`.
#' @return `list(ems = list(ms, k, df), coefs)`.
#' @noRd
ems_unbalanced <- function(desc, components, target, counts, design) {
  syn <- ems_synthesis(counts, design)
  grp <- if (design == "crossed") {
    c(A = desc$grp_A, B = desc$grp_B, AB = desc$grp_AB, Residual = "Residual")
  } else {
    c(A = desc$grp_alpha, B = desc$grp_beta, Residual = "Residual")
  }
  sigma2 <- stats::setNames(components[grp[syn$terms]], syn$terms)
  list(ems = ems_list_from_synthesis(syn, sigma2),
       coefs = .target_coefs(components, target))
}
