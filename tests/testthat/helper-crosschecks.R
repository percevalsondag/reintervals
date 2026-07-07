# Cross-check machinery for the unbalanced EMS synthesis (tests only).

# Francq Section 4.4 cell-count table (alpha=4 rows, beta=5 cols), from the
# paper's replicate vectors ((1,2,1,9)',(2,9,8,1)',(1,1,1,3)',(7,1,6,2)',(1,3,1,1)').
build_44_counts <- function() {
  sapply(list(c(1, 2, 1, 9), c(2, 9, 8, 1), c(1, 1, 1, 3),
              c(7, 1, 6, 2), c(1, 3, 1, 1)), identity)
}

# §4.4 data: expand the cell-count table to rows with an arbitrary deterministic
# response (the EMS coefficients are count-based, so y is irrelevant to the gate).
build_44_data <- function() {
  nij <- build_44_counts()
  rows <- do.call(rbind, lapply(1:4, function(i) do.call(rbind, lapply(1:5, function(j)
            data.frame(a = i, b = j)[rep(1, nij[i, j]), ]))))
  rows$y <- 2 * rows$a + 1.5 * rows$b +
    rep_len(c(-1.3, 0.4, 1.1, -0.7, 0.5, -0.2, 0.9), nrow(rows))
  rows$a <- factor(rows$a)
  rows$b <- factor(rows$b)
  rows
}

# --- n0-approximation arm (ti_approx coefficients) ----------
# This is NOT Type III; it substitutes an effective replicate into the balanced
# formula. Used ONLY to characterise divergence, never as an engine.
.n_eff <- function(counts, method = "n0") {
  nn <- sum(counts)
  kk <- length(counts)
  if (method == "flat") return(nn / kk)
  if (kk < 2) return(mean(counts))
  (nn - sum(counts^2) / nn) / (kk - 1)
}

n0_approx_k <- function(counts, design, method = "n0") {
  if (design == "crossed") {
    a <- nrow(counts); b <- ncol(counts)
    n <- .n_eff(as.vector(counts), method)
    setNames(c(1 / (n * b), 1 / (n * a), 1 / n - 1 / (n * a) - 1 / (n * b), 1 - 1 / n),
             c("A", "B", "AB", "Residual"))
  } else {
    a <- length(counts)
    bbar <- sum(lengths(counts)) / a
    n <- .n_eff(unlist(counts), method)
    setNames(c(1 / (n * bbar), 1 / n - 1 / (n * bbar), 1 - 1 / n),
             c("A", "B", "Residual"))
  }
}

# --- independent closed-form nested EMS (derived from the block-projection
# trace algebra, coded separately from the matrix-projection synthesis) -------
nested_ems_closed <- function(counts) {
  a <- length(counts)
  ni <- vapply(counts, sum, numeric(1))
  nn <- sum(ni)
  kf <- sum(lengths(counts))
  s1 <- sum(vapply(counts, function(v) sum(v^2) / sum(v), numeric(1)))
  s2 <- sum(unlist(counts)^2) / nn
  e_aa <- (nn - sum(ni^2) / nn) / (a - 1)        # sigma2_alpha coeff in EMS_A
  e_ab <- (s1 - s2) / (a - 1)                    # sigma2_beta coeff in EMS_A
  e_bb <- (nn - s1) / (kf - a)                   # sigma2_beta coeff in EMS_B(A)
  rbind(A = c(A = e_aa, B = e_ab, Residual = 1),
        B = c(A = 0, B = e_bb, Residual = 1),
        Residual = c(A = 0, B = 0, Residual = 1))
}

# full n0-approx ems list (ms, k, df) so its TI width can be compared
n0_approx_ems <- function(counts, design, sigma2, method = "n0") {
  k <- n0_approx_k(counts, design, method)
  if (design == "crossed") {
    a <- nrow(counts); b <- ncol(counts); nn <- sum(counts)
    n <- .n_eff(as.vector(counts), method)
    s2a <- sigma2[["A"]]; s2b <- sigma2[["B"]]; s2ab <- sigma2[["AB"]]; s2e <- sigma2[["Residual"]]
    ms <- c(A = s2e + n * s2ab + n * b * s2a, B = s2e + n * s2ab + n * a * s2b,
            AB = s2e + n * s2ab, Residual = s2e)
    df <- c(A = a - 1, B = b - 1, AB = (a - 1) * (b - 1), Residual = nn - a * b)
  } else {
    a <- length(counts); kf <- sum(lengths(counts)); nn <- sum(unlist(counts))
    bbar <- kf / a; n <- .n_eff(unlist(counts), method)
    s2a <- sigma2[["A"]]; s2b <- sigma2[["B"]]; s2e <- sigma2[["Residual"]]
    ms <- c(A = s2e + n * s2b + n * bbar * s2a, B = s2e + n * s2b, Residual = s2e)
    df <- c(A = a - 1, B = kf - a, Residual = nn - kf)
  }
  list(ms = ms, k = k, df = df)
}

# TI width from an ems list, holding everything else fixed (for divergence dir.)
ti_width_from_ems <- function(ems, sigma2, var_mean = 0.01, content = 0.95, conf = 0.90) {
  comp <- re_components(
    components = sigma2, dfs = c(pi = Inf), mean = 0, var_mean = var_mean,
    ems = ems, design = "x"
  )
  ti <- ti_francq(comp, content = content, conf = conf)
  ti$upper - ti$lower
}
