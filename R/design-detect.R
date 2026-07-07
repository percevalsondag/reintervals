#' Classify the random structure of a mixed model (pure)
#'
#' Operates on plain data --- the
#' grouping-factor list and the sample size --- so it carries no `lme4`
#' dependency and is unit-testable with hand-built factors. The extraction layer
#' (`vc_extract()`) pulls `flist`/`N` from the fit and hands them here.
#'
#' @param flist A list of grouping factors (as from `lme4::getME(model,
#'   "flist")`).
#' @param N Integer sample size.
#' @return A list describing the design: `type` is one of `"oneway"`,
#'   `"nested"`, `"crossed"`, or an `"unsupported-*"` / `"crossed-no-interaction"`
#'   marker, with the level counts, effective replicate(s), grouping names, and
#'   balance flag that the EMS adapters need.
#' @noRd
classify_design <- function(flist, N) {
  fl <- flist
  nf <- length(fl)
  refines <- function(i, j) all(rowSums(table(fl[[i]], fl[[j]]) > 0) == 1L)

  if (nf == 1L) {
    ni <- as.integer(table(fl[[1]]))
    A <- length(ni)
    return(list(type = "oneway", A = A, N = N,
                n0 = (N - sum(ni^2) / N) / (A - 1),
                balanced = length(unique(ni)) == 1L && all(ni >= 2),
                grp = names(fl)[1]))
  }
  if (nf == 2L) {
    r12 <- refines(1, 2)
    r21 <- refines(2, 1)
    if (xor(r12, r21)) {                                   # one nested in the other
      fine <- if (r12) 1L else 2L
      coarse <- if (r12) 2L else 1L
      A <- nlevels(fl[[coarse]])
      AB <- nlevels(fl[[fine]])
      n <- N / AB
      # Balanced nesting needs equal fine-per-coarse AND equal reps; otherwise it
      # is a valid UNBALANCED nested design (handled by the EMS synthesis).
      bal <- AB %% A == 0 &&
        length(unique(table(fl[[coarse]]))) == 1L &&
        length(unique(table(fl[[fine]]))) == 1L &&
        isTRUE(all.equal(n, round(n))) && n >= 2
      return(list(type = "nested", A = A,
                  B = if (bal) AB / A else NA_real_,
                  n = if (bal) n else NA_real_,
                  N = N, balanced = bal,
                  grp_alpha = names(fl)[coarse], grp_beta = names(fl)[fine]))
    }
    return(list(type = "crossed-no-interaction"))
  }
  if (nf == 3L) {
    is_int <- function(i) {
      o <- setdiff(1:3, i)
      refines(i, o[1]) && refines(i, o[2])
    }
    intc <- which(vapply(1:3, is_int, logical(1)))
    if (length(intc) == 1L) {
      m <- setdiff(1:3, intc)
      if (refines(m[1], m[2]) || refines(m[2], m[1])) {
        return(list(type = "unsupported-3factor"))         # mains not mutually crossed
      }
      A <- nlevels(fl[[m[1]]])
      B <- nlevels(fl[[m[2]]])
      n <- N / (A * B)
      tab <- table(fl[[m[1]]], fl[[m[2]]])
      bal <- length(unique(as.vector(tab))) == 1L &&
        isTRUE(all.equal(n, round(n))) && n >= 2 &&
        nlevels(fl[[intc]]) == A * B
      return(list(type = "crossed", A = A, B = B, n = n, N = N, balanced = bal,
                  grp_A = names(fl)[m[1]], grp_B = names(fl)[m[2]],
                  grp_AB = names(fl)[intc]))
    }
    return(list(type = "unsupported-3factor"))
  }
  list(type = paste0("unsupported-", nf, "factor"))
}

## Target selector a_i: observable = all ones; true_value zeroes the residual.
.target_coefs <- function(components, target) {
  a <- stats::setNames(rep(1, length(components)), names(components))
  if (identical(target, "true_value")) a[["Residual"]] <- 0
  a
}

## Dispatch to the per-design EMS adapter. Returns list(ems, coefs, note);
## ems is NULL (with an explanatory note) outside the closed-form catalog
## (unbalanced nested/crossed -> deferred v2; off-catalog designs).
.design_components <- function(desc, components, target, counts = NULL) {
  unbal <- function(design) {
    if (is.null(counts)) return(NULL)                     # need the cell-count table
    ems_unbalanced(desc, components, target, counts, design)
  }
  built <- switch(
    desc$type,
    oneway  = ems_oneway(desc, components, target),        # balanced AND unbalanced (n0)
    nested  = if (isTRUE(desc$balanced)) ems_nested(desc, components, target) else unbal("nested"),
    crossed = if (isTRUE(desc$balanced)) ems_crossed(desc, components, target) else unbal("crossed"),
    NULL
  )
  if (is.null(built)) {
    note <- paste0("Design '", desc$type,
                   "' has no closed-form EMS decomposition; ems = NULL.")
    return(list(ems = NULL, coefs = .target_coefs(components, target), note = note))
  }
  c(built, list(note = NULL))
}
