# B3: vectorised P / level / sides. One call computes every content,
# confidence and side; the rows are the full cross product with the evaluation
# points, t0 fastest, then sides, then P, then level. The identity contract:
# with an explicit seed, each block of a vectorised call is identical() to the
# corresponding scalar call, per-row diagnostics included (the GPQ engines
# share one draw set and one tau2 solve per t0 across the combinations; every
# other engine loops its scalar code). Plumbing only: the scalar outputs are
# pinned by the existing tests, which pass unchanged.

## Two tiers. The small identity set (one fit per engine, default arms,
## P = c(0.90, 0.99), level = c(0.90, 0.95), all three sides, two t0) runs
## everywhere, CRAN included. The all-arms loop and the full grids
## (P = c(0.90, 0.95, 0.99), three t0, 11 TI route cases) are skip_on_cran():
## they run in local and evidence runs (NOT_CRAN = "true").
vec_P     <- c(0.90, 0.95, 0.99)
vec_level <- c(0.90, 0.95)
vec_sides <- c("two.sided", "lower", "upper")
sm_P      <- c(0.90, 0.99)
sm_level  <- c(0.90, 0.95)

## the fits are built once per test file run (memoised)
vec_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- vec_fits_build()
    cache
  }
})

vec_fits_build <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  fit <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, control = lc)))
  set.seed(20261110)
  ow <- data.frame(g = factor(rep(1:6, each = 4)))
  ow$y <- 10 + rnorm(6, sd = 2)[ow$g] + rnorm(24)
  ne <- expand.grid(r = 1:2, b = 1:3, a = factor(1:5))
  ne$b <- factor(paste(ne$a, ne$b))
  ne$y <- 5 + rnorm(5, sd = 1.5)[ne$a] + rnorm(15)[ne$b] + rnorm(nrow(ne))
  cr <- expand.grid(r = 1:2, a = factor(1:4), b = factor(1:6))
  cr$y <- 3 + rnorm(4)[cr$a] + rnorm(6, sd = 1.2)[cr$b] +
    rnorm(24, sd = 0.5)[interaction(cr$a, cr$b)] + rnorm(nrow(cr))
  ## fixed slope with single-observation lots (the anova-mls TI route)
  so <- data.frame(
    y    = c(70, 71, 69, 68, 66, 65, 71, 69, 68, 66, 66, 71, 69,
             73, 73, 71, 71, 70, 72, 73, 70),
    time = c(0, 1, 3, 6, 9, 12, 0, 1, 3, 6, 9, 0, 3, 0, 1, 3, 6, 9, 0, 1, 0),
    lot  = factor(rep(c("A", "B", "C", "D", "E", "F"), c(6, 5, 2, 5, 2, 1))))
  ss <- lme4::sleepstudy
  d3 <- data.frame(Days = c(0, 5, 9))
  list(
    oneway      = list(fit = fit(y ~ 1 + (1 | g), ow), nd = NULL),
    nested      = list(fit = fit(y ~ 1 + (1 | a/b), ne), nd = NULL),
    crossed     = list(fit = fit(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), cr), nd = NULL),
    fixed_slope = list(fit = fit(Reaction ~ Days + (1 | Subject), ss), nd = d3),
    fs_single   = list(fit = fit(y ~ time + (1 | lot), so),
                       nd = data.frame(time = c(0, 6, 12))),
    rs_un       = list(fit = fit(Reaction ~ Days + (1 + Days | Subject), ss), nd = d3),
    rs_ind      = list(fit = fit(Reaction ~ Days + (1 + Days || Subject), ss), nd = d3))
}

## the result without its attributes (data columns only)
vec_strip <- function(x) {
  y <- x
  attr(y, "diagnostics") <- NULL; attr(y, "grouping") <- NULL
  class(y) <- "data.frame"; rownames(y) <- NULL
  y
}

## per-call diagnostics (one value per call, never row-aligned)
vec_per_call <- c("gpq_seed", "rng_kind", "mc_sections", "M", "beta_hat", "Vhat_Y",
                  "pi_df_fallback", "singular_check", "grouping_reason")

## Assert the identity contract for one verb call. `f(P, level, sides)` runs
## the verb; P = NULL for the coverage verbs. Returns the number of blocks.
expect_vec_identity <- function(f, P = NULL, level = vec_level, sides = vec_sides,
                                label = "") {
  vec <- suppressWarnings(f(P, level, sides))
  grid <- expand.grid(sides = sides, P = if (is.null(P)) NA_real_ else P,
                      level = level, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vd <- attr(vec, "diagnostics"); vs <- vec_strip(vec)
  n_tot <- 0L
  for (j in seq_len(nrow(grid))) {
    p <- if (is.null(P)) NULL else grid$P[j]
    sc <- suppressWarnings(f(p, grid$level[j], grid$sides[j]))
    n <- nrow(sc)
    idx <- n_tot + seq_len(n)
    n_tot <- n_tot + n
    lab <- sprintf("%s [sides=%s P=%s level=%s]", label, grid$sides[j],
                   format(grid$P[j]), format(grid$level[j]))
    blk <- vs[idx, , drop = FALSE]; rownames(blk) <- NULL
    expect_identical(blk, vec_strip(sc), label = lab)
    sd <- attr(sc, "diagnostics")
    expect_setequal(names(vd), names(sd))
    for (fld in names(sd)) {
      if (fld %in% vec_per_call) {
        expect_identical(vd[[fld]], sd[[fld]], label = paste(lab, fld))
      } else {
        expect_identical(vd[[fld]][idx], rep_len(sd[[fld]], n),
                         label = paste(lab, fld))
      }
    }
    expect_identical(attr(vec, "grouping"), attr(sc, "grouping"))
  }
  expect_identical(nrow(vec), n_tot, label = paste(label, "row count"))
  nrow(grid)
}

## every value of every arm (the defaults included)
vec_arms <- list(
  ci_df            = c("satterthwaite", "kenward-roger", "containment"),
  slope_ci_df      = c("hessian-satterthwaite", "kenward-roger", "k-minus-1"),
  pi_df            = c("francq-satterthwaite", "containment"),
  ngm_method       = c("package", "broad-inference", "hts", "classical"),
  slope_ngm_method = c("reml", "broad-inference", "gpq", "fv-anova"))
vec_arm_verb <- c(ci_df = "ci", slope_ci_df = "ci", pi_df = "pi",
                  ngm_method = "ngm", slope_ngm_method = "ngm")

## two evaluation points per fit with a time variable (RI fits: intercept only)
sm_nd <- function(nm) {
  switch(nm,
         fs_single = data.frame(time = c(0, 12)),
         rs_un = , rs_ind = data.frame(Days = c(0, 9)),
         NULL)
}

test_that("identity (small set): one fit per engine, default arms", {
  skip_if_not_installed("lme4")
  fits <- vec_fits()
  ctl <- lmm_interval_control(seed = 11L, M = 2000L)
  ## CI / PI / NGM, default arms: EMS (one-way, nested, crossed), the fixed-slope
  ## single-obs fit, and the random-slope REML engine (rs_un; the || fit runs
  ## the same engine and is in the full grid below)
  for (nm in c("oneway", "nested", "crossed", "fs_single", "rs_un")) {
    m <- fits[[nm]]$fit; nd <- sm_nd(nm)
    for (verb_nm in c("ci", "pi", "ngm")) {
      verb <- switch(verb_nm, ci = ci_lmm, pi = pi_lmm, ngm = new_group_mean_lmm)
      nb <- expect_vec_identity(function(P, level, sides)
        verb(m, newdata = nd, level = level, sides = sides, control = ctl),
        level = sm_level, label = sprintf("small %s %s", nm, verb_nm))
      expect_identical(nb, 6L)
    }
  }
  ## TI: Hoffman-Kringle (one-way), EMS (nested, crossed), anova-mls (fixed
  ## slope, single-obs lots), GPQ (RS UN both targets, RS ||)
  cases <- list(
    list(nm = "oneway",    over = "observation"),
    list(nm = "nested",    over = "observation"),
    list(nm = "crossed",   over = "observation"),
    list(nm = "fs_single", over = "observation"),
    list(nm = "rs_un",     over = "observation"),
    list(nm = "rs_un",     over = "group_mean"),
    list(nm = "rs_ind",    over = "observation"))
  for (cs in cases) {
    m <- fits[[cs$nm]]$fit; nd <- sm_nd(cs$nm)
    nb <- expect_vec_identity(function(P, level, sides)
      ti_lmm(m, newdata = nd, P = P, over = cs$over, level = level, sides = sides,
             control = ctl),
      P = sm_P, level = sm_level, label = sprintf("small ti %s %s", cs$nm, cs$over))
    expect_identical(nb, 12L)
  }
})

test_that("identity: ci / pi / new_group_mean, every fit x every arm value", {
  skip_on_cran()   # full grid; runs in local and evidence runs
  skip_if_not_installed("lme4")
  skip_if_not_installed("pbkrtest")
  fits <- vec_fits()
  verbs <- list(ci = ci_lmm, pi = pi_lmm, ngm = new_group_mean_lmm)
  n_cells <- 0L
  for (nm in names(fits)) {
    m <- fits[[nm]]$fit; nd <- fits[[nm]]$nd
    for (opt in names(vec_arms)) for (val in vec_arms[[opt]]) {
      args <- list(seed = 11L, M = 2000L); args[[opt]] <- val
      ctl <- do.call(lmm_interval_control, args)
      verb <- verbs[[vec_arm_verb[[opt]]]]
      ## the fixed-slope single-obs fit has no NGM GPQ / KR issues; every
      ## design the verb rejects is skipped (none here)
      expect_vec_identity(function(P, level, sides)
        verb(m, newdata = nd, level = level, sides = sides, control = ctl),
        label = sprintf("%s %s %s=%s", nm, vec_arm_verb[[opt]], opt, val))
      n_cells <- n_cells + 1L
    }
  }
  expect_identical(n_cells, length(fits) * sum(lengths(vec_arms)))
})

test_that("identity: ti_lmm on every TI engine (GPQ both targets, EMS, anova-mls)", {
  skip_on_cran()   # full grid; runs in local and evidence runs
  skip_if_not_installed("lme4")
  fits <- vec_fits()
  ctl <- lmm_interval_control(seed = 5L, M = 2000L)
  cases <- list(
    list(nm = "oneway",      method = "auto",    over = "observation"),  # Hoffman-Kringle
    list(nm = "oneway",      method = "ems-mls", over = "observation"),  # Francq EMS
    list(nm = "nested",      method = "auto",    over = "observation"),
    list(nm = "crossed",     method = "auto",    over = "observation"),
    list(nm = "fixed_slope", method = "auto",    over = "observation"),  # EMS
    list(nm = "fixed_slope", method = "anova-mls", over = "observation"),
    list(nm = "fs_single",   method = "auto",    over = "observation"),  # Montes
    list(nm = "rs_un",       method = "auto",    over = "observation"),  # GPQ
    list(nm = "rs_un",       method = "auto",    over = "group_mean"),
    list(nm = "rs_ind",      method = "auto",    over = "observation"),
    list(nm = "rs_ind",      method = "gpq",     over = "group_mean"))
  for (cs in cases) {
    m <- fits[[cs$nm]]$fit; nd <- fits[[cs$nm]]$nd
    nb <- expect_vec_identity(function(P, level, sides)
      ti_lmm(m, newdata = nd, P = P, over = cs$over, level = level, sides = sides,
             method = cs$method, control = ctl),
      P = vec_P, label = sprintf("ti %s %s %s", cs$nm, cs$method, cs$over))
    expect_identical(nb, 18L)
  }
})

test_that("identity: lmm_predict forwards vector level / sides / P", {
  skip_on_cran()   # full grid; runs in local and evidence runs
  skip_if_not_installed("lme4")
  fits <- vec_fits()
  ctl <- lmm_interval_control(seed = 9L, M = 2000L)
  m <- fits$rs_un$fit; nd <- fits$rs_un$nd
  q <- function(x) vec_strip(suppressWarnings(x))
  expect_identical(
    q(lmm_predict(m, "tolerance", nd, P = vec_P, level = vec_level, sides = vec_sides,
                  control = ctl)),
    q(ti_lmm(m, nd, P = vec_P, level = vec_level, sides = vec_sides, control = ctl)))
  for (iv in c("confidence", "prediction", "new_group_mean")) {
    verb <- switch(iv, confidence = ci_lmm, prediction = pi_lmm,
                   new_group_mean = new_group_mean_lmm)
    expect_identical(
      q(lmm_predict(m, iv, nd, level = vec_level, sides = vec_sides, control = ctl)),
      q(verb(m, nd, level = vec_level, sides = vec_sides, control = ctl)))
    expect_vec_identity(function(P, level, sides)
      lmm_predict(m, iv, nd, level = level, sides = sides, control = ctl),
      label = paste("lmm_predict", iv))
  }
  ## alpha is the same pair as level
  expect_identical(
    q(ci_lmm(m, nd, alpha = 1 - vec_level, sides = "lower")),
    q(ci_lmm(m, nd, level = vec_level, sides = "lower")))
})

test_that("row order: t0 fastest, then sides, then P, then level", {
  skip_if_not_installed("lme4")
  m <- vec_fits()$rs_un$fit
  nd <- data.frame(Days = c(0, 5, 9))
  v <- ti_lmm(m, nd, P = vec_P, level = vec_level, sides = vec_sides,
              control = lmm_interval_control(seed = 1L, M = 2000L))
  g <- expand.grid(Days = nd$Days, sides = vec_sides, P = vec_P, level = vec_level,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  expect_identical(nrow(v), 54L)
  expect_identical(v$Days, g$Days)
  expect_identical(v$sides, g$sides)
  expect_identical(v$P, g$P)
  expect_identical(v$level, g$level)
  ## the order follows the order the values are given in
  w <- ci_lmm(m, nd, level = c(0.95, 0.90), sides = c("upper", "lower"))
  expect_identical(w$sides, rep(rep(c("upper", "lower"), each = 3), 2))
  expect_identical(w$level, rep(c(0.95, 0.90), each = 6))
})

test_that("sides: the result depends only on the value passed", {
  skip_if_not_installed("lme4")
  m <- vec_fits()$fixed_slope$fit
  f <- function(m, sides = c("two.sided", "lower", "upper")) ti_lmm(m, sides = sides)
  r3 <- f(m)
  r1 <- ti_lmm(m)
  expect_identical(r3$sides, rep(vec_sides, each = nrow(r1)))
  expect_identical(r1$sides, rep("two.sided", nrow(r1)))
  blk <- vec_strip(r3)[seq_len(nrow(r1)), , drop = FALSE]; rownames(blk) <- NULL
  expect_identical(blk, vec_strip(r1))
  ## the same for the coverage verbs and lmm_predict
  g <- function(m, sides = c("two.sided", "lower", "upper")) ci_lmm(m, sides = sides)
  expect_identical(nrow(g(m)), 3L * nrow(ci_lmm(m)))
  h <- function(m, sides = c("two.sided", "lower", "upper"))
    lmm_predict(m, "prediction", sides = sides)
  expect_identical(nrow(h(m)), 3L * nrow(lmm_predict(m, "prediction")))
})

test_that("scalar defaults: sides default is two-sided; all three is explicit", {
  skip_if_not_installed("lme4")
  m <- vec_fits()$fixed_slope$fit
  nd <- data.frame(Days = c(0, 5))
  expect_identical(unique(ci_lmm(m, nd)$sides), "two.sided")
  expect_identical(nrow(ci_lmm(m, nd)), 2L)
  expect_identical(nrow(ci_lmm(m, nd, sides = c("two.sided", "lower", "upper"))), 6L)
  expect_identical(unique(lmm_predict(m, "confidence", nd)$sides), "two.sided")
  expect_identical(nrow(lmm_predict(m, "confidence", nd, sides = vec_sides)), 6L)
})

test_that("GPQ: one draw set per call; per-call seed, rng kind and MC sections", {
  skip_if_not_installed("lme4")
  fits <- vec_fits()
  m <- fits$rs_un$fit; nd <- fits$rs_un$nd
  ## seed = NULL: a vectorised call draws ONE seed (one draw set) for every
  ## combination; separate scalar calls each draw their own seed from the
  ## caller's stream, so equality holds only with an explicit seed.
  ctl0 <- lmm_interval_control(M = 2000L)
  set.seed(42); v <- ti_lmm(m, nd, P = c(0.90, 0.99), sides = c("lower", "upper"),
                            control = ctl0)
  dv <- attr(v, "diagnostics")
  expect_length(dv$gpq_seed, 1L)
  expect_identical(dv$mc_sections, 20L)
  expect_length(dv$rng_kind, 3L)
  expect_length(dv$mc_se_lower, nrow(v))
  set.seed(42)
  s1 <- ti_lmm(m, nd, P = 0.90, sides = "lower", control = ctl0)
  s2 <- ti_lmm(m, nd, P = 0.90, sides = "upper", control = ctl0)
  ## the first scalar call draws the same seed as the vectorised call ...
  expect_identical(attr(s1, "diagnostics")$gpq_seed, dv$gpq_seed)
  expect_identical(vec_strip(v)[1:3, ], vec_strip(s1))
  ## ... the second draws another, so its rows differ from the vectorised block
  expect_false(identical(attr(s2, "diagnostics")$gpq_seed, dv$gpq_seed))
  blk2 <- vec_strip(v)[4:6, ]; rownames(blk2) <- NULL
  expect_false(identical(blk2, vec_strip(s2)))
  ## with the vectorised call's seed made explicit, the scalar call matches
  s2e <- ti_lmm(m, nd, P = 0.90, sides = "upper",
                control = lmm_interval_control(seed = dv$gpq_seed, M = 2000L))
  expect_identical(blk2, vec_strip(s2e))
  ## the vectorised call advances the caller's stream by exactly one draw
  set.seed(42); invisible(ti_lmm(m, nd, P = vec_P, control = ctl0)); a <- runif(1)
  set.seed(42); invisible(sample.int(.Machine$integer.max, 1L)); b <- runif(1)
  expect_identical(a, b)
  ## the GPQ new-group-mean arm: one seed per call, too
  cg <- lmm_interval_control(M = 2000L, slope_ngm_method = "gpq")
  set.seed(3)
  g <- suppressWarnings(new_group_mean_lmm(m, nd, level = vec_level, sides = vec_sides,
                                           control = cg))
  expect_length(attr(g, "diagnostics")$gpq_seed, 1L)
  expect_length(attr(g, "diagnostics")$mc_se_upper, 18L)
})

test_that("warnings are signalled once per call, not once per combination", {
  skip_if_not_installed("lme4")
  fits <- vec_fits()
  count <- function(expr, class) {
    n <- 0L
    withCallingHandlers(expr, warning = function(w) {
      if (inherits(w, class)) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    n
  }
  rs <- fits$rs_un$fit
  far <- data.frame(Days = c(0, 15))                     # 15 is outside [0, 9]
  ctl <- lmm_interval_control(seed = 1L, M = 2000L)
  expect_identical(count(ti_lmm(rs, far, P = vec_P, level = vec_level,
                                sides = vec_sides, control = ctl),
                         "reintervals_warning_extrapolation"), 1L)
  expect_identical(count(ci_lmm(rs, far, level = vec_level, sides = vec_sides),
                         "reintervals_warning_extrapolation"), 1L)
  ## unvalidated arm: once per call
  expect_identical(count(new_group_mean_lmm(rs, fits$rs_un$nd, level = vec_level,
                                            sides = vec_sides,
                                            control = lmm_interval_control(
                                              slope_ngm_method = "fv-anova")),
                         "reintervals_warning_unvalidated"), 1L)
  expect_identical(count(new_group_mean_lmm(rs, fits$rs_un$nd, level = vec_level,
                                            sides = vec_sides,
                                            control = lmm_interval_control(
                                              M = 2000L, slope_ngm_method = "gpq")),
                         "reintervals_warning_unvalidated"), 1L)
  ## limitation (fixed slope with single-observation lots): once per call
  fs1 <- fits$fs_single
  expect_identical(count(ci_lmm(fs1$fit, fs1$nd, level = vec_level, sides = vec_sides),
                         "reintervals_warning_limitation"), 1L)
  ## singular fit: the EMS TI singular warning and the default df fallback,
  ## once per call (they are raised inside the looped scalar engine)
  dsg <- sqrt(0.001253)
  y <- as.vector(sapply(rep(0.981, 6), function(mm) mm + c(-dsg, 0, dsg)))
  sg <- suppressMessages(suppressWarnings(lme4::lmer(
    y ~ 1 + (1 | run), data = data.frame(y = y, run = factor(rep(1:6, each = 3))))))
  expect_identical(count(ti_lmm(sg, P = vec_P, level = vec_level, sides = vec_sides,
                                method = "ems-mls"),
                         "reintervals_warning_singular"), 1L)
  expect_identical(count(ci_lmm(sg, level = vec_level, sides = vec_sides),
                         "reintervals_warning_df_fallback"), 1L)
})

test_that("input validation for vector P / level / alpha / sides", {
  skip_if_not_installed("lme4")
  m <- vec_fits()$fixed_slope$fit
  err <- "reintervals_error_input"
  expect_error(ti_lmm(m, P = c(0.9, 0.9)), class = err)
  expect_error(ti_lmm(m, P = c(0.9, NA)), class = err)
  expect_error(ti_lmm(m, P = c(0.9, 1)), class = err)
  expect_error(ti_lmm(m, level = c(0.9, 0.9)), class = err)
  expect_error(ci_lmm(m, alpha = c(0.05, 0)), class = err)
  expect_error(ci_lmm(m, level = c(0.9, NA)), class = err)
  expect_error(ci_lmm(m, alpha = 0.05, level = c(0.9, 0.95)), "only one")
  expect_error(ci_lmm(m, sides = c("lower", "lower")), class = err)
  expect_error(ci_lmm(m, sides = c("lower", NA)), class = err)
  expect_error(ci_lmm(m, sides = "both"))
  ## partial matching, as match.arg
  expect_identical(ci_lmm(m, sides = c("low", "up"))$sides, c("lower", "upper"))
  ## over stays scalar (the full choice vector is the default, as match.arg)
  expect_error(ti_lmm(m, over = c("group_mean", "observation")))
  ## lmm_predict: P is still tolerance-only
  expect_error(lmm_predict(m, "confidence", P = c(0.9, 0.95)), class = err)
})
