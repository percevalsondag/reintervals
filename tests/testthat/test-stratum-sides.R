# Part 2d: (A) stratum df adjusted for the fixed effects; (B) one-sided intervals.
lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
fit_q <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, control = lc)))
dg <- function(x) attr(x, "diagnostics")

## ---- (A) stratum df ---------------------------------------------------------
v5b_data <- function() {          # the 0.2.0 validation's V5(b) design, A = 5
  set.seed(5e8 + 1e7 * 2 + 1e5 * 2 + 1e4 * 2 + 1)
  A <- 5; g <- factor(rep(seq_len(A), each = 4)); a <- rnorm(A); z <- runif(A, 0, 10)
  data.frame(g = g, z = z[g], y = 100 + 0.5 * z[g] + a[g] + rnorm(length(g)))
}

test_that("intercept-only designs keep the classical stratum df (oracles unchanged)", {
  skip_if_not_installed("lme4")
  expect_equal(unname(.stratum_df(fit_q(conc ~ 1 + (1 | run), hk_assay))), c(5, 12))
  set.seed(3); A <- 4; B <- 3; n <- 2
  fa <- factor(rep(1:A, each = B * n)); fc <- factor(rep(1:(A * B), each = n))
  d <- data.frame(fa, fc, y = rnorm(A)[fa] + rnorm(A * B)[fc] + rnorm(A * B * n))
  expect_equal(unname(.stratum_df(fit_q(y ~ 1 + (1 | fa) + (1 | fc), d))), c(A - 1, A * (B - 1), A * B * (n - 1)))
})

test_that("V5(b) between-group covariate: top stratum df = A - 2, used by Gm and containment", {
  skip_if_not_installed("lme4")
  fm <- fit_q(y ~ z + (1 | g), v5b_data())
  sd <- .stratum_df(fm)
  expect_equal(sd[["g"]], 3)                                  # A - 2, A = 5
  expect_equal(sd[["Residual"]], 20 - 5)                      # within df untouched
  x <- new_group_mean_lmm(fm, newdata = data.frame(z = 5))
  expect_equal(dg(x)$df, 3)
  expect_identical(dg(x)$df_rule, "top-factor")
  expect_equal(.ddf_containment(fm), 3)
  dc <- .design_components_fit(fm, .design_of(fm), .var_components(fm), "observable")
  expect_equal(unname(dc$ems$df), c(3, 15))                   # the EMS TI uses the same df
})

test_that("one-way with a within-group 2-level factor: residual df = N - A - 1, top df A - 1", {
  skip_if_not_installed("lme4")
  d <- v5b_data(); d$grp <- factor(rep(c("a", "b"), length.out = nrow(d)))
  fm <- fit_q(y ~ grp + (1 | g), d)
  N <- nrow(d); A <- nlevels(d$g)
  expect_equal(unname(.stratum_df(fm)), c(A - 1, N - A - 1))  # 4, 14
  dc <- .design_components_fit(fm, .design_of(fm), .var_components(fm), "observable")
  expect_equal(unname(dc$ems$df), c(A - 1, N - A - 1))
})

test_that("nested with a between-cell covariate: the cell stratum loses the df", {
  skip_if_not_installed("lme4")
  set.seed(3); A <- 4; B <- 3; n <- 2
  fa <- factor(rep(1:A, each = B * n)); fc <- factor(rep(1:(A * B), each = n)); w <- rnorm(A * B)
  d <- data.frame(fa, fc, w = w[fc], y = rnorm(A)[fa] + rnorm(A * B)[fc] + rnorm(A * B * n))
  fm <- fit_q(y ~ w + (1 | fa) + (1 | fc), d)
  expect_equal(unname(.stratum_df(fm)), c(A - 1, A * (B - 1) - 1, A * B * (n - 1)))
  ## a covariate constant within the TOP factor instead costs the top stratum
  d$v <- rnorm(A)[d$fa]
  fv <- fit_q(y ~ v + (1 | fa) + (1 | fc), d)
  expect_equal(unname(.stratum_df(fv)), c(A - 2, A * (B - 1), A * B * (n - 1)))
})

## ---- (B) one-sided intervals ------------------------------------------------
balanced_slope <- function() {
  set.seed(11)
  d <- expand.grid(time = c(0, 3, 6, 9, 12), lot = factor(1:6))
  d$y <- 100 - 0.5 * d$time + rnorm(6, 0, 1.5)[d$lot] + rnorm(nrow(d), 0, 0.7)
  d
}
balanced_nested <- function() {
  set.seed(4); A <- 5; B <- 3; n <- 2
  fa <- factor(rep(1:A, each = B * n)); fc <- factor(rep(1:(A * B), each = n))
  data.frame(fa, fc, y = 10 + rnorm(A, 0, 2)[fa] + rnorm(A * B)[fc] + rnorm(A * B * n))
}

test_that("one-sided 0.95 bound == two-sided 0.90 bound for CI / PI / new-group mean (EMS and anova paths)", {
  skip_if_not_installed("lme4")
  cases <- list(
    hk     = list(fit_q(conc ~ 1 + (1 | run), hk_assay), NULL),                       # NGM anova (HK)
    slope  = list(fit_q(y ~ time + (1 | lot), balanced_slope()), data.frame(time = 6)), # NGM anova (Fa)
    nested = list(fit_q(y ~ 1 + (1 | fa) + (1 | fc), balanced_nested()), NULL))        # NGM ems (Gm)
  for (nm in names(cases)) for (v in c("ci_lmm", "pi_lmm", "new_group_mean_lmm")) {
    fm <- cases[[nm]][[1]]; nd <- cases[[nm]][[2]]
    two <- get(v)(fm, newdata = nd, level = 0.90)
    lo <- get(v)(fm, newdata = nd, level = 0.95, sides = "lower")
    up <- get(v)(fm, newdata = nd, level = 0.95, sides = "upper")
    info <- paste(nm, v)
    expect_equal(lo$lower, two$lower, tolerance = 1e-12, info = info)
    expect_equal(up$upper, two$upper, tolerance = 1e-12, info = info)
    expect_identical(c(lo$upper, up$lower), c(Inf, -Inf), info = info)
    expect_identical(c(lo$sides, up$sides), c("lower", "upper"), info = info)
    expect_identical(lo$method, two$method, info = info)
    expect_true(all(lo$lower < lo$estimate & lo$estimate < lo$upper), info = info)
    expect_true(all(up$lower < up$estimate & up$estimate < up$upper), info = info)
  }
})

test_that("one-sided TI routing: random slope -> gpq (route_reason one_sided); reml-mls when forced", {
  skip_if_not_installed("lme4")
  fs <- fit_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
  nd <- data.frame(Days = c(2, 5)); ctl <- lmm_interval_control(seed = 1, M = 2000)
  for (s in c("lower", "upper")) {
    g <- ti_lmm(fs, nd, P = 0.95, sides = s, control = ctl)
    expect_identical(g$method, c("gpq", "gpq"))
    expect_identical(dg(g)$route_reason, c("one_sided", "one_sided"))
    expect_true(all(g$lower < g$estimate & g$estimate < g$upper))
    r <- ti_lmm(fs, nd, P = 0.95, sides = s, method = "reml-mls")
    expect_identical(r$method, c("reml-mls", "reml-mls"))
    expect_identical(dg(r)$route_reason, c("method", "method"))
  }
  expect_identical(dg(ti_lmm(fs, nd, P = 0.95))$route_reason, c("none", "none"))   # two-sided unchanged
})

test_that("one-sided TI on the ems-mls route (provisional): Hoffman's limit on the Francq bound", {
  skip_if_not_installed("lme4")
  fm <- fit_q(y ~ 1 + (1 | fa) + (1 | fc), balanced_nested())           # EMS (Francq)
  lo <- ti_lmm(fm, P = 0.95, level = 0.90, sides = "lower")
  expect_identical(lo$method, "ems-mls")
  e <- .design_components_fit(fm, .design_of(fm), .var_components(fm), "observable")$ems
  h <- e$df / qchisq(0.10, e$df) - 1
  U <- sum(e$k * e$ms) + sqrt(sum(h^2 * e$k^2 * e$ms^2))
  VF <- unname(.fixed(fm, NULL)$var_fix)
  expect_equal(lo$estimate - lo$lower, qnorm(0.95) * sqrt(U) + qnorm(0.90) * sqrt(VF), tolerance = 1e-12)
  fs <- fit_q(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)      # reml-mls, forced
  cp <- vfun_extract(fs, t0 = 5)
  r <- ti_lmm(fs, data.frame(Days = 5), P = 0.95, level = 0.90, sides = "upper", method = "reml-mls")
  HG <- cp$dfs[["G"]] / qchisq(0.10, cp$dfs[["G"]]) - 1; He <- cp$dfs[["e"]] / qchisq(0.10, cp$dfs[["e"]]) - 1
  Uv <- cp$V_T + sqrt((HG * cp$V_G)^2 + (He * cp$phi[["s2e"]])^2)
  expect_equal(r$upper - r$estimate, qnorm(0.95) * sqrt(Uv) + qnorm(0.90) * sqrt(cp$var_mean), tolerance = 1e-12)
  fb <- fit_q(Y ~ Month + (1 | Lot), mbl_appendix_a)                     # Montes (anova-mls)
  up <- ti_lmm(fb, data.frame(Month = 12), P = 0.99, sides = "upper")
  expect_identical(up$method, "anova-mls")
  expect_true(up$lower == -Inf && up$estimate < up$upper)
})
