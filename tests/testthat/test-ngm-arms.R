# B2: the step-13 new-group-mean (NGM) arms as lmm_interval_control() options
# (prereg H4 ngm_method on random-intercept models, H3 slope_ngm_method on
# random-slope models). Not oracles, except the SAS PROC MIXED comparison in (f)
# and the lot-line reference pins in (e). The default pins in (a) were computed
# on the B1 commit 14d3e4d (plus the B2 item-0 KR change, which does not touch
# the NGM) BEFORE any B2 edit, so they prove that no default output changed.
# GPQ tests use small M: Monte-Carlo checks, not oracles.

ngm_arm_fits <- function() {
  lc <- lme4::lmerControl(check.conv.singular = "ignore", calc.derivs = FALSE)
  fit <- function(f, d) suppressWarnings(suppressMessages(lme4::lmer(f, data = d, control = lc)))
  set.seed(20261108)
  ow <- data.frame(g = factor(rep(1:6, each = 4)))
  ow$y <- 10 + rnorm(6, sd = 2)[ow$g] + rnorm(24)
  ou <- data.frame(g = factor(rep(1:7, times = c(2, 5, 3, 4, 1, 6, 3))))
  ou$y <- 20 + rnorm(7, sd = 1.5)[ou$g] + rnorm(nrow(ou))
  ne <- expand.grid(r = 1:2, b = 1:3, a = factor(1:5))
  ne$b <- factor(paste(ne$a, ne$b))
  ne$y <- 5 + rnorm(5, sd = 1.5)[ne$a] + rnorm(15)[ne$b] + rnorm(nrow(ne))
  cr <- expand.grid(r = 1:2, a = factor(1:4), b = factor(1:6))
  cr$y <- 3 + rnorm(4)[cr$a] + rnorm(6, sd = 1.2)[cr$b] + rnorm(24, sd = 0.5)[interaction(cr$a, cr$b)] + rnorm(nrow(cr))
  oc <- data.frame(g = factor(rep(1:8, each = 3)), x = rep(c(-1, 0, 1), 8))
  oc$y <- 4 + 0.7 * oc$x + rnorm(8, sd = 1.3)[oc$g] + rnorm(24, sd = 0.8)
  ss <- lme4::sleepstudy
  d3 <- data.frame(Days = c(0, 5, 12))
  list(
    oneway_bal   = list(fit = fit(y ~ 1 + (1 | g), ow), nd = NULL),
    oneway_unbal = list(fit = fit(y ~ 1 + (1 | g), ou), nd = NULL),
    nested       = list(fit = fit(y ~ 1 + (1 | a/b), ne), nd = NULL),
    crossed      = list(fit = fit(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), cr), nd = NULL),
    fixed_slope  = list(fit = fit(Reaction ~ Days + (1 | Subject), ss), nd = d3),
    oneway_cov   = list(fit = fit(y ~ x + (1 | g), oc), nd = data.frame(x = c(-1, 0.5))),
    rs_un        = list(fit = fit(Reaction ~ Days + (1 + Days | Subject), ss), nd = d3),
    rs_ind       = list(fit = fit(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), ss), nd = d3))
}

ngm_default_pins <- 
list(oneway_bal.two.sided = list(estimate = 9.2427550562483418, 
    lower = 3.6186058835710924, upper = 14.866904228925591, df = 5, 
    method = "anova-mls", df_rule = "anova-single"), oneway_bal.lower = list(
    estimate = 9.2427550562483418, lower = 4.834051748430781, 
    upper = Inf, df = 5, method = "anova-mls", df_rule = "anova-single"), 
    oneway_bal.upper = list(estimate = 9.2427550562483418, lower = -Inf, 
        upper = 13.651458364065903, df = 5, method = "anova-mls", 
        df_rule = "anova-single"), oneway_unbal.two.sided = list(
        estimate = 19.504128856611779, lower = 16.12022477522828, 
        upper = 22.888032937995277, df = 6, method = "anova-mls", 
        df_rule = "anova-single"), oneway_unbal.lower = list(
        estimate = 19.504128856611779, lower = 16.816849509371373, 
        upper = Inf, df = 6, method = "anova-mls", df_rule = "anova-single"), 
    oneway_unbal.upper = list(estimate = 19.504128856611779, 
        lower = -Inf, upper = 22.191408203852184, df = 6, method = "anova-mls", 
        df_rule = "anova-single"), nested.two.sided = list(estimate = 4.4094252906530693, 
        lower = -1.2498296610563822, upper = 10.068680242362522, 
        df = 4, method = "ems-mls", df_rule = "top-factor"), 
    nested.lower = list(estimate = 4.4094252906530693, lower = 0.06406133614207743, 
        upper = Inf, df = 4, method = "ems-mls", df_rule = "top-factor"), 
    nested.upper = list(estimate = 4.4094252906530693, lower = -Inf, 
        upper = 8.7547892451640621, df = 4, method = "ems-mls", 
        df_rule = "top-factor"), crossed.two.sided = list(estimate = 2.9478954746320909, 
        lower = -0.71168309231212756, upper = 6.6074740415763094, 
        df = 10.940425177986535, method = "ems-mls", df_rule = "satterthwaite-vf"), 
    crossed.lower = list(estimate = 2.9478954746320909, lower = -0.03762660090435066, 
        upper = Inf, df = 10.940425177986535, method = "ems-mls", 
        df_rule = "satterthwaite-vf"), crossed.upper = list(estimate = 2.9478954746320909, 
        lower = -Inf, upper = 5.9334175501685325, df = 10.940425177986535, 
        method = "ems-mls", df_rule = "satterthwaite-vf"), fixed_slope.two.sided = list(
        estimate = c(251.40510484848488, 303.74153464646469, 
        377.0125363636364), lower = c(167.93376370382401, 220.51103655203306, 
        293.10618249633643), upper = c(334.87644599314575, 386.97203274089634, 
        460.91889023093637), df = c(17.287236039257884, 17.003533017312794, 
        17.803156426577548), method = "anova-mls", df_rule = c("anova-satterthwaite", 
        "anova-satterthwaite", "anova-satterthwaite")), fixed_slope.lower = list(
        estimate = c(251.40510484848488, 303.74153464646469, 
        377.0125363636364), lower = c(182.55915584594965, 235.11520696927161, 
        307.77113704326501), upper = c(Inf, Inf, Inf), df = c(17.287236039257884, 
        17.003533017312794, 17.803156426577548), method = "anova-mls", 
        df_rule = c("anova-satterthwaite", "anova-satterthwaite", 
        "anova-satterthwaite")), fixed_slope.upper = list(estimate = c(251.40510484848488, 
    303.74153464646469, 377.0125363636364), lower = c(-Inf, -Inf, 
    -Inf), upper = c(320.25105385102012, 372.36786232365773, 
    446.25393568400779), df = c(17.287236039257884, 17.003533017312794, 
    17.803156426577548), method = "anova-mls", df_rule = c("anova-satterthwaite", 
    "anova-satterthwaite", "anova-satterthwaite")), oneway_cov.two.sided = list(
        estimate = c(3.4123256486403002, 4.8559687518979837), 
        lower = c(-0.88368358481097209, 0.56435428602630644), 
        upper = c(7.708334882091572, 9.147583217769661), df = c(7.1864148984516181, 
        7.0464840597856444), method = "anova-mls", df_rule = c("anova-satterthwaite", 
        "anova-satterthwaite")), oneway_cov.lower = list(estimate = c(3.4123256486403002, 
    4.8559687518979837), lower = c(-0.034368778857155835, 1.4162667111042642
    ), upper = c(Inf, Inf), df = c(7.1864148984516181, 7.0464840597856444
    ), method = "anova-mls", df_rule = c("anova-satterthwaite", 
    "anova-satterthwaite")), oneway_cov.upper = list(estimate = c(3.4123256486403002, 
    4.8559687518979837), lower = c(-Inf, -Inf), upper = c(6.8590200761377567, 
    8.2956707926917037), df = c(7.1864148984516181, 7.0464840597856444
    ), method = "anova-mls", df_rule = c("anova-satterthwaite", 
    "anova-satterthwaite")), rs_un.two.sided = list(estimate = c(251.40510484848207, 
    303.74153464646298, 377.01253636363629), lower = c(193.33269655923857, 
    216.77142383198913, 207.88794299556997), upper = c(309.47751313772557, 
    390.71164546093684, 546.13712973170254), df = c(8.9848849686968428, 
    15.636602553060809, 14.379594631555205), method = "reml", 
        df_rule = c("satterthwaite", "satterthwaite", "satterthwaite"
        )), rs_un.lower = list(estimate = c(251.40510484848207, 
    303.74153464646298, 377.01253636363629), lower = c(204.34985429716522, 
    232.14906559841745, 238.04174658726356), upper = c(Inf, Inf, 
    Inf), df = c(8.9848849686968428, 15.636602553060809, 14.379594631555205
    ), method = "reml", df_rule = c("satterthwaite", "satterthwaite", 
    "satterthwaite")), rs_un.upper = list(estimate = c(251.40510484848207, 
    303.74153464646298, 377.01253636363629), lower = c(-Inf, 
    -Inf, -Inf), upper = c(298.46035539979891, 375.33400369450851, 
    515.98332614000901), df = c(8.9848849686968428, 15.636602553060809, 
    14.379594631555205), method = "reml", df_rule = c("satterthwaite", 
    "satterthwaite", "satterthwaite")), rs_ind.two.sided = list(
        estimate = c(251.40510484848457, 303.74153464646389, 
        377.01253636363492), lower = c(193.34557166664109, 220.66237091594064, 
        210.20200224320755), upper = c(309.46463803032805, 386.82069837698714, 
        543.82307048406233), df = c(9.7863609045254556, 22.940807639459496, 
        15.234369588067219), method = "reml", df_rule = c("satterthwaite", 
        "satterthwaite", "satterthwaite")), rs_ind.lower = list(
        estimate = c(251.40510484848457, 303.74153464646389, 
        377.01253636363492), lower = c(204.21224291794297, 234.91328146461825, 
        239.7709104568718), upper = c(Inf, Inf, Inf), df = c(9.7863609045254556, 
        22.940807639459496, 15.234369588067219), method = "reml", 
        df_rule = c("satterthwaite", "satterthwaite", "satterthwaite"
        )), rs_ind.upper = list(estimate = c(251.40510484848457, 
    303.74153464646389, 377.01253636363492), lower = c(-Inf, 
    -Inf, -Inf), upper = c(298.59796677902619, 372.56978782830953, 
    514.25416227039807), df = c(9.7863609045254556, 22.940807639459496, 
    15.234369588067219), method = "reml", df_rule = c("satterthwaite", 
    "satterthwaite", "satterthwaite")))

F_b2 <- ngm_arm_fits()

## muffle the documented limitation warnings only
quiet <- function(expr) withCallingHandlers(expr,
  reintervals_warning_limitation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))
## ... and the unvalidated-arm warning (tested separately in (g))
quiet_arm <- function(expr) withCallingHandlers(expr,
  reintervals_warning_unvalidated = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_limitation  = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
  reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))
ngm <- function(fx, ..., sides = "two.sided")
  quiet_arm(new_group_mean_lmm(fx$fit, newdata = fx$nd, sides = sides,
                               control = lmm_interval_control(...)))
n_unvalidated <- function(expr) {
  n <- 0L
  withCallingHandlers(expr,
    reintervals_warning_unvalidated = function(w) { n <<- n + 1L; invokeRestart("muffleWarning") },
    reintervals_warning_limitation  = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_extrapolation = function(w) invokeRestart("muffleWarning"),
    reintervals_warning_df_fallback   = function(w) invokeRestart("muffleWarning"))
  n
}
## an interval without its arm fields (which record the option, not the numbers)
drop_arm <- function(r) {
  d <- attr(r, "diagnostics"); d$ngm_arm <- NULL; d$slope_ngm_arm <- NULL
  out <- as.data.frame(r); attr(out, "diagnostics") <- d; out
}
ow_parts <- function(fm) {                         # one-way pieces off a fit
  g <- lme4::getME(fm, "flist")[[1]]; y <- lme4::getME(fm, "y")
  vc <- as.data.frame(lme4::VarCorr(fm))
  list(I = nlevels(g), n_i = as.numeric(table(g)), gm = as.numeric(tapply(y, g, mean)),
       mu = unname(lme4::fixef(fm)[1]), VF = as.numeric(stats::vcov(fm)[1, 1]),
       s2g = vc$vcov[1], s2e = vc$vcov[2])
}

test_that("B2 (a): default NGM outputs are unchanged (pinned before the B2 edits)", {
  for (key in names(ngm_default_pins)) {
    nm <- sub("\\..*$", "", key); sd <- sub("^[^.]*\\.", "", key)
    p <- ngm_default_pins[[key]]
    r <- quiet(new_group_mean_lmm(F_b2[[nm]]$fit, newdata = F_b2[[nm]]$nd, sides = sd))
    d <- attr(r, "diagnostics")
    expect_equal(r$estimate, p$estimate, tolerance = 1e-12, info = key)
    expect_equal(r$lower, p$lower, tolerance = 1e-12, info = key)
    expect_equal(r$upper, p$upper, tolerance = 1e-12, info = key)
    expect_equal(d$df, p$df, tolerance = 1e-12, info = key)
    expect_identical(unique(r$method), p$method, info = key)
    expect_identical(d$df_rule, p$df_rule, info = key)
    rs <- grepl("^rs_", nm)
    expect_identical(unique(d$ngm_arm), if (rs) "not_applicable" else "package", info = key)
    expect_identical(unique(d$slope_ngm_arm), if (rs) "reml" else "not_applicable", info = key)
  }
})

test_that("B2 (b): broad-inference recomputes from W and its Satterthwaite df", {
  ## one-way, unbalanced: Var(mu_hat) = 1 / sum n_i / (s2e + n_i s2g), whose
  ## gradient in (s2g, s2e) is written here in closed form (independent of the
  ## package's numerical gradient); A = the observed-information covariance.
  fm <- F_b2$oneway_unbal$fit; p <- ow_parts(fm)
  d_i <- p$s2e + p$n_i * p$s2g
  VF <- 1 / sum(p$n_i / d_i)
  g <- VF^2 * c(sum(p$n_i^2 / d_i^2), sum(p$n_i / d_i^2)) + c(1, 0)
  A <- .theta_cov(fm)$cov_theta
  W <- p$VF + p$s2g
  nu <- 2 * W^2 / as.numeric(t(g) %*% A %*% g)
  r <- ngm(F_b2$oneway_unbal, ngm_method = "broad-inference")
  d <- attr(r, "diagnostics")
  expect_equal(VF, p$VF, tolerance = 1e-6)               # lme4's vcov at the optimum
  expect_equal(d$W, W, tolerance = 1e-12)
  expect_equal(d$df, nu, tolerance = 1e-6)
  expect_false(d$ngm_df_failed)
  expect_identical(d$df_rule, "satterthwaite-broad")
  hw <- stats::qt(0.975, d$df) * sqrt(W)
  expect_equal(c(r$lower, r$upper), p$mu + c(-1, 1) * hw, tolerance = 1e-12)
  ## one-sided: t_{1 - alpha} on the same W and df
  lo <- ngm(F_b2$oneway_unbal, ngm_method = "broad-inference", sides = "lower")
  expect_equal(lo$lower, p$mu - stats::qt(0.95, d$df) * sqrt(W), tolerance = 1e-12)
  expect_identical(lo$upper, Inf)
  ## nested and crossed: W = V_F + every between-group variance component
  for (nm in c("nested", "crossed", "fixed_slope", "oneway_cov")) {
    fx <- F_b2[[nm]]
    r <- ngm(fx, ngm_method = "broad-inference")
    vc <- as.data.frame(lme4::VarCorr(fx$fit))
    nd <- if (is.null(fx$nd)) stats::model.frame(fx$fit)[1, , drop = FALSE] else fx$nd
    L <- stats::model.matrix(stats::delete.response(stats::terms(fx$fit, fixed.only = TRUE)), nd)
    VF <- diag(L %*% as.matrix(stats::vcov(fx$fit)) %*% t(L))
    expect_equal(attr(r, "diagnostics")$W, unname(VF + sum(vc$vcov[vc$grp != "Residual"])),
                 tolerance = 1e-12, info = nm)
    expect_true(all(is.finite(attr(r, "diagnostics")$df)), info = nm)
  }
})

test_that("B2 (b): broad-inference gives NA limits, no fallback, when its df fails", {
  ## a singular one-way fit: the observed-information covariance is unavailable
  set.seed(3)
  d <- data.frame(g = factor(rep(1:5, each = 4)))
  d$y <- rnorm(20)
  fm <- suppressMessages(suppressWarnings(lme4::lmer(y ~ 1 + (1 | g), data = d)))
  skip_if_not(lme4::isSingular(fm), "fixture fit is not singular")
  r <- quiet_arm(new_group_mean_lmm(fm, control = lmm_interval_control(ngm_method = "broad-inference")))
  expect_true(is.na(r$lower) && is.na(r$upper))
  expect_true(is.finite(r$estimate))
  expect_true(attr(r, "diagnostics")$ngm_df_failed)
})

test_that("B2 (b): hts and classical recompute from their formulas; classical == 2b on balanced data", {
  fx <- F_b2$oneway_bal; p <- ow_parts(fx$fit)
  h <- ngm(fx, ngm_method = "hts")
  expect_equal(c(h$lower, h$upper),
               p$mu + c(-1, 1) * stats::qt(0.975, p$I - 2) * sqrt(p$s2g + p$VF), tolerance = 1e-12)
  expect_identical(attr(h, "diagnostics")$df, p$I - 2)
  hu <- ngm(fx, ngm_method = "hts", sides = "upper")
  expect_equal(hu$upper, p$mu + stats::qt(0.95, p$I - 2) * sqrt(p$s2g + p$VF), tolerance = 1e-12)
  ## classical (Hahn & Meeker): Ybar +/- t_{I-1} S_L sqrt(1 + 1/I) on the group means
  cl <- ngm(fx, ngm_method = "classical")
  hw <- stats::qt(0.975, p$I - 1) * stats::sd(p$gm) * sqrt(1 + 1 / p$I)
  expect_equal(c(cl$lower, cl$upper), mean(p$gm) + c(-1, 1) * hw, tolerance = 1e-12)
  ## identity with the package's 2b (anova-single) on balanced one-way data
  for (sd in c("two.sided", "lower", "upper")) {
    a <- ngm(fx, ngm_method = "classical", sides = sd)
    b <- ngm(fx, sides = sd)
    expect_equal(c(a$estimate, a$lower, a$upper), c(b$estimate, b$lower, b$upper),
                 tolerance = 1e-12, info = sd)
    expect_identical(attr(a, "diagnostics")$df, attr(b, "diagnostics")$df)
  }
})

test_that("B2 (c): an arm not defined for the model or design is a no-op, with no warning", {
  cases <- list(
    list("nested",       list(ngm_method = "hts")),
    list("crossed",      list(ngm_method = "classical")),
    list("oneway_unbal", list(ngm_method = "classical")),
    list("oneway_cov",   list(ngm_method = "hts")),
    list("fixed_slope",  list(ngm_method = "classical")),
    list("oneway_bal",   list(slope_ngm_method = "gpq")),
    list("oneway_unbal", list(slope_ngm_method = "broad-inference")),
    list("nested",       list(slope_ngm_method = "broad-inference")),
    list("fixed_slope",  list(slope_ngm_method = "broad-inference")),
    list("fixed_slope",  list(slope_ngm_method = "fv-anova")),
    list("rs_un",        list(ngm_method = "broad-inference")),
    list("rs_ind",       list(ngm_method = "hts")))
  for (cs in cases) {
    fx <- F_b2[[cs[[1]]]]; info <- paste(cs[[1]], names(cs[[2]]), cs[[2]][[1]])
    ctl <- do.call(lmm_interval_control, cs[[2]])
    expect_identical(n_unvalidated(r <- new_group_mean_lmm(fx$fit, newdata = fx$nd, control = ctl)), 0L,
                     info = info)
    b <- quiet(new_group_mean_lmm(fx$fit, newdata = fx$nd))
    d <- attr(r, "diagnostics")
    expect_equal(drop_arm(r), drop_arm(b), tolerance = 0, info = info)
    opt <- names(cs[[2]]); field <- if (opt == "ngm_method") "ngm_arm" else "slope_ngm_arm"
    expect_identical(unique(d[[field]]), "not_applicable", info = info)
  }
  ## and the NGM options are no-ops on the other verbs
  ctl <- lmm_interval_control(ngm_method = "broad-inference", slope_ngm_method = "gpq")
  for (nm in c("oneway_bal", "rs_un")) {
    fx <- F_b2[[nm]]
    expect_identical(n_unvalidated(a <- ci_lmm(fx$fit, newdata = fx$nd, control = ctl)), 0L)
    expect_equal(as.data.frame(a), as.data.frame(quiet(ci_lmm(fx$fit, newdata = fx$nd))), tolerance = 0)
    ## (CI rows record the NGM fields as not_applicable under any control)
    expect_identical(unique(attr(a, "diagnostics")$ngm_arm), "not_applicable")
    expect_identical(unique(attr(a, "diagnostics")$slope_ngm_arm), "not_applicable")
  }
})

test_that("B2 (d): GPQ NGM leaves the TI pivots unchanged and uses the A2 machinery", {
  ## Z_new is drawn after U, Q, Z: the TI draws are identical with or without it
  a <- .gpq_pivots(df_resid = 140, K = 18, M = 2000, seed = 11L)
  b <- .gpq_pivots(df_resid = 140, K = 18, M = 2000, seed = 11L, z_new = TRUE)
  expect_identical(a[c("U", "Q", "Z", "seed")], b[c("U", "Q", "Z", "seed")])
  expect_null(a$Z_new); expect_length(b$Z_new, 2000L)

  fit <- F_b2$rs_un$fit
  g <- function(t0, seed = 5L, sides = "two.sided", M = 2000L)
    quiet_arm(new_group_mean_lmm(fit, newdata = data.frame(Days = t0), sides = sides,
              control = lmm_interval_control(slope_ngm_method = "gpq", seed = seed, M = M)))
  ## common random numbers: a t0's limits do not depend on the other t0
  all3 <- g(c(0, 5, 12)); one <- g(5)
  expect_identical(c(all3$lower[2], all3$upper[2]), c(one$lower, one$upper))
  ## seed contract: the recorded gpq_seed reproduces the call
  set.seed(99); r1 <- g(5, seed = NULL)
  d1 <- attr(r1, "diagnostics")
  r2 <- g(5, seed = d1$gpq_seed)
  expect_identical(c(r1$lower, r1$upper), c(r2$lower, r2$upper))
  expect_identical(d1$rng_kind, RNGkind())
  expect_identical(d1$mc_sections, 20L)
  expect_true(all(is.finite(c(d1$mc_se_lower, d1$mc_se_upper)) & c(d1$mc_se_lower, d1$mc_se_upper) > 0))
  expect_identical(d1$n_na_draws, 0L)
  expect_identical(d1$slope_ngm_arm, "gpq")
  ## one-sided: the alpha / 1 - alpha quantiles, the other limit infinite
  lo <- g(5, sides = "lower"); up <- g(5, sides = "upper")
  expect_identical(lo$upper, Inf); expect_identical(up$lower, -Inf)
  expect_true(lo$lower > one$lower && up$upper < one$upper)
  expect_true(is.na(attr(lo, "diagnostics")$mc_se_upper))
  ## M >= 2000 rule (the control) and the per-lot time-point guard
  expect_error(lmm_interval_control(slope_ngm_method = "gpq", M = 1000), class = "reintervals_error_input")
  ss <- lme4::sleepstudy; ss <- ss[!(ss$Subject == "308" & ss$Days > 0), ]
  f1 <- suppressMessages(suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = ss)))
  for (arm in c("gpq", "fv-anova")) {
    expect_error(quiet_arm(new_group_mean_lmm(f1, newdata = data.frame(Days = 5),
                 control = lmm_interval_control(slope_ngm_method = arm, seed = 1L, M = 2000))),
                 class = "reintervals_error_time_points", info = arm)
  }
})

test_that("B2 (d): GPQ NGM is within 15% (half-width) of the reml NGM on sleepstudy UN at M = 20000", {
  ## large-K asymptotics: both target the same new-lot mean; a sanity check,
  ## not an oracle (K = 18)
  skip_on_cran()
  nd <- data.frame(Days = 5)
  a <- quiet_arm(new_group_mean_lmm(F_b2$rs_un$fit, newdata = nd,
                 control = lmm_interval_control(slope_ngm_method = "gpq", seed = 20261108L)))
  b <- quiet(new_group_mean_lmm(F_b2$rs_un$fit, newdata = nd))
  hw <- (b$upper - b$lower) / 2
  expect_lt(abs(a$lower - b$lower), 0.15 * hw)
  expect_lt(abs(a$upper - b$upper), 0.15 * hw)
})

## Reference pins for fv-anova: lotline() of the simstudy pilot,
## ~/dev/reintervals-simstudy pilot/methods.R at commit f4ce68e (only the
## lotline() definition evaluated; the file is not sourced here), run on
## 2026-11-08. Fixtures: balanced sleepstudy; sleepstudy with Days > 5 dropped
## for subjects 308, 309, 310, 330, 331, 332; stab_unbal_v2_a.csv (simstudy
## oracles/sas_ngm/part2, copied to fixtures/sas_ngm_part2).
lotline_pins <- list(
  list(fixture = "ss_bal", t0 = 0, Ybar = 251.4051048484848, S2L = 838.34229346931807),
  list(fixture = "ss_bal", t0 = 5, Ybar = 303.74153464646457, S2L = 1652.4035567741391),
  list(fixture = "ss_bal", t0 = 15, Ybar = 408.41439424242418, S2L = 9732.0770255355601),
  list(fixture = "ss_unbal", t0 = 0, Ybar = 250.29596049062047, S2L = 970.07180249464454),
  list(fixture = "ss_unbal", t0 = 5, Ybar = 305.10737590668589, S2L = 2126.3709244837728),
  list(fixture = "ss_unbal", t0 = 15, Ybar = 414.73020673881672, S2L = 17304.888176338911),
  list(fixture = "v2a", t0 = 0, Ybar = 100.08906479362707, S2L = 1.7139804552481497),
  list(fixture = "v2a", t0 = 12, Ybar = 97.032951671300268, S2L = 1.3852426512589382),
  list(fixture = "v2a", t0 = 48, Ybar = 87.864612304319834, S2L = 12.741702923289127))

lotline_data <- function() {
  ss <- lme4::sleepstudy
  ub <- ss[!(ss$Subject %in% c("308", "309", "310", "330", "331", "332") & ss$Days > 5), ]
  v2 <- utils::read.csv(test_path("fixtures", "sas_ngm_part2", "stab_unbal_v2_a.csv"),
                        stringsAsFactors = TRUE)
  list(ss_bal = list(y = ss$Reaction, t = ss$Days, g = ss$Subject, d = ss, f = Reaction ~ Days + (1 + Days | Subject), tv = "Days"),
       ss_unbal = list(y = ub$Reaction, t = ub$Days, g = ub$Subject, d = ub, f = Reaction ~ Days + (1 + Days | Subject), tv = "Days"),
       v2a = list(y = v2$y, t = v2$month, g = v2$lot, d = v2, f = y ~ month + (1 + month | lot), tv = "month"))
}

# oracle: simstudy-lotline-f4ce68e (fixtures/oracles.csv)
test_that("B2 (e): fv-anova reproduces the pilot lotline() pins and recomputes from its formula", {
  D <- lotline_data()
  for (p in lotline_pins) {
    x <- D[[p$fixture]]; info <- paste(p$fixture, p$t0)
    L <- .lotline(x$y, x$t, x$g, p$t0)
    expect_equal(L$Ybar, p$Ybar, tolerance = 1e-10, info = info)
    expect_equal(L$S2L, p$S2L, tolerance = 1e-10, info = info)
  }
  for (nm in names(D)) {
    x <- D[[nm]]
    fm <- suppressMessages(suppressWarnings(lme4::lmer(x$f, data = x$d)))
    t0 <- vapply(Filter(function(p) p$fixture == nm, lotline_pins), function(p) p$t0, 0)
    nd <- stats::setNames(data.frame(t0), x$tv)
    pin <- Filter(function(p) p$fixture == nm, lotline_pins)
    K <- nlevels(factor(x$g))
    for (sd in c("two.sided", "lower", "upper")) {
      r <- quiet_arm(new_group_mean_lmm(fm, newdata = nd, sides = sd,
                     control = lmm_interval_control(slope_ngm_method = "fv-anova")))
      Yb <- vapply(pin, function(p) p$Ybar, 0); W <- vapply(pin, function(p) p$S2L, 0) * (1 + 1 / K)
      tq <- if (sd == "two.sided") stats::qt(0.975, K - 1) else stats::qt(0.95, K - 1)
      expect_equal(r$estimate, Yb, tolerance = 1e-10, info = paste(nm, sd))
      if (sd != "upper") expect_equal(r$lower, Yb - tq * sqrt(W), tolerance = 1e-10, info = paste(nm, sd))
      if (sd != "lower") expect_equal(r$upper, Yb + tq * sqrt(W), tolerance = 1e-10, info = paste(nm, sd))
      d <- attr(r, "diagnostics")
      expect_identical(unique(d$df), K - 1); expect_identical(unique(d$df_rule), "k-minus-1")
    }
  }
})

## ---- (f) SAS PROC MIXED DDFM=SATTERTH, ##NGM rows (prereg 7 item 9) ----------
## Tolerances: SE 1e-4 relative, df 0.05 absolute.
## broad-inference: all 13 random-intercept rows pass (OW_BAL, OW_UNBAL,
##   NEST_BAL, FS_BAL x 5, FS_UNBAL x 5).
## reml (random-slope default, NOT tuned): no row passes.
##   RS_UN_BAL (5 rows), RS_VC_BAL (5 rows): SE within 1e-4 except RS_VC_BAL at
##   month 36 and 48 (1.1e-4, 1.2e-4; lme4's REML slope variance differs from
##   SAS's by 3e-4 relative on a flat surface); df 0.13-0.40 below SAS. Cause:
##   the package's nu_G is the Satterthwaite df of V_G alone, while SAS's is the
##   df of the whole prediction variance V_F + V_G (gradient of V_F plus that of
##   V_G); with that gradient the package's own Hessian reproduces SAS to 1e-3
##   (checked below).
##   RS_UN_UNBAL (5 rows): singular fit (SAS UN(2,2) = 0); SAS reports no SE at
##   months 12-48 and the package uses its containment df; boundary case only.
sas_ngm_rows <- function() {
  res <- readLines(test_path("fixtures", "sas_ngm", "ngm_sas_check_satterth.results.txt"))
  m <- grep("^##NGM +[A-Z_]+ +SATTERTH ", res, value = TRUE)
  num <- function(key, x) {
    hit <- grepl(paste0(key, "= *[-0-9]"), x); out <- rep(NA_real_, length(x))
    out[hit] <- as.numeric(sub(paste0(".*", key, "= *([-0-9.eE]+).*"), "\\1", x[hit])); out
  }
  data.frame(fit = sub("^##NGM +([A-Z_]+) .*", "\\1", m), t0 = num("month", m),
             se = num("SE", m), df = num("DF", m), stringsAsFactors = FALSE)
}
sas_ngm_specs <- list(
  OW_BAL      = list(y ~ 1 + (1 | grp),                    "oneway_bal.csv"),
  OW_UNBAL    = list(y ~ 1 + (1 | grp),                    "oneway_unbal.csv"),
  NEST_BAL    = list(y ~ 1 + (1 | a) + (1 | a:b),          "nested_bal.csv"),
  FS_BAL      = list(y ~ month + (1 | lot),                "stab_bal.csv"),
  FS_UNBAL    = list(y ~ month + (1 | lot),                "stab_unbal.csv"),
  RS_UN_BAL   = list(y ~ month + (1 + month | lot),        "stab_bal.csv"),
  RS_VC_BAL   = list(y ~ month + (1 | lot) + (0 + month | lot), "stab_bal.csv"))
sas_ngm_fit <- function(tag) {
  d <- utils::read.csv(test_path("fixtures", "sas_ngm", sas_ngm_specs[[tag]][[2]]), stringsAsFactors = TRUE)
  suppressMessages(lme4::lmer(sas_ngm_specs[[tag]][[1]], data = d, REML = TRUE))
}
## package values (pins) for the reml rows, which fail the SAS tolerance
reml_rs_pins <- list(
  RS_UN_BAL = list(se = c(0.98796391, 0.71759462, 0.95068522, 1.46378065, 2.05706008),
                   df = c(4.3894657, 4.9485572, 3.4016469, 2.9568377, 2.8698415)),
  RS_VC_BAL = list(se = c(0.79751077, 0.91554433, 1.23727375, 1.64739133, 2.09460664),
                   df = c(3.7359849, 4.3707037, 3.3617061, 2.5888613, 2.2480654)))

# oracle: sas-part1-ngm-ri-broad (fixtures/oracles.csv)
test_that("B2 (f): broad-inference reproduces SAS DDFM=SATTERTH on every random-intercept ##NGM row", {
  sas <- sas_ngm_rows()
  for (tag in c("OW_BAL", "OW_UNBAL", "NEST_BAL", "FS_BAL", "FS_UNBAL")) {
    s <- sas[sas$fit == tag, ]; fm <- sas_ngm_fit(tag)
    nd <- if (all(is.na(s$t0))) NULL else data.frame(month = s$t0)
    r <- quiet_arm(new_group_mean_lmm(fm, newdata = nd,
                   control = lmm_interval_control(ngm_method = "broad-inference")))
    d <- attr(r, "diagnostics")
    expect_true(all(abs(d$se / s$se - 1) < 1e-4), info = tag)
    expect_true(all(abs(d$df - s$df) < 0.05), info = tag)
  }
})

# oracle: sas-part1-ngm-rs-reml (fixtures/oracles.csv)
test_that("B2 (f): reml random-slope NGM vs SAS: pinned, SE within tolerance, df gap diagnosed", {
  sas <- sas_ngm_rows()
  for (tag in names(reml_rs_pins)) {
    s <- sas[sas$fit == tag, ]; fm <- sas_ngm_fit(tag)
    r <- quiet(new_group_mean_lmm(fm, newdata = data.frame(month = s$t0)))
    d <- attr(r, "diagnostics")
    expect_equal(d$se, reml_rs_pins[[tag]]$se, tolerance = 1e-6, info = tag)
    expect_equal(d$df, reml_rs_pins[[tag]]$df, tolerance = 1e-6, info = tag)
    ok_se <- if (tag == "RS_VC_BAL") 1:3 else 1:5      # SE rows inside 1e-4
    expect_true(all(abs(d$se[ok_se] / s$se[ok_se] - 1) < 1e-4), info = tag)
    ## diagnosis: the Satterthwaite df of V_F + V_G (not of V_G alone) is SAS's df
    parsed <- .vfun_parse(fm); cv <- .vfun_vc_cov(parsed)$cov; phi <- parsed$phi
    nu_W <- vapply(s$t0, function(t0) {
      x0 <- stats::setNames(c(1, t0), c("(Intercept)", parsed$fe_time))
      VF <- as.numeric(t(x0) %*% parsed$vcov_beta %*% x0)
      VG <- .vfun_value_VG(parsed$type, phi, t0)
      gF <- numDeriv::grad(function(p) parsed$vF_fun(p, x0), phi)
      gG <- .vfun_grad_VG(parsed$type, t0)[names(phi)]
      .vfun_df_satterthwaite(VF + VG, gF + gG, cv)
    }, 0)
    expect_true(all(abs(nu_W - s$df) < 0.05), info = tag)
  }
})

## ---- B2b: slope_ngm_method = "broad-inference" (random slope) ---------------
## SAS DDFM=SATTERTH broad inference of a new lot: the default "reml" W =
## V_F + V_G(t0), with the Satterthwaite df of the whole W (the B2 (f)
## diagnosis, made an arm). Package vs SAS on RS_UN_BAL, RS_VC_BAL (all t0):
##   RS_UN_BAL: SE within 5.2e-6 relative, df within 5.4e-5: pass x 5.
##   RS_VC_BAL: df within 4.2e-4: pass x 5; SE 1.1e-5, 4.7e-5, 8.8e-5 (pass)
##     and 1.08e-4, 1.18e-4 at months 36 and 48 (fail, > 1e-4).
## The same W and df formulas evaluated AT the SAS covariance-parameter
## estimates (##COVP lines; A = inverse Hessian of the REML NLL at those
## estimates, same X and Z): SE within 1.4e-12 relative, df within 1.1e-4 on all
## ten rows. The RS_VC_BAL SE gap is therefore lme4's optimum (slope variance
## 0.00135877 vs SAS 0.00135838), not the method.
broad_rs_pins <- list(
  RS_UN_BAL = list(se = c(0.987963915, 0.717594616, 0.950685223, 1.46378065, 2.05706008),
                   df = c(4.6520663, 5.1607579, 3.741319, 3.324586, 3.2424758)),
  RS_VC_BAL = list(se = c(0.797510767, 0.915544326, 1.23727375, 1.64739133, 2.09460664),
                   df = c(4.131501, 4.5031842, 3.5829781, 2.8717316, 2.5643233)))
sas_covp <- function(tag) {
  res <- readLines(test_path("fixtures", "sas_ngm", "ngm_sas_check_satterth.results.txt"))
  m <- grep(paste0("^##COVP +", tag, " +SATTERTH "), res, value = TRUE)
  est <- as.numeric(sub(".*Estimate= *([-0-9.eE]+).*", "\\1", m))
  if (tag == "RS_UN_BAL") c(s2_0 = est[1], s01 = est[2], s2_1 = est[3], s2e = est[4])
  else c(s2_0 = est[1], s2_1 = est[2], s2e = est[3])
}

# oracle: sas-part1-ngm-rs-broad (fixtures/oracles.csv)
test_that("B2b: random-slope broad-inference vs SAS DDFM=SATTERTH: pinned, passing rows, and at the SAS estimates", {
  sas <- sas_ngm_rows()
  for (tag in names(broad_rs_pins)) {
    s <- sas[sas$fit == tag, ]; fm <- sas_ngm_fit(tag)
    r <- quiet_arm(new_group_mean_lmm(fm, newdata = data.frame(month = s$t0),
                   control = lmm_interval_control(slope_ngm_method = "broad-inference")))
    d <- attr(r, "diagnostics")
    expect_equal(d$se, broad_rs_pins[[tag]]$se, tolerance = 1e-6, info = tag)
    expect_equal(d$df, broad_rs_pins[[tag]]$df, tolerance = 1e-6, info = tag)
    expect_true(all(abs(d$df - s$df) < 0.05), info = tag)
    ok_se <- if (tag == "RS_VC_BAL") 1:3 else 1:5      # SE rows inside 1e-4
    expect_true(all(abs(d$se[ok_se] / s$se[ok_se] - 1) < 1e-4), info = tag)
    ## the same formulas at SAS's covariance-parameter estimates reproduce SAS
    parsed <- .vfun_parse(fm); phiS <- sas_covp(tag)
    covS <- solve(numDeriv::hessian(parsed$reml_nll, phiS))
    atS <- vapply(s$t0, function(t0) .ngm_rs_broad_df(parsed, phiS, covS, t0), c(W = 0, df = 0))
    expect_true(all(abs(sqrt(atS["W", ]) / s$se - 1) < 1e-8), info = tag)
    expect_true(all(abs(atS["df", ] - s$df) < 1e-3), info = tag)
  }
})

test_that("B2b: random-slope broad-inference recomputes from its formula; SE == reml, df >= reml", {
  for (nm in c("rs_un", "rs_ind")) {
    fx <- F_b2[[nm]]; fm <- fx$fit
    for (sd in c("two.sided", "lower", "upper")) {
      r <- ngm(fx, slope_ngm_method = "broad-inference", sides = sd)
      b <- quiet(new_group_mean_lmm(fm, newdata = fx$nd, sides = sd))
      d <- attr(r, "diagnostics"); db <- attr(b, "diagnostics")
      expect_identical(d$se, db$se, info = paste(nm, sd))
      expect_identical(r$estimate, b$estimate, info = paste(nm, sd))
      expect_identical(unique(d$df_rule), "satterthwaite-broad")
      expect_identical(unique(d$slope_ngm_arm), "broad-inference")
      expect_false(any(d$ngm_df_failed))
      ## formula: 2 W^2 / (g' A g), g = grad V_F + grad V_G, A = .vfun_vc_cov
      parsed <- .vfun_parse(fm); A <- .vfun_vc_cov(parsed)$cov; phi <- parsed$phi
      nu <- vapply(fx$nd$Days, function(t0) {
        x0 <- c(1, t0)
        W <- as.numeric(t(x0) %*% parsed$vcov_beta %*% x0) + .vfun_value_VG(parsed$type, phi, t0)
        g <- numDeriv::grad(function(p) parsed$vF_fun(p, x0), phi) +
          .vfun_grad_VG(parsed$type, t0)[names(phi)]
        2 * W^2 / as.numeric(t(g) %*% A %*% g)
      }, 0)
      expect_equal(d$df, nu, tolerance = 1e-10, info = paste(nm, sd))
      tq <- if (sd == "two.sided") stats::qt(0.975, nu) else stats::qt(0.95, nu)
      if (sd != "upper") expect_equal(r$lower, r$estimate - tq * d$se, tolerance = 1e-12)
      if (sd != "lower") expect_equal(r$upper, r$estimate + tq * d$se, tolerance = 1e-12)
    }
  }
  ## df >= reml df on the SAS fixtures (the whole-W df adds the V_F information)
  for (tag in names(broad_rs_pins))
    expect_true(all(broad_rs_pins[[tag]]$df >= reml_rs_pins[[tag]]$df), info = tag)
})

test_that("B2b: random-slope broad-inference: NA limits when the df fails", {
  fx <- F_b2$rs_un
  local_mocked_bindings(.vfun_df_satterthwaite = function(value, grad_vec, vc_cov) NA_real_)
  r <- ngm(fx, slope_ngm_method = "broad-inference")
  expect_true(all(is.na(r$lower) & is.na(r$upper)))
  expect_true(all(is.finite(r$estimate)))
  expect_true(all(attr(r, "diagnostics")$ngm_df_failed))
})

test_that("B2b: random-slope broad-inference accepts a lot with one time point, as the default does", {
  ss <- lme4::sleepstudy
  ss <- ss[!(ss$Subject == "308" & ss$Days > 0), ]           # lot 308: Day 0 only
  fm <- suppressMessages(lme4::lmer(Reaction ~ Days + (1 + Days | Subject), data = ss))
  nd <- data.frame(Days = c(0, 5))
  b <- quiet(new_group_mean_lmm(fm, newdata = nd))
  r <- quiet_arm(new_group_mean_lmm(fm, newdata = nd,
                 control = lmm_interval_control(slope_ngm_method = "broad-inference")))
  expect_identical(attr(r, "diagnostics")$se, attr(b, "diagnostics")$se)
  expect_true(all(is.finite(c(r$lower, r$upper))))
  expect_error(quiet_arm(new_group_mean_lmm(fm, newdata = nd,
                 control = lmm_interval_control(slope_ngm_method = "fv-anova"))),
               class = "reintervals_error_time_points")
})

test_that("B2 (g): the unvalidated warning fires once per call for each applied non-default arm only", {
  for (a in c("broad-inference", "hts", "classical"))
    expect_identical(n_unvalidated(new_group_mean_lmm(F_b2$oneway_bal$fit,
      control = lmm_interval_control(ngm_method = a))), 1L, info = a)
  for (a in c("broad-inference", "gpq", "fv-anova"))
    expect_identical(n_unvalidated(new_group_mean_lmm(F_b2$rs_un$fit, newdata = data.frame(Days = c(0, 5)),
      control = lmm_interval_control(slope_ngm_method = a, seed = 1L, M = 2000))), 1L, info = a)
  ## defaults never warn
  expect_identical(n_unvalidated(new_group_mean_lmm(F_b2$oneway_bal$fit)), 0L)
  expect_identical(n_unvalidated(new_group_mean_lmm(F_b2$rs_un$fit, newdata = data.frame(Days = 5))), 0L)
  ## the warning names the option and the default
  w <- tryCatch(new_group_mean_lmm(F_b2$oneway_bal$fit, control = lmm_interval_control(ngm_method = "hts")),
                reintervals_warning_unvalidated = function(w) w)
  expect_identical(w$option, "ngm_method"); expect_identical(w$value, "hts")
  expect_match(conditionMessage(w), "ngm_method = \"package\"", fixed = TRUE)
})

## ---- PENDING oracle: chow-shao1991-t23 (fixtures/oracles.csv) -----------------
## Chow & Shao (1991), Biometrics 47(3):1071-1079: the stability data and the
## exact t_23 limits. To collect: transcribe the data and the printed limits to
## fixtures/chow_shao1991/ (data.csv: batch, time, y; limits.csv: t0, side,
## estimate, limit, dp; `#` header: reference, table, transcriber, date), then
## compare slope_ngm_method = "fv-anova" (mean of the per-lot OLS lines; df
## K - 1, which is 23 only if the data have 24 batches: check, and record any
## other df in the notes) at each printed t0: estimate and limit at the printed
## precision, and the df exactly (23).
## test_that("chow-shao1991-t23: fv-anova reproduces the Chow & Shao (1991) limits", {
##   fx <- oracle_fixture("chow-shao1991-t23")
##   fm <- lme4::lmer(y ~ time + (1 + time | batch), data = fx$data)
##   for (sd in unique(fx$limits$side)) {
##     L <- fx$limits[fx$limits$side == sd, ]
##     r <- quiet_arm(new_group_mean_lmm(fm, newdata = data.frame(time = L$t0), sides = sd,
##                    control = lmm_interval_control(slope_ngm_method = "fv-anova")))
##     expect_identical(unique(attr(r, "diagnostics")$df), 23)
##     expect_equal(round(r$estimate, L$dp), L$estimate)
##     expect_equal(round(if (sd == "lower") r$lower else r$upper, L$dp), L$limit)
##   }
## })

## ---- PENDING oracle: sas-part2-ngm-rs-unbal (fixtures/oracles.csv) ------------
## SAS part 2 (simstudy oracles/sas_ngm/part2, ngm_sas_check_part2.sas; data
## stab_unbal_v2_a/b/c.csv): unbalanced random-slope fits, once the results file
## exists. To collect: copy ngm_sas_check_part2.results.txt and the three CSVs to
## fixtures/sas_ngm_part2/ unchanged (add their hashes to the README and
## SHA256SUMS), then compare on every ##NGM SATTERTH row, as part 1 (B2 (f) /
## B2b above): slope_ngm_method = "broad-inference" SE 1e-4 relative and df 0.05
## absolute; the default "reml" SE 1e-4 relative, its df pinned with the gap
## diagnosed as in part 1; ##COVP vs lme4::VarCorr as a fit-agreement gate
## (1e-4 relative). Rows on a singular SAS fit (SE ".") are reported, not compared.
## test_that("sas-part2-ngm-rs-unbal: random-slope NGM vs SAS part 2", {
##   sas <- oracle_fixture("sas-part2-ngm-rs-unbal")
##   s <- sas$ngm[sas$ngm$ddfm == "SATTERTH" & !is.na(sas$ngm$se), ]
##   for (tag in unique(s$fit)) {
##     st <- s[s$fit == tag, ]
##     fm <- sas_part2_fit(tag)                   # tag -> csv + y ~ month + (1 + month | lot)
##     r <- quiet_arm(new_group_mean_lmm(fm, newdata = data.frame(month = st$t0),
##                    control = lmm_interval_control(slope_ngm_method = "broad-inference")))
##     d <- attr(r, "diagnostics")
##     expect_true(all(abs(d$se / st$se - 1) < 1e-4), info = tag)
##     expect_true(all(abs(d$df - st$df) < 0.05), info = tag)
##   }
## })

## ---- PENDING oracle: sas-part3-ngm-designs (fixtures/oracles.csv) -------------
## SAS part 3 (simstudy bcf4661, oracles/sas_ngm/part3, ngm_sas_check_part3.sas;
## cr_bal, cr_unbal, nest_unbal, ow_covw, ow_covb, nest_covw, cr_covw), once the
## results file exists. To collect: copy the results file and the CSVs to
## fixtures/sas_ngm_part3/ unchanged (README with hashes; SHA256SUMS), then
## compare ngm_method = "broad-inference" on every ##NGM SATTERTH row: SE 1e-4
## relative, df 0.05 absolute (as part 1, B2 (f)); ##MEAN SATTERTH vs ci_lmm()
## with the tolerances of test-oracle-sas-ci.R; ##COVP vs lme4::VarCorr (1e-4
## relative) as the fit-agreement gate. The lmer formula of each tag is taken from
## the PROC MIXED statements in the program and listed in the fixture README.
## test_that("sas-part3-ngm-designs: broad-inference NGM vs SAS part 3", {
##   sas <- oracle_fixture("sas-part3-ngm-designs")
##   s <- sas$ngm[sas$ngm$ddfm == "SATTERTH", ]
##   for (tag in unique(s$fit)) {
##     st <- s[s$fit == tag, ]
##     fm <- sas_part3_fit(tag)                   # tag -> csv + lmer formula (fixture README)
##     r <- quiet_arm(new_group_mean_lmm(fm, newdata = sas_part3_newdata(tag, st),
##                    control = lmm_interval_control(ngm_method = "broad-inference")))
##     d <- attr(r, "diagnostics")
##     expect_true(all(abs(d$se / st$se - 1) < 1e-4), info = tag)
##     expect_true(all(abs(d$df - st$df) < 0.05), info = tag)
##   }
## })
