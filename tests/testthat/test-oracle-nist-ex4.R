# oracle: nist-recipe-ex4 (fixtures/oracles.csv)
# oracle: hoffman2010-s52-nist-ex4 (fixtures/oracles.csv)
#
# NIST RECIPE Example 4 (Vangel; RECIPE 1.0, fixtures/nist_ex4/ex4.dat and
# ex4.out): 72 strength values on 8 batches at two fixed temperatures (75 F and
# -67 F; unbalanced: batch 6 is at 75 F only, batches 7 and 8 at -67 F only),
# lower one-sided tolerance limits (P = 0.90, confidence 0.95) at 7 points. The
# same data are reanalysed in Hoffman D (2010), Technometrics 52(3):303-312,
# doi:10.1198/TECH.2010.09171, section 5.2.
#
# Identification (read from the package source; see the registry notes):
#   - RECIPE's "Regression" column is the OLS line of strength on temperature
#     (checked below to 1e-5). Its limit uses that OLS center.
#   - Hoffman's "NEW" limit is his eq. 1 on a WLS line through the cell
#     (batch x temperature) means.
#   - ti_lmm(method = "anova-mls"), the default route here (the design is not
#     clean balanced, .m1_not_clean_balanced), is the Montes et al. (2019)
#     unweighted lot-mean construction (.montes_core, R/ti-montes.R): center
#     ybar_star + beta_w (t0 - tbar_star), beta_w the pooled within-lot slope,
#     with Hoffman's one-sided H limit on that center.
#   - ti_lmm(method = "ems-mls") (.ti_ems, R/ti-lmm.R; ti_francq) centers on the
#     REML GLS line (.fixed) with the package's extension of Hoffman's limit.
# Neither package center is the OLS line or a WLS cell-mean line, so no package
# limit is the RECIPE or the Hoffman NEW estimator: nothing is asserted equal.
# The comparison table is recorded, and only a sanity screen is asserted:
# lower < center, finite, monotone in t0 in the direction of the fitted slope,
# and within 10% of the RECIPE limit. The 10% band is a gross-error screen, NOT
# an oracle tolerance.

nist_t0 <- c(-67, -50, -25, 0, 25, 50, 75)

## RECIPE 1.0 output (ex4.out), P = 0.90, confidence 0.95
nist_out <- function() {
  l <- readLines(test_path("fixtures", "nist_ex4", "ex4.out"), warn = FALSE)
  l <- trimws(l)
  rows <- l[grepl("^0\\.90\\s+0\\.95\\s", l)]
  x <- utils::read.table(text = rows, col.names = c("prob", "conf", "regression", "limit"))
  cbind(t0 = nist_t0, x)
}

## RECIPE input (ex4.dat): control line, 2 level rows (1 temp), 72 data rows
## (level batch strength), 7 point rows (1 t0)
nist_dat <- function() {
  l <- trimws(readLines(test_path("fixtures", "nist_ex4", "ex4.dat"), warn = FALSE))
  l <- l[nzchar(l) & !startsWith(l, "#")]
  ctl <- strsplit(l[1], "\\s+")[[1]]
  n <- as.integer(ctl[1]); nlvl <- as.integer(ctl[2]); npts <- as.integer(ctl[5])
  lev <- utils::read.table(text = l[2:(1 + nlvl)])[[2]]
  d <- utils::read.table(text = l[(2 + nlvl):(1 + nlvl + n)],
                         col.names = c("level", "batch", "strength"))
  pts <- utils::read.table(text = l[(2 + nlvl + n):(1 + nlvl + n + npts)])[[2]]
  d$temp <- lev[d$level]
  d$batch <- factor(d$batch)
  list(data = d, levels = lev, t0 = pts, n = n, prob = as.numeric(sub("d0$", "", ctl[6])),
       conf = as.numeric(sub("d0$", "", ctl[7])))
}

## Hoffman (2010) section 5.2, lower limits (beta 0.90, gamma 0.95) at nist_t0
hoffman_new  <- c(284.908, 284.437, 283.315, 281.493, 278.906, 275.713, 272.132)  # analytical "NEW"
hoffman_boot <- c(287.853, 287.381, 286.260, 284.437, 281.851, 278.658, 275.076)  # bootstrap-adjusted

nist_fit <- function(d) suppressWarnings(suppressMessages(
  lme4::lmer(strength ~ temp + (1 | batch), data = d, REML = TRUE)))

test_that("nist-recipe-ex4: the parsed input reproduces RECIPE's regression column (OLS, 1e-5)", {
  x <- nist_dat(); o <- nist_out()
  expect_identical(x$n, 72L); expect_identical(nrow(x$data), 72L)
  expect_identical(nlevels(x$data$batch), 8L)
  expect_equal(x$levels, c(75, -67))
  expect_equal(x$t0, nist_t0)
  expect_equal(c(x$prob, x$conf), c(0.90, 0.95))
  expect_identical(nrow(o), 7L)
  f <- stats::lm(strength ~ temp, data = x$data)
  expect_equal(unname(stats::predict(f, data.frame(temp = nist_t0))), o$regression,
               tolerance = 1e-5 / 330)                         # 1e-5 absolute at ~330
  expect_true(all(abs(stats::predict(f, data.frame(temp = nist_t0)) - o$regression) < 1e-5))
  expect_true(all(abs(unname(stats::coef(f)) - c(322.098785, -0.0811720)) < 1e-5))
})

test_that("nist-recipe-ex4 / hoffman2010-s52-nist-ex4: package one-sided TI vs RECIPE and Hoffman (comparison; sanity screen only)", {
  skip_if_not_installed("lme4")
  x <- nist_dat(); o <- nist_out(); fm <- nist_fit(x$data)
  expect_true(.m1_not_clean_balanced(fm))                 # -> auto routes to anova-mls
  nd <- data.frame(temp = nist_t0)
  ti <- function(m) suppressWarnings(ti_lmm(fm, nd, P = 0.90, level = 0.95, sides = "lower",
                                           method = m))
  res <- list(auto = ti("auto"), `anova-mls` = ti("anova-mls"), `ems-mls` = ti("ems-mls"))
  ## method = "gpq" (step 13b) is compared with Liao, Lin & Iyer (2005) Example 3
  ## in test-oracle-lli2005.R
  expect_identical(res$auto$method, rep("anova-mls", 7))
  expect_identical(res$`ems-mls`$method, rep("ems-mls", 7))
  expect_identical(res$auto$lower, res$`anova-mls`$lower)

  ## the identification, checked: anova-mls centers on the Montes lot-mean line,
  ## ems-mls on the REML GLS line; neither is the OLS (RECIPE) line
  d <- x$data
  lm_i <- tapply(d$strength, d$batch, mean); tm_i <- tapply(d$temp, d$batch, mean)
  tc <- d$temp - tm_i[d$batch]; yc <- d$strength - lm_i[d$batch]
  bw <- sum(tc * yc) / sum(tc^2)
  expect_equal(res$`anova-mls`$estimate, mean(lm_i) + bw * (nist_t0 - mean(tm_i)), tolerance = 1e-10)
  fe <- lme4::fixef(fm)
  expect_equal(res$`ems-mls`$estimate, unname(fe[1] + fe[2] * nist_t0), tolerance = 1e-10)
  for (m in names(res)) expect_gt(max(abs(res[[m]]$estimate - o$regression)), 0.1)

  slope_sign <- sign(fe[["temp"]])
  tab <- NULL
  for (m in c("anova-mls", "ems-mls")) {
    r <- res[[m]]
    ## sanity screen (not an oracle)
    expect_true(all(is.finite(c(r$estimate, r$lower))), info = m)
    expect_true(all(r$lower < r$estimate), info = m)
    expect_true(all(sign(diff(r$lower)) == slope_sign), info = m)
    expect_true(all(abs(r$lower - o$limit) / o$limit < 0.10), info = m)
    tab <- rbind(tab, data.frame(method = m, t0 = nist_t0,
      center = r$estimate, lower = r$lower, recipe_reg = o$regression,
      recipe_lim = o$limit, hoffman_new = hoffman_new, hoffman_boot = hoffman_boot,
      d_center_recipe = r$estimate - o$regression, d_recipe = r$lower - o$limit,
      d_hoffman_new = r$lower - hoffman_new, d_hoffman_boot = r$lower - hoffman_boot))
  }
  cat("\nnist-recipe-ex4 / hoffman2010-s52-nist-ex4 comparison (lower, P 0.90, conf 0.95):\n")
  print(format(tab, digits = 7), row.names = FALSE)
})
