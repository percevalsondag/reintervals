# A4: the MLS upper bound on sum_j k_j E[MS_j] (.mls_upper, R/mls.R).
# Graybill-Wang when every k_j >= 0; Ting-Burdick-Graybill (1990) when a k_j < 0.

test_that("TBG bound reproduces a hand computation for one c2 < 0 case", {
  ## theta = c1 th1 + c2 th2 with c1 = 1, MS1 = 10, df1 = 5; c2 = -0.5, MS2 = 4,
  ## df2 = 20; confidence 0.95 (a = 0.05). Burdick & Graybill (1992) notation.
  ##   estimate  = 1 * 10 - 0.5 * 4                                   =  8
  ##   chi2_{.05}(5)  = 1.1454762261      H1  = 5 / 1.1454762261 - 1  =  3.3649967465
  ##   chi2_{.95}(20) = 31.4104328442     G2  = 1 - 20 / 31.4104328442 = 0.3632688827
  ##   F_{.95:5,20} = qf(.05, 5, 20)      = 0.2193881420
  ##   H12 = ((1 - F)^2 - H1^2 F^2 - G2^2) / F                        = -0.3081673970
  ##   VU  = (H1 * 1 * 10)^2 + (G2 * 0.5 * 4)^2 + H12 * 1 * 0.5 * 10 * 4
  ##       = 1132.3203 + 0.5278516 - 6.1633479                       = 1126.6848196073
  ##   U   = 8 + sqrt(1126.6848196073) = 8 + 33.5661260739          = 41.5661260739
  ## (identical to reintervals-simstudy pilot/methods.R mls_upper(1, 10, 5, -0.5, 4, 20, 0.95))
  u <- .mls_upper(ms = c(10, 4), k = c(1, -0.5), df = c(5, 20), conf = 0.95)
  expect_identical(u$bound, "ting-burdick-graybill")
  expect_equal(u$excess, 33.5661260739, tolerance = 1e-10)
  expect_equal(u$U, 41.5661260739, tolerance = 1e-10)
})

test_that("two or more negative coefficients are unsupported (no validated bound)", {
  expect_error(.mls_upper(ms = c(a = 10, b = 4, e = 2), k = c(a = 1.4, b = -0.2, e = -0.2),
                          df = c(a = 5, b = 20, e = 30), conf = 0.95),
               class = "reintervals_error_unsupported")
  expect_error(.mls_upper(ms = c(a = 10, b = 4, e = 2), k = c(a = 1.4, b = -0.2, e = -0.2),
                          df = c(a = 5, b = 20, e = 30), conf = 0.95),
               "more than one negative EMS coefficient is not implemented or validated")
  ## through ti_francq as well
  comp <- re_components(
    components = c(a = 1, b = 1, ab = 0.5, Residual = 1),
    dfs = c(pi = Inf), mean = 10, var_mean = 0.2,
    ems = list(ms = c(a = 9, b = 8, ab = 3, Residual = 1),
               k  = c(a = 0.9, b = 0.5, ab = -0.2, Residual = -0.2),
               df = c(a = 3, b = 4, ab = 12, Residual = 40)),
    design = "crossed")
  expect_error(ti_francq(comp), class = "reintervals_error_unsupported")
})

test_that("Graybill-Wang branch is the 0.2.0 arithmetic bit for bit", {
  ms <- c(run = 0.001253 + 3 * 0.000681, Residual = 0.001253)     # Francq 5.2
  k  <- c(run = 1 / 3, Residual = 2 / 3); df <- c(run = 5, Residual = 12)
  u  <- .mls_upper(ms, k, df, conf = 0.90)
  a  <- 1 - 0.90                       # as in ti_francq (0.2.0): not exactly 0.10
  h  <- df / stats::qchisq(a, df = df) - 1
  expect_identical(u$bound, "graybill-wang")
  expect_identical(u$excess, sqrt(sum(h^2 * k^2 * ms^2)))
})

test_that("round-off negative k (exact zero) does not switch to TBG", {
  ## the 2 x 2 crossed table below has exact k_AB = 0; the synthesis returns
  ## -3.6e-17. A bare `k < 0` test would switch the bound on round-off.
  syn <- ems_synthesis(matrix(c(3, 4, 1, 1), 2, 2), "crossed")
  expect_lt(syn$k[["AB"]], 0)
  expect_gt(syn$k[["AB"]], -1e-15)
  u <- .mls_upper(ms = c(5, 4, 3, 1), k = syn$k, df = syn$df, conf = 0.95)
  expect_identical(u$bound, "graybill-wang")
})

test_that("a negative k_j makes ti_francq use TBG, report it, and stay finite", {
  ## No unbalanced design reachable through the synthesis gave a negative k_j in
  ## the searches run for A4 (625 2x2 crossed tables with n_ij in 1..5; 3000+
  ## random crossed tables; 3000 random nested designs: smallest k was
  ## -2.8e-16 round-off). So the negative coefficient is constructed at the
  ## container level: an EMS combination with k_AB < 0 and sum(k) = 1.
  comp <- re_components(
    components = c(a = 1, b = 1, ab = 0.5, Residual = 1),
    dfs = c(pi = Inf), mean = 10, var_mean = 0.2,
    ems = list(ms = c(a = 9, b = 8, ab = 3, Residual = 1),
               k  = c(a = 0.6, b = 0.5, ab = -0.3, Residual = 0.2),
               df = c(a = 3, b = 4, ab = 12, Residual = 40)),
    design = "crossed")
  for (s in c("two", "lower", "upper")) {
    ti <- ti_francq(comp, content = 0.95, conf = 0.95, sides = s)
    expect_identical(attr(ti, "mls_bound"), "ting-burdick-graybill", info = s)
    expect_true(is.finite(if (s == "upper") ti$upper else ti$lower), info = s)
  }
  gw <- comp; gw$ems$k[] <- c(0.3, 0.3, 0.2, 0.2)
  expect_identical(attr(ti_francq(gw), "mls_bound"), "graybill-wang")
})

test_that("mls_bound reaches the diagnostics of every MLS tolerance interval", {
  skip_if_not_installed("lme4")
  fe <- ti_lmm(fit_52(), method = "ems-mls")                          # Francq EMS
  expect_identical(attr(fe, "diagnostics")$mls_bound, "graybill-wang")
  fh <- ti_lmm(lme4::lmer(conc ~ 1 + (1 | run), data = hk_assay))     # HK anova-mls
  expect_identical(attr(fh, "diagnostics")$mls_bound, "graybill-wang")
  ## the random-slope TI is GPQ (no MLS bound): mls_bound is NA there
  fs <- suppressWarnings(lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
                                    data = lme4::sleepstudy))
  fr <- ti_lmm(fs, newdata = data.frame(Days = 5),
               control = lmm_interval_control(seed = 1, M = 2000))  # not an oracle
  expect_null(attr(fr, "diagnostics")$mls_bound)
})

test_that("crossed NGM df: an undefined Satterthwaite term falls back to containment, never NaN", {
  skip_if_not_installed("lme4")
  set.seed(1); A <- 3; B <- 4; n <- 2
  g <- expand.grid(r = 1:n, a = factor(1:A), b = factor(1:B)); g$ab <- interaction(g$a, g$b)
  g$y <- 10 + rnorm(A, 0, 2)[g$a] + rnorm(B, 0, 2)[g$b] + rnorm(A * B, 0, 1)[g$ab] + rnorm(nrow(g))
  fm <- suppressMessages(lme4::lmer(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), data = g))
  ok <- new_group_mean_lmm(fm)
  expect_identical(attr(ok, "diagnostics")$df_rule, "satterthwaite-vf")
  ## force k_A EMS_A + k_B EMS_B = 0 -> pA = 0/0 (the undefined Satterthwaite term)
  real <- .design_components_fit
  local_mocked_bindings(.design_components_fit = function(...) {
    dc <- real(...); dc$ems$k[1:2] <- 0; dc
  })
  r <- new_group_mean_lmm(fm)
  expect_identical(attr(r, "diagnostics")$df_rule, "containment")
  expect_equal(attr(r, "diagnostics")$df, .ddf_containment(fm))      # min(A, B) - 1 = 2
  expect_true(all(is.finite(c(r$lower, r$upper))))
})
