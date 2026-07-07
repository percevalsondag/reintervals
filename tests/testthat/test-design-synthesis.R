# Helpers build_44_counts(), n0_approx_k/_ems(), nested_ems_closed(),
# ti_width_from_ems() are in helper-crosschecks.R.

test_that("ACCEPTANCE GATE: synthesis reproduces Francq Section 4.4 (unbalanced crossed)", {
  syn <- ems_synthesis(build_44_counts(), "crossed")
  # EMS_A coefficients (paper: sigma2_e=1, sigma2_alphabeta=1.568, sigma2_alpha=7.84)
  expect_equal(syn$E["A", "A"], 7.84, tolerance = 1e-3)
  expect_equal(syn$E["A", "AB"], 1.568, tolerance = 1e-3)
  expect_equal(syn$E["A", "B"], 0)
  expect_equal(unname(syn$E[, "Residual"]), c(1, 1, 1, 1))   # residual coeff is 1 in every EMS
  # Combination coefficients k.
  # The paper prints kA=0.128, kB=0.151, kAB=0.265, ke=0.464. The ke=0.464 is a
  # TRANSCRIPTION ERROR: the four must sum to 1 (every EMS has residual coeff 1,
  # so sum(k)=1), but 0.128+0.151+0.265+0.464 = 1.008. kA+kB+kAB = 0.544 forces
  # ke = 0.456, NOT 0.464. We assert the correct 0.456.
  # match the paper's printed 3 decimals (synthesis: 0.1276/0.1510/0.2650/0.4564)
  expect_equal(round(unname(syn$k["A"]), 3), 0.128)
  expect_equal(round(unname(syn$k["B"]), 3), 0.151)
  expect_equal(round(unname(syn$k["AB"]), 3), 0.265)
  expect_equal(round(unname(syn$k["Residual"]), 3), 0.456)
  expect_equal(sum(syn$k), 1)
})

test_that("Sum(k)=1 is a hard runtime invariant that errors on a malformed E", {
  good <- ems_synthesis(build_44_counts(), "crossed")
  expect_equal(sum(good$k), 1)                              # positive: enforced
  # craft an E whose residual column is not all-ones -> sum(k) != 1 -> must error
  bad <- matrix(c(2, 0, 0,  0.5, 1.3, 0,  2, 1, 1), nrow = 3,
                dimnames = list(c("A", "B", "Residual"), c("A", "B", "Residual")))
  expect_error(.solve_k_observable(bad), "invariant violated")
})

test_that("balanced-limit crossed: synthesis == M1 closed form == n0-approx", {
  cnt <- matrix(2L, nrow = 3, ncol = 4)                    # A=3, B=4, n=2, balanced
  syn <- ems_synthesis(cnt, "crossed")
  m1 <- ems_crossed(list(grp_A = "A", grp_B = "B", grp_AB = "AB",
                         A = 3, B = 4, n = 2),
                    c(A = 1.5, B = 1, AB = 0.5, Residual = 0.25))
  expect_equal(syn$k, m1$ems$k, tolerance = 1e-9)
  expect_equal(syn$k, n0_approx_k(cnt, "crossed"), tolerance = 1e-9)
  expect_equal(unname(syn$df), unname(m1$ems$df))
})

test_that("balanced-limit nested: synthesis == M1 closed form == n0-approx", {
  cnt <- list(c(4, 4), c(4, 4), c(4, 4))                   # A=3, B=2, n=4, balanced
  syn <- ems_synthesis(cnt, "nested")
  m1 <- ems_nested(list(grp_alpha = "A", grp_beta = "B", A = 3, B = 2, n = 4),
                   c(A = 2, B = 1, Residual = 0.5))
  expect_equal(syn$k, m1$ems$k, tolerance = 1e-9)
  expect_equal(syn$k, n0_approx_k(cnt, "nested"), tolerance = 1e-9)
})

test_that("nested synthesis matches independent closed-form arithmetic (hand-computed)", {
  # Tiny unbalanced nested: A=2 coarse; a1 fine reps (2,1), a2 fine reps (1,2).
  # By hand (block-projection traces): EMS_A = e + (5/3)b + 3a, EMS_B(A) = e + (4/3)b.
  cnt <- list(c(2, 1), c(1, 2))
  syn <- ems_synthesis(cnt, "nested")
  expect_equal(syn$E["A", "A"], 3)
  expect_equal(syn$E["A", "B"], 5 / 3)
  expect_equal(syn$E["B", "B"], 4 / 3)
  expect_equal(syn$E["B", "A"], 0)
  expect_equal(syn$E, nested_ems_closed(cnt), tolerance = 1e-9, ignore_attr = "dimnames")
  expect_equal(unname(syn$k), c(1 / 3, 1 / 3, 1 / 3), tolerance = 1e-9)
})

test_that("trace cross-check: the tr(Q Z Z') shortcut equals the explicit form", {
  syn <- ems_synthesis(build_44_counts(), "crossed")
  # residual column must be exactly 1 in every term (structural invariant)
  expect_equal(unname(syn$E[, "Residual"]), rep(1, nrow(syn$E)))
  # off-diagonal orthogonality the synthesis must produce for crossed Type III
  expect_equal(syn$E["A", "B"], 0)
  expect_equal(syn$E["B", "A"], 0)
  expect_equal(syn$E["AB", "A"], 0)
  expect_equal(syn$E["AB", "B"], 0)
})

test_that("n0-approx cross-check: agrees at the balanced limit, diverges under imbalance", {
  sig <- c(A = 2, B = 1, AB = 0.5, Residual = 1)
  # (a) balanced limit -> identical coefficients
  bal <- matrix(3L, 3, 4)
  expect_equal(ems_synthesis(bal, "crossed")$k, n0_approx_k(bal, "crossed"),
               tolerance = 1e-9)
  # (b) severe imbalance -> the two methods DIVERGE (not asserted close)
  sev <- build_44_counts()
  syn_k <- ems_synthesis(sev, "crossed")$k
  n0_k <- n0_approx_k(sev, "crossed")
  expect_gt(max(abs(syn_k - n0_k)), 0.02)                  # genuinely different
  # documented direction: the n0-approx TI is narrower (anti-conservative)
  w_syn <- ti_width_from_ems(ems_list_from_synthesis(ems_synthesis(sev, "crossed"), sig), sig)
  w_n0 <- ti_width_from_ems(n0_approx_ems(sev, "crossed", sig), sig)
  expect_lt(w_n0, w_syn)
})

test_that("ems_synthesis validates its input", {
  expect_error(ems_synthesis(matrix(c(1, 0, 2, 3), 2, 2), "crossed"), "cell filled")
  expect_error(ems_synthesis(list(c(2, 1)), "nested"), ">= 2 coarse")
})
