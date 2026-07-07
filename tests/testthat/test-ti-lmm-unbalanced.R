# M2: unbalanced nested/crossed TI end-to-end through ti_lmm / vc_extract.
# Helpers build_44_counts/_data() are in helper-crosschecks.R.

fit_44 <- function() {
  suppressWarnings(suppressMessages(
    lme4::lmer(y ~ 1 + (1 | a) + (1 | b) + (1 | a:b), data = build_44_data(),
               control = lme4::lmerControl(check.conv.singular = "ignore"))
  ))
}

test_that("Section 4.4 END-TO-END: full pipeline yields kε=0.456 (plumbing, not just math)", {
  skip_if_not_installed("lme4")
  m <- fit_44()
  desc <- reintervals:::.design_of(m)
  expect_identical(desc$type, "crossed")
  expect_false(isTRUE(desc$balanced))
  # the k coefficients carried on the extracted re_components (model -> vc_extract
  # -> count extraction -> synthesis). kε is orientation-invariant.
  comp <- vc_extract(m, target = "observable")
  expect_false(is.null(comp$ems))
  expect_equal(round(unname(comp$ems$k[["Residual"]]), 3), 0.456)
  expect_equal(sort(round(unname(comp$ems$k), 3)), c(0.128, 0.151, 0.265, 0.456))
  # and ti_lmm returns a finite interval with the refusal note gone (the §4.4
  # fit is singular -- single-rep interaction cells -- so suppress that expected
  # warning; the EMS coefficients are count-based and unaffected)
  ti <- suppressWarnings(.ti_ems(m, level = 0.95, conf = 0.90))
  expect_null(attr(ti, "note"))
  expect_true(all(is.finite(c(ti$lower, ti$upper))))
})

test_that("cell-count extraction equals the known Section 4.4 table (catches transposition)", {
  skip_if_not_installed("lme4")
  m <- fit_44()
  desc <- reintervals:::.design_of(m)
  extracted <- reintervals:::.cell_counts(m, desc)
  # orient the known table to whichever factor lme4 labeled grp_A. A transposed
  # extraction would fail this even though Sigma(k)=1 would still pass.
  fl <- lme4::getME(m, "flist")
  expected <- if (nlevels(fl[[desc$grp_A]]) == 4) build_44_counts() else t(build_44_counts())
  expect_equal(extracted, expected, ignore_attr = TRUE)
  # marginals must match the true alpha/beta level counts (a second transposition trap)
  expect_setequal(rowSums(extracted), as.integer(table(fl[[desc$grp_A]])))
  expect_setequal(colSums(extracted), as.integer(table(fl[[desc$grp_B]])))
})

test_that("arbitrary unbalanced NESTED fit returns a finite TI (no paper oracle; refusal gone)", {
  skip_if_not_installed("lme4")
  # arbitrary unbalanced nested with UNEQUAL fine-levels-per-coarse (3, 2, 4) and
  # unequal reps -- the genuinely irregular case (not just unequal reps).
  dat <- data.frame(
    a = factor(c(rep("c1", 7), rep("c2", 5), rep("c3", 6))),
    b = factor(c("f1", "f1", "f1", "f2", "f2", "f3", "f3",
                 "f4", "f4", "f4", "f5", "f5",
                 "f6", "f6", "f7", "f8", "f8", "f9")),
    y = c(5.1, 4.8, 5.4, 6.2, 5.9, 4.4, 4.7,
          7.1, 6.8, 8.0, 5.5, 5.9,
          3.2, 3.5, 2.9, 3.8, 4.6, 6.1)
  )
  m <- suppressWarnings(suppressMessages(
    lme4::lmer(y ~ 1 + (1 | a / b), data = dat,
               control = lme4::lmerControl(check.conv.singular = "ignore"))
  ))
  desc <- reintervals:::.design_of(m)
  expect_identical(desc$type, "nested")
  expect_false(isTRUE(desc$balanced))

  ti <- .ti_ems(m, level = 0.95, conf = 0.90)
  expect_identical(attr(ti, "design"), "nested")
  expect_null(attr(ti, "note"))                          # NA-refusal removed
  expect_true(all(is.finite(c(ti$lower, ti$upper))))
  expect_lt(ti$lower, ti$estimate)
  expect_gt(ti$upper, ti$estimate)
})
