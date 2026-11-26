# lmm_interval_control(): every choice argument (validated with match.arg)
# signals an error of class reintervals_error_input for an invalid value, with
# a message naming the argument and listing its allowed values. Valid values,
# partial matching and the defaults are unchanged.

control_choices <- list(
  ci_df            = c("satterthwaite", "kenward-roger", "containment", "capped-kr"),
  slope_ci_df      = c("hessian-satterthwaite", "kenward-roger", "k-minus-1"),
  pi_df            = c("francq-satterthwaite", "containment", "capped", "ratio-guard"),
  ngm_method       = c("package", "broad-inference", "hts", "classical"),
  slope_ngm_method = c("reml", "broad-inference", "gpq", "fv-anova"),
  slope_pi_method  = c("reml", "fv-anova", "hybrid"),
  gpq_two_sided    = c("approx", "content")
)

test_that("the choice arguments of lmm_interval_control() are the ones listed here", {
  f <- formals(lmm_interval_control)
  is_choice <- vapply(f, function(x) is.call(x) && identical(x[[1]], as.name("c")), logical(1))
  expect_setequal(names(f)[is_choice], names(control_choices))
  for (nm in names(control_choices))
    expect_identical(eval(f[[nm]]), control_choices[[nm]], info = nm)
})

for (nm in names(control_choices)) {
  local({
    arg <- nm; ch <- control_choices[[nm]]
    test_that(sprintf("lmm_interval_control(%s = <invalid>) is a reintervals_error_input naming it", arg), {
      allowed <- paste0("\"", ch, "\"", collapse = ", ")
      for (bad in list("no-such-value", c(ch[1], "x"), 1, NA_character_)) {
        e <- tryCatch(do.call(lmm_interval_control, stats::setNames(list(bad), arg)),
                      error = identity)
        expect_s3_class(e, "reintervals_error_input")
        expect_s3_class(e, "reintervals_error")
        expect_identical(conditionMessage(e), sprintf("`%s` must be one of %s.", arg, allowed))
        expect_identical(e$argument, arg)
        expect_identical(e$choices, ch)
      }
      ## valid values, partial matching and the default are unchanged
      for (v in ch) expect_identical(do.call(lmm_interval_control,
                                             stats::setNames(list(v), arg))[[arg]], v)
      expect_identical(do.call(lmm_interval_control,
                               stats::setNames(list(substr(ch[2], 1, 5)), arg))[[arg]], ch[2])
      expect_identical(lmm_interval_control()[[arg]], ch[1])
    })
  })
}
