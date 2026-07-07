#' Interval result object
#'
#' The shallow S3 value returned by the EMS method-pure constructors. It records
#' the point estimate, the two-sided bounds, and the metadata identifying which
#' interval was computed. It is internal plumbing: the user-facing verbs return
#' the unified `lmm_interval` data frame, not this object.
#'
#' Provides `print` and `as.data.frame` methods.
#'
#' @param estimate Numeric scalar point estimate (`lb-hat`).
#' @param lower,upper Numeric scalar interval bounds.
#' @param type Character: `"CI"`, `"PI"`, or `"TI"`.
#' @param conf Confidence level (CI: `1 - alpha`; TI: `gamma`). `NA` if not used.
#' @param level Prediction level (PI: `1 - psi`). `NA` if not used.
#' @param content Content fraction (TI: `beta`). `NA` if not used.
#' @param df Degrees of freedom used for the bound (`NA` for the normal-based TI).
#' @param sides Character: `"two"` (two-sided intervals only).
#' @param target Character target distribution (`"observable"` / `"true_value"`),
#'   or `NA` for the CI (which concerns the mean, not a future value).
#' @param design Character naming the random structure, or `""`.
#' @param call The originating call, or `NULL`.
#'
#' @return An object of class `re_interval`.
#' @aliases re_interval
#' @keywords internal
new_re_interval <- function(estimate, lower, upper, type,
                            conf = NA_real_, level = NA_real_,
                            content = NA_real_, df = NA_real_,
                            sides = "two", target = NA_character_,
                            design = "", call = NULL) {
  structure(
    list(estimate = estimate, lower = lower, upper = upper, type = type,
         content = content, conf = conf, level = level, df = df,
         sides = sides, method = "francq", target = target,
         design = design, call = call),
    class = "re_interval"
  )
}

#' @export
print.re_interval <- function(x, digits = 4, ...) {
  lvl <- if (!is.na(x$conf)) x$conf else x$level
  hdr <- paste0(x$type, if (!is.na(x$content)) paste0(" (", x$content, " content)"),
                if (!is.na(lvl)) paste0(" at ", lvl, " level"))
  cat("<re_interval>", hdr, "\n")
  cat("  estimate:", format(x$estimate, digits = digits), "\n")
  cat("  ", x$sides, "-sided: [", format(x$lower, digits = digits), ", ",
      format(x$upper, digits = digits), "]\n", sep = "")
  if (!is.na(x$df)) cat("  df      :", format(x$df, digits = digits), "\n")
  invisible(x)
}

#' @export
as.data.frame.re_interval <- function(x, ...) {
  data.frame(
    type = x$type, estimate = x$estimate, lower = x$lower, upper = x$upper,
    df = x$df, content = x$content, conf = x$conf, level = x$level,
    sides = x$sides, target = x$target, design = x$design,
    stringsAsFactors = FALSE, row.names = NULL
  )
}
