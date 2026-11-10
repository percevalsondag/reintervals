# Fixture integrity: every fixture data file matches the SHA-256 recorded in
# fixtures/SHA256SUMS (copies of external files must not be edited; see
# fixtures/README.md). tools::sha256sum() exists from R 4.5.0.

test_that("fixture files match fixtures/SHA256SUMS", {
  skip_if_not(exists("sha256sum", envir = asNamespace("tools")), "tools::sha256sum() needs R >= 4.5")
  sums <- readLines(test_path("fixtures", "SHA256SUMS"))
  sums <- sums[nzchar(sums) & !startsWith(sums, "#")]
  hash <- sub(" .*", "", sums)
  path <- sub("^[0-9a-f]{64} [ *]", "", sums)
  expect_gt(length(path), 0)
  got <- unname(tools::sha256sum(test_path("fixtures", path)))
  expect_identical(got, hash, info = paste(path[got != hash], collapse = ", "))
})

test_that("every fixture data file is listed in fixtures/SHA256SUMS", {
  sums <- readLines(test_path("fixtures", "SHA256SUMS"))
  listed <- sub("^[0-9a-f]{64} [ *]", "", sums[nzchar(sums) & !startsWith(sums, "#")])
  all_files <- list.files(test_path("fixtures"), recursive = TRUE)
  data_files <- all_files[!grepl("\\.md$", all_files) & !all_files %in% c("oracles.csv", "SHA256SUMS")]
  expect_setequal(listed, data_files)
})
