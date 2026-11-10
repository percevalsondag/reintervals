# The oracle registry (fixtures/oracles.csv; columns and statuses in
# fixtures/README.md) stays in step with the fixtures and the tests. Pending
# oracles are registry status, never skipped tests: this file has no skip().

reg <- oracle_registry()

## id appears in a test_that() description or in a comment of the test file
id_traced <- function(id, test_file) {
  lines <- readLines(test_path(test_file), warn = FALSE)
  pat <- paste0("(^|[^a-z0-9-])", id, "($|[^a-z0-9-])")
  hit <- grepl(pat, lines)
  any(hit & (grepl("^\\s*#", lines) | grepl("test_that\\(", lines)))
}

test_that("the registry has the documented columns, statuses and source kinds", {
  expect_identical(names(reg), c("id", "verb_engine", "source_citation", "source_kind",
                                 "fixture_path", "test_file", "status", "tolerance", "notes"))
  expect_gt(nrow(reg), 0)
  expect_true(all(reg$status %in% c("collected", "pending")))
  expect_true(all(reg$source_kind %in%
    c("published_table", "sas_run", "published_formula", "independent_code")))
})

test_that("ids are unique and use [a-z0-9-] only", {
  expect_false(anyDuplicated(reg$id) > 0, info = paste(reg$id[duplicated(reg$id)], collapse = ", "))
  expect_true(all(grepl("^[a-z0-9][a-z0-9-]*$", reg$id)),
              info = paste(reg$id[!grepl("^[a-z0-9][a-z0-9-]*$", reg$id)], collapse = ", "))
})

test_that("source_citation, verb_engine, tolerance and test_file are non-empty everywhere", {
  for (col in c("source_citation", "verb_engine", "tolerance", "test_file"))
    expect_true(all(nzchar(trimws(reg[[col]]))), info = col)
})

test_that("every collected oracle has an existing fixture and test file that names its id", {
  col <- reg[reg$status == "collected", ]
  expect_gt(nrow(col), 0)
  for (i in seq_len(nrow(col))) {
    id <- col$id[i]; fp <- col$fixture_path[i]; tf <- col$test_file[i]
    expect_true(nzchar(fp), info = id)
    if (startsWith(fp, "data:")) {
      nm <- sub("^data:", "", fp)
      expect_true(nm %in% utils::data(package = "reintervals")$results[, "Item"], info = id)
    } else if (startsWith(fp, "inline:")) {
      expect_true(file.exists(test_path(sub("^inline:", "", fp))), info = id)
    } else {
      expect_true(file.exists(test_path("fixtures", fp)), info = paste(id, fp))
    }
    expect_true(grepl("^test-.*\\.R$", tf) && file.exists(test_path(tf)), info = paste(id, tf))
    expect_true(id_traced(id, tf), info = paste(id, "not named in", tf))
  }
})

test_that("every pending oracle has an empty fixture_path and a template naming its id", {
  pen <- reg[reg$status == "pending", ]
  expect_true(all(pen$fixture_path == ""), info = paste(pen$id[pen$fixture_path != ""], collapse = ", "))
  for (i in seq_len(nrow(pen))) {
    expect_true(file.exists(test_path(pen$test_file[i])), info = pen$id[i])
    expect_true(id_traced(pen$id[i], pen$test_file[i]), info = pen$id[i])
  }
})

test_that("oracle_fixture() reads a collected fixture and refuses pending and inline rows", {
  sas <- oracle_fixture("sas-part1-ci-satterth")
  expect_named(sas, c("ngm", "mean", "covp", "fixd", "conv"))
  expect_identical(nrow(sas$mean), 56L)
  expect_s3_class(oracle_fixture("montes2019-appa-ti"), "data.frame")
  expect_error(oracle_fixture("chow-shao1991-t23"), "pending")
  expect_error(oracle_fixture("francq2019-s52-ci"), "inline")
  expect_error(oracle_fixture("no-such-oracle"), "no oracle")
})
