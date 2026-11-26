# Oracle registry and fixture readers (see fixtures/README.md).

oracle_registry <- function() {
  utils::read.csv(test_path("fixtures", "oracles.csv"), stringsAsFactors = FALSE,
                  na.strings = character(), colClasses = "character")
}

## The parsed fixture of a registry id:
##   *.csv -> read.csv(comment.char = "#");  *.rds -> readRDS();
##   data:<name> -> the package dataset.
## inline: and pending rows have no fixture to read and are an error.
oracle_fixture <- function(id) {
  reg <- oracle_registry()
  row <- reg[reg$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("no oracle with id '", id, "' in fixtures/oracles.csv")
  fp <- row$fixture_path
  if (row$status != "collected" || !nzchar(fp))
    stop("oracle '", id, "' is ", row$status, " and has no fixture")
  if (startsWith(fp, "inline:"))
    stop("oracle '", id, "' is inline in ", sub("^inline:", "", fp), "; there is no fixture file")
  if (startsWith(fp, "data:")) {
    e <- new.env()
    utils::data(list = sub("^data:", "", fp), package = "reintervals", envir = e)
    return(get(sub("^data:", "", fp), envir = e))
  }
  path <- test_path("fixtures", fp)
  switch(tolower(tools::file_ext(fp)),
    csv = utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE),
    rds = readRDS(path),
    stop("no reader for fixture '", fp, "'"))
}
