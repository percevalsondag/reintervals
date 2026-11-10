# Oracle registry and fixture readers (see fixtures/README.md).

oracle_registry <- function() {
  utils::read.csv(test_path("fixtures", "oracles.csv"), stringsAsFactors = FALSE,
                  na.strings = character(), colClasses = "character")
}

## The parsed fixture of a registry id:
##   *.results.txt -> parse_sas_results();  *.csv -> read.csv(comment.char = "#");
##   *.rds -> readRDS();  data:<name> -> the package dataset.
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
  if (grepl("\\.results\\.txt$", fp)) return(parse_sas_results(path))
  switch(tolower(tools::file_ext(fp)),
    csv = utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE),
    rds = readRDS(path),
    stop("no reader for fixture '", fp, "'"))
}

## Parse the `##` result lines of a SAS NGM oracle results file into one data
## frame per line kind. Line layout (fields separated by runs of blanks):
##   ##CONV <fit> <ddfm> Status=<n> Reason=<text>
##   ##COVP <fit> <ddfm> <covparm> [<subject>] Estimate=<x>
##   ##FIXD <fit> <ddfm> <effect> Est=<x> SE=<x> DF=<x>
##   ##NGM / ##MEAN <fit> <ddfm> <var>= <x> Pred=<x> SE=<x> DF=<x> Lower=<x> Upper=<x>
## A SAS missing value (".") is NA. For ##NGM / ##MEAN, `x_var` is the evaluation
## variable ("one" for an intercept-only fit, else e.g. "month"), `x` its value,
## and `t0` = x, or NA when x_var == "one". ##BEGIN / ##END / ##SKIP lines and the
## `#` header are ignored.
parse_sas_results <- function(path) {
  res <- readLines(path, warn = FALSE)
  lines_of <- function(kind) grep(paste0("^##", kind, " "), res, value = TRUE)
  num <- function(key, x) {
    out <- rep(NA_real_, length(x))
    hit <- grepl(paste0("(^| )", key, "= *[-0-9.]"), x)
    v <- sub(paste0(".*(^| )", key, "= *([-0-9.eE]+|\\.).*"), "\\2", x[hit])
    out[hit] <- ifelse(v == ".", NA_real_, suppressWarnings(as.numeric(v)))
    out
  }
  head3 <- function(x) {
    f <- strsplit(trimws(x), " +")
    data.frame(fit = vapply(f, `[`, "", 2L), ddfm = vapply(f, `[`, "", 3L),
               stringsAsFactors = FALSE)
  }
  rest <- function(x) sub("^##[A-Z]+ +[^ ]+ +[^ ]+ +", "", x)

  pred <- function(kind) {
    m <- lines_of(kind)
    r <- rest(m)
    x_var <- sub("^([A-Za-z_][A-Za-z0-9_]*)=.*", "\\1", r)
    x <- as.numeric(sub("^[A-Za-z_][A-Za-z0-9_]*= *([-0-9.eE]+).*", "\\1", r))
    cbind(head3(m), x_var = x_var, x = x, t0 = ifelse(x_var == "one", NA_real_, x),
          pred = num("Pred", m), se = num("SE", m), df = num("DF", m),
          lower = num("Lower", m), upper = num("Upper", m),
          stringsAsFactors = FALSE)
  }

  covp <- lines_of("COVP")
  cp <- strsplit(trimws(sub(" *Estimate=.*", "", rest(covp))), " +")
  fixd <- lines_of("FIXD")
  conv <- lines_of("CONV")
  list(
    ngm  = pred("NGM"),
    mean = pred("MEAN"),
    covp = cbind(head3(covp),
                 covparm = vapply(cp, `[`, "", 1L),
                 subject = vapply(cp, function(z) if (length(z) > 1L && z[2] != ".") z[2] else NA_character_, ""),
                 estimate = num("Estimate", covp), stringsAsFactors = FALSE),
    fixd = cbind(head3(fixd),
                 effect = sub(" .*", "", rest(fixd)),
                 estimate = num("Est", fixd), se = num("SE", fixd), df = num("DF", fixd),
                 stringsAsFactors = FALSE),
    conv = cbind(head3(conv),
                 status = as.integer(num("Status", conv)),
                 reason = sub(".*Reason=", "", conv), stringsAsFactors = FALSE)
  )
}
