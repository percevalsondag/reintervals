# Test fixtures and the oracle registry

`oracles.csv` lists every **external oracle** the test suite checks the package
against: a number the package did not produce itself (a published table or
formula value, a SAS run, or independent code). Self-consistency, cross-method
and regression-pin tests are not oracles and are not listed (see plan section
T.1). `test-oracle-registry.R` keeps the registry and the tests in step.

## Columns

| Column | Content |
|---|---|
| `id` | Unique, lower case, `[a-z0-9-]` only. The same id must appear in `test_file`, in a `test_that()` description or a comment (`# oracle: <id>`), so `grep -rn <id> tests/` finds the check. |
| `verb_engine` | What is checked: `<interval>/<engine or arm>` (e.g. `ti/anova-mls`, `ngm/random intercept, ngm_method = broad-inference`). |
| `source_citation` | Full reference: authors, year, title, journal, volume(issue):pages, DOI, and the table / section / example. For a SAS run: software and version, run date, program and the source repository commit. |
| `source_kind` | `published_table` (numbers printed in a paper), `sas_run` (output of a SAS program run outside this repo), `published_formula` (a closed form stated in a paper, evaluated here), `independent_code` (code that shares nothing with the package: another R package, or a pinned external script). |
| `fixture_path` | Where the reference values (or the data they are computed on) live; see below. Empty while `pending`. |
| `test_file` | The test file under `tests/testthat/` that compares the package to the oracle (for a pending row: the file holding its commented-out template). |
| `status` | `collected` or `pending`. |
| `tolerance` | The tolerance the test uses, next to the printed precision it reflects. |
| `notes` | Anything a reviewer needs: data inputs, known deviations, corrections of misprints. |

### `fixture_path` forms

- `dir/file.ext`: a file under `tests/testthat/fixtures/`.
- `data:<name>`: a package dataset (`data/<name>.rda`, built by a script in
  `data-raw/`; the dataset's `@source` in `R/data.R` gives the reference).
- `inline:<file>`: the printed values are typed into the test code itself (in
  `<file>` under `tests/testthat/`, a test or a helper), with the reference in a
  comment next to them. Plan step T.2 will move these into fixture files; until
  then the registry records where they are.

## Statuses

- `collected`: the fixture exists and a passing test compares the package to it.
- `pending`: the oracle is identified but its values are not in the repo yet. Its
  `fixture_path` is empty and its `test_file` carries a commented-out template
  naming the id, what will be compared and at what tolerance. Pending oracles are
  registry status, never `skip()`ped tests.

A pending row becomes `collected` only when **both** its fixture file exists and
its test (the uncommented template) passes; change `status` and fill
`fixture_path` in the same commit. `test-oracle-registry.R` enforces the shape
(collected rows point at existing fixtures and test files containing the id;
pending rows have no fixture path; ids are unique; citations non-empty).

## Provenance and integrity

- Each fixture directory has a `README.md` with the source, how the files were
  produced (copied from which repository and commit, SAS version and run date, or
  generator), and the SHA-256 at the time of copying.
- **CSV convention.** A CSV copied from a source of record (the simstudy SAS
  fixtures) carries **no** in-file header: a header would change its bytes and
  break the hash match with the source; its provenance is the directory README.
  A CSV transcribed in this repo (a paper table) starts with `#` comment lines
  giving the reference, the table, who transcribed it and when, and who checked
  it; read it with `utils::read.csv(path, comment.char = "#")` (as
  `oracle_fixture()` does). Do not put `#` in data fields.
- A SAS results file (`*.results.txt`) starts with `#` header lines (program,
  software, run date) and is parsed by `parse_sas_results()` (helper-oracles.R).
  It holds only the `##` result lines, never the raw SAS log (host details).
- An `.rds` cannot carry a header: its README documents it (`slope_oracle.md`).
- `SHA256SUMS` lists the SHA-256 of every fixture data file;
  `test-fixture-sha256.R` checks it (R >= 4.5, `tools::sha256sum()`). Copies of
  external files must not be edited: re-copy them and update the hashes.

## Adding an oracle

1. Put the values in `fixtures/<id or source>/` with a README (or a `#` header
   for a transcribed CSV), and add its hash to `SHA256SUMS`.
2. Uncomment the template in `test_file` (or write the test), naming the id.
3. Set `status = collected` and `fixture_path` in `oracles.csv`.
4. Read it in the test with `oracle_fixture("<id>")`.
