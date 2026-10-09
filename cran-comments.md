## Submission type

New submission (first release of `reintervals`).

## Test environments

- Local: macOS (aarch64-apple-darwin20), R 4.5.3.
- GitHub Actions: ubuntu-latest (R devel, R release, R oldrel-1),
  macOS-latest (R release), windows-latest (R release) --- dependencies installed
  from current CRAN.

## R CMD check results

`R CMD check --as-cran`: 0 errors | 0 warnings | 1 note.

The single NOTE is the expected one for a first submission:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Perceval Sondag <reintervals@gmail.com>'

New submission
```

## Reverse dependencies

None --- this is a new submission.

## Notes

The package's interval arithmetic is validated against published and
external-source references in the test suite. The suite has exactly one skipped
test, on every platform: a documented placeholder in `test-vc-extract.R` for
reproducing one published degrees-of-freedom value (Francq, Lin and Hoyer 2019,
Section 5.3) exactly, which needs the original Hoffman-Kringle raw assay
values rather than the constructed data used in the suite. The same code path
is exercised by a range check in the preceding test.
