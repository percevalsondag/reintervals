## Submission type

New submission (first release of `reintervals`).

## Test environments

- Local: macOS (aarch64-apple-darwin20), R 4.5.3.
- GitHub Actions: ubuntu-latest (R devel, R release, R oldrel-1),
  macOS-latest (R release), windows-latest (R release) --- dependencies installed
  from current CRAN.

## R CMD check results

`R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes.

(On a local run without network access, a single NOTE
`checking for future file timestamps ... unable to verify current time` can
appear; it is an environment artifact from the check host being unable to reach
the time-verification service and does not appear on the time-synchronized CI
runners. It is not related to the package.)

## Reverse dependencies

None --- this is a new submission.

## Notes

The package's interval arithmetic is validated against published and
external-source references in the test suite. A few end-to-end, fit-dependent
checks are skipped on CI because they depend on `lme4`'s platform-specific
numerical behavior; the platform-portable checks (comparing committed extracted
components against the references) run on every platform.
