## Resubmission / submission history

Version 0.1.0 was submitted and then withdrawn by the maintainer on 2026-10-27
(email to CRAN-submissions) before acceptance, because post-submission
validation found coverage defects in two interval functions. 0.2.0 corrects
them (the new-group-mean interval; a confidence-interval df that no longer
depends on which packages are installed), adds the Hoffman & Kringle (2005)
one-way engine and one-sided intervals, and is validated by simulation and
against published worked examples. The earlier CRAN feedback on 0.1.0
(restoring `options()`, not modifying `.GlobalEnv`) remains addressed.

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 3 notes.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Perceval Sondag <reintervals@gmail.com>'

New submission
```

The other two local NOTEs come from the local machine, not the package:
"checking for future file timestamps ... unable to verify current time" (no
time-server access) and "checking HTML version of manual ... Skipping checking
HTML validation" (local HTML Tidy too old).

- win-builder devel: <pending>
- win-builder release: <pending>
- mac-builder: <pending>

## Test environments

- Local: macOS 27.0 (aarch64-apple-darwin20), R 4.5.3.
- win-builder: R-devel and R-release.
- mac-builder: R-release.

## Downstream dependencies

None (new package).
