# Contributing to reintervals

Thanks for your interest in `reintervals`. The package computes confidence,
prediction, tolerance, and new-group-mean intervals for linear mixed models, and
every shipped interval is validated against a published or external-source
oracle. Contributions that preserve that validation discipline are welcome.

## Reporting issues

Please open an issue at
<https://github.com/percevalsondag/reintervals/issues>. A good bug report
includes:

- a minimal reproducible example (a small `lme4::lmer()` fit and the verb call),
- what you expected and what you got (paste the `lmm_interval` output),
- your `sessionInfo()` (R, `lme4`, and `reintervals` versions).

For a numerical discrepancy, say which engine (`method`) produced the interval
and, if you can, what you believe the correct value is and its source.

## Proposing changes

1. Open an issue first to discuss anything beyond a small fix.
2. Fork, branch from `main`, and keep changes focused.
3. Add or update tests (see below) — a change that moves a number must show why,
   and a new capability must come with an oracle to test it against.
4. Open a pull request describing the change and the validation.

## Development workflow

```r
# from the package root
devtools::load_all()        # load for interactive work
devtools::document()        # regenerate man/ + NAMESPACE after roxygen edits
devtools::test()            # run the test suite
devtools::check()           # full R CMD check (use --as-cran before submitting)
```

The test suite checks the interval **arithmetic** against the reference
implementations to ~`1e-9` (on committed fixed components, so the checks are
platform-portable), and reproduces the published point numbers (the
Montes-Burdick-Leblond Appendix A interval and the Oliva-Aviles & Hauser content
variance). Please keep these gates green.

## Continuous integration

Every push and pull request is checked by `R CMD check` on Linux (devel /
release / oldrel-1), macOS, and Windows via GitHub Actions. A contribution is
expected to be green on all five before merge.

## Code of conduct

By participating in this project you agree to abide by its
[Code of Conduct](CODE_OF_CONDUCT.md).
