# nist_ex4

Oracles `nist-recipe-ex4` and `hoffman2010-s52-nist-ex4` (`../oracles.csv`),
checked by `test-oracle-nist-ex4.R`.

- **Source.** NIST RECIPE (Vangel MG; RECIPE 1.0, April 1995: one-sided
  random-effect regression tolerance limits), NIST software, public domain.
  Archived homepage:
  <https://webarchive.library.unt.edu/web/20120924222228/http://www.itl.nist.gov/div898/software/recipe/homepage.html>
- **Files.** Retrieved by Perceval Sondag, 2026-11-22. Unmodified copies of the
  distribution files (no header was added; their bytes match the source).
  - `ex4.dat`: the RECIPE input of Example 4. A `#`-commented text file: the
    control line (ntot 72, nlvl 2, nbch 8, npar 2, npts 7, prob 0.90,
    conf 0.95), the two fixed levels (1 = 75 F, 2 = -67 F, as design rows
    `1 temp`), 72 rows `level batch strength`, then the 7 design rows of the
    points t0 = -67, -50, -25, 0, 25, 50, 75.
  - `ex4.out`: the RECIPE 1.0 output for it. The simulated pivot critical value
    file `ex4.crt` was absent, so RECIPE used its Satterthwaite approximation.
    Columns: probability, confidence, regression (the fitted line at t0), lower
    tolerance limit.
- **Hoffman (2010) section 5.2** reanalyses the same data. Its printed limits
  (analytical "NEW" and bootstrap-adjusted) are typed into the test, with the
  reference next to them: Hoffman D (2010). *Technometrics* 52(3):303-312.
  doi:10.1198/TECH.2010.09171.
- **SHA-256** at retrieval:
  - `ex4.dat`: `02f390eaf382dbe4de5f36151d95f260ece960c100868f2bac050ba7327e6e0b`
  - `ex4.out`: `6272b8b79000d1327c0dfc50ab6937bf6447898a4bb88b0c70a90e286313d441`
