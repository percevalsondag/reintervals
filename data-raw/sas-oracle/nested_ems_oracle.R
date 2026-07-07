# =============================================================================
# SAS EMS oracle hand-off for the nested tolerance-interval synthesis.
#
# Not shipped (data-raw/ is .Rbuildignored). Self-contained + reproducible:
# builds the oracle datasets, prints what reintervals' synthesis produces for
# each (computed from the data, NOT pre-loaded from SAS), and emits ready-to-run
# SAS programs whose `datalines` are generated from the SAME data frames (so the
# R and SAS datasets are provably identical, row-for-row).
#
# WHY this exists (see SPEC 6.1): nested uses the SEQUENTIAL (Type 1) EMS
# decomposition, crossed uses Type III. Type III is degenerate for nested
# (equal sigma2_a/sigma2_b coeffs in EMS_a force k_b(a)=0; MLS Eq.26 needs
# independent mean squares). The SAS Type 1 result is the external oracle; the
# synthesis must not certify itself. Tier-4 nested oracle test is added only
# AFTER the SAS Type 1 numbers are pasted back and match (by what each
# coefficient multiplies, label-reconciled). If Type 1 does not match -> STOP
# and flag a real bug.
#
# Usage:  Rscript data-raw/sas-oracle/nested_ems_oracle.R   (run from package root)
# =============================================================================

suppressMessages(pkgload::load_all(".", quiet = TRUE))

emit_sas <- function(dat, dataname, procs) {
  cat(sprintf("\n--- SAS program (%s) -----------------------------------\n", dataname))
  cat(sprintf("data %s;\n  input a $ b $ y;\n  datalines;\n", dataname))
  for (i in seq_len(nrow(dat))) cat(sprintf("%s %s %s\n", dat$a[i], dat$b[i], format(dat$y[i], nsmall = 1)))
  cat(";\nrun;\n")
  cat(procs)
}

ems_line <- function(syn) sprintf(
  "EMS_a = e + %.4f b(a) + %.4f a ;  EMS_b(a) = e + %.4f b(a) ;  df a/b(a)/resid = %g/%g/%g",
  syn$E["A", "B"], syn$E["A", "A"], syn$E["B", "B"],
  syn$df["A"], syn$df["B"], syn$df["Residual"])

# ---- CASE 1: PRIMARY oracle -- small UNBALANCED nested ----------------------
# a1: b1(3),b2(2) | a2: b3(2),b4(1),b5(2) | a3: b6(1),b7(3)   (A=3, K=7, N=14)
unbal <- data.frame(
  a = c("a1","a1","a1","a1","a1","a2","a2","a2","a2","a2","a3","a3","a3","a3"),
  b = c("b1","b1","b1","b2","b2","b3","b3","b4","b5","b5","b6","b7","b7","b7"),
  y = c(9.8,10.2,10.0,11.0,10.6,12.1,11.9,12.5,11.5,11.7,13.5,14.2,14.0,13.8))
cat("CASE 1 (UNBALANCED, primary oracle) -- reintervals sequential synthesis:\n  ",
    ems_line(ems_synthesis(list(c(3,2), c(2,1,2), c(1,3)), "nested")), "\n")
cat("  EXPECT SAS METHOD=TYPE1 a-row to match: 2.3071 b(a), 4.6429 a\n")
cat("  EXPECT SAS METHOD=TYPE3 a-row degenerate: 2.625 b(a), 2.625 a -> k_b(a)=0 (why type3 is rejected)\n")
emit_sas(unbal, "nested_unbal", paste0(
  "\n/* Type 1 (sequential) = the package's nested EMS; the ORACLE */\n",
  "proc mixed data=nested_unbal method=type1; class a b; model y = ; random a b(a); run;\n",
  "/* Type 3 = degenerate for nested (expect 2.625/2.625, k_b(a)=0) */\n",
  "proc mixed data=nested_unbal method=type3; class a b; model y = ; random a b(a); run;\n"))

# ---- CASE 2: OPTIONAL anchor -- balanced nested (Type1 = Type3) -------------
# A=3 coarse, B=2 fine each, n=3 reps (N=18). Anchors the balanced nested/crossed
# row (currently formula-only) to an external number at near-zero cost.
bal <- data.frame(
  a = rep(c("a1","a2","a3"), each = 6),
  b = rep(c("b1","b2","b3","b4","b5","b6"), each = 3),
  y = c(10.1,9.8,10.3,11.0,10.7,11.2,12.2,11.9,12.4,13.1,12.8,13.3,9.1,8.8,9.4,10.0,9.7,10.2))
cat("\nCASE 2 (BALANCED, optional anchor) -- reintervals synthesis:\n  ",
    ems_line(ems_synthesis(list(c(3,3), c(3,3), c(3,3)), "nested")), "\n")
cat("  EXPECT SAS METHOD=TYPE1 (= TYPE3 for balanced) a-row: 3 b(a), 6 a\n")
emit_sas(bal, "nested_bal", paste0(
  "\n/* balanced -> Type1 = Type3; expect a-row 3 b(a), 6 a */\n",
  "proc mixed data=nested_bal method=type1; class a b; model y = ; random a b(a); run;\n"))

cat("\n--- Copy back: Source / DF / Expected Mean Square for a, b(a), Residual from each table. ---\n")
