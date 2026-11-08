# =============================================================================
# Modified large-sample (MLS) upper confidence bound on a linear combination of
# expected mean squares,  theta = sum_j k_j E[MS_j],  MS_j independent with
# MS_j * df_j / E[MS_j] ~ chi^2(df_j). Method-pure (base R).
#
# * All k_j >= 0 -> Graybill & Wang (1980):
#     U = theta_hat + sqrt( sum_j H_j^2 k_j^2 MS_j^2 ),  H_j = df_j / chi2_{a}(df_j) - 1.
# * Any k_j < 0  -> Ting, Burdick & Graybill (1990), Burdick & Graybill (1992)
#   notation (F_{d:a,b} has area d to the RIGHT). With q over the positive and r
#   over the negative coefficients, d_r = -k_r:
#     U = theta_hat + sqrt( sum_q H_q^2 k_q^2 MS_q^2 + sum_r G_r^2 d_r^2 MS_r^2
#                           + sum_q sum_r H_qr k_q d_r MS_q MS_r ),
#     G_r  = 1 - df_r / chi2_{1-a}(df_r)          (= 1 - 1/F_{a:df_r,inf})
#     F_qr = F_{1-a:df_q,df_r}  (= qf(a, df_q, df_r))
#     H_qr = ((1 - F_qr)^2 - H_q^2 F_qr^2 - G_r^2) / F_qr,
#   the sqrt argument floored at 0. For one positive and one negative term this
#   is exactly the reference implementation (reintervals-simstudy
#   pilot/methods.R, mls_upper()); with several terms it is the sum over all
#   (positive, negative) pairs.
#
# "Negative" means k_j < -tol * max|k| (tol = 1e-12). The unbalanced EMS
# synthesis returns k_AB = -2.8e-16 for 2 x 2 crossed tables whose exact k_AB is
# 0; a bare k < 0 test would switch those to the TBG bound on round-off.
# =============================================================================

## Returns list(U, excess, bound): U = estimate + excess, excess = the sqrt term
## (kept separate so the Graybill-Wang arithmetic of ti_francq() is unchanged
## bit for bit), bound = "graybill-wang" or "ting-burdick-graybill".
.mls_upper <- function(ms, k, df, conf, tol = 1e-12) {
  a <- 1 - conf
  est <- sum(k * ms)
  neg <- k < -tol * max(abs(k))
  if (!any(neg)) {
    h <- df / stats::qchisq(a, df = df) - 1
    excess <- sqrt(sum(h^2 * k^2 * ms^2))
    return(list(U = est + excess, excess = excess, bound = "graybill-wang"))
  }
  pos <- !neg
  Hq <- df[pos] / stats::qchisq(a, df = df[pos]) - 1
  cq <- k[pos]; Sq <- ms[pos]; nq <- df[pos]
  Gr <- 1 - df[neg] / stats::qchisq(1 - a, df = df[neg])
  dr <- -k[neg]; Sr <- ms[neg]; nr <- df[neg]
  cross <- 0
  for (i in seq_along(cq)) for (j in seq_along(dr)) {
    Fqr <- stats::qf(a, nq[i], nr[j])
    Hqr <- ((1 - Fqr)^2 - Hq[i]^2 * Fqr^2 - Gr[j]^2) / Fqr
    cross <- cross + Hqr * cq[i] * dr[j] * Sq[i] * Sr[j]
  }
  VU <- sum((Hq * cq * Sq)^2) + sum((Gr * dr * Sr)^2) + cross
  excess <- sqrt(max(VU, 0))
  list(U = est + excess, excess = excess, bound = "ting-burdick-graybill")
}
