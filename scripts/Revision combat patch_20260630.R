###############################################################################
## PATCH: finish the ComBat analysis using objects already in the session
## Requires Mcb, Msub, meta from revision_ComBat.R (same session).
## Fixes the summarization step (robust coefficient extraction by name).
###############################################################################
suppressPackageStartupMessages(library(data.table))
stopifnot(exists("Mcb"), exists("Msub"), exists("meta"))

out_dir <- here::here("results", "TP53", "20260630", "R1-4_ComBat")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

APM_genes   <- c("B2M","TAP1","TAP2","TAPBP","HLA-A","HLA-B","HLA-C","NLRC5")
Tcell_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMA","GZMB","PRF1")
score <- function(m, g) colMeans(m[intersect(g, rownames(m)), , drop = FALSE])

make_df <- function(mat) data.frame(
  meta, LAG3 = mat["LAG3", ],
  APM = score(mat, APM_genes), Tcell = score(mat, Tcell_genes)
)
df_raw <- make_df(Msub)     # same gene subset, uncorrected
df_cb  <- make_df(Mcb)      # ComBat-corrected

# robust extractor: pull the TP53 row by name, build a clean one-row data.frame
get_tp53 <- function(m) {
  co <- summary(m)$coefficients
  rn <- grep("tp53", rownames(co), value = TRUE)[1]
  ci <- confint(m)[rn, ]
  data.frame(beta = unname(co[rn, "Estimate"]),
             se   = unname(co[rn, "Std. Error"]),
             lo   = unname(ci[1]),
             hi   = unname(ci[2]),
             p    = unname(co[rn, "Pr(>|t|)"]))
}

fit_tp53 <- function(d, label) {
  m0 <- lm(LAG3 ~ tp53 + source + idh, data = d)
  m1 <- lm(LAG3 ~ tp53 + source + idh + Tcell + APM, data = d)
  rbind(
    data.frame(data = label, model = "M0 (TP53+source+IDH)", get_tp53(m0)),
    data.frame(data = label, model = "M1 (+Tcell+APM)",      get_tp53(m1))
  )
}

res <- rbind(fit_tp53(df_raw, "raw"), fit_tp53(df_cb, "ComBat"))
rownames(res) <- NULL

# Delta beta (ComBat vs raw) per model
wide <- merge(
  res[res$data == "raw",    c("model","beta")],
  res[res$data == "ComBat", c("model","beta")],
  by = "model", suffixes = c("_raw","_ComBat"))
wide$delta_beta_pct <- 100 * (wide$beta_ComBat - wide$beta_raw) / wide$beta_raw

cat("=== TP53 coefficient for LAG3: raw vs ComBat-corrected ===\n")
print(res, row.names = FALSE, digits = 4)
cat("\n=== Change in TP53 beta (ComBat vs raw) ===\n")
print(wide, row.names = FALSE, digits = 4)

fwrite(res,  file.path(out_dir, "ComBat_TP53_coefficients.csv"))
fwrite(wide, file.path(out_dir, "ComBat_delta_beta.csv"))
cat("\nDone. Outputs in:\n", out_dir, "\n")
