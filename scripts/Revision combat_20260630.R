###############################################################################
## Revision analysis : ComBat batch-correction sensitivity  (Reviewer 1-4)
## Question: was any batch correction / normalization harmonization attempted
## across the GDC sub-cohorts (TCGA vs CPTAC/HCMI)?
##
## Design:
##   - GDC grade-4 expression (rows = samples[pair_id], cols = genes) -> subset
##     442 samples -> transpose to genes x samples
##   - Apply ComBat (sva) with batch = source (TCGA / CPTAC_HCMI),
##     PRESERVING biology via mod = ~ tp53_status + idh_status
##   - Refit the TP53 model on ComBat-corrected LAG3 and compare the TP53
##     coefficient to the original (uncorrected) model: M0 and M1
##   - Immune scores (Tcell, APM) recomputed from ComBat-corrected genes so the
##     adjustment is internally consistent
## Output -> 20260630/R1-4_ComBat
###############################################################################
suppressPackageStartupMessages({
  inst <- function(p, bioc = FALSE) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) { if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(p, update = FALSE, ask = FALSE) }
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  inst("data.table"); inst("sva", bioc = TRUE)
  library(data.table); library(sva)
})

orig_dir <- here::here("results", "TP53", "20260221")
out_dir  <- here::here("results", "TP53", "20260630", "R1-4_ComBat")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_HVG <- 5000   # number of high-variance genes used for ComBat (plus analysis genes)

## ------------------------------------------------------------------ ##
## 1. Load metadata (grade-4 frame) and full expression matrix
## ------------------------------------------------------------------ ##
gdc <- readRDS(file.path(orig_dir, "30_screening/analysis_dataset_grade4.rds"))
meta <- data.frame(
  pair_id = gdc$pair_id,
  tp53    = factor(ifelse(grepl("^mut", tolower(gdc$tp53_status)), "mutant", "wildtype"),
                   levels = c("wildtype","mutant")),
  idh     = factor(ifelse(grepl("^mut", tolower(gdc$idh_status)),  "mutant", "wildtype"),
                   levels = c("wildtype","mutant")),
  source  = factor(gdc$source)
)
cat("Samples:", nrow(meta), "| batches:\n"); print(table(meta$source))

expr_file <- file.path(orig_dir, "27_expression_matrix/expression_full_log2tpm_wide.csv")
E <- fread(expr_file, check.names = FALSE, data.table = FALSE)   # rows = samples, col1 = pair_id
E <- E[E$pair_id %in% meta$pair_id, ]
E <- E[match(meta$pair_id, E$pair_id), ]                          # align to meta order
stopifnot(all(E$pair_id == meta$pair_id))
genes_all <- setdiff(names(E), "pair_id")
M <- t(as.matrix(E[, genes_all]))                                 # genes x samples
colnames(M) <- meta$pair_id
cat("Expression matrix (genes x samples):", nrow(M), "x", ncol(M), "\n")

## ------------------------------------------------------------------ ##
## 2. Choose genes for ComBat: top-variance HVGs + all analysis genes
## ------------------------------------------------------------------ ##
analysis_genes <- c("LAG3",
                    "B2M","TAP1","TAP2","TAPBP","HLA-A","HLA-B","HLA-C","NLRC5",            # APM
                    "STAT1","IRF1","IRF9","CXCL9","CXCL10","CXCL11","GBP1","GBP2","GBP4","GBP5","IDO1", # IFNg
                    "CD3D","CD3E","CD3G","CD8A","CD8B","GZMA","GZMB","PRF1")                # Tcell
analysis_genes <- analysis_genes[analysis_genes %in% rownames(M)]

# drop all-constant genes, then rank by variance
v   <- apply(M, 1, var)
keep_hvg <- names(sort(v[v > 0], decreasing = TRUE))[seq_len(min(N_HVG, sum(v > 0)))]
keep <- union(keep_hvg, analysis_genes)
Msub <- M[keep, , drop = FALSE]
cat("Genes used for ComBat:", nrow(Msub), "(HVG", length(keep_hvg),
    "+ analysis", length(analysis_genes), ")\n")

## ------------------------------------------------------------------ ##
## 3. Run ComBat (preserve TP53 + IDH)
## ------------------------------------------------------------------ ##
mod <- model.matrix(~ tp53 + idh, data = meta)
set.seed(42)
Mcb <- ComBat(dat = Msub, batch = meta$source, mod = mod,
              par.prior = TRUE, prior.plots = FALSE)

## ------------------------------------------------------------------ ##
## 4. Build analysis frames: BEFORE (raw) vs AFTER (ComBat) and refit
## ------------------------------------------------------------------ ##
score <- function(mat, genes) colMeans(mat[intersect(genes, rownames(mat)), , drop = FALSE])
APM_genes   <- c("B2M","TAP1","TAP2","TAPBP","HLA-A","HLA-B","HLA-C","NLRC5")
Tcell_genes <- c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMA","GZMB","PRF1")

make_df <- function(mat) data.frame(
  meta,
  LAG3  = mat["LAG3", ],
  APM   = score(mat, APM_genes),
  Tcell = score(mat, Tcell_genes)
)
df_raw <- make_df(Msub)     # raw (subset, same genes) for an apples-to-apples comparison
df_cb  <- make_df(Mcb)      # ComBat-corrected

fit_tp53 <- function(d, label) {
  m0 <- lm(LAG3 ~ tp53 + source + idh, data = d)                       # M0
  m1 <- lm(LAG3 ~ tp53 + source + idh + Tcell + APM, data = d)         # M1
  g  <- function(m) { s <- summary(m)$coefficients; ci <- confint(m)["tp53mutant", ]
  c(beta = s["tp53mutant","Estimate"], se = s["tp53mutant","Std. Error"],
    p = s["tp53mutant","Pr(>|t|)"], lo = ci[1], hi = ci[2]) }
  rbind(
    data.frame(model = "M0 (TP53+source+IDH)",            t(g(m0))),
    data.frame(model = "M1 (+Tcell+APM)",                 t(g(m1)))
  ) -> r
  r$data <- label; r
}

res <- rbind(fit_tp53(df_raw, "raw"), fit_tp53(df_cb, "ComBat"))
res <- res[, c("data","model","beta","se","lo","hi","p")]

# Δβ relative to raw, per model
wide <- reshape(res[, c("data","model","beta")], idvar = "model",
                timevar = "data", direction = "wide")
wide$delta_beta_pct <- 100 * (wide$beta.ComBat - wide$beta.raw) / wide$beta.raw

cat("\n=== TP53 coefficient for LAG3: raw vs ComBat-corrected ===\n")
print(res, row.names = FALSE, digits = 4)
cat("\n=== Change in TP53 beta (ComBat vs raw) ===\n")
print(wide, row.names = FALSE, digits = 4)

## ------------------------------------------------------------------ ##
## 5. Save outputs
## ------------------------------------------------------------------ ##
fwrite(res,  file.path(out_dir, "ComBat_TP53_coefficients.csv"))
fwrite(wide, file.path(out_dir, "ComBat_delta_beta.csv"))

cat("\nInterpretation: if the ComBat TP53 beta is close to the raw beta (small |delta_beta_pct|)\n",
    "and remains significant, the TP53-LAG3 association is robust to explicit cross-cohort\n",
    "batch correction, supporting that cohort-as-covariate already harmonized the data.\n", sep = "")
cat("\nDone. Outputs in:\n", out_dir, "\n")
