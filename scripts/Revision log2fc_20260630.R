###############################################################################
## Revision analysis : log2 fold-change of 7 checkpoint genes  (Reviewer 2-1)
## Structure-aware version:
##   - checkpoint genes are taken from the FULL expression matrix
##     expression_full_log2tpm_wide.csv  (rows = samples[pair_id], cols = genes)
##   - joined to the grade-4 frame by pair_id to attach tp53_status / source
##   - GLASS: LAG3 from glass_analysis_dataset.rds; other 6 from GLASS matrix
## Output -> 20260630/R2-1_log2FC
###############################################################################
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  library(data.table)
  have_ggplot <- requireNamespace("ggplot2", quietly = TRUE); if (have_ggplot) library(ggplot2)
})

orig_dir <- here::here("results", "TP53", "20260221")
out_dir  <- here::here("results", "TP53", "20260630", "R2-1_log2FC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

checkpoint <- c("LAG3","PDCD1","CD274","PDCD1LG2","CTLA4","TIGIT","HAVCR2")
cliffs_delta <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))

## ------------------------------------------------------------------ ##
## 1. GDC grade 4 : pull 7 genes (log2tpm) from full matrix, attach TP53
## ------------------------------------------------------------------ ##
gdc  <- readRDS(file.path(orig_dir, "30_screening/analysis_dataset_grade4.rds"))
expr_file <- file.path(orig_dir, "27_expression_matrix/expression_full_log2tpm_wide.csv")

# read only pair_id + the 7 checkpoint columns (avoid loading all 59k cols)
hdr   <- names(fread(expr_file, nrows = 0, check.names = FALSE))
found <- checkpoint[checkpoint %in% hdr]
miss  <- setdiff(checkpoint, hdr)
if (length(miss)) cat("Genes NOT in expression matrix:", paste(miss, collapse = ", "), "\n")
expr  <- fread(expr_file, select = c("pair_id", found), check.names = FALSE, data.table = FALSE)

# values are log2(TPM+1); recover TPM for an interpretable median-based log2FC
key <- "pair_id"
meta <- gdc[, c(key, "tp53_status", "source")]
dat  <- merge(meta, expr, by = key)
dat$grp <- ifelse(grepl("^mut", tolower(dat$tp53_status)), "Mutant", "Wild-type")
cat(sprintf("GDC grade4 matched: %d (Mutant %d / WT %d)\n",
            nrow(dat), sum(dat$grp == "Mutant"), sum(dat$grp == "Wild-type")))

per_gene <- function(g, d) {
  l <- as.numeric(d[[g]]); mu <- l[d$grp == "Mutant"]; wt <- l[d$grp == "Wild-type"]
  tpm_mu <- 2^mu - 1; tpm_wt <- 2^wt - 1            # back-transform from log2(TPM+1)
  data.frame(
    gene             = g,
    n_mut            = sum(!is.na(mu)), n_wt = sum(!is.na(wt)),
    log2FC_mean      = mean(mu, na.rm = TRUE) - mean(wt, na.rm = TRUE),   # diff of mean log2(TPM+1)
    median_diff_log2 = median(mu, na.rm = TRUE) - median(wt, na.rm = TRUE), # matches manuscript
    log2FC_medianTPM = log2((median(tpm_mu, na.rm = TRUE) + 1) /
                              (median(tpm_wt, na.rm = TRUE) + 1)),          # interpretable FC
    cliffs_delta     = cliffs_delta(mu[!is.na(mu)], wt[!is.na(wt)]),
    wilcox_p         = suppressWarnings(wilcox.test(mu, wt)$p.value)
  )
}
res <- do.call(rbind, lapply(found, per_gene, d = dat))
res$wilcox_p_BH <- p.adjust(res$wilcox_p, method = "BH")
res <- res[order(-res$log2FC_mean), ]
cat("\n=== log2FC of checkpoint genes (GDC grade 4, TP53-mut vs WT) ===\n")
print(res, row.names = FALSE)
fwrite(res, file.path(out_dir, "checkpoint_log2FC_GDC_grade4.csv"))

## ------------------------------------------------------------------ ##
## 2. GLASS : LAG3 from frame; other 6 from GLASS expression matrix
## ------------------------------------------------------------------ ##
glass <- readRDS(file.path(orig_dir, "32_glass_validation_screening/glass_analysis_dataset.rds"))
glass$grp <- ifelse(grepl("wt|wild", tolower(glass$tp53_status)), "Wild-type", "Mutant")

gexpr_file <- file.path(orig_dir, "27_expression_matrix/glass_expression_full_log2tpm_wide.csv")
ghdr  <- names(fread(gexpr_file, nrows = 0, check.names = FALSE))
gfound<- checkpoint[checkpoint %in% ghdr]
gexpr <- fread(gexpr_file, select = c("pair_id", gfound), check.names = FALSE, data.table = FALSE)
gdat  <- merge(glass[, c("pair_id", "tp53_status", "grp")], gexpr, by = "pair_id")

# LAG3 in the frame is already log2tpm; prefer matrix for uniformity if present
res_glass <- do.call(rbind, lapply(gfound, per_gene, d = gdat))
res_glass$wilcox_p_BH <- p.adjust(res_glass$wilcox_p, method = "BH")
res_glass <- res_glass[order(-res_glass$log2FC_mean), ]
cat("\n=== log2FC of checkpoint genes (GLASS, TP53-mut vs WT) ===\n")
print(res_glass, row.names = FALSE)
fwrite(res_glass, file.path(out_dir, "checkpoint_log2FC_GLASS.csv"))

## ------------------------------------------------------------------ ##
## 3. Figure : log2FC barplot (both cohorts), LAG3 highlighted
## ------------------------------------------------------------------ ##
if (have_ggplot) {
  res$cohort <- "GDC"; res_glass$cohort <- "GLASS"
  both <- rbind(res[, c("gene","log2FC_mean","wilcox_p_BH","cohort")],
                res_glass[, c("gene","log2FC_mean","wilcox_p_BH","cohort")])
  both$gene <- factor(both$gene, levels = checkpoint)
  p <- ggplot(both, aes(gene, log2FC_mean, fill = wilcox_p_BH < 0.05)) +
    geom_col(width = 0.7) + coord_flip() + facet_wrap(~cohort) +
    scale_fill_manual(values = c("FALSE"="grey70","TRUE"="#fc9272"), name="BH p < 0.05") +
    labs(x = NULL, y = "log2 fold-change (TP53-mut vs WT)",
         title = "Checkpoint gene log2FC by cohort") +
    theme_classic(base_size = 13) +
    theme(axis.text = element_text(color = "black", size = 12),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "checkpoint_log2FC_both_cohorts.pdf"), p, width = 7, height = 4, dpi = 450)
  ggsave(file.path(out_dir, "checkpoint_log2FC_both_cohorts.png"), p, width = 7, height = 4, dpi = 450)
}

cat("\nDone. Outputs in:\n", out_dir, "\n")
