###############################################################################
## CPTAC-GBM: LAG3 protein abundance by TP53 mutation status
## NOA-D-26-00229 revision  (Reviewer #1 Comment 1 / Reviewer #2 Comment 2)
##
## Workflow (agreed with Dr. Ishikawa):
##   STEP 1  Check whether LAG3 is detected at an analyzable level (missingness)
##   STEP 2  If detected -> compare LAG3 protein between TP53-mut vs WT
##           (same framework as the manuscript RNA analysis:
##            Wilcoxon, Cliff's delta, Hodges-Lehmann, median difference)
##   STEP 3  Output supplementary-ready figure + results table
##
## Input : cBioPortal study gbm_cptac_2021 (downloaded tarball)
## Run   : just press Source in RStudio. Only edit `tar_path` if needed.
###############################################################################

## ------------------------------------------------------------------ ##
## 0. Settings & packages
## ------------------------------------------------------------------ ##
# Script location (for reference):
#   D:/Projects/GBM_Analysis/scripts/TP53/20260630
# Results root (input tarball lives here; output folder is created here):
results_root <- here::here("results", "TP53", "20260630")
tar_path     <- file.path(results_root, "gbm_cptac_2021.tar.gz")

# Count a sample as TP53-mutant only for non-silent variants (set FALSE to count any)
EXCLUDE_SILENT <- TRUE

# Output directory = NEW folder inside the results root
out_dir <- file.path(results_root, "CPTAC_LAG3_protein_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Packages (installed only if missing)
need <- c("data.table", "ggplot2")
for (p in need) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

set.seed(42)  # match manuscript

## ------------------------------------------------------------------ ##
## 1. Extract tarball & locate files
## ------------------------------------------------------------------ ##
extract_root <- file.path(dirname(tar_path), "gbm_cptac_2021_extracted")
if (!dir.exists(extract_root)) {
  dir.create(extract_root, showWarnings = FALSE)
  untar(tar_path, exdir = extract_root)
}
# Find the directory that actually holds the data files
all_files  <- list.files(extract_root, recursive = TRUE, full.names = TRUE)
study_dir  <- dirname(all_files[grepl("meta_study\\.txt$", all_files)][1])
if (is.na(study_dir)) study_dir <- dirname(
  all_files[grepl("data_protein_quantification", all_files)][1])

# Select the GLOBAL proteome file (exclude PTM-specific & z-score files)
prot_all <- list.files(study_dir, pattern = "protein_quantification",
                       full.names = TRUE)
cat("Proteomics files found in study:\n"); print(basename(prot_all))
ptm_or_z <- "zscore|phospho|acetyl|ubiquit|glyco|sumoyl|methyl|redox"
prot_candidates <- prot_all[!grepl(ptm_or_z, basename(prot_all), ignore.case = TRUE)]
exact <- prot_candidates[basename(prot_candidates) == "data_protein_quantification.txt"]
prot_file <- if (length(exact)) exact[1] else prot_candidates[1]
if (is.na(prot_file) || length(prot_file) == 0)
  stop("No global proteome file found among: ",
       paste(basename(prot_all), collapse = ", "))
using_zscore <- isTRUE(grepl("zscore", basename(prot_file), ignore.case = TRUE))

mut_file <- list.files(study_dir, pattern = "data_mutations",
                       full.names = TRUE)[1]

cat("Study dir :", study_dir, "\n")
cat("Proteomics:", basename(prot_file),
    if (using_zscore) "(z-scores — units differ)" else "", "\n")
cat("Mutations :", basename(mut_file), "\n\n")

## ------------------------------------------------------------------ ##
## 2. Read proteomics, isolate LAG3, sample columns
## ------------------------------------------------------------------ ##
prot <- fread(prot_file, check.names = FALSE, data.table = FALSE)

# Identifier columns vs sample columns (drop known ID columns; rest = samples)
id_like <- c("Hugo_Symbol","Entrez_Gene_Id","Composite.Element.REF",
             "ENTITY_STABLE_ID","NAME","GENE_SYMBOL","Gene","gene","Description")
id_cols     <- intersect(id_like, names(prot))
if (length(id_cols) == 0) id_cols <- names(prot)[1]   # fallback: 1st column = identifier
sample_cols <- setdiff(names(prot), id_cols)

# Locate the gene-symbol column that contains LAG3 (allow synonyms)
gene_syn <- c("LAG3", "LAG-3", "CD223")
gene_col <- id_cols[sapply(id_cols, function(cc)
  any(prot[[cc]] %in% gene_syn, na.rm = TRUE))][1]
lag3_rows <- if (!is.na(gene_col))
  prot[prot[[gene_col]] %in% gene_syn, , drop = FALSE] else prot[0, , drop = FALSE]

# Graceful exit if LAG3 is simply not in the global proteome (= the STEP 1 answer)
if (nrow(lag3_rows) == 0) {
  cat("\n==================== STEP 1: LAG3 detection ====================\n")
  cat("LAG3 was NOT found in the global proteome file (",
      basename(prot_file), ").\n", sep = "")
  cat("0 of ", length(sample_cols), " samples have LAG3 protein quantified.\n\n", sep = "")
  cat("Interpretation: LAG3 is a low-abundance immune cell-surface protein and is\n",
      "below the detection limit of bulk TMT mass spectrometry in this cohort.\n",
      "Recommended response: address Reviewer #1 C1 / #2 C2 by TEXT -- report that\n",
      "LAG3 protein was not reliably quantified and cite the MS detection limit as a\n",
      "Discussion limitation, rather than adding an underpowered supplementary figure.\n", sep = "")
  fwrite(data.frame(
    metric = c("global_proteome_file", "samples_in_matrix", "LAG3_detected", "detection_rate"),
    value  = c(basename(prot_file), length(sample_cols), 0, 0)),
    file.path(out_dir, "LAG3_detection_summary.csv"))
  cat("\nWrote LAG3_detection_summary.csv to:\n", out_dir, "\n")
  stop("LAG3 not present in the global proteome (this null result IS the finding; see message above).")
}

# Coerce sample values to numeric
lag3_mat <- suppressWarnings(
  apply(lag3_rows[, sample_cols, drop = FALSE], 2, as.numeric))
if (is.null(dim(lag3_mat))) lag3_mat <- matrix(lag3_mat, nrow = nrow(lag3_rows))

# If multiple LAG3 rows (isoforms), keep the one with the most measurements
keep <- which.max(rowSums(!is.na(lag3_mat)))
lag3 <- setNames(as.numeric(lag3_mat[keep, ]), sample_cols)

## ------------------------------------------------------------------ ##
## STEP 1 :  LAG3 detection / missingness
## ------------------------------------------------------------------ ##
n_total    <- length(lag3)
n_detected <- sum(!is.na(lag3))
det_rate   <- n_detected / n_total
cat("==================== STEP 1: LAG3 detection ====================\n")
cat(sprintf("Samples in proteomics matrix : %d\n", n_total))
cat(sprintf("Samples with LAG3 quantified : %d (%.1f%%)\n",
            n_detected, 100 * det_rate))
cat(sprintf("Missing LAG3                  : %d (%.1f%%)\n\n",
            n_total - n_detected, 100 * (1 - det_rate)))

## ------------------------------------------------------------------ ##
## 3. Read mutations -> TP53 status per sample
## ------------------------------------------------------------------ ##
# 'skip="Hugo_Symbol"' makes fread start at the header, skipping any #version line
mut <- fread(mut_file, skip = "Hugo_Symbol", check.names = FALSE, data.table = FALSE)

silent_classes <- c("Silent","Intron","3'UTR","5'UTR","3'Flank","5'Flank",
                    "IGR","RNA","lincRNA")
tp53 <- mut[mut$Hugo_Symbol == "TP53", , drop = FALSE]
if (EXCLUDE_SILENT && "Variant_Classification" %in% names(mut))
  tp53 <- tp53[!tp53$Variant_Classification %in% silent_classes, , drop = FALSE]

tp53_mut_samples <- unique(tp53$Tumor_Sample_Barcode)
sequenced        <- unique(mut$Tumor_Sample_Barcode)   # profiled samples

# Optional: use the official sequenced case list if present (more complete WT set)
seq_list <- list.files(file.path(study_dir, "case_lists"),
                       pattern = "sequenced", full.names = TRUE)
if (length(seq_list) >= 1) {
  ln <- readLines(seq_list[1])
  ids <- ln[grepl("^case_list_ids", ln)]
  if (length(ids))
    sequenced <- union(sequenced,
                       strsplit(sub("^case_list_ids:\\s*", "", ids), "\t")[[1]])
}

## ------------------------------------------------------------------ ##
## 4. Match proteomics <-> mutation samples, build analysis frame
## ------------------------------------------------------------------ ##
prot_samples <- names(lag3)
common <- intersect(prot_samples, sequenced)

# Fallback: if direct ID match is poor, try trimming common suffixes
if (length(common) < 0.3 * length(prot_samples)) {
  trim <- function(x) sub("(-|\\.)(0\\d|T|Tumor)$", "", x)
  m1 <- setNames(prot_samples, trim(prot_samples))
  m2 <- setNames(sequenced,    trim(sequenced))
  hit <- intersect(names(m1), names(m2))
  common <- unname(m1[hit])
  tp53_mut_samples <- trim(tp53_mut_samples)
  names(lag3) <- trim(names(lag3))
}

df <- data.frame(
  sample = common,
  LAG3   = lag3[common],
  TP53   = ifelse(common %in% tp53_mut_samples |
                    sub("(-|\\.)(0\\d|T|Tumor)$","",common) %in% tp53_mut_samples,
                  "Mutant", "Wild-type"),
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$LAG3), ]                       # need LAG3 value for comparison
df$TP53 <- factor(df$TP53, levels = c("Wild-type", "Mutant"))

n_mut <- sum(df$TP53 == "Mutant"); n_wt <- sum(df$TP53 == "Wild-type")
cat("============ Matched samples with LAG3 + TP53 status ============\n")
cat(sprintf("TP53-mutant : n = %d\nTP53 wild-type : n = %d\n\n", n_mut, n_wt))

## ------------------------------------------------------------------ ##
## STEP 2 :  TP53-mut vs WT comparison (only if both groups are usable)
## ------------------------------------------------------------------ ##
cliffs_delta <- function(a, b) {
  m <- outer(a, b, function(x, y) sign(x - y))
  mean(m)                                          # P(a>b) - P(a<b)
}

results <- NULL
if (n_mut >= 3 && n_wt >= 3) {
  x <- df$LAG3[df$TP53 == "Mutant"]; y <- df$LAG3[df$TP53 == "Wild-type"]
  wt <- wilcox.test(x, y, conf.int = TRUE)         # HL estimate via conf.int
  results <- data.frame(
    n_mut          = n_mut,
    n_wt           = n_wt,
    median_mut     = median(x),
    median_wt      = median(y),
    median_diff    = median(x) - median(y),        # mut - wt
    mean_diff      = mean(x) - mean(y),
    HodgesLehmann  = unname(wt$estimate),
    cliffs_delta   = cliffs_delta(x, y),
    wilcox_p       = wt$p.value,
    units          = ifelse(using_zscore, "z-score", "log2-ratio (CPTAC)")
  )
  cat("==================== STEP 2: Group comparison ====================\n")
  print(t(results))
  cat(sprintf("\nDirection: LAG3 is %s in TP53-mutant tumors.\n",
              ifelse(results$median_diff > 0, "HIGHER", "lower")))
} else {
  cat("==================== STEP 2: SKIPPED ====================\n")
  cat("Too few samples per group after requiring LAG3 detection.\n",
      "-> Protein-level comparison is underpowered; recommend addressing\n",
      "   Reviewer comments by text (state MS detection limit for low-abundance\n",
      "   immune surface proteins as a limitation in the Discussion).\n\n")
}

## ------------------------------------------------------------------ ##
## STEP 3 :  Figure + saved outputs
## ------------------------------------------------------------------ ##
if (!is.null(results)) {
  ann <- sprintf("median diff = %+.3f\nCliff's \u03b4 = %+.3f\nWilcoxon p = %s\nn = %d vs %d",
                 results$median_diff, results$cliffs_delta,
                 format.pval(results$wilcox_p, digits = 2), n_mut, n_wt)
  ylab <- ifelse(using_zscore, "LAG3 protein (z-score)",
                 "LAG3 protein abundance (log2-ratio)")
  
  p <- ggplot(df, aes(TP53, LAG3, fill = TP53)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.6) +
    scale_fill_manual(values = c("Wild-type" = "#9ecae1", "Mutant" = "#fc9272")) +
    annotate("text", x = 0.62, y = max(df$LAG3, na.rm = TRUE),
             label = ann, hjust = 0, vjust = 1, size = 4.2) +
    labs(x = NULL, y = ylab,
         title = "CPTAC-GBM: LAG3 protein by TP53 status") +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14),
          axis.text  = element_text(size = 13, color = "black"))
  
  ggsave(file.path(out_dir, "FigS13_LAG3_protein_TP53.pdf"), p,
         width = 5, height = 5, dpi = 450)
  ggsave(file.path(out_dir, "FigS13_LAG3_protein_TP53.png"), p,
         width = 5, height = 5, dpi = 450)
  fwrite(results, file.path(out_dir, "LAG3_protein_TP53_results.csv"))
}

# Always save the per-sample table and a detection summary (useful either way)
fwrite(df, file.path(out_dir, "LAG3_protein_per_sample.csv"))
fwrite(data.frame(
  metric = c("samples_in_matrix","LAG3_detected","detection_rate",
             "matched_with_TP53","n_TP53_mut","n_TP53_wt","exclude_silent"),
  value  = c(n_total, n_detected, round(det_rate, 3),
             nrow(df), n_mut, n_wt, EXCLUDE_SILENT)),
  file.path(out_dir, "LAG3_detection_summary.csv"))

cat("\nDone. Outputs written to:\n", out_dir, "\n")
cat("Files: LAG3_detection_summary.csv, LAG3_protein_per_sample.csv",
    if (!is.null(results)) ", LAG3_protein_TP53_results.csv, FigS13_*.pdf/.png" else "", "\n")
