# =============================================================================
# NOA-D-26-00229 Revision (R1)
# Reviewer comment R1-3b: validate tumor-purity adjustment with an orthogonal,
#   DNA-based purity estimate (ABSOLUTE) in addition to the expression-based
#   ESTIMATE used in the manuscript.
#
# STEP 2 (THIS SCRIPT): join ABSOLUTE purity to the TCGA grade-4 samples and
#   refit the immune-adjusted base model, swapping ONLY the purity variable
#   (ESTIMATE vs ABSOLUTE) on the SAME samples, then compare Delta-beta.
#
# DESIGN (TCGA-only sensitivity analysis; ABSOLUTE is TCGA-only)
#   base   : LAG3_log2tpm ~ TP53 + IDH + Tcell + APM         (immune-adjusted)
#   m_est  : base + ESTIMATE purity   -> Delta-beta vs base
#   m_abs  : base + ABSOLUTE purity   -> Delta-beta vs base
#   All three fitted on the SAME rows (those with non-missing ABSOLUTE purity)
#   so the TP53 coefficients are directly comparable.
#   NOTE: This Delta-beta is computed within the TCGA subset, so it is NOT
#   numerically identical to the manuscript's full-GDC ESTIMATE Delta-beta
#   (-9.7%). The comparison of interest is ESTIMATE-vs-ABSOLUTE *within this
#   same subset* -> if both attenuations are similar, the purity-adjustment
#   conclusion is robust to the purity-estimation method.
#
# WORKFLOW NOTES
#   - Run in RStudio on the local machine, AFTER STEP 1. No downloads here.
#   - Requires the ABSOLUTE table already saved locally (see STEP 1 download).
# =============================================================================

set.seed(42)

## ---- 0. Paths (edit if your directory nesting differs) ----------------------
RESULTS_IN  <- here::here("results", "TP53", "20260221")
RESULTS_OUT <- here::here("results", "TP53", "20260630")
OUT_DIR     <- file.path(RESULTS_OUT, "R1-3b_ABSOLUTE")
RAW_DIR     <- file.path(OUT_DIR, "raw")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ABSOLUTE_RAW <- file.path(RAW_DIR, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
GRADE4_RDS   <- file.path(RESULTS_IN, "30_screening", "analysis_dataset_grade4.rds")

LOG_FILE     <- file.path(OUT_DIR, "step2_join_refit_log.txt")
sink(LOG_FILE, split = TRUE); on.exit(sink(), add = TRUE)

## ---- 1. Load grade-4 frame, subset TCGA ------------------------------------
stopifnot(file.exists(GRADE4_RDS))
g4 <- readRDS(GRADE4_RDS)
g4_tcga <- g4[g4$source == "TCGA", ]
cat("TCGA grade-4 n:", nrow(g4_tcga), "\n")

# LAG3 outcome column is the bare "LAG3" column, which already holds the
# log2(TPM+1) value (verified: LAG3 == log2(LAG3_tpm + 1)). The 27 immune-score
# genes carry explicit *_tpm / *_log2tpm columns; LAG3 (the outcome) does not.
lag3_col <- "LAG3"
stopifnot(lag3_col %in% names(g4_tcga), "LAG3_tpm" %in% names(g4_tcga))
.chk <- max(abs(g4_tcga[[lag3_col]] - log2(g4_tcga[["LAG3_tpm"]] + 1)), na.rm = TRUE)
cat("LAG3 outcome column:", lag3_col,
    "| max|LAG3 - log2(LAG3_tpm+1)| =", signif(.chk, 3), "(expect ~0)\n")

## ---- 2. Load ABSOLUTE table, detect columns --------------------------------
# Resolve the file even if saved under a slightly different name or as .gz.
if (!file.exists(ABSOLUTE_RAW)) {
  cand <- list.files(RAW_DIR, pattern = "abs_tables.*\\.(txt|tsv|gz)$",
                     full.names = TRUE, ignore.case = TRUE)
  if (!length(cand))
    cand <- list.files(RAW_DIR, pattern = "\\.(txt|tsv|gz)$",
                       full.names = TRUE, ignore.case = TRUE)
  if (length(cand)) ABSOLUTE_RAW <- cand[1]
}
if (!file.exists(ABSOLUTE_RAW)) {
  cat("\n[!] ABSOLUTE table not found in:\n   ", RAW_DIR, "\n")
  cat("    Download 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt' (open access) from\n")
  cat("    https://gdc.cancer.gov/about-data/publications/pancanatlas , save it there, re-run.\n")
  stop("ABSOLUTE table missing -- see message above.")
}
cat("Using ABSOLUTE file:", ABSOLUTE_RAW, "\n")
.con <- if (grepl("\\.gz$", ABSOLUTE_RAW, ignore.case = TRUE)) gzfile(ABSOLUTE_RAW) else ABSOLUTE_RAW
abs_tab <- read.delim(.con, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, check.names = FALSE)
sample_col <- grep("sample|array|barcode", colnames(abs_tab), ignore.case = TRUE, value = TRUE)[1]
purity_col <- grep("purity", colnames(abs_tab), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(sample_col), !is.na(purity_col))
cat("ABSOLUTE sample column:", sample_col, "| purity column:", purity_col, "\n")

abs_tab$.purity  <- suppressWarnings(as.numeric(abs_tab[[purity_col]]))
abs_tab$.patient <- substr(abs_tab[[sample_col]], 1, 12)            # TCGA-XX-XXXX
abs_tab$.stype   <- substr(abs_tab[[sample_col]], 14, 15)           # 01=primary, etc.

## ---- 3. Collapse ABSOLUTE to one purity per patient -------------------------
# Preference order within a patient: (1) primary solid tumor (01) > recurrent
# (02) > others; (2) ABSOLUTE call status "called" (high-confidence solution)
# over legacy/maf/snp calls; (3) a non-missing purity value. (Patient-level
# join; the RNA-seq sample is matched at the patient level.)
stype_rank <- function(s) ifelse(s == "01", 1L, ifelse(s == "02", 2L, 3L))
abs_tab$.rank <- stype_rank(abs_tab$.stype)
cs_col <- grep("call.?status|^status$", colnames(abs_tab), ignore.case = TRUE, value = TRUE)[1]
abs_tab$.csrank <- if (!is.na(cs_col))
  ifelse(tolower(trimws(abs_tab[[cs_col]])) == "called", 1L, 2L) else 1L
abs_ord <- abs_tab[order(abs_tab$.patient, abs_tab$.rank, abs_tab$.csrank,
                         is.na(abs_tab$.purity)), ]
abs_one <- abs_ord[!duplicated(abs_ord$.patient),
                   c(".patient", ".purity", ".stype")]
names(abs_one) <- c("patient", "purity_absolute", "abs_sampletype")
cat("ABSOLUTE: unique patients after collapse:", nrow(abs_one), "\n")

## ---- 4. Join to TCGA frame --------------------------------------------------
g4_tcga$patient <- substr(g4_tcga$case_barcode, 1, 12)
df <- merge(g4_tcga, abs_one, by = "patient", all.x = TRUE)
cat("matched (non-missing ABSOLUTE purity):",
    sum(!is.na(df$purity_absolute)), "of", nrow(df), "\n")

# Sanity: correlation between the two purity estimates (matched samples)
ok <- !is.na(df$purity_absolute) & !is.na(df$TumorPurity)
r_p <- cor(df$TumorPurity[ok], df$purity_absolute[ok], method = "pearson")
r_s <- cor(df$TumorPurity[ok], df$purity_absolute[ok], method = "spearman")
cat(sprintf("ESTIMATE vs ABSOLUTE purity: Pearson r = %.3f, Spearman rho = %.3f (n = %d)\n",
            r_p, r_s, sum(ok)))

## ---- 5. Build the common modeling frame ------------------------------------
# Same rows for all models -> comparable TP53 coefficients.
df$y    <- df[[lag3_col]]
df$TP53 <- factor(df$tp53_status)
# Put wild-type as reference so the TP53 term is the mutant effect (expect +).
wt_lvl  <- grep("wild|wt", levels(df$TP53), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(wt_lvl)) df$TP53 <- relevel(df$TP53, ref = wt_lvl)
df$IDH  <- factor(df$idh_status)

need <- c("y", "TP53", "IDH", "Tcell", "APM", "TumorPurity", "purity_absolute")
cc   <- complete.cases(df[, need])
m    <- df[cc, ]
cat("\nCommon modeling sample n:", nrow(m), "\n")
cat("TP53 levels (ref first):", paste(levels(m$TP53), collapse = " / "), "\n")
cat("IDH  levels (ref first):", paste(levels(m$IDH),  collapse = " / "), "\n")

# Drop IDH if constant within this subset (avoids rank-deficiency)
idh_ok <- nlevels(droplevels(m$IDH)) > 1
base_rhs <- if (idh_ok) "TP53 + IDH + Tcell + APM" else "TP53 + Tcell + APM"
if (!idh_ok) cat("NOTE: IDH constant in TCGA subset -> dropped from models.\n")

## ---- 6. Fit the three models -----------------------------------------------
f_base <- as.formula(paste("y ~", base_rhs))
f_est  <- as.formula(paste("y ~", base_rhs, "+ TumorPurity"))
f_abs  <- as.formula(paste("y ~", base_rhs, "+ purity_absolute"))

m_base <- lm(f_base, data = m)
m_est  <- lm(f_est,  data = m)
m_abs  <- lm(f_abs,  data = m)

# Helper: pull the TP53 (mutant) coefficient, SE, 95% CI, p
tp53_row <- function(fit) {
  cf  <- summary(fit)$coefficients
  trm <- grep("^TP53", rownames(cf), value = TRUE)[1]
  ci  <- confint(fit)[trm, ]
  data.frame(term = trm,
             beta = unname(cf[trm, "Estimate"]),
             se   = unname(cf[trm, "Std. Error"]),
             ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
             p    = unname(cf[trm, "Pr(>|t|)"]),
             stringsAsFactors = FALSE)
}

b_base <- tp53_row(m_base)
b_est  <- tp53_row(m_est)
b_abs  <- tp53_row(m_abs)

coef_tab <- rbind(
  cbind(model = "base (immune-adjusted)",     b_base),
  cbind(model = "base + ESTIMATE purity",     b_est),
  cbind(model = "base + ABSOLUTE purity",     b_abs)
)
coef_tab$delta_beta_pct <- (coef_tab$beta - b_base$beta) / b_base$beta * 100
rownames(coef_tab) <- NULL

cat("\n==== TP53 coefficient for LAG3 (TCGA grade-4, n =", nrow(m), ") ====\n")
print(format(coef_tab, digits = 4))

delta_tab <- data.frame(
  comparison        = c("ESTIMATE purity adj.", "ABSOLUTE purity adj."),
  beta_base         = b_base$beta,
  beta_adjusted     = c(b_est$beta, b_abs$beta),
  delta_beta_pct    = c((b_est$beta - b_base$beta)/b_base$beta*100,
                        (b_abs$beta - b_base$beta)/b_base$beta*100),
  p_adjusted        = c(b_est$p, b_abs$p),
  n                 = nrow(m),
  purity_corr_pearson = r_p
)
cat("\n==== Delta-beta (vs immune-adjusted base) ====\n")
print(format(delta_tab, digits = 4))

## ---- 7. Write outputs -------------------------------------------------------
write.csv(coef_tab,  file.path(OUT_DIR, "ABSOLUTE_TP53_coefficients.csv"), row.names = FALSE)
write.csv(delta_tab, file.path(OUT_DIR, "ABSOLUTE_delta_beta.csv"),        row.names = FALSE)

# Purity scatter (ESTIMATE vs ABSOLUTE), 450 dpi to match figure standard
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- ggplot(m, aes(TumorPurity, purity_absolute)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "ESTIMATE purity", y = "ABSOLUTE purity",
         title = sprintf("TCGA grade-4: ESTIMATE vs ABSOLUTE purity (r = %.2f, n = %d)",
                         r_p, nrow(m))) +
    theme_bw()
  ggsave(file.path(OUT_DIR, "purity_ESTIMATE_vs_ABSOLUTE_scatter.png"),
         p, width = 5, height = 5, dpi = 450)
}

cat("\nDONE. Outputs written to:\n  ", OUT_DIR, "\n")
cat("Report back coef_tab + delta_tab so we can finalize the R1-3b response text.\n")
