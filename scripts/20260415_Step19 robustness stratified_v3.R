# =============================================================================
# step19_robustness_stratified_v3.R
# GBM/Glioma TP53×LAG3 解析 - Step 19: 免疫コンテクスト層別 Robustness
# v3: BH補正の適用範囲を修正（v2からの変更点）
#
# 【v2からの変更点】
#   - BH補正の適用範囲を修正:
#       v2（誤）: All (ref)を含む5層×7遺伝子=35件にBH適用
#       v3（正）: All (ref)を除いた4層×7遺伝子=28件にBH適用
#     → 論文Methods/図キャプション記載の "BH correction (28 tests)" と一致
#   - All (ref)行のp_adj_BHはNAとして保存（28件補正の対象外）
#   - 出力ファイルを別名保存:
#       step19_stratified_results_v3.csv
#       figS6_robustness_v3.pdf / figS6_robustness_v3_600dpi.png
#
# 出力先: 19_robustness/
# 作成日: 2026-03-06 → v3修正: 2026-04-15
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "19_robustness")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

CP_GENES <- c("LAG3", "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "TIGIT", "HAVCR2")
CP_LABEL <- c(
  LAG3     = "LAG3",
  PDCD1    = "PDCD1 (PD-1)",
  CD274    = "CD274 (PD-L1)",
  PDCD1LG2 = "PDCD1LG2 (PD-L2)",
  CTLA4    = "CTLA4",
  TIGIT    = "TIGIT",
  HAVCR2   = "HAVCR2 (TIM-3)"
)

MU_TCELL <- 1.2569; SD_TCELL <- 0.6840
MU_APM   <- 6.8355; SD_APM   <- 0.7742

COL_MUT  <- "#E64B35"
COL_WT   <- "#AAAAAA"
COL_LAG3 <- "#E64B35"

GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step19_v3_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 19 v3: 免疫コンテクスト層別 Robustness 開始 ===")
log_msg("  [v3変更点] BH補正: 35件(v2) -> 28件(4層x7遺伝子, All(ref)除外)")

# =============================================================================
# 2. データ読み込み・スコア計算（Step17/18と同一ロジック）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH,  show_col_types = FALSE) %>% filter(include_flag == TRUE)
wide     <- read_csv(WIDE_PATH, show_col_types = FALSE)

get_log2_col <- function(gene, df) {
  c1 <- paste0(gene, "_log2tpm")
  c2 <- paste0(gsub("-", ".", gene), "_log2tpm")
  if (c1 %in% names(df)) return(c1)
  if (c2 %in% names(df)) return(c2)
  NA_character_
}

calc_score <- function(df, genes, label) {
  cols  <- sapply(genes, get_log2_col, df = df)
  found <- cols[!is.na(cols)]
  log_msg(sprintf("  %s: %d/%d遺伝子", label, length(found), length(genes)))
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

all_genes        <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols_found <- intersect(
  c(paste0(all_genes, "_log2tpm"), paste0(gsub("-",".",all_genes), "_log2tpm")),
  names(wide)
)
wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")), any_of(score_cols_found))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_score %>% filter(!is.na(case_barcode)),  by = "case_barcode",  suffix = c("",".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_score %>% filter(!is.na(wxs_sample_id)), by = "wxs_sample_id", suffix = c("",".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)
for (col in score_cols_found) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) { gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]]); gdc[[wcol]] <- NULL }
}

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM")
  )

gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT")
  )

log_msg(sprintf("GDC Grade4: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4), sum(gdc_g4$tp53_bin==1), sum(gdc_g4$tp53_bin==0)))

# =============================================================================
# 3. 層別変数の作成
# =============================================================================

log_msg("--- 層別変数作成 ---")
log_msg(sprintf("  T-cell カットオフ: Low<%.4f, High>%.4f",
                MU_TCELL - SD_TCELL, MU_TCELL + SD_TCELL))
log_msg(sprintf("  APM    カットオフ: Low<%.4f, High>%.4f",
                MU_APM - SD_APM, MU_APM + SD_APM))

gdc_g4 <- gdc_g4 %>%
  mutate(
    tcell_group = case_when(
      tcell_score < (MU_TCELL - SD_TCELL) ~ "Tcell_Low",
      tcell_score > (MU_TCELL + SD_TCELL) ~ "Tcell_High",
      TRUE ~ "Tcell_Mid"
    ),
    apm_group = case_when(
      apm_score < (MU_APM - SD_APM) ~ "APM_Low",
      apm_score > (MU_APM + SD_APM) ~ "APM_High",
      TRUE ~ "APM_Mid"
    )
  )

for (grp in c("Tcell_Low", "Tcell_Mid", "Tcell_High")) {
  sub <- gdc_g4 %>% filter(tcell_group == grp)
  log_msg(sprintf("  %s: n=%d (Mut=%d, WT=%d)",
                  grp, nrow(sub), sum(sub$tp53_bin==1), sum(sub$tp53_bin==0)))
}
for (grp in c("APM_Low", "APM_Mid", "APM_High")) {
  sub <- gdc_g4 %>% filter(apm_group == grp)
  log_msg(sprintf("  %s: n=%d (Mut=%d, WT=%d)",
                  grp, nrow(sub), sum(sub$tp53_bin==1), sum(sub$tp53_bin==0)))
}

# =============================================================================
# 4. 層別回帰（Low/High のみ使用）
# =============================================================================

log_msg("--- 層別回帰 ---")

fit_all   <- lm(LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin, data = gdc_g4)
ci_all    <- as.numeric(confint(fit_all, "tp53_bin"))
beta_all  <- round(as.numeric(coef(fit_all)["tp53_bin"]), 4)
ci_all_lo <- round(ci_all[1], 4)
ci_all_hi <- round(ci_all[2], 4)
log_msg(sprintf("  全体（参照）: beta=%+.4f [%+.4f, %+.4f]", beta_all, ci_all_lo, ci_all_hi))

make_row <- function(stratum_label, gene_name, n, n_mut, n_wt,
                     beta, se, ci_lo, ci_hi, pval, rsq) {
  data.frame(
    stratum     = as.character(stratum_label),
    gene        = as.character(gene_name),
    n           = as.integer(n),
    n_mut       = as.integer(n_mut),
    n_wt        = as.integer(n_wt),
    beta        = round(as.double(beta),  4),
    se          = round(as.double(se),    4),
    ci_lower    = round(as.double(ci_lo), 4),
    ci_upper    = round(as.double(ci_hi), 4),
    p_value     = as.double(pval),
    r_squared   = round(as.double(rsq),   4),
    significant = FALSE,
    stringsAsFactors = FALSE
  )
}

run_strat <- function(df, gene_col, stratum_label, n_min = 15) {
  df_use <- df %>% select(all_of(c(gene_col, "tp53_bin", "source_bin", "idh_bin"))) %>%
    na.omit()
  if (nrow(df_use) < n_min || length(unique(df_use$tp53_bin)) < 2) {
    log_msg(sprintf("    SKIP %s (%s): n=%d", stratum_label, gene_col, nrow(df_use)))
    return(NULL)
  }
  fit   <- lm(as.formula(paste(gene_col, "~ tp53_bin + source_bin + idh_bin")), data = df_use)
  coefs <- summary(fit)$coefficients
  ci_v  <- tryCatch(as.double(confint(fit, "tp53_bin")),
                    error = function(e) c(NA_real_, NA_real_))
  make_row(
    stratum_label,
    gsub("_log2tpm", "", gene_col),
    nrow(df_use),
    sum(df_use$tp53_bin == 1),
    sum(df_use$tp53_bin == 0),
    coefs["tp53_bin", "Estimate"],
    coefs["tp53_bin", "Std. Error"],
    ci_v[1], ci_v[2],
    coefs["tp53_bin", "Pr(>|t|)"],
    summary(fit)$r.squared
  )
}

log_msg("--- チェックポイント遺伝子列の確認 ---")
for (g in CP_GENES) {
  col_found <- get_log2_col(g, gdc_g4)
  log_msg(sprintf("  %s -> %s", g,
                  ifelse(is.na(col_found), "NOT FOUND", col_found)))
}

gene_cols <- sapply(CP_GENES, get_log2_col, df = gdc_g4)
gene_cols <- gene_cols[!is.na(gene_cols)]
log_msg(sprintf("gene_cols final (%d): %s",
                length(gene_cols), paste(gene_cols, collapse = ", ")))

if (length(gene_cols) < 2) {
  log_msg("  !! gene_colsが不足 -> wideからCP遺伝子列を直接gdc_g4に追加")
  cp_cols_needed <- unlist(lapply(CP_GENES, function(g) {
    c(paste0(g, "_log2tpm"), paste0(gsub("-", ".", g), "_log2tpm"))
  }))
  cp_cols_in_wide <- intersect(cp_cols_needed, names(wide))
  log_msg(sprintf("  wide内のCP列 (%d): %s",
                  length(cp_cols_in_wide), paste(cp_cols_in_wide, collapse = ", ")))
  
  wide_cp <- wide %>%
    select(any_of(c("case_barcode", "wxs_sample_id")), any_of(cp_cols_in_wide))
  
  gdc_g4_tcga  <- gdc_g4 %>% filter(source == "TCGA") %>%
    left_join(wide_cp %>% filter(!is.na(case_barcode)),
              by = "case_barcode", suffix = c("", ".cp"))
  gdc_g4_cptac <- gdc_g4 %>% filter(source == "CPTAC_HCMI") %>%
    left_join(wide_cp %>% filter(!is.na(wxs_sample_id)),
              by = "wxs_sample_id", suffix = c("", ".cp"))
  gdc_g4 <- bind_rows(gdc_g4_tcga, gdc_g4_cptac)
  for (col in cp_cols_in_wide) {
    wcol <- paste0(col, ".cp")
    if (wcol %in% names(gdc_g4)) {
      gdc_g4[[col]] <- coalesce(gdc_g4[[wcol]], gdc_g4[[col]])
      gdc_g4[[wcol]] <- NULL
    }
  }
  
  gene_cols <- sapply(CP_GENES, get_log2_col, df = gdc_g4)
  gene_cols <- gene_cols[!is.na(gene_cols)]
  log_msg(sprintf("  補完後 gene_cols (%d): %s",
                  length(gene_cols), paste(gene_cols, collapse = ", ")))
}

strata <- list(
  Tcell_Low  = gdc_g4 %>% filter(tcell_group == "Tcell_Low"),
  Tcell_High = gdc_g4 %>% filter(tcell_group == "Tcell_High"),
  APM_Low    = gdc_g4 %>% filter(apm_group   == "APM_Low"),
  APM_High   = gdc_g4 %>% filter(apm_group   == "APM_High")
)

strat_results <- do.call(rbind, Filter(Negate(is.null),
                                       lapply(names(strata), function(st) {
                                         rows <- Filter(Negate(is.null),
                                                        lapply(gene_cols, function(gc) run_strat(strata[[st]], gc, st)))
                                         if (length(rows) == 0) return(NULL)
                                         do.call(rbind, rows)
                                       })
))

ref_rows <- lapply(gene_cols, function(gc) {
  gv <- gsub("_log2tpm", "", gc)
  if (gv == "LAG3") {
    pval_all <- as.double(summary(fit_all)$coefficients["tp53_bin", "Pr(>|t|)"])
    make_row("All (ref)", gv, nrow(gdc_g4),
             sum(gdc_g4$tp53_bin==1), sum(gdc_g4$tp53_bin==0),
             beta_all, NA_real_, ci_all_lo, ci_all_hi,
             pval_all, summary(fit_all)$r.squared)
  } else {
    fit_g   <- lm(as.formula(paste(gc, "~ tp53_bin + source_bin + idh_bin")),
                  data = gdc_g4 %>%
                    select(all_of(c(gc, "tp53_bin", "source_bin", "idh_bin"))) %>%
                    na.omit())
    coefs_g <- summary(fit_g)$coefficients
    ci_g    <- as.double(confint(fit_g, "tp53_bin"))
    make_row("All (ref)", gv,
             nrow(fit_g$model),
             sum(fit_g$model$tp53_bin == 1),
             sum(fit_g$model$tp53_bin == 0),
             coefs_g["tp53_bin", "Estimate"],
             coefs_g["tp53_bin", "Std. Error"],
             ci_g[1], ci_g[2],
             coefs_g["tp53_bin", "Pr(>|t|)"],
             summary(fit_g)$r.squared)
  }
})
ref_df        <- do.call(rbind, Filter(Negate(is.null), ref_rows))
strat_results <- do.call(rbind, list(ref_df, strat_results))

strat_results$gene_label <- dplyr::recode(as.character(strat_results$gene),
                                          "LAG3"     = "LAG3",
                                          "PDCD1"    = "PDCD1 (PD-1)",
                                          "CD274"    = "CD274 (PD-L1)",
                                          "PDCD1LG2" = "PDCD1LG2 (PD-L2)",
                                          "CTLA4"    = "CTLA4",
                                          "TIGIT"    = "TIGIT",
                                          "HAVCR2"   = "HAVCR2 (TIM-3)"
)

# ===========================================================================
# 【v3修正箇所】BH補正: All (ref)を除いた28件のみに適用
# ---------------------------------------------------------------------------
# v2（誤）: p.adjust(strat_results$p_value, method="BH")  # 35件
# v3（正）: All (ref)行を除いた28件（4層×7遺伝子）にのみBH適用
#           All (ref)行のp_adj_BHはNAとして保存
# ===========================================================================
log_msg("--- BH補正 (v3修正: 28件 = 4層 x 7遺伝子, All(ref)除外) ---")

is_ref_row <- strat_results$stratum == "All (ref)"
n_tests    <- sum(!is_ref_row)
log_msg(sprintf("  補正対象: %d件 (All(ref)の%d行を除外)", n_tests, sum(is_ref_row)))

strat_results$p_adj_BH <- NA_real_
strat_results$p_adj_BH[!is_ref_row] <- p.adjust(
  strat_results$p_value[!is_ref_row],
  method = "BH"
)
strat_results$significant <- ifelse(
  is_ref_row,
  NA,  # All(ref)行は有意判定しない
  strat_results$p_adj_BH < 0.05
)

strat_results$is_lag3      <- strat_results$gene == "LAG3"
strat_results$stratum_type <- ifelse(
  grepl("Tcell", strat_results$stratum) | strat_results$stratum == "All (ref)",
  "T-cell", "APM"
)

write_csv(strat_results, file.path(OUT_DIR, "step19_stratified_results_v3.csv"))
log_msg("保存: step19_stratified_results_v3.csv")

log_msg("--- LAG3 層別結果 (p_adj_BH: 28件補正) ---")
strat_results %>% filter(gene == "LAG3") %>%
  arrange(stratum_type, stratum) %>%
  { for (i in seq_len(nrow(.))) {
    r <- .[i,]
    sig_mark <- if (is.na(r$significant)) "  " else if (r$significant) "★" else "  "
    padj_str <- if (is.na(r$p_adj_BH)) "p_BH=N/A (ref)" else sprintf("p_BH=%.4f", r$p_adj_BH)
    log_msg(sprintf("  %s %-15s n=%3d  beta=%+.4f [%+.4f, %+.4f] p_raw=%.4f  %s",
                    sig_mark, r$stratum, r$n, r$beta,
                    r$ci_lower, r$ci_upper, r$p_value, padj_str))
  }}

# =============================================================================
# 5. 図：figS6_robustness_v3（Forest風）
# =============================================================================

log_msg("--- figS6_robustness_v3 作成 ---")

stratum_order <- c("All (ref)", "Tcell_Low", "Tcell_High", "APM_Low", "APM_High")

# Panel A: LAG3のみ
lag3_data <- strat_results %>%
  filter(gene == "LAG3", stratum %in% stratum_order) %>%
  mutate(
    stratum_chr = as.character(stratum),
    stratum_lab = factor(
      dplyr::recode(stratum_chr,
                    "All (ref)"  = "All (ref, n=442)",
                    "Tcell_Low"  = "T-cell Low (<mean-1SD)",
                    "Tcell_High" = "T-cell High (>mean+1SD)",
                    "APM_Low"    = "APM Low (<mean-1SD)",
                    "APM_High"   = "APM High (>mean+1SD)"
      ),
      levels = rev(c("All (ref, n=442)","T-cell Low (<mean-1SD)",
                     "T-cell High (>mean+1SD)","APM Low (<mean-1SD)",
                     "APM High (>mean+1SD)"))
    ),
    is_ref  = (stratum_chr == "All (ref)"),
    n_label = sprintf("n=%d", n),
    # All(ref)はp生値のみ表示、層別はBH補正後のp値を表示
    p_label = dplyr::case_when(
      stratum_chr == "All (ref)" ~
        sprintf("p=%.4f (ref)", p_value),
      p_adj_BH < 0.001 ~
        "p_BH<0.001",
      p_adj_BH < 0.05 ~
        sprintf("p_BH=%.3f", p_adj_BH),
      TRUE ~
        sprintf("p_BH=%.3f (ns)", p_adj_BH)
    )
  )

pa <- ggplot(lag3_data, aes(x = beta, y = stratum_lab,
                            color = is_ref, fill = is_ref)) +
  geom_vline(xintercept = 0,        linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_vline(xintercept = beta_all, linetype = "dotted", color = COL_MUT,   linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.2, linewidth = 0.8) +
  geom_point(shape = 21, size = 4, stroke = 1.0) +
  geom_text(aes(label = p_label, x = ci_upper),
            hjust = -0.15, size = 2.8, color = "grey30") +
  scale_color_manual(values = c("TRUE" = COL_MUT, "FALSE" = COL_LAG3), guide = "none") +
  scale_fill_manual(values  = c("TRUE" = COL_MUT, "FALSE" = COL_LAG3), guide = "none") +
  scale_x_continuous(
    name   = expression(italic(TP53) ~ "coefficient (\u03b2) with 95% CI"),
    breaks = seq(0, 0.8, by = 0.2)
  ) +
  coord_cartesian(xlim = c(
    min(lag3_data$ci_lower, na.rm = TRUE) - 0.05,
    max(lag3_data$ci_upper, na.rm = TRUE) * 1.35
  )) +
  labs(y = NULL, title = expression(bold("A") ~ ~ italic(LAG3) ~ "across immune strata")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 8.5),
        plot.title         = element_text(size = 10))

# Panel B: 全7遺伝子 × 4層（All(ref)除外）
heat_data <- strat_results %>%
  filter(stratum != "All (ref)") %>%
  mutate(
    stratum_chr = as.character(stratum),
    stratum_lab = factor(
      dplyr::recode(stratum_chr,
                    "Tcell_Low"  = "T-cell Low (<mean-1SD)",
                    "Tcell_High" = "T-cell High (>mean+1SD)",
                    "APM_Low"    = "APM Low (<mean-1SD)",
                    "APM_High"   = "APM High (>mean+1SD)"
      ),
      levels = c("T-cell Low (<mean-1SD)","T-cell High (>mean+1SD)",
                 "APM Low (<mean-1SD)","APM High (>mean+1SD)")
    ),
    gene_chr = as.character(gene),
    gene_lab = factor(
      dplyr::recode(gene_chr,
                    "LAG3"     = "LAG3",
                    "PDCD1"    = "PDCD1 (PD-1)",
                    "CD274"    = "CD274 (PD-L1)",
                    "PDCD1LG2" = "PDCD1LG2 (PD-L2)",
                    "CTLA4"    = "CTLA4",
                    "TIGIT"    = "TIGIT",
                    "HAVCR2"   = "HAVCR2 (TIM-3)"
      ),
      levels = c("HAVCR2 (TIM-3)","TIGIT","CTLA4",
                 "PDCD1LG2 (PD-L2)","CD274 (PD-L1)","PDCD1 (PD-1)","LAG3")
    ),
    sig_cat = case_when(
      significant & is_lag3  ~ "lag3_sig",
      significant & !is_lag3 ~ "other_sig",
      TRUE                   ~ "ns"
    ),
    beta_text = ifelse(
      significant,
      sprintf("%+.2f", beta),
      sprintf("(%+.2f)", beta)
    )
  )

pb <- ggplot(heat_data, aes(x = stratum_lab, y = gene_lab)) +
  geom_tile(aes(fill = sig_cat),
            color = "white", linewidth = 0.4) +
  geom_text(aes(label = beta_text, color = sig_cat),
            size = 3.0) +
  scale_fill_manual(
    values = c(
      "lag3_sig"  = "#FFCCBB",
      "other_sig" = "#CCCCCC",
      "ns"        = "#F5F5F5"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "lag3_sig"  = COL_LAG3,
      "other_sig" = "#333333",
      "ns"        = "#999999"
    ),
    guide = "none"
  ) +
  scale_x_discrete(labels = function(x) gsub(" \\(.*\\)", "", x)) +
  labs(
    x     = NULL,
    y     = NULL,
    title = expression(
      bold("B") ~ "All 7 checkpoint genes   " *
        "(significant: colored; ns: grey in parentheses)"
    )
  ) +
  scale_y_discrete(labels = function(x) {
    sapply(x, function(lab) {
      gene_part  <- sub(" \\(.*\\)", "", lab)
      alias_part <- regmatches(lab, regexpr("\\(.*\\)", lab))
      if (length(alias_part) > 0) {
        bquote(italic(.(gene_part)) ~ .(alias_part))
      } else {
        bquote(italic(.(gene_part)))
      }
    })
  }) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x  = element_text(size = 8.5, angle = 20, hjust = 1),
    axis.text.y  = element_text(size = 8.5),
    panel.grid   = element_blank(),
    plot.title   = element_text(size = 9.5)
  )

fig_s6 <- pa / pb +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    caption = paste0(
      "GDC Grade4 (n=442). Strata: T-cell/APM score < mean-1SD (Low) or > mean+1SD (High).\n",
      "Panel A: LAG3 TP53 \u03b2 coefficient in each stratum. Red dotted line: overall \u03b2.",
      " All (ref) shows raw p value; strata show BH-adjusted p.\n",
      "Panel B: All 7 checkpoint genes; \u03b2 value shown in each cell; ",
      "colored = BH-adjusted p<0.05 (28 tests = 4 strata x 7 genes); grey in parentheses = ns.\n",
      "Model: gene ~ tp53_bin + source_bin + idh_bin within each stratum. ",
      "T-cell/APM score cutoffs from Step18 (mean+/-1SD of GDC Grade4 distribution)."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

# =============================================================================
# 6. 出力（別名・600dpi PNG + PDF）
# =============================================================================

log_msg("--- ファイル出力 ---")

pdf_path <- file.path(OUT_DIR, "figS6_robustness_v3.pdf")
pdf(pdf_path, width = 9, height = 9)
print(fig_s6)
dev.off()
log_msg(sprintf("保存: %s", pdf_path))

png_path <- file.path(OUT_DIR, "figS6_robustness_v3_600dpi.png")
agg_png(png_path, width = 9, height = 9, units = "in", res = 600)
print(fig_s6)
dev.off()
log_msg(sprintf("保存: %s", png_path))

# =============================================================================
# 7. 完了
# =============================================================================

log_msg("=== Step 19 v3: 完了 ===")

cat("\n============================\n")
cat("Step 19 v3 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s\n", file.path(OUT_DIR, "step19_stratified_results_v3.csv")))
cat(sprintf("  %s\n", pdf_path))
cat(sprintf("  %s\n", png_path))
cat("============================\n")
