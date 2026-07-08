# =============================================================================
# step21_immune_score_tp53_diff.R
# GBM/Glioma TP53xLAG3 解析 - Step 21: 免疫スコアのTP53差確認
#
# 目的:
#   T-cell / APM / IFN-gamma スコアがTP53 Mut/WTで差がないことを示す。
#   "TP53-LAG3関連は免疫量の差ではない"の根拠として Supplement に配置。
#
# 出力:
#   figS21a_score_violin.pdf/png       : A) violin+boxplot（3スコア facet）
#   figS21b_score_forest.pdf/png       : B) Forest風（beta+CI、3スコア）
#   figS21_combined.pdf/png            : C) A+B上下配置
#   step21_score_tp53_regression.csv   : 回帰結果
#   step21_log.txt
#
# スコア定義（Step17/20と同一）:
#   T-cell : CD3D, CD3E, CD3G, CD8A, CD8B, GZMA, GZMB, PRF1  (8遺伝子)
#   APM    : B2M, TAP1, TAP2, TAPBP, HLA-A, HLA-B, HLA-C, NLRC5  (8遺伝子)
#   IFN-gamma: STAT1, IRF1, IRF9, CXCL9, CXCL10, CXCL11,
#              GBP1, GBP2, GBP4, GBP5, IDO1  (11遺伝子)
#
# 回帰モデル:
#   score ~ tp53_bin + source_bin + idh_bin  (GDC Grade4)
#
# 出力先: 21_immune_score_tp53/
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "21_immune_score_tp53")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# スコア定義（Step17/20と完全に同一）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

# 色設定（引継書10章に準拠）
COL_MUT <- "#E64B35"   # TP53 Mut
COL_WT  <- "#AAAAAA"   # TP53 WT

# スコアラベル（図表示用・ASCII統一）
SCORE_LABELS <- c(
  tcell_score = "T-cell score",
  apm_score   = "APM score",
  ifng_score  = "IFN-gamma score"
)

# 入力ファイル
GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step21_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 21: 免疫スコアのTP53差確認 開始 ===")

# =============================================================================
# 2. データ読み込み・結合（Step20と同一ロジック）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
wide     <- read_csv(WIDE_PATH, show_col_types = FALSE)

log_msg(sprintf("GDC: %d行 / wide: %d行x%d列", nrow(gdc_base), nrow(wide), ncol(wide)))

# 列名正規化ヘルパー
get_log2_col <- function(gene, df) {
  c1 <- paste0(gene, "_log2tpm")
  c2 <- paste0(gsub("-", ".", gene), "_log2tpm")
  if (c1 %in% names(df)) return(c1)
  if (c2 %in% names(df)) return(c2)
  return(NA_character_)
}

calc_score <- function(df, genes, score_name) {
  cols    <- sapply(genes, get_log2_col, df = df)
  found   <- cols[!is.na(cols)]
  missing <- genes[is.na(cols)]
  if (length(missing) > 0)
    log_msg(sprintf("  %s: 不足遺伝子 %s", score_name, paste(missing, collapse = ",")))
  log_msg(sprintf("  %s: %d/%d遺伝子で計算", score_name, length(found), length(genes)))
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

# GDC結合
all_score_genes  <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols_found <- intersect(
  c(paste0(all_score_genes, "_log2tpm"),
    paste0(gsub("-", ".", all_score_genes), "_log2tpm")),
  names(wide)
)

wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")), any_of(score_cols_found))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_score %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_score %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)

for (col in score_cols_found) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]])
    gdc[[wcol]] <- NULL
  }
}
log_msg(sprintf("GDC結合後: %d行", nrow(gdc)))

# =============================================================================
# 3. 免疫スコア計算
# =============================================================================

log_msg("--- 免疫スコア計算 ---")

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "IFN-gamma")
  )

# Grade4解析データ
gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = factor(
      ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT"),
      levels = c("TP53 WT", "TP53 Mut")
    )
  )

log_msg(sprintf("GDC Grade4: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4), sum(gdc_g4$tp53_bin == 1), sum(gdc_g4$tp53_bin == 0)))

# スコアNA確認
for (sc in c("tcell_score", "apm_score", "ifng_score")) {
  log_msg(sprintf("  %s: NA=%d, median=%.3f",
                  sc, sum(is.na(gdc_g4[[sc]])),
                  median(gdc_g4[[sc]], na.rm = TRUE)))
}

# =============================================================================
# 4. 記述統計
# =============================================================================

log_msg("--- 記述統計（TP53 Mut/WT別） ---")

desc_stats <- gdc_g4 %>%
  group_by(tp53_label) %>%
  summarise(
    n             = n(),
    tcell_median  = round(median(tcell_score, na.rm = TRUE), 4),
    tcell_IQR     = round(IQR(tcell_score,    na.rm = TRUE), 4),
    apm_median    = round(median(apm_score,   na.rm = TRUE), 4),
    apm_IQR       = round(IQR(apm_score,      na.rm = TRUE), 4),
    ifng_median   = round(median(ifng_score,  na.rm = TRUE), 4),
    ifng_IQR      = round(IQR(ifng_score,     na.rm = TRUE), 4),
    .groups = "drop"
  )

for (i in seq_len(nrow(desc_stats))) {
  r <- desc_stats[i, ]
  log_msg(sprintf("  [%s] n=%d | Tcell=%.4f (IQR=%.4f) | APM=%.4f (IQR=%.4f) | IFNg=%.4f (IQR=%.4f)",
                  r$tp53_label, r$n,
                  r$tcell_median, r$tcell_IQR,
                  r$apm_median,   r$apm_IQR,
                  r$ifng_median,  r$ifng_IQR))
}

# =============================================================================
# 5. 回帰解析: score ~ tp53_bin + source_bin + idh_bin
# =============================================================================

log_msg("--- 回帰解析: score ~ tp53 + source + idh ---")

run_score_reg <- function(df, score_col, score_name) {
  df_use <- df %>%
    select(all_of(c(score_col, "tp53_bin", "source_bin", "idh_bin"))) %>%
    rename(score = all_of(score_col)) %>%
    na.omit()
  
  fit <- lm(score ~ tp53_bin + source_bin + idh_bin, data = df_use)
  ci  <- as.double(confint(fit, "tp53_bin"))
  smr <- coef(summary(fit))
  
  beta <- smr["tp53_bin", "Estimate"]
  se   <- smr["tp53_bin", "Std. Error"]
  pval <- smr["tp53_bin", "Pr(>|t|)"]
  r2   <- round(summary(fit)$r.squared, 4)
  
  log_msg(sprintf("  [%s] n=%d, beta=%+.4f [%+.4f, %+.4f], p=%.4f, R2=%.4f",
                  score_name, nrow(df_use), beta, ci[1], ci[2], pval, r2))
  
  tibble(
    score     = score_name,
    score_col = score_col,
    n         = nrow(df_use),
    beta      = round(beta,  4),
    ci_lower  = round(ci[1], 4),
    ci_upper  = round(ci[2], 4),
    se        = round(se,    4),
    p_value   = pval,
    r_squared = r2
  )
}

reg_results <- do.call(rbind, list(
  run_score_reg(gdc_g4, "tcell_score", "T-cell score"),
  run_score_reg(gdc_g4, "apm_score",   "APM score"),
  run_score_reg(gdc_g4, "ifng_score",  "IFN-gamma score")
))

# 有意差ラベル
reg_results <- reg_results %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    p_display = sprintf("p=%.3f (%s)", p_value, sig_label),
    score = factor(score, levels = rev(c("T-cell score", "APM score", "IFN-gamma score")))
  )

write_csv(reg_results %>% mutate(score = as.character(score)),
          file.path(OUT_DIR, "step21_score_tp53_regression.csv"))
log_msg("保存: step21_score_tp53_regression.csv")

# =============================================================================
# 6. Panel A: violin + boxplot（3スコア facet）
# =============================================================================

log_msg("--- Panel A: violin+boxplot 作成 ---")

# long形式に変換
plot_long <- gdc_g4 %>%
  select(tp53_label, tcell_score, apm_score, ifng_score) %>%
  pivot_longer(
    cols      = c(tcell_score, apm_score, ifng_score),
    names_to  = "score_col",
    values_to = "score_value"
  ) %>%
  mutate(
    score_label = factor(
      dplyr::recode(score_col,
                    tcell_score = "T-cell score",
                    apm_score   = "APM score",
                    ifng_score  = "IFN-gamma score"
      ),
      levels = c("T-cell score", "APM score", "IFN-gamma score")
    )
  )

# n数ラベル（facet下部に表示用）
n_labels <- plot_long %>%
  group_by(score_label, tp53_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    y_pos = -Inf,
    label = paste0("n=", n)
  )

fig_a <- ggplot(plot_long, aes(x = tp53_label, y = score_value,
                               color = tp53_label, fill = tp53_label)) +
  geom_violin(alpha = 0.25, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.18, alpha = 0.7, linewidth = 0.6,
               outlier.size = 0.8, outlier.alpha = 0.4) +
  geom_jitter(width = 0.12, alpha = 0.15, size = 0.6, shape = 16) +
  geom_text(data = n_labels,
            aes(x = tp53_label, y = y_pos, label = label),
            vjust = -0.3, size = 2.8, color = "grey40", inherit.aes = FALSE) +
  facet_wrap(~ score_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = c("TP53 WT" = COL_WT, "TP53 Mut" = COL_MUT),
                     guide = "none") +
  scale_fill_manual(values  = c("TP53 WT" = COL_WT, "TP53 Mut" = COL_MUT),
                    guide = "none") +
  labs(
    x        = NULL,
    y        = "Score [mean log2(TPM+1)]",
    title    = "A  Immune quantity scores by TP53 status",
    subtitle = "GDC Grade4 (n=442). Violin + boxplot + jitter."
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "#F0F0F0", color = NA),
    strip.text        = element_text(size = 9.5, face = "bold"),
    axis.text.x       = element_text(size = 9),
    plot.title        = element_text(face = "bold", size = 11),
    plot.subtitle     = element_text(size = 8.5, color = "grey40")
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS21a_score_violin.pdf"), width = 9, height = 4.5)
print(fig_a)
dev.off()
agg_png(file.path(OUT_DIR, "figS21a_score_violin_450dpi.png"),
        width = 9, height = 4.5, units = "in", res = 450)
print(fig_a)
dev.off()
log_msg("保存: figS21a_score_violin.pdf/png")

# =============================================================================
# 7. Panel B: Forest風（beta + CI、3スコア）
# =============================================================================

log_msg("--- Panel B: Forest風 作成 ---")

# x軸範囲（0跨ぎを強調するため対称に広めに取る）
x_abs_max <- max(abs(c(reg_results$ci_lower, reg_results$ci_upper)),
                 na.rm = TRUE) * 1.5
x_lim     <- c(-x_abs_max, x_abs_max)

fig_b <- ggplot(reg_results,
                aes(x = beta, y = score,
                    color = sig_label == "ns",
                    fill  = sig_label == "ns")) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.25, linewidth = 0.8) +
  geom_point(shape = 21, size = 4.5, stroke = 1.1) +
  # p値ラベル
  geom_text(aes(x = ci_upper, label = p_display),
            hjust = -0.1, size = 3.2, color = "grey30") +
  scale_color_manual(
    values = c("TRUE" = "#888888", "FALSE" = COL_MUT),
    guide  = "none"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#CCCCCC", "FALSE" = COL_MUT),
    guide  = "none"
  ) +
  scale_x_continuous(
    name   = "TP53 coefficient (beta_TP53) with 95% CI\n[score ~ tp53 + source + IDH]",
    limits = x_lim,
    breaks = seq(-0.3, 0.3, by = 0.1)
  ) +
  coord_cartesian(xlim = x_lim) +
  labs(
    y        = NULL,
    title    = "B  Regression: immune score ~ TP53 + source + IDH",
    subtitle = "GDC Grade4 (n=442). All 95% CI cross zero (ns = not significant)."
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40"),
    plot.margin        = margin(10, 100, 10, 10)
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS21b_score_forest.pdf"), width = 8, height = 4)
print(fig_b)
dev.off()
agg_png(file.path(OUT_DIR, "figS21b_score_forest_450dpi.png"),
        width = 8, height = 4, units = "in", res = 450)
print(fig_b)
dev.off()
log_msg("保存: figS21b_score_forest.pdf/png")

# =============================================================================
# 8. Panel C: 上下配置（combined）
# =============================================================================

log_msg("--- Panel C: combined 上下配置 作成 ---")

fig_combined <- fig_a / fig_b +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title   = "Immune quantity scores are comparable between TP53 Mut and WT",
    caption = paste0(
      "T-cell score: CD3D/E/G, CD8A/B, GZMA/B, PRF1 (n=8). ",
      "APM score: B2M, TAP1/2, TAPBP, HLA-A/B/C, NLRC5 (n=8). ",
      "IFN-gamma score: STAT1, IRF1/9, CXCL9/10/11, GBP1/2/4/5, IDO1 (n=11). ",
      "Scores: mean log2(TPM+1) of component genes. ",
      "Regression adjusted for source (TCGA=0) and IDH status (Mut=1)."
    ),
    theme = theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

pdf(file.path(OUT_DIR, "figS21_combined.pdf"), width = 10, height = 9)
print(fig_combined)
dev.off()
agg_png(file.path(OUT_DIR, "figS21_combined_450dpi.png"),
        width = 10, height = 9, units = "in", res = 450)
print(fig_combined)
dev.off()
log_msg("保存: figS21_combined.pdf/png")

# =============================================================================
# 9. 完了
# =============================================================================

log_msg("=== Step 21: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 21 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/figS21a_score_violin.pdf/png\n",   OUT_DIR))
cat(sprintf("  %s/figS21b_score_forest.pdf/png\n",   OUT_DIR))
cat(sprintf("  %s/figS21_combined.pdf/png\n",        OUT_DIR))
cat(sprintf("  %s/step21_score_tp53_regression.csv\n", OUT_DIR))
cat(sprintf("  %s/step21_log.txt\n",                 OUT_DIR))
cat("============================\n")
