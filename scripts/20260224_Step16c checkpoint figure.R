# =============================================================================
# step16c_checkpoint_figure.R
# GBM/Glioma TP53×LAG3 解析 - Step 16c: チェックポイント特異性 可視化
# （修正版3：遺伝子順序バグ修正 - LAG3最上段固定）
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

RESULT_DIR <- here::here("results", "TP53", "20260221")
OUT_DIR    <- file.path(RESULT_DIR, "16_checkpoint_specificity")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 遺伝子表示順（上から）: LAG3 > PDCD1 > CD274 > PDCD1LG2 > CTLA4 > TIGIT > HAVCR2
# ggplotのy軸は下から上に並ぶため、rev()で反転することを前提に「逆順」で定義する
# rev()後: LAG3が最大factor値 = 最上段、HAVCR2が最小 = 最下段
GENE_ORDER_DISPLAY <- c("LAG3", "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "TIGIT", "HAVCR2")
GENE_ORDER_FACTOR  <- rev(GENE_ORDER_DISPLAY)   # factor levels: HAVCR2...LAG3（下→上）

# y軸表示ラベル（遺伝子記号 → 読者向け表示名）
GENE_LABEL <- c(
  LAG3     = "LAG3",
  PDCD1    = "PDCD1 (PD-1)",
  CD274    = "CD274 (PD-L1)",
  PDCD1LG2 = "PDCD1LG2 (PD-L2)",
  CTLA4    = "CTLA4",
  TIGIT    = "TIGIT",
  HAVCR2   = "HAVCR2 (TIM-3)"
)

COHORT_ORDER <- c("GDC_Grade4_all", "GDC_Grade4_TCGA",
                  "GDC_Grade4_CPTAC", "GLASS_WXS")
COHORT_LABEL <- c(
  GDC_Grade4_all   = "GDC Grade4\n(All, n=442)",
  GDC_Grade4_TCGA  = "GDC Grade4\n(TCGA, n=245)",
  GDC_Grade4_CPTAC = "GDC Grade4\n(CPTAC/HCMI, n=197)",
  GLASS_WXS        = "GLASS\n(WXS, n=79)"
)

COL_LAG3  <- "#E64B35"
COL_OTHER <- "#AAAAAA"

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step16c_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 16c: チェックポイント特異性 可視化 開始（修正版3） ===")
log_msg(sprintf("表示順（上から）: %s", paste(GENE_ORDER_DISPLAY, collapse = " > ")))
log_msg(sprintf("factor順（下から）: %s", paste(GENE_ORDER_FACTOR, collapse = " > ")))

# =============================================================================
# 2. データ読み込み
# =============================================================================

log_msg("--- データ読み込み ---")

reg_path <- file.path(OUT_DIR, "step16b_regression_results.csv")
wil_path <- file.path(OUT_DIR, "step16a_v2_results.csv")

if (!file.exists(reg_path)) {
  log_msg("ERROR: step16b_regression_results.csv が見つかりません")
  close(log_con); stop("Step16b未実行")
}
if (!file.exists(wil_path)) {
  log_msg("ERROR: step16a_v2_results.csv が見つかりません")
  close(log_con); stop("Step16a未実行")
}

reg <- read_csv(reg_path, show_col_types = FALSE)
wil <- read_csv(wil_path, show_col_types = FALSE)
log_msg(sprintf("回帰結果: %d行 / Wilcoxon結果: %d行", nrow(reg), nrow(wil)))

# =============================================================================
# 3. プロット用データ整形
# =============================================================================

log_msg("--- データ整形 ---")

plot_data <- reg %>%
  filter(cohort %in% COHORT_ORDER, note == "ok",
         gene %in% GENE_ORDER_DISPLAY) %>%
  mutate(
    # GENE_ORDER_FACTOR順でfactor化（HAVCR2=1, ..., LAG3=7）
    gene        = factor(gene, levels = GENE_ORDER_FACTOR),
    cohort      = factor(cohort, levels = COHORT_ORDER),
    cohort_lab  = factor(COHORT_LABEL[as.character(cohort)],
                         levels = COHORT_LABEL[COHORT_ORDER]),
    is_lag3     = (as.character(gene) == "LAG3"),
    point_color = ifelse(is_lag3, COL_LAG3, COL_OTHER),
    point_fill  = case_when(
      significant & is_lag3  ~ COL_LAG3,
      significant & !is_lag3 ~ COL_OTHER,
      TRUE                   ~ "white"
    )
  )

# Wilcoxon結果をmerge
wil_flag <- wil %>%
  filter(cohort %in% COHORT_ORDER, gene %in% GENE_ORDER_DISPLAY, note == "ok") %>%
  mutate(
    cohort = factor(cohort, levels = COHORT_ORDER),
    gene   = factor(gene,   levels = GENE_ORDER_FACTOR)
  ) %>%
  select(cohort, gene, cliffs_delta, p_BH_wilcox = p_BH, sig_wilcox = significant)

plot_data <- plot_data %>%
  left_join(wil_flag, by = c("cohort", "gene")) %>%
  mutate(gene = factor(as.character(gene), levels = GENE_ORDER_FACTOR))

log_msg(sprintf("プロットデータ: %d行", nrow(plot_data)))

# factor levelの確認（デバッグ用）
log_msg("gene factor levels（下から上）:")
log_msg(sprintf("  %s", paste(levels(plot_data$gene), collapse = " < ")))
log_msg("期待: HAVCR2 < TIGIT < CTLA4 < PDCD1LG2 < CD274 < PDCD1 < LAG3")

x_min <- floor(min(plot_data$ci_lower, na.rm = TRUE) * 10) / 10 - 0.05
x_max <- ceiling(max(plot_data$ci_upper, na.rm = TRUE) * 10) / 10 + 0.05
x_lim <- c(max(x_min, -1.2), min(x_max, 1.2))
log_msg(sprintf("x軸範囲: [%.2f, %.2f]", x_lim[1], x_lim[2]))

log_msg("LAG3 beta確認:")
plot_data %>% filter(as.character(gene) == "LAG3") %>%
  select(cohort, beta, ci_lower, ci_upper, p_BH, significant) %>%
  { for (i in seq_len(nrow(.))) {
    r <- .[i,]
    sig <- if (!is.na(r$significant) && r$significant) "★" else "  "
    log_msg(sprintf("  %s %-25s beta=%+.3f [%+.3f, %+.3f] p_BH=%.4f",
                    sig, as.character(r$cohort),
                    r$beta, r$ci_lower, r$ci_upper,
                    ifelse(is.na(r$p_BH), NA, r$p_BH)))
  }}

# =============================================================================
# 4. y軸テキストスタイル
# GENE_ORDER_FACTORの順（下→上）でfaceとcolorを設定
# LAG3はGENE_ORDER_FACTORの最後の要素（= 最上段）
# =============================================================================

y_face  <- ifelse(GENE_ORDER_FACTOR == "LAG3", "bold", "plain")
y_color <- ifelse(GENE_ORDER_FACTOR == "LAG3", COL_LAG3, "black")

# scale_y_discrete用ラベル（GENE_ORDER_FACTOR順）
y_labels_ordered <- GENE_LABEL[GENE_ORDER_FACTOR]

# LAG3のy位置（factor内での数値 = length = 7 = 最大値 = 最上段）
lag3_y <- length(GENE_ORDER_FACTOR)
log_msg(sprintf("LAG3のy位置: %d（factor最大値 = 最上段）", lag3_y))

# =============================================================================
# 5. プロット作成
# =============================================================================

log_msg("--- プロット作成 ---")

p <- ggplot(plot_data,
            aes(x = beta, y = gene,
                color = point_color, fill = point_fill)) +
  
  # LAG3行背景ハイライト（y = lag3_y = 7）
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = lag3_y - 0.5, ymax = lag3_y + 0.5,
           fill = "#FFF0EE", alpha = 0.6) +
  
  facet_wrap(~ cohort_lab, nrow = 1, scales = "free_x") +
  
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#444444", linewidth = 0.4) +
  
  # CIバー（ggplot2 4.0対応）
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y",
                width = 0.25, linewidth = 0.5) +
  
  # 点推定
  geom_point(shape = 21, size = 3.2, stroke = 0.9) +
  
  scale_color_identity() +
  scale_fill_identity() +
  
  scale_x_continuous(
    name   = expression(paste("Regression coefficient (", beta[TP53], ")")),
    breaks = seq(-1.0, 1.0, by = 0.5),
    labels = function(x) sprintf("%+.1f", x),
    limits = x_lim
  ) +
  
  # y軸ラベルを読者向け表示名に変換
  scale_y_discrete(
    name   = NULL,
    labels = y_labels_ordered
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    strip.background   = element_rect(fill = "#F5F5F5", color = "grey60"),
    strip.text         = element_text(size = 8.5, face = "bold"),
    panel.grid.major.y = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.text.y  = element_text(size = 9, face = y_face, color = y_color),
    axis.text.x  = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    legend.position = "none",
    plot.title    = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 8,  color = "grey40"),
    plot.caption  = element_text(size = 7,  color = "grey50", hjust = 0)
  ) +
  
  labs(
    title    = "Checkpoint gene specificity: TP53 mutation effect on immune checkpoint expression",
    subtitle = "Adjusted regression coefficient with 95% CI | Filled circle: BH-adjusted p < 0.05 | Red: LAG3",
    caption  = paste0(
      "Models: GDC Grade4 (All/IDH-WT): gene ~ tp53_bin + source_bin + idh_bin; ",
      "GDC (TCGA/CPTAC): gene ~ tp53_bin + idh_bin; GLASS: gene ~ tp53_bin\n",
      "BH correction: 7 genes per cohort. ",
      "tp53_bin: Mut=1/WT=0; source_bin: CPTAC_HCMI=1/TCGA=0; idh_bin: IDH Mut=1/IDH WT=0"
    )
  )

# =============================================================================
# 6. PDF出力
# =============================================================================

pdf_path <- file.path(OUT_DIR, "figS_checkpoint_specificity.pdf")
pdf(pdf_path, width = 10, height = 5.5)
print(p)
dev.off()
log_msg("保存: figS_checkpoint_specificity.pdf")

# =============================================================================
# 7. PNG出力（450dpi）
# =============================================================================

png_path <- file.path(OUT_DIR, "figS_checkpoint_specificity_450dpi.png")
agg_png(png_path, width = 10, height = 5.5, units = "in", res = 450)
print(p)
dev.off()
log_msg("保存: figS_checkpoint_specificity_450dpi.png")

# =============================================================================
# 8. Supplement Table
# =============================================================================

log_msg("--- Supplement Table作成 ---")

supp_wil <- wil %>%
  filter(cohort %in% COHORT_ORDER, gene %in% GENE_ORDER_DISPLAY, note == "ok") %>%
  select(cohort, gene, n_mut, n_wt,
         median_mut, median_wt, median_diff,
         cliffs_delta, p_wilcox, p_BH_wilcox = p_BH, sig_wilcox = significant)

supp_reg <- reg %>%
  filter(cohort %in% COHORT_ORDER, gene %in% GENE_ORDER_DISPLAY, note == "ok") %>%
  select(cohort, gene, beta, ci_lower, ci_upper,
         p_value_reg = p_value, p_BH_reg = p_BH, sig_reg = significant, r_squared)

supp_table <- supp_wil %>%
  left_join(supp_reg, by = c("cohort", "gene")) %>%
  mutate(
    gene_label   = GENE_LABEL[gene],
    gene         = factor(gene,   levels = GENE_ORDER_DISPLAY),
    cohort       = factor(cohort, levels = COHORT_ORDER),
    median_diff  = round(median_diff,  3),
    cliffs_delta = round(cliffs_delta, 3),
    p_wilcox     = signif(p_wilcox,    3),
    p_BH_wilcox  = signif(p_BH_wilcox, 3),
    beta         = round(beta,    4),
    ci_lower     = round(ci_lower, 4),
    ci_upper     = round(ci_upper, 4),
    p_value_reg  = signif(p_value_reg, 3),
    p_BH_reg     = signif(p_BH_reg,   3),
    r_squared    = round(r_squared,    4),
    ci_95        = sprintf("[%+.3f, %+.3f]", ci_lower, ci_upper)
  ) %>%
  arrange(cohort, gene) %>%
  select(cohort, gene, gene_label, n_mut, n_wt,
         median_mut, median_wt, median_diff,
         cliffs_delta, p_wilcox, p_BH_wilcox, sig_wilcox,
         beta, ci_95, p_value_reg, p_BH_reg, sig_reg, r_squared)

write_csv(supp_table, file.path(OUT_DIR, "step16c_supplement_table.csv"))
log_msg(sprintf("保存: step16c_supplement_table.csv (%d行)", nrow(supp_table)))

# =============================================================================
# 9. 完了サマリー
# =============================================================================

log_msg("--- 有意遺伝子サマリー ---")
sig_summary <- reg %>%
  filter(cohort %in% COHORT_ORDER, note == "ok",
         !is.na(significant), significant == TRUE) %>%
  mutate(gene_label = GENE_LABEL[gene]) %>%
  group_by(cohort) %>%
  summarise(sig_genes = paste(gene_label, collapse = ", "), .groups = "drop")

if (nrow(sig_summary) == 0) {
  log_msg("  （なし）")
} else {
  for (i in seq_len(nrow(sig_summary))) {
    log_msg(sprintf("  %-25s: %s", sig_summary$cohort[i], sig_summary$sig_genes[i]))
  }
}

log_msg("=== Step 16c: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 16c 完了\n")
cat("遺伝子順序（上から）:\n")
cat("  LAG3 > PDCD1 (PD-1) > CD274 (PD-L1) > PDCD1LG2 (PD-L2) > CTLA4 > TIGIT > HAVCR2 (TIM-3)\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/figS_checkpoint_specificity.pdf\n",        OUT_DIR))
cat(sprintf("  %s/figS_checkpoint_specificity_450dpi.png\n", OUT_DIR))
cat(sprintf("  %s/step16c_supplement_table.csv\n",           OUT_DIR))
cat(sprintf("  %s/step16c_log.txt\n",                        OUT_DIR))
cat("============================\n")
