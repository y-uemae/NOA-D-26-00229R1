# =============================================================================
# step31_visualization.R 【v1.4 write_excel_csv使用版】
# 修正: BOM付きCSVをwrite_excel_csv()で一発保存（writeBin廃止）
# 作成日: 2026-02-27
# バージョン履歴:
#   v1.0 - 初版
#   v1.1 - 文字化け修正、NA対処
#   v1.2 - BOM付きCSV保存方法を修正
#   v1.3 - CSVヘッダー消失修正、Fig B日本語ラベル英語化
#   v1.4 - write_excel_csv()でBOM付き保存に統一（バイナリ書き込みエラー解消）
# =============================================================================

library(tidyverse)

# =============================================================================
# 0. パス設定
# =============================================================================
base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "31_visualization")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("フォルダ作成:", out_dir, "\n")
} else {
  cat("出力フォルダ確認済み:", out_dir, "\n")
}

stage2   <- readRDS(file.path(base_dir, "30_screening", "stage2_results.rds"))
check_bc <- readRDS(file.path(base_dir, "30_screening", "check_bc_delta_beta_extended.rds"))

# =============================================================================
# 1. NA確認
# =============================================================================
cat("NA確認:\n")
na_rows <- stage2 %>% filter(is.na(p_score) | is.na(p_score_adj))
if (nrow(na_rows) > 0) {
  cat("p_scoreにNAがある行:\n")
  print(na_rows %>% select(score, beta_score, p_score, p_score_adj))
} else {
  cat("NAなし\n")
}

# =============================================================================
# 2. ラベル整形用ヘルパー
# =============================================================================
clean_label <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

# =============================================================================
# 3. Fig A: Stage2 Delta-beta% bar plot（上位10）
# =============================================================================
cat("\nFig A: Delta-beta% bar plot作成中...\n")

top10 <- stage2 %>%
  arrange(desc(abs(delta_beta_pct))) %>%
  head(10) %>%
  mutate(
    label     = clean_label(score),
    sig       = case_when(
      is.na(p_score_adj)  ~ "NA",
      p_score_adj < 0.001 ~ "***",
      p_score_adj < 0.01  ~ "**",
      p_score_adj < 0.05  ~ "*",
      TRUE                ~ "ns"
    ),
    direction = if_else(delta_beta_pct < 0, "Attenuation", "Amplification"),
    label     = fct_reorder(label, abs(delta_beta_pct))
  )

fig_a <- ggplot(top10, aes(x = delta_beta_pct, y = label, fill = direction)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(
    label = paste0(sig, "  ", sprintf("%.1f%%", delta_beta_pct)),
    x     = delta_beta_pct + ifelse(delta_beta_pct < 0, -0.3, 0.3)
  ),
  hjust  = ifelse(top10$delta_beta_pct < 0, 1, 0),
  size   = 3.2, color = "gray20") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "gray40") +
  scale_fill_manual(
    values = c("Attenuation" = "#4C72B0", "Amplification" = "#DD8452")
  ) +
  scale_x_continuous(
    limits = c(-22, 8),
    breaks = seq(-20, 5, by = 5),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = expression("Change in " * beta[TP53] * " after adding Hallmark scores"),
    subtitle = "Base model: LAG3 ~ TP53 + source + IDH + Tcell + APM + IFNg + Purity",
    x        = expression(Delta * beta[TP53] * " (%)"),
    y        = NULL,
    fill     = NULL,
    caption  = "Score significance (BH-adjusted): *** p<0.001, ** p<0.01, * p<0.05, ns"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "gray40"),
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.caption       = element_text(size = 8, color = "gray50")
  )

ggsave(file.path(out_dir, "figS_delta_beta_barplot.pdf"),
       fig_a, width = 7, height = 5.5)
ggsave(file.path(out_dir, "figS_delta_beta_barplot.png"),
       fig_a, width = 7, height = 5.5, dpi = 300)
cat("Fig A保存完了\n")

# =============================================================================
# 4. Fig B: モデル比較 dot plot（Check BC）
# =============================================================================
cat("Fig B: モデル比較 dot plot作成中...\n")

model_df <- check_bc %>%
  mutate(
    model = recode(model,
                   "Base（Purityなし）"             = "Base (w/o Purity)",
                   "Base（Purityあり）"             = "Base (w/ Purity)",
                   "Base + HALLMARK_G2M_CHECKPOINT" = "Base + G2M Checkpoint",
                   "Base + top2"                    = "Base + G2M + MYC_V1",
                   "Base + top3"                    = "Base + G2M + MYC_V1 + E2F"
    ),
    model = factor(model, levels = rev(c(
      "Base (w/o Purity)",
      "Base (w/ Purity)",
      "Base + G2M Checkpoint",
      "Base + G2M + MYC_V1",
      "Base + G2M + MYC_V1 + E2F"
    ))),
    sig_label = formatC(p_tp53, format = "e", digits = 1)
  )

beta_base_nopurity <- model_df %>%
  filter(model == "Base (w/o Purity)") %>%
  pull(beta_tp53_test)

fig_b <- ggplot(model_df, aes(x = beta_tp53_test, y = model)) +
  geom_point(size = 3.5, color = "#4C72B0") +
  geom_vline(xintercept = beta_base_nopurity,
             linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_text(aes(label = paste0(
    "b=", sprintf("%.3f", beta_tp53_test),
    "  p=", sig_label,
    "  D=", sprintf("%.1f%%", delta_beta_pct)
  )),
  hjust = -0.08, size = 3, color = "gray20") +
  scale_x_continuous(
    limits = c(0.15, 0.42),
    breaks = seq(0.15, 0.40, by = 0.05)
  ) +
  labs(
    title    = expression("Robustness of " * beta[TP53] * " across adjusted models"),
    subtitle = "TP53-LAG3 association after progressive covariate adjustment",
    x        = expression(beta[TP53]),
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "gray40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(file.path(out_dir, "figS_model_comparison.pdf"),
       fig_b, width = 8.5, height = 3.5)
ggsave(file.path(out_dir, "figS_model_comparison.png"),
       fig_b, width = 8.5, height = 3.5, dpi = 300)
cat("Fig B保存完了\n")

# =============================================================================
# 5. Supplement Table
# ★ write_excel_csv() = UTF-8 BOM付き、ヘッダーあり、Excelで開いても化けない
# =============================================================================
cat("Supplement Table作成中...\n")

suppl_table <- stage2 %>%
  arrange(desc(abs(delta_beta_pct))) %>%
  mutate(
    Gene_Set        = clean_label(score),
    Beta_TP53_Base  = round(beta_tp53_base, 4),
    Beta_TP53_Test  = round(beta_tp53_test, 4),
    Delta_Beta_Pct  = round(delta_beta_pct, 1),
    P_TP53_Test     = formatC(p_tp53_test,  format = "e", digits = 2),
    Beta_Score      = ifelse(is.na(beta_score), "NA",
                             as.character(round(beta_score, 4))),
    P_Score_Nominal = ifelse(is.na(p_score), "NA",
                             formatC(p_score, format = "e", digits = 2)),
    P_Score_BH_Adj  = ifelse(is.na(p_score_adj), "NA",
                             formatC(p_score_adj, format = "e", digits = 2)),
    AdjR2_Test      = round(r2_adj_test, 3)
  ) %>%
  select(
    Gene_Set,
    Beta_TP53_Base, Beta_TP53_Test, Delta_Beta_Pct,
    P_TP53_Test,
    Beta_Score, P_Score_Nominal, P_Score_BH_Adj,
    AdjR2_Test
  )

out_suppl <- file.path(out_dir, "suppltable_screening.csv")

# ★write_excel_csv = BOM付きUTF-8、col_names=TRUE（デフォルト）
write_excel_csv(suppl_table, out_suppl)
cat("Supplement Table保存完了:", out_suppl, "\n")
cat("行数:", nrow(suppl_table), "/ 列数:", ncol(suppl_table), "\n")

# =============================================================================
# 6. 完了サマリー
# =============================================================================
cat("\n========================================\n")
cat("Step 31 完了\n")
cat("========================================\n")
cat("出力先:", out_dir, "\n")
cat("  figS_delta_beta_barplot.pdf / .png\n")
cat("  figS_model_comparison.pdf / .png\n")
cat("  suppltable_screening.csv (UTF-8 BOM / Excel対応)\n")
