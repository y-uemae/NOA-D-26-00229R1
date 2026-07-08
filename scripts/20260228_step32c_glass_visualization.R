# ============================================================
# ファイル名: step32c_glass_visualization.R
# 目的: GDC vs GLASS の Δβ比較図（Supplementary Figure S12）
# 出力:
#   - 32_glass_validation_screening/figS12_gdc_glass_delta_beta.pdf
#   - 32_glass_validation_screening/figS12_gdc_glass_delta_beta.png
# ============================================================

library(tidyverse)
library(ggplot2)

base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "32_glass_validation_screening")

# ============================================================
# 1. データ準備
# ============================================================

# GLASS結果
glass_res <- readRDS(file.path(out_dir, "glass_stage2_results.rds"))

# GDC結果（確定値）
gdc_res <- tibble(
  model         = c("M1_Base", "M2_Base_Purity",
                    "M3_Purity_G2M", "M4_Purity_top3"),
  beta_TP53     = c(0.2665, 0.2406, 0.2007, 0.1973),
  se_TP53       = c(NA, NA, NA, NA),   # 図にはΔβ%のみ使用
  p_TP53        = c(9.30e-7, 3.24e-6, 1.11e-4, 1.44e-4),
  delta_beta_pct= c(0, -9.7, -16.6, -18.0),
  adj_r2        = c(0.303, 0.371, 0.389, 0.390),
  n             = 442
)

# モデルラベルと順序を統一
model_labels <- c(
  "M1_Base"       = "Base\n(TP53 + immune scores)",
  "M2_Base_Purity"= "Base\n+ Tumor purity",
  "M3_Purity_G2M" = "Base + Purity\n+ G2M",
  "M4_Purity_top3"= "Base + Purity\n+ G2M + MYC + E2F"
)

# 両コホートを縦結合
plot_df <- bind_rows(
  gdc_res   %>% mutate(cohort = "GDC (n = 442)"),
  glass_res %>% select(model, beta_TP53, se_TP53, p_TP53,
                       delta_beta_pct, adj_r2, n) %>%
    mutate(cohort = "GLASS (n = 79)")
) %>%
  mutate(
    model_label = factor(model_labels[model],
                         levels = model_labels),
    cohort      = factor(cohort,
                         levels = c("GDC (n = 442)", "GLASS (n = 79)")),
    # 有意水準ラベル
    sig_label   = case_when(
      p_TP53 < 0.001 ~ "***",
      p_TP53 < 0.01  ~ "**",
      p_TP53 < 0.05  ~ "*",
      TRUE           ~ "ns"
    ),
    # Δβ%テキスト（基準モデルは空白）
    delta_label = ifelse(delta_beta_pct == 0, "",
                         paste0(delta_beta_pct, "%"))
  )

# ============================================================
# 2. Panel A: β_TP53の推移（dot + line plot）
# ============================================================

# seが欠損のGDCはエラーバーなし
# GLASSはseあり → 95%CI = ±1.96×se
plot_df <- plot_df %>%
  mutate(
    ci_lo = ifelse(!is.na(se_TP53), beta_TP53 - 1.96 * se_TP53, NA),
    ci_hi = ifelse(!is.na(se_TP53), beta_TP53 + 1.96 * se_TP53, NA)
  )

col_gdc   <- "#2166AC"   # 青
col_glass <- "#D6604D"   # 橙赤

p_beta <- ggplot(plot_df,
                 aes(x = model_label, y = beta_TP53,
                     color = cohort, group = cohort)) +
  # ゼロ参照線
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  # 折れ線
  geom_line(linewidth = 0.7, alpha = 0.8) +
  # GLASSのみ95%CI
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.12, linewidth = 0.6,
                data = filter(plot_df, cohort == "GLASS (n = 79)")) +
  # 点
  geom_point(size = 3.2, shape = 21, fill = "white",
             stroke = 1.5) +
  # 有意水準ラベル（点の上）
  geom_text(aes(label = sig_label, y = beta_TP53 + 0.04),
            size = 3.5, vjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c("GDC (n = 442)" = col_gdc,
                                "GLASS (n = 79)" = col_glass)) +
  scale_y_continuous(limits = c(-0.05, 0.75),
                     breaks = seq(0, 0.7, 0.1)) +
  labs(
    title = "A  TP53 regression coefficient (β) across adjustment models",
    x     = NULL,
    y     = expression(beta[TP53]~(log[2](TPM+1))),
    color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.4, "cm"),
    axis.text.x     = element_text(size = 9, lineheight = 1.2),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# ============================================================
# 3. Panel B: Δβ%の並列棒グラフ（基準モデルを除く3モデル）
# ============================================================

plot_df_delta <- plot_df %>%
  filter(model != "M1_Base") %>%
  mutate(
    model_label_short = factor(
      c("M2_Base_Purity" = "+ Purity",
        "M3_Purity_G2M"  = "+ Purity\n+ G2M",
        "M4_Purity_top3" = "+ Purity\n+ G2M\n+ MYC + E2F")[as.character(model)],
      levels = c("+ Purity", "+ Purity\n+ G2M",
                 "+ Purity\n+ G2M\n+ MYC + E2F")
    )
  )

p_delta <- ggplot(plot_df_delta,
                  aes(x = model_label_short,
                      y = delta_beta_pct,
                      fill = cohort)) +
  geom_col(position = position_dodge(width = 0.6),
           width = 0.5, alpha = 0.85) +
  geom_text(aes(label = paste0(delta_beta_pct, "%"),
                y = delta_beta_pct - 0.8),
            position = position_dodge(width = 0.6),
            size = 3.2, vjust = 1, color = "white", fontface = "bold") +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c("GDC (n = 442)" = col_gdc,
                               "GLASS (n = 79)" = col_glass)) +
  scale_y_continuous(limits = c(-35, 2),
                     breaks = seq(-30, 0, 10),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "B  Attenuation of TP53 coefficient (Δβ%) relative to base model",
    x     = "Additional covariate(s)",
    y     = "Δβ (% change from base model)",
    fill  = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11),
    legend.position = c(0.82, 0.12),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.4, "cm"),
    axis.text.x     = element_text(size = 9, lineheight = 1.2),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# ============================================================
# 4. 2パネルを縦に結合して保存
# ============================================================

# patchworkで結合
library(patchwork)

fig_s12 <- p_beta / p_delta +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    caption = paste0(
      "Sig. labels: *** p < 0.001, ** p < 0.01, * p < 0.05, ns p ≥ 0.05 (unadjusted).\n",
      "Error bars for GLASS: 95% CI (±1.96 × SE). GDC SE not available for error bars.\n",
      "Base model — GDC: LAG3 ~ TP53 + Tcell + APM + IFNγ; ",
      "GLASS: same (source and IDH excluded).\n",
      "Top-3 scores (G2M, MYC targets v1, E2F targets) were entered simultaneously."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, color = "grey40",
                                  lineheight = 1.3, hjust = 0)
    )
  )

# PDF（論文投稿用）
ggsave(
  file.path(out_dir, "figS12_gdc_glass_delta_beta.pdf"),
  fig_s12, width = 7, height = 9, device = cairo_pdf
)

# PNG（確認用・450 dpi）
ggsave(
  file.path(out_dir, "figS12_gdc_glass_delta_beta.png"),
  fig_s12, width = 7, height = 9, dpi = 450
)

cat("✅ Step 32c 完了\n")
cat("出力:\n")
cat("  figS12_gdc_glass_delta_beta.pdf\n")
cat("  figS12_gdc_glass_delta_beta.png\n")
