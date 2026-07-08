# =============================================================================
# step11c_forest_plot.R
# Forest plot 可視化専用スクリプト（計算なし・CSV読み込みのみ）
#
# 入力:
#   11_meta_analysis/step11_subgroup_regression.csv  （個別コホートβ・SE）
#   11_meta_analysis/step11_meta_results.csv          （pooled推定値）
# 出力:
#   11_meta_analysis/step11c_forest_plot.pdf
#   11_meta_analysis/step11c_forest_plot_450dpi.png
#
# 修正内容:
#   - x軸上限を拡張（Fixed Effectのバーが切れる問題を解消）
#   - I²・τ²・Q p値をfootnoteに明示
#   - 統計注記はfootnote行に集約（列追加なし）
#   - Unicode不使用（文字化け対策）
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ragg)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR <- here::here("results", "TP53", "20260221")
OUT_DIR  <- file.path(BASE_DIR, "11_meta_analysis")

DPI        <- 450
WIDTH_MM   <- 174
HEIGHT_MM  <- 110
MM_TO_INCH <- 1 / 25.4

# ── 1. CSV読み込み ────────────────────────────────────────────────────────────
subgroup <- read_csv(file.path(OUT_DIR, "step11_subgroup_regression.csv"),
                     show_col_types = FALSE)
meta_res <- read_csv(file.path(OUT_DIR, "step11_meta_results.csv"),
                     show_col_types = FALSE)

cat("=== 個別コホート ===\n")
subgroup %>%
  select(study, n_mut, n_wt, beta_TP53, ci_lo, ci_hi, p_val) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

cat("\n=== Pooled推定 ===\n")
meta_res %>%
  select(model, beta_pooled, ci_lo, ci_hi, p_val, I2, tau2, Q_p) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

# ── 2. プロット用データ整形 ───────────────────────────────────────────────────

# 個別コホートの表示ラベルとn数
individual_df <- subgroup %>%
  transmute(
    label    = case_when(
      study == "GDC_TCGA"       ~ "GDC - TCGA",
      study == "GDC_CPTAC_HCMI" ~ "GDC - CPTAC/HCMI",
      study == "GLASS_WXS"      ~ "GLASS - WXS (non-TCGA)",
      TRUE ~ study
    ),
    sublabel = sprintf("Mut n=%d, WT n=%d", n_mut, n_wt),
    beta     = beta_TP53,
    ci_lo,
    ci_hi,
    p_val,
    type     = "individual",
    shape    = 15L,   # 塗りつぶし正方形
    size     = 3.0,
    color    = "gray30"
  )

# Pooled推定
fe_row <- meta_res %>% filter(model == "Fixed Effect")
re_row <- meta_res %>% filter(model == "Random Effect")

pooled_df <- tibble(
  label    = c("Pooled (Fixed Effect)", "Pooled (Random Effect)"),
  sublabel = c("", ""),
  beta     = c(fe_row$beta_pooled, re_row$beta_pooled),
  ci_lo    = c(fe_row$ci_lo,       re_row$ci_lo),
  ci_hi    = c(fe_row$ci_hi,       re_row$ci_hi),
  p_val    = c(fe_row$p_val,       re_row$p_val),
  type     = "pooled",
  shape    = c(18L, 23L),   # 18=菱形塗, 23=菱形枠
  size     = c(4.0, 5.0),
  color    = c("gray50", "#E64B35")
)

# y軸順序（上から下）
# 個別3行 → 空白行 → Fixed → Random
plot_df <- bind_rows(
  individual_df,
  tibble(label = "separator", sublabel = "", beta = NA,
         ci_lo = NA, ci_hi = NA, p_val = NA,
         type = "sep", shape = NA_integer_,
         size = NA_real_, color = NA_character_),
  pooled_df
) %>%
  mutate(
    y_pos = rev(row_number()),
    label = factor(label, levels = rev(label))
  )

# 表示用p値テキスト（ASCII）
plot_df <- plot_df %>%
  mutate(
    p_text = case_when(
      is.na(p_val)  ~ "",
      p_val < 0.001 ~ "p < 0.001",
      p_val < 0.01  ~ sprintf("p = %.3f", p_val),
      TRUE          ~ sprintf("p = %.3f", p_val)
    ),
    beta_text = case_when(
      is.na(beta) ~ "",
      TRUE ~ sprintf("%.3f (%.3f, %.3f)", beta, ci_lo, ci_hi)
    )
  )

cat("\n=== プロット用データ ===\n")
plot_df %>%
  select(label, y_pos, beta, ci_lo, ci_hi, p_val, type) %>%
  print()

# ── 3. x軸範囲（バーが切れない設定）─────────────────────────────────────────
x_max_data <- max(plot_df$ci_hi, na.rm = TRUE)
x_min_data <- min(plot_df$ci_lo, na.rm = TRUE)
x_lo   <- floor(x_min_data * 10) / 10 - 0.1
x_hi   <- ceiling(x_max_data * 10) / 10 + 0.15   # 右に余白
x_lo   <- min(x_lo, -0.25)   # 最低でも-0.25まで表示

cat(sprintf("\nx軸範囲: %.2f to %.2f\n", x_lo, x_hi))

# ── 4. footnote テキスト（I²等）─────────────────────────────────────────────
# ASCII表記で統一
footnote <- sprintf(
  "Heterogeneity: I2 = %.1f%%, tau2 = %.4f, Q p = %.3f  |  Random Effect model (REML) is the primary pooled estimate",
  re_row$I2, re_row$tau2, re_row$Q_p
)
cat("\nfootnote:", footnote, "\n")

# ── 5. Forest plot 作成 ───────────────────────────────────────────────────────

# separatorのy位置
sep_y <- filter(plot_df, type == "sep")$y_pos

theme_forest <- theme_classic(base_size = 9) +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.line.y      = element_blank(),
    axis.title       = element_text(size = 9),
    axis.text.x      = element_text(size = 8),
    plot.title       = element_text(size = 10, face = "bold"),
    plot.caption     = element_text(size = 7, color = "gray40",
                                    hjust = 0, margin = margin(t = 6)),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

p_forest <- ggplot(filter(plot_df, type != "sep")) +
  # 参照線（β=0）
  geom_vline(xintercept = 0, linetype = "solid",
             color = "black", linewidth = 0.5) +
  # CI バー（個別）
  geom_errorbarh(
    data = filter(plot_df, type == "individual"),
    aes(xmin = ci_lo, xmax = ci_hi, y = y_pos),
    height = 0.25, linewidth = 0.6, color = "gray30"
  ) +
  # CI バー（pooled）
  geom_errorbarh(
    data = filter(plot_df, type == "pooled"),
    aes(xmin = ci_lo, xmax = ci_hi, y = y_pos,
        color = label),
    height = 0.3, linewidth = 0.8
  ) +
  # ポイント（個別: 正方形）
  geom_point(
    data = filter(plot_df, type == "individual"),
    aes(x = beta, y = y_pos,
        size = 1 / (beta - ci_lo)^2),   # SEに基づくサイズ（重みの視覚化）
    shape = 15, color = "gray30"
  ) +
  # ポイント（Fixed Effect: 小菱形・灰色）
  geom_point(
    data = filter(plot_df, type == "pooled", color == "gray50"),
    aes(x = beta, y = y_pos),
    shape = 18, size = 4, color = "gray50"
  ) +
  # ポイント（Random Effect: 大菱形・赤）
  geom_point(
    data = filter(plot_df, type == "pooled", color == "#E64B35"),
    aes(x = beta, y = y_pos),
    shape = 23, size = 5,
    fill = "#E64B35", color = "#E64B35"
  ) +
  # 区切り線
  geom_hline(yintercept = sep_y + 0.5,
             linetype = "dashed", color = "gray60", linewidth = 0.4) +
  # 左列: study ラベル
  geom_text(
    data = filter(plot_df, type == "individual"),
    aes(x = x_lo, y = y_pos + 0.28, label = label),
    hjust = 0, size = 3.0, color = "gray20", fontface = "bold"
  ) +
  geom_text(
    data = filter(plot_df, type == "individual"),
    aes(x = x_lo, y = y_pos - 0.22, label = sublabel),
    hjust = 0, size = 2.4, color = "gray50"
  ) +
  geom_text(
    data = filter(plot_df, type == "pooled"),
    aes(x = x_lo, y = y_pos, label = label, color = label),
    hjust = 0, size = 3.0, fontface = "bold"
  ) +
  # 右列: β [95% CI]
  geom_text(
    data = filter(plot_df, type != "sep"),
    aes(x = x_hi + 0.01, y = y_pos, label = beta_text),
    hjust = 0, size = 2.6, color = "gray20"
  ) +
  # 右列: p値
  geom_text(
    data = filter(plot_df, type != "sep"),
    aes(x = x_hi + 0.16, y = y_pos, label = p_text),
    hjust = 0, size = 2.6, color = "gray20"
  ) +
  # ヘッダー
  annotate("text", x = x_lo,   y = max(plot_df$y_pos, na.rm=TRUE) + 0.7,
           label = "Study", hjust = 0, size = 3.0, fontface = "bold") +
  annotate("text", x = x_hi + 0.01, y = max(plot_df$y_pos, na.rm=TRUE) + 0.7,
           label = "beta [95% CI]", hjust = 0, size = 2.8, fontface = "bold") +
  annotate("text", x = x_hi + 0.16, y = max(plot_df$y_pos, na.rm=TRUE) + 0.7,
           label = "p value", hjust = 0, size = 2.8, fontface = "bold") +
  # スケール
  scale_x_continuous(
    limits = c(x_lo, x_hi + 0.32),
    breaks = seq(round(x_lo, 1), round(x_hi, 1), by = 0.2),
    expand = c(0, 0)
  ) +
  scale_size_continuous(range = c(2.5, 5.0), guide = "none") +
  scale_color_manual(
    values = c(
      "Pooled (Fixed Effect)"  = "gray50",
      "Pooled (Random Effect)" = "#E64B35"
    ),
    guide = "none"
  ) +
  labs(
    title   = "C   Meta-analysis (GDC + GLASS)",
    x       = "Regression coefficient (beta_TP53) with 95% CI",
    y       = NULL,
    caption = footnote
  ) +
  theme_forest

# ── 6. PDF出力 ────────────────────────────────────────────────────────────────
pdf_path <- file.path(OUT_DIR, "step11c_forest_plot.pdf")
pdf(pdf_path,
    width  = WIDTH_MM * MM_TO_INCH,
    height = HEIGHT_MM * MM_TO_INCH)
print(p_forest)
dev.off()
cat("PDF:", pdf_path, "\n")

# ── 7. PNG出力（ragg）────────────────────────────────────────────────────────
png_path <- file.path(OUT_DIR,
                      sprintf("step11c_forest_plot_%ddpi.png", DPI))
agg_png(png_path,
        width   = WIDTH_MM * MM_TO_INCH,
        height  = HEIGHT_MM * MM_TO_INCH,
        units   = "in",
        res     = DPI,
        scaling = 1.0)
print(p_forest)
dev.off()
cat("PNG:", png_path, "\n")

cat("\n=== step11c 完了 ===\n")
cat("入力: step11_subgroup_regression.csv, step11_meta_results.csv\n")
cat("出力: step11c_forest_plot.pdf / .png\n")
cat("I2 =", round(re_row$I2, 1), "% | tau2 =", round(re_row$tau2, 4),
    "| Q p =", round(re_row$Q_p, 3), "\n")
