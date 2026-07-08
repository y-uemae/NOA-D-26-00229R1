# =============================================================================
# step13b_forest_plot_v3.R
# Forest plot 可視化専用（図のみ・全テキスト列なし・n数なし）
#
# キャップ対応:
#   - 個別コホート3行にも端点キャップを追加（FEとのバランス）
#   - FE: 破線 + lineend="butt" + キャップ
#   - 個別: 実線 + lineend="round" + キャップ
#   - RE: 太実線のみ（ダイヤが大きいのでキャップ不要）
#   - coord_cartesian(clip="off") 維持
#
# 入力:
#   11_meta_analysis/step11_meta_results.csv
# 出力（上書き）:
#   13_visualization/fig13b_forest_v3.pdf
#   13_visualization/fig13b_forest_v3_450dpi.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ragg)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR   <- here::here("results", "TP53", "20260221")
RESULT_CSV <- file.path(BASE_DIR, "11_meta_analysis/step11_meta_results.csv")
OUT_DIR    <- file.path(BASE_DIR, "13_visualization")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DPI       <- 450
WIDTH_IN  <- 6.0
HEIGHT_IN <- 5.0

COL_GDC_TCGA  <- "#3C5488"
COL_GDC_CPTAC <- "#E07B54"
COL_GLASS     <- "#8491B4"
COL_RE        <- "#E64B35"
COL_FE        <- "#777777"

# ── 1. データ読み込み ─────────────────────────────────────────────────────────
result_raw <- read_csv(RESULT_CSV, show_col_types = FALSE)
re_row <- filter(result_raw, model == "Random Effect")
fe_row <- filter(result_raw, model == "Fixed Effect")

cat("=== Pooled estimates ===\n")
result_raw %>%
  select(model, beta_pooled, ci_lo, ci_hi, p_val, I2, tau2, Q_p) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

# ── 2. 描画データ ─────────────────────────────────────────────────────────────
plot_df <- tibble(
  study    = c("GDC_TCGA", "GDC_CPTAC_HCMI", "GLASS_WXS",
               "Pooled_RE", "Pooled_FE"),
  y        = c(6, 5, 4, 2, 1),
  beta     = c(0.375, 0.135, 0.497,
               re_row$beta_pooled, fe_row$beta_pooled),
  ci_lo    = c(0.188, -0.009, 0.192,
               re_row$ci_lo, fe_row$ci_lo),
  ci_hi    = c(0.563,  0.279, 0.802,
               re_row$ci_hi, fe_row$ci_hi),
  col      = c(COL_GDC_TCGA, COL_GDC_CPTAC, COL_GLASS, COL_RE, COL_FE),
  row_type = c("ind", "ind", "ind", "re", "fe")
)

# キャップを付ける行（ind + fe）
cap_df <- filter(plot_df, row_type %in% c("ind", "fe"))

# キャップ縦幅（ind と fe で同じ高さに統一）
CAP_H_IND <- 0.15
CAP_H_FE  <- 0.15

cat("\n=== CI範囲確認 ===\n")
cat(sprintf("ci_lo min: %.3f\n", min(plot_df$ci_lo)))
cat(sprintf("ci_hi max: %.3f\n", max(plot_df$ci_hi)))

# ── 3. 軸範囲 ─────────────────────────────────────────────────────────────────
XMAX   <- max(plot_df$ci_hi, na.rm = TRUE) + 0.12
XMIN   <- min(plot_df$ci_lo, na.rm = TRUE) - 0.08
XF_MIN <- -0.25
XF_MAX <- XMAX - 0.05

Y_MIN <- 0.2
Y_MAX <- 6.8
Y_SEP <- 3.0

cat(sprintf("\ncoord_cartesian xlim: %.3f to %.3f\n", XMIN, XMAX))

# ── 4. footer ─────────────────────────────────────────────────────────────────
footer <- sprintf(
  "Heterogeneity: I2 = %.1f%%, tau2 = %.4f, Q p = %.3f  |  Random Effect model (REML) is the primary pooled estimate",
  re_row$I2, re_row$tau2, re_row$Q_p
)

# ── 5. 描画 ───────────────────────────────────────────────────────────────────
SZ_FOOT <- 2.6
SZ_TICK <- 3.0
SZ_XTTL <- 3.4

fig <- ggplot() +
  
  # ゼロライン
  annotate("segment",
           x = 0, xend = 0,
           y = Y_MIN + 0.05, yend = Y_MAX - 0.1,
           color = "black", linewidth = 0.5, lineend = "round") +
  
  # x軸ライン
  annotate("segment",
           x = XF_MIN, xend = XF_MAX,
           y = Y_MIN, yend = Y_MIN,
           color = "black", linewidth = 0.6, lineend = "round") +
  
  # x軸目盛り
  annotate("segment",
           x    = seq(-0.2, 0.8, by = 0.2),
           xend = seq(-0.2, 0.8, by = 0.2),
           y = Y_MIN, yend = Y_MIN - 0.12,
           color = "black", linewidth = 0.45) +
  annotate("text",
           x = seq(-0.2, 0.8, by = 0.2),
           y = Y_MIN - 0.30,
           label = sprintf("%.1f", seq(-0.2, 0.8, by = 0.2)),
           size = SZ_TICK, color = "black", hjust = 0.5) +
  
  # x軸タイトル
  annotate("text",
           x = (XF_MIN + XF_MAX) / 2,
           y = Y_MIN - 0.65,
           label = "beta coefficient (TP53 Mutant vs Wild-type)",
           size = SZ_XTTL, color = "black", hjust = 0.5) +
  
  # 区切り破線
  annotate("segment",
           x = XMIN, xend = XMAX,
           y = Y_SEP, yend = Y_SEP,
           color = "gray70", linewidth = 0.4, linetype = "dashed") +
  
  # ── CI線: 個別コホート（実線・round）
  geom_segment(
    data = filter(plot_df, row_type == "ind"),
    aes(x = ci_lo, xend = ci_hi, y = y, yend = y, color = col),
    linewidth = 1.3, lineend = "round"
  ) +
  
  # ── CI端点キャップ: 個別コホート（左）
  geom_segment(
    data = filter(plot_df, row_type == "ind"),
    aes(x = ci_lo, xend = ci_lo,
        y = y - CAP_H_IND, yend = y + CAP_H_IND,
        color = col),
    linewidth = 1.0, lineend = "butt"
  ) +
  
  # ── CI端点キャップ: 個別コホート（右）
  geom_segment(
    data = filter(plot_df, row_type == "ind"),
    aes(x = ci_hi, xend = ci_hi,
        y = y - CAP_H_IND, yend = y + CAP_H_IND,
        color = col),
    linewidth = 1.0, lineend = "butt"
  ) +
  
  # ── CI線: Random Effect（太・実線・round・キャップなし）
  geom_segment(
    data = filter(plot_df, row_type == "re"),
    aes(x = ci_lo, xend = ci_hi, y = y, yend = y),
    color = COL_RE, linewidth = 2.8, lineend = "round"
  ) +
  
  # ── CI線: Fixed Effect（破線・butt）
  geom_segment(
    data = filter(plot_df, row_type == "fe"),
    aes(x = ci_lo, xend = ci_hi, y = y, yend = y),
    color = COL_FE, linewidth = 1.4, linetype = "dashed",
    lineend = "butt"
  ) +
  
  # ── CI端点キャップ: Fixed Effect（左）
  geom_segment(
    data = filter(plot_df, row_type == "fe"),
    aes(x = ci_lo, xend = ci_lo,
        y = y - CAP_H_FE, yend = y + CAP_H_FE),
    color = COL_FE, linewidth = 1.0, lineend = "butt"
  ) +
  
  # ── CI端点キャップ: Fixed Effect（右）
  geom_segment(
    data = filter(plot_df, row_type == "fe"),
    aes(x = ci_hi, xend = ci_hi,
        y = y - CAP_H_FE, yend = y + CAP_H_FE),
    color = COL_FE, linewidth = 1.0, lineend = "butt"
  ) +
  
  # 点: 個別コホート（正方形）
  geom_point(
    data = filter(plot_df, row_type == "ind"),
    aes(x = beta, y = y, color = col),
    shape = 15, size = 4.2
  ) +
  
  # 点: Random Effect（大ダイヤ・赤）
  geom_point(
    data = filter(plot_df, row_type == "re"),
    aes(x = beta, y = y),
    shape = 23, size = 6.5,
    fill = COL_RE, color = COL_RE
  ) +
  
  # 点: Fixed Effect（小ダイヤ・灰）
  geom_point(
    data = filter(plot_df, row_type == "fe"),
    aes(x = beta, y = y),
    shape = 23, size = 4.5,
    fill = COL_FE, color = COL_FE
  ) +
  
  # footer
  annotate("text",
           x = XMIN, y = Y_MIN - 1.05,
           label = footer,
           hjust = 0, size = SZ_FOOT, color = "gray40") +
  
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(Y_MIN - 1.35, Y_MAX + 0.2),
                     expand = c(0, 0)) +
  scale_color_identity() +
  coord_cartesian(xlim = c(XMIN, XMAX), clip = "off") +
  
  theme_void() +
  theme(plot.margin = margin(8, 18, 8, 8, unit = "mm"))

# ── 6. PDF出力（上書き）──────────────────────────────────────────────────────
pdf_path <- file.path(OUT_DIR, "fig13b_forest_v3.pdf")
pdf(pdf_path, width = WIDTH_IN, height = HEIGHT_IN)
print(fig)
dev.off()
cat("PDF:", pdf_path, "\n")

# ── 7. PNG出力（ragg・上書き）────────────────────────────────────────────────
png_path <- file.path(OUT_DIR,
                      sprintf("fig13b_forest_v3_%ddpi.png", DPI))
agg_png(png_path,
        width = WIDTH_IN, height = HEIGHT_IN,
        units = "in", res = DPI, scaling = 1.0)
print(fig)
dev.off()
cat("PNG:", png_path, "\n")

cat("\n=== Step 13b v3 完了（上書き）===\n")
cat("キャップ: 個別コホート3行 + FE 全行に統一\n")
cat("RE: 太実線のみ（ダイヤで端点明示のためキャップ不要）\n")
cat(sprintf("XMIN=%.3f, XMAX=%.3f\n", XMIN, XMAX))
cat(sprintf("I2=%.1f%% | tau2=%.4f | Q p=%.3f\n",
            re_row$I2, re_row$tau2, re_row$Q_p))
