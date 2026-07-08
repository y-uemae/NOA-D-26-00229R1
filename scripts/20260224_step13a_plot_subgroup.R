# =============================================================================
# 20260224_step13a_plot_subgroup.R
# Step 13a: TP53 4群サブグループ boxplot（Supplement Figure用）
#
# 内容: Grade4全体（TCGA+CPTAC/HCMI）の4群比較
#   WT / Hotspot / Truncating / Other_missense
#   jitter + boxplot、事前比較2本のみp値注記
#
# 入力: 12_subgroup/step12b_subgroup_classified.csv
# 出力:
#   13_visualization/figS_subgroup_4group.pdf
#   13_visualization/figS_subgroup_4group_450dpi.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ragg)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR <- here::here("results", "TP53", "20260221")
IN_CSV   <- file.path(BASE_DIR, "12_subgroup/step12b_subgroup_classified.csv")
OUT_DIR  <- file.path(BASE_DIR, "13_visualization")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DPI       <- 450
WIDTH_IN  <-  7.0
HEIGHT_IN <-  5.5

# 引継書固定色
COL_4 <- c(
  "WT"             = "#AAAAAA",
  "Hotspot"        = "#E64B35",
  "Truncating"     = "#4DBBD5",
  "Other_missense" = "#00A087"
)

# ── 1. データ読み込み・整形 ───────────────────────────────────────────────────
df <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(
    tp53_class4 = factor(tp53_class4,
                         levels = c("WT","Hotspot","Truncating","Other_missense"))
  )

# x軸ラベル（n付き）
n_per_group <- df %>% count(tp53_class4)
x_labels <- n_per_group %>%
  mutate(label = sprintf("%s\n(n=%d)", tp53_class4, n)) %>%
  pull(label)
names(x_labels) <- n_per_group$tp53_class4

cat("=== 群別n ===\n")
print(n_per_group)

# ── 2. 統計注記用p値計算 ──────────────────────────────────────────────────────
wt_vals <- df$LAG3_log2tpm[df$tp53_class4 == "WT"]

calc_p <- function(group) {
  vals <- df$LAG3_log2tpm[df$tp53_class4 == group]
  wilcox.test(vals, wt_vals, exact = FALSE)$p.value
}

p_hot   <- calc_p("Hotspot")
p_trunc <- calc_p("Truncating")

# BH補正（3本まとめて・事前比較2本+Other）
p_other <- calc_p("Other_missense")
p_bh    <- p.adjust(c(p_hot, p_trunc, p_other), method = "BH")

fmt_p <- function(p) {
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

cat(sprintf("\nHotspot vs WT:    p_BH = %.4f\n", p_bh[1]))
cat(sprintf("Truncating vs WT: p_BH = %.4f\n", p_bh[2]))

# ── 3. ブラケット・p値注記の座標 ─────────────────────────────────────────────
y_max   <- max(df$LAG3_log2tpm, na.rm = TRUE)
y_top   <- ceiling(y_max * 10) / 10 + 0.2
y_lim   <- c(-0.15, y_top + 0.85)

# ブラケット高さ（2段）
brk1_y  <- y_top + 0.15   # Hotspot vs WT
brk2_y  <- y_top + 0.50   # Truncating vs WT
p1_y    <- brk1_y + 0.10
p2_y    <- brk2_y + 0.10

# x位置（factor levelの整数対応: WT=1, Hot=2, Trunc=3, Other=4）
x_wt    <- 1
x_hot   <- 2
x_trunc <- 3

# ブラケット用tibble
brk_df <- tibble(
  x1    = c(x_wt,    x_wt),
  x2    = c(x_hot,   x_trunc),
  y     = c(brk1_y,  brk2_y),
  py    = c(p1_y,    p2_y),
  label = c(fmt_p(p_bh[1]), fmt_p(p_bh[2]))
)

# ── 4. テーマ ────────────────────────────────────────────────────────────────
theme_paper <- theme_classic(base_size = 11) +
  theme(
    axis.title   = element_text(size = 11),
    axis.text.x  = element_text(size = 10, lineheight = 0.9),
    axis.text.y  = element_text(size = 10),
    legend.position = "none",
    plot.title   = element_text(size = 12, face = "bold"),
    plot.margin  = margin(8, 8, 8, 8, unit = "mm")
  )

# ── 5. 描画 ───────────────────────────────────────────────────────────────────
fig <- ggplot(df, aes(x = tp53_class4, y = LAG3_log2tpm,
                      color = tp53_class4, fill = tp53_class4)) +
  
  # jitter
  geom_jitter(width = 0.18, alpha = 0.30, size = 0.9, shape = 16) +
  
  # boxplot
  geom_boxplot(width = 0.45, outlier.shape = NA,
               alpha = 0.15, color = "black", linewidth = 0.5) +
  
  # median横線
  stat_summary(fun = median, geom = "crossbar",
               width = 0.40, linewidth = 0.75,
               color = "black", fatten = 1) +
  
  # ブラケット横線
  geom_segment(data = brk_df,
               aes(x = x1 + 0.08, xend = x2 - 0.08,
                   y = y, yend = y),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.4) +
  
  # ブラケット縦線（左）
  geom_segment(data = brk_df,
               aes(x = x1 + 0.08, xend = x1 + 0.08,
                   y = y - 0.06, yend = y),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.4) +
  
  # ブラケット縦線（右）
  geom_segment(data = brk_df,
               aes(x = x2 - 0.08, xend = x2 - 0.08,
                   y = y - 0.06, yend = y),
               inherit.aes = FALSE,
               color = "black", linewidth = 0.4) +
  
  # p値テキスト
  geom_text(data = brk_df,
            aes(x = (x1 + x2) / 2, y = py, label = label),
            inherit.aes = FALSE,
            size = 3.3, color = "black", fontface = "italic") +
  
  scale_color_manual(values = COL_4) +
  scale_fill_manual( values = COL_4) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0),
                     breaks = seq(0, floor(y_top) + 1, by = 1)) +
  
  labs(
    title = "TP53 subgroup analysis  –  GDC Grade 4",
    x     = NULL,
    y     = "LAG3 expression  [log2(TPM+1)]"
  ) +
  theme_paper

# ── 6. 出力 ──────────────────────────────────────────────────────────────────
pdf_path <- file.path(OUT_DIR, "figS_subgroup_4group.pdf")
pdf(pdf_path, width = WIDTH_IN, height = HEIGHT_IN)
print(fig)
dev.off()
cat("PDF:", pdf_path, "\n")

png_path <- file.path(OUT_DIR,
                      sprintf("figS_subgroup_4group_%ddpi.png", DPI))
agg_png(png_path,
        width = WIDTH_IN, height = HEIGHT_IN,
        units = "in", res = DPI, scaling = 1.0)
print(fig)
dev.off()
cat("PNG:", png_path, "\n")

cat("\n=== Step 13a 完了 ===\n")
cat("  figS_subgroup_4group.pdf\n")
cat("  figS_subgroup_4group_450dpi.png\n")
