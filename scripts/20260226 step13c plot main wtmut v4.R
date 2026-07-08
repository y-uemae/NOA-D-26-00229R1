# =============================================================================
# 20260226_step13c_plot_main_wtmut_v4.R
# Step 13c v4
#
# v4での変更点:
#   - Panel A: GDC All + TCGA + CPTAC/HCMI の3パネルfacet
#   - 図内の統計注記（p値・Δmedian・δ）を全て削除（手書き対応）
#   - Unicode特殊文字を全てASCIIに変更（文字化け対策）
#   - fatten → middle.linewidth に変更（ggplot2 4.0.0対応）
#   - y軸上限を拡張（高い点が隠れないよう余裕を追加）
#
# 入力:
#   08_final_cohort/final_cohort.csv
#   05c_glass/glass_final_cohort_wxs_notcga.csv
# 出力:
#   13_visualization/fig13c_main_wtmut_v4.pdf
#   13_visualization/fig13c_main_wtmut_v4_450dpi.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ragg)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR  <- here::here("results", "TP53", "20260221")
GDC_CSV   <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
GLASS_CSV <- file.path(BASE_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")
OUT_DIR   <- file.path(BASE_DIR, "13_visualization")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 出力仕様
DPI        <- 450
WIDTH_MM   <- 174          # 2カラム幅
HEIGHT_MM  <- 90
MM_TO_INCH <- 1 / 25.4

# 共通色
COL_WT  <- "#AAAAAA"
COL_MUT <- "#E64B35"

# 共通テーマ
theme_paper <- theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 9, face = "bold"),
    axis.title       = element_text(size = 9),
    axis.text        = element_text(size = 8),
    axis.text.x      = element_text(size = 9, face = "bold"),
    legend.position  = "none",
    plot.title       = element_text(size = 10, face = "bold"),
    panel.spacing    = unit(4, "mm")
  )

# ── 1. GDCデータ整形 ──────────────────────────────────────────────────────────
gdc_raw <- read_csv(GDC_CSV, show_col_types = FALSE) %>%
  filter(include_flag == TRUE, grade == "Grade4") %>%
  mutate(
    tp53_label   = factor(if_else(tp53_status == "mutant", "Mut", "WT"),
                          levels = c("WT", "Mut")),
    source_label = factor(
      case_when(
        source == "TCGA"      ~ "TCGA\n(n=245)",
        source == "CPTAC_HCMI" ~ "CPTAC/HCMI\n(n=197)",
        TRUE                  ~ source
      ),
      levels = c("TCGA\n(n=245)", "CPTAC/HCMI\n(n=197)")
    )
  )

# GDC All 用データ（source_labelをAllに統一）
gdc_all <- gdc_raw %>%
  mutate(source_label = factor("All GDC\n(n=442)",
                               levels = "All GDC\n(n=442)"))

# 結合（All + TCGA + CPTAC/HCMI）
gdc_combined <- bind_rows(gdc_all, gdc_raw) %>%
  mutate(source_label = factor(source_label,
                               levels = c("All GDC\n(n=442)",
                                          "TCGA\n(n=245)",
                                          "CPTAC/HCMI\n(n=197)")))

# n注記用（facet別）
gdc_n_annot <- gdc_combined %>%
  group_by(source_label, tp53_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    x = if_else(tp53_label == "WT", 1, 2),
    label = sprintf("n=%d", n)
  )

cat("=== GDC n per facet ===\n")
print(gdc_n_annot)

# ── 2. GLASSデータ整形 ────────────────────────────────────────────────────────
glass_raw <- read_csv(GLASS_CSV, show_col_types = FALSE) %>%
  mutate(
    tp53_label = factor(if_else(tp53_status == "Mut", "Mut", "WT"),
                        levels = c("WT", "Mut"))
  )

glass_n_annot <- tibble(
  x     = c(1, 2),
  label = c(sprintf("n=%d", sum(glass_raw$tp53_label == "WT")),
            sprintf("n=%d", sum(glass_raw$tp53_label == "Mut")))
)

cat("\n=== GLASS n ===\n")
print(glass_n_annot)

# ── 3. y軸範囲の統一（余裕を大きめに）────────────────────────────────────────
y_max  <- max(c(gdc_combined$LAG3_log2tpm,
                glass_raw$LAG3_log2tpm), na.rm = TRUE)
y_top  <- ceiling(y_max) + 0.5    # 最大値の切り上げ + 0.5（点が隠れない余裕）
y_lim  <- c(-0.2, y_top)
n_y    <- y_lim[1] + 0.08         # n注記y位置（下端）

cat(sprintf("\ny_max=%.2f, y_top=%.2f\n", y_max, y_top))

# ── 4. Panel A: GDC 3パネル ───────────────────────────────────────────────────
panel_A <- ggplot(gdc_combined,
                  aes(x = tp53_label, y = LAG3_log2tpm,
                      color = tp53_label)) +
  geom_jitter(width = 0.18, alpha = 0.30, size = 0.7, shape = 16) +
  geom_boxplot(aes(fill = tp53_label),
               width = 0.42, outlier.shape = NA,
               alpha = 0.15, color = "black", linewidth = 0.45) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.38, linewidth = 0.7,
               middle.linewidth = 0.7,
               color = "black") +
  # n注記（x軸下）
  geom_text(data = gdc_n_annot,
            aes(x = x, y = n_y, label = label),
            inherit.aes = FALSE,
            size = 2.3, color = "gray40") +
  facet_wrap(~ source_label, nrow = 1) +
  scale_color_manual(values = c("WT" = COL_WT, "Mut" = COL_MUT)) +
  scale_fill_manual( values = c("WT" = COL_WT, "Mut" = COL_MUT)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0),
                     breaks = seq(0, floor(y_top), by = 1)) +
  scale_x_discrete(labels = c("WT", "Mut")) +
  labs(title = "A   GDC Grade 4",
       x = NULL,
       y = "LAG3 expression [log2(TPM+1)]") +
  theme_paper

# ── 5. Panel B: GLASS ─────────────────────────────────────────────────────────
panel_B <- ggplot(glass_raw,
                  aes(x = tp53_label, y = LAG3_log2tpm,
                      color = tp53_label)) +
  geom_jitter(width = 0.18, alpha = 0.40, size = 0.9, shape = 16) +
  geom_boxplot(aes(fill = tp53_label),
               width = 0.42, outlier.shape = NA,
               alpha = 0.15, color = "black", linewidth = 0.45) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.38, linewidth = 0.7,
               middle.linewidth = 0.7,
               color = "black") +
  # n注記
  geom_text(data = glass_n_annot,
            aes(x = x, y = n_y, label = label),
            inherit.aes = FALSE,
            size = 2.3, color = "gray40") +
  scale_color_manual(values = c("WT" = COL_WT, "Mut" = COL_MUT)) +
  scale_fill_manual( values = c("WT" = COL_WT, "Mut" = COL_MUT)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0),
                     breaks = seq(0, floor(y_top), by = 1)) +
  scale_x_discrete(labels = c("WT", "Mut")) +
  labs(title = "B   GLASS (WXS, non-TCGA)",
       x = NULL, y = NULL) +
  theme_paper

# ── 6. 結合（Panel A : Panel B = 3 : 1）─────────────────────────────────────
fig_combined <- panel_A + panel_B +
  plot_layout(widths = c(3, 1))

# ── 7. PDF出力 ────────────────────────────────────────────────────────────────
pdf_path <- file.path(OUT_DIR, "fig13c_main_wtmut_v4.pdf")
pdf(pdf_path,
    width  = WIDTH_MM * MM_TO_INCH,
    height = HEIGHT_MM * MM_TO_INCH)
print(fig_combined)
dev.off()
cat("PDF:", pdf_path, "\n")

# ── 8. PNG出力（ragg・高画質）────────────────────────────────────────────────
png_path <- file.path(OUT_DIR,
                      sprintf("fig13c_main_wtmut_v4_%ddpi.png", DPI))
agg_png(png_path,
        width   = WIDTH_MM * MM_TO_INCH,
        height  = HEIGHT_MM * MM_TO_INCH,
        units   = "in",
        res     = DPI,
        scaling = 1.0)
print(fig_combined)
dev.off()
cat("PNG:", png_path, "\n")

cat("\n=== Step 13c v4 完了 ===\n")
cat("DPI:", DPI, "/ Size:", WIDTH_MM, "x", HEIGHT_MM, "mm\n")
cat("Panel A: All GDC + TCGA + CPTAC/HCMI (3 facets)\n")
cat("Panel B: GLASS\n")
cat("Notes: 統計注記なし（手書き対応）\n")
cat("ragg version:", as.character(packageVersion("ragg")), "\n")
