# ============================================================
# step_S1_grade23_violin.R
# Fig. S1: LAG3 expression by TP53 status in WHO Grade 2/3
# Input : 08_final_cohort/final_cohort.csv
# Output: 09_statistics/figS1_grade23.pdf / .png (450 dpi)
# ============================================================

library(tidyverse)
library(ggplot2)

# ---- 0. パス設定 ----
BASE_DIR  <- here::here("results", "TP53", "20260221")
IN_FILE   <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
OUT_DIR   <- file.path(BASE_DIR, "09_statistics")
OUT_PDF   <- file.path(OUT_DIR, "figS1_grade23.pdf")
OUT_PNG   <- file.path(OUT_DIR, "figS1_grade23.png")

# ---- 1. 色定数（引継書共通設定）----
COL_MUT <- "#E64B35"
COL_WT  <- "#AAAAAA"

# ---- 2. データ読み込み・整形 ----
df_raw <- read_csv(IN_FILE, show_col_types = FALSE) %>%
  filter(include_flag == TRUE,
         grade %in% c("Grade2", "Grade3"),
         tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_label = factor(
      ifelse(tp53_status == "mutant", "TP53 Mutant", "TP53 Wild-type"),
      levels = c("TP53 Wild-type", "TP53 Mutant")
    ),
    grade_label = factor(
      ifelse(grade == "Grade2", "WHO Grade 2", "WHO Grade 3"),
      levels = c("WHO Grade 2", "WHO Grade 3")
    )
  ) %>%
  select(grade_label, tp53_label, LAG3 = LAG3_log2tpm) %>%
  drop_na(LAG3)

cat("Grade2/3 読込完了\n")
print(df_raw %>% count(grade_label, tp53_label))

# ---- 3. アノテーション用統計量を summary CSVから読み込む ----
# （step09_lag3_summary.csv の SubAnalysis 行を使用）
df_sum <- read_csv(
  file.path(OUT_DIR, "step09_lag3_summary.csv"),
  show_col_types = FALSE
) %>%
  filter(analysis == "SubAnalysis") %>%
  mutate(
    grade_label = factor(
      ifelse(subset == "Grade2", "WHO Grade 2", "WHO Grade 3"),
      levels = c("WHO Grade 2", "WHO Grade 3")
    ),
    p_label = case_when(
      p_BH < 0.001 ~ formatC(p_BH, format = "e", digits = 2),
      TRUE         ~ formatC(p_BH, format = "f", digits = 4)
    ),
    annot = paste0(
      "p_BH = ", p_label, "\n",
      "delta = ", sprintf("%.3f", cliffs_delta), "\n",
      "n(Mut) = ", n_mut, ", n(WT) = ", n_wt
    )
  )

# アノテーションのy位置（各グレードの最大値 + 余白）
df_ymax <- df_raw %>%
  group_by(grade_label) %>%
  summarise(y_max = max(LAG3, na.rm = TRUE), .groups = "drop")

df_annot <- df_sum %>%
  left_join(df_ymax, by = "grade_label") %>%
  mutate(y_pos = y_max + 0.15)

# ---- 4. プロット作成 ----
p <- ggplot(df_raw, aes(x = tp53_label, y = LAG3, fill = tp53_label)) +
  
  # バイオリン
  geom_violin(
    trim      = TRUE,
    scale     = "width",
    alpha     = 0.6,
    linewidth = 0.3,
    color     = "white"
  ) +
  
  # ボックスプロット（中央値・IQR）
  geom_boxplot(
    width         = 0.18,
    outlier.shape = NA,
    alpha         = 0.85,
    linewidth     = 0.45,
    color         = "grey20"
  ) +
  
  # アノテーション（p_BH / delta / n）
  geom_text(
    data        = df_annot,
    aes(x       = 1.5, y = y_pos, label = annot),
    inherit.aes = FALSE,
    size        = 2.8,
    hjust       = 0.5,
    vjust       = 0,
    color       = "grey25",
    lineheight  = 1.35,
    family      = "sans"
  ) +
  
  scale_fill_manual(
    values = c("TP53 Wild-type" = COL_WT,
               "TP53 Mutant"   = COL_MUT)
  ) +
  
  scale_y_continuous(
    name   = "LAG3 expression [log2(TPM+1)]",
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  
  facet_wrap(~ grade_label, nrow = 1) +
  
  labs(
    x       = NULL,
    title   = NULL,
    caption = paste0(
      "Two-sided Wilcoxon rank-sum test with BH correction. ",
      "Boxes show median and IQR; whiskers extend to 1.5x IQR."
    )
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.background = element_blank(),
    strip.text       = element_text(size = 11, face = "bold"),
    axis.text.x      = element_text(size = 10, color = "grey20"),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 10),
    plot.caption     = element_text(size = 7.5, color = "grey50", hjust = 0),
    panel.spacing    = unit(1.5, "lines")
  )

# ---- 5. 出力 ----
ggsave(OUT_PDF, plot = p, width = 6.5, height = 5.5, device = "pdf")
ggsave(OUT_PNG, plot = p, width = 6.5, height = 5.5, dpi = 450)

cat("Done.\n")
cat("PDF:", OUT_PDF, "\n")
cat("PNG:", OUT_PNG, "\n")
