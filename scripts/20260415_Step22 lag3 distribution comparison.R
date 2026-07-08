# =============================================================================
# step22_lag3_distribution_comparison_v3.R
# GBM/Glioma TP53xLAG3 解析 - Step 22: TCGA vs CPTAC/HCMI LAG3分布比較
# v3: mutant SDをログ・CSV・図に追加（v2からの変更点）
#
# 【v2からの変更点】
#   - mutant SDをログに出力（WT SDと同じ形式）
#   - step22_descriptive_stats.csvはv2と同一（元から全groupのSDを含む）
#   - step22_source_regression.csvにmut_sd列を追加
#   - Panel Aのfacet内注記をWT SDのみ → WT SD + Mut SDの2行表示に変更
#   - 出力ファイルを別名保存:
#       figS8_lag3_distribution_v3.pdf
#       figS8_lag3_distribution_v3_600dpi.png
#       step22_source_regression_v3.csv（mut_sd列追加版）
#
# 出力先: 22_lag3_distribution/
# 作成日: 2026-03-06 → v3修正: 2026-04-15
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
OUT_DIR    <- file.path(RESULT_DIR, "22_lag3_distribution")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

COL_MUT   <- "#E64B35"
COL_WT    <- "#AAAAAA"
COL_TCGA  <- "#3C5488"
COL_CPTAC <- "#E07B54"

GDC_PATH <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step22_v3_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 22 v3: 2パネル構成 開始 ===")
log_msg("  [v3変更点] mutant SDをログ・CSV・図に追加")

# =============================================================================
# 2. データ読み込み
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

gdc_g4 <- gdc_base %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = factor(
      ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT"),
      levels = c("TP53 WT", "TP53 Mut")
    ),
    source_label = factor(source,
                          levels = c("TCGA", "CPTAC_HCMI"),
                          labels = c("TCGA", "CPTAC/HCMI"))
  )

log_msg(sprintf("GDC Grade4: n=%d (TCGA=%d, CPTAC/HCMI=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$source == "TCGA"),
                sum(gdc_g4$source == "CPTAC_HCMI")))

# =============================================================================
# 3. 記述統計（source x tp53別）
# =============================================================================

log_msg("--- 記述統計 ---")

desc_stats <- gdc_g4 %>%
  group_by(source_label, tp53_label) %>%
  summarise(
    n         = n(),
    mean      = round(mean(LAG3_log2tpm,            na.rm = TRUE), 4),
    sd        = round(sd(LAG3_log2tpm,              na.rm = TRUE), 4),
    median    = round(median(LAG3_log2tpm,          na.rm = TRUE), 4),
    Q1        = round(quantile(LAG3_log2tpm, 0.25,  na.rm = TRUE), 4),
    Q3        = round(quantile(LAG3_log2tpm, 0.75,  na.rm = TRUE), 4),
    IQR       = round(IQR(LAG3_log2tpm,             na.rm = TRUE), 4),
    range_min = round(min(LAG3_log2tpm,             na.rm = TRUE), 4),
    range_max = round(max(LAG3_log2tpm,             na.rm = TRUE), 4),
    .groups   = "drop"
  ) %>%
  mutate(range_width = round(range_max - range_min, 4))

write_csv(desc_stats, file.path(OUT_DIR, "step22_descriptive_stats.csv"))
log_msg("保存: step22_descriptive_stats.csv")

# ---------------------------------------------------------------------------
# 【v3追加】WT群・Mut群それぞれのSDをログ出力
# ---------------------------------------------------------------------------
wt_sd <- desc_stats %>%
  filter(tp53_label == "TP53 WT") %>%
  select(source_label, sd) %>%
  mutate(sd_label = sprintf("WT SD = %.3f", sd))

mut_sd <- desc_stats %>%
  filter(tp53_label == "TP53 Mut") %>%
  select(source_label, sd) %>%
  rename(mut_sd = sd) %>%
  mutate(mut_sd_label = sprintf("Mut SD = %.3f", mut_sd))

log_msg("--- WT群SD ---")
for (i in seq_len(nrow(wt_sd))) {
  log_msg(sprintf("  [%s] WT SD = %.4f", wt_sd$source_label[i], wt_sd$sd[i]))
}

log_msg("--- Mut群SD ---")  # 【v3追加】
for (i in seq_len(nrow(mut_sd))) {
  log_msg(sprintf("  [%s] Mut SD = %.4f", mut_sd$source_label[i], mut_sd$mut_sd[i]))
}

# =============================================================================
# 4. source別回帰（IDH調整あり）
# =============================================================================

log_msg("--- source別回帰: LAG3 ~ tp53_bin + idh_bin ---")

run_source_reg <- function(df, source_name) {
  df_use <- df %>%
    filter(source_label == source_name) %>%
    select(LAG3_log2tpm, tp53_bin, idh_bin) %>%
    na.omit()
  
  fit <- lm(LAG3_log2tpm ~ tp53_bin + idh_bin, data = df_use)
  ci  <- as.double(confint(fit, "tp53_bin"))
  smr <- coef(summary(fit))
  
  beta <- smr["tp53_bin", "Estimate"]
  se   <- smr["tp53_bin", "Std. Error"]
  pval <- smr["tp53_bin", "Pr(>|t|)"]
  r2   <- round(summary(fit)$r.squared, 4)
  n_use <- nrow(df_use)
  
  log_msg(sprintf("  [%s] n=%d | beta=%+.4f [%+.4f, %+.4f] | p=%.4f | SE=%.4f | R2=%.4f",
                  source_name, n_use, beta, ci[1], ci[2], pval, se, r2))
  
  tibble(
    source    = as.character(source_name),
    n         = n_use,
    beta      = round(beta,  4),
    ci_lower  = round(ci[1], 4),
    ci_upper  = round(ci[2], 4),
    se        = round(se,    4),
    p_value   = pval,
    r_squared = r2
  )
}

reg_results <- do.call(rbind, list(
  run_source_reg(gdc_g4, "TCGA"),
  run_source_reg(gdc_g4, "CPTAC/HCMI")
)) %>%
  left_join(wt_sd  %>% rename(source = source_label), by = "source") %>%
  left_join(mut_sd %>% rename(source = source_label), by = "source") %>%  # 【v3追加】
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    source = factor(source, levels = c("CPTAC/HCMI", "TCGA"))
  )

# 【v3変更】出力CSVをv3別名で保存（mut_sd列追加）
write_csv(reg_results %>% mutate(source = as.character(source)),
          file.path(OUT_DIR, "step22_source_regression_v3.csv"))
log_msg("保存: step22_source_regression_v3.csv")

# =============================================================================
# 5. Panel A: density plot（WT SD + Mut SD注記）
# =============================================================================

log_msg("--- Panel A: density plot 作成 ---")

# 【v3変更】注記をWT SD + Mut SDの2行に変更
sd_annot <- wt_sd %>%
  rename(source_label = source_label) %>%
  left_join(mut_sd %>% rename(source_label = source_label), by = "source_label") %>%
  mutate(sd_label2 = paste0(sd_label, "\n", mut_sd_label))

fig_a <- ggplot(gdc_g4,
                aes(x = LAG3_log2tpm, color = tp53_label, fill = tp53_label)) +
  geom_density(alpha = 0.2, linewidth = 0.7) +
  geom_text(
    data = sd_annot,
    aes(x = Inf, y = Inf, label = sd_label2),  # 【v3変更】2行ラベル
    hjust = 1.1, vjust = 1.5,
    size = 3.0, color = "#555555", fontface = "italic",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ source_label, nrow = 1) +
  scale_color_manual(
    values = c("TP53 WT" = COL_WT, "TP53 Mut" = COL_MUT),
    name   = NULL
  ) +
  scale_fill_manual(
    values = c("TP53 WT" = COL_WT, "TP53 Mut" = COL_MUT),
    name   = NULL
  ) +
  labs(
    x        = expression(italic(LAG3) ~ "expression [log"[2] * "(TPM+1)]"),
    y        = "Density",
    title    = expression(bold("A") ~ ~ italic(LAG3) ~
                            "distribution by source and" ~ italic(TP53) ~ "status"),
    subtitle = "GDC Grade 4. CPTAC/HCMI shows narrower range (smaller wild-type SD)."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 10),
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "#F0F0F0", color = NA),
    strip.text        = element_text(size = 11, face = "bold"),
    plot.title        = element_text(size = 12),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    axis.title        = element_text(size = 10)
  )

# =============================================================================
# 6. Panel B: Forest風回帰係数プロット（WT SD注記あり）
# =============================================================================

log_msg("--- Panel B: Forest風プロット 作成 ---")

fig_b <- ggplot(reg_results,
                aes(x = beta, y = source, color = source, fill = source)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.18, linewidth = 1.0) +
  geom_point(shape = 21, size = 5, stroke = 1.2) +
  geom_text(
    aes(x = ci_upper,
        label = sprintf("b = %+.3f, p = %.3f (%s)", beta, p_value, sig_label)),
    hjust = -0.08, size = 3.4, color = "grey30"
  ) +
  geom_text(
    aes(x = ci_lower, label = sd_label),
    hjust = 1.1, size = 3.2, color = "grey45"
  ) +
  scale_color_manual(
    values = c("TCGA" = COL_TCGA, "CPTAC/HCMI" = COL_CPTAC),
    guide  = "none"
  ) +
  scale_fill_manual(
    values = c("TCGA" = COL_TCGA, "CPTAC/HCMI" = COL_CPTAC),
    guide  = "none"
  ) +
  scale_x_continuous(
    name   = expression(italic(TP53) ~
                          "coefficient (\u03b2) with 95% CI   [LAG3 ~ tp53 + IDH]"),
    breaks = seq(-0.1, 0.6, by = 0.1)
  ) +
  coord_cartesian(xlim = c(
    min(reg_results$ci_lower) - 0.22,
    max(reg_results$ci_upper) + 0.30
  )) +
  labs(
    y        = NULL,
    title    = expression(bold("B") ~ ~ italic(TP53) ~
                            "effect size and precision by source"),
    subtitle = paste0(
      "Narrower LAG3 distribution in CPTAC/HCMI (smaller wild-type SD) ",
      "reduces precision (wider CI), not effect direction."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.title.x       = element_text(size = 10),
    plot.title         = element_text(size = 12),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    plot.margin        = margin(10, 130, 10, 10)
  )

# =============================================================================
# 7. 結合（A上・B下）
# =============================================================================

log_msg("--- 結合図 作成 ---")

fig_combined <- fig_a / fig_b +
  plot_layout(heights = c(1.1, 1.0)) +
  plot_annotation(
    caption = paste0(
      "GDC Grade 4 (TCGA n = 245, CPTAC/HCMI n = 197). ",
      "Regression: LAG3 ~ tp53 + IDH (source-stratified). ",
      "Both sources show positive TP53 effect direction; ",
      "narrower LAG3 range in CPTAC/HCMI reduces statistical precision (larger SE/CI). ",
      "Wild-type and mutant SD shown in Panel A."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0)
    )
  )

# =============================================================================
# 8. 出力（別名・600dpi PNG + PDF）
# =============================================================================

log_msg("--- ファイル出力 ---")

pdf_path <- file.path(OUT_DIR, "figS8_lag3_distribution_v3.pdf")
pdf(pdf_path, width = 10, height = 10)
print(fig_combined)
dev.off()
log_msg(sprintf("保存: %s", pdf_path))

png_path <- file.path(OUT_DIR, "figS8_lag3_distribution_v3_600dpi.png")
agg_png(png_path, width = 10, height = 10, units = "in", res = 600)
print(fig_combined)
dev.off()
log_msg(sprintf("保存: %s", png_path))

# =============================================================================
# 9. 完了
# =============================================================================

log_msg("=== Step 22 v3: 完了 ===")

cat("\n============================\n")
cat("Step 22 v3 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s\n", file.path(OUT_DIR, "step22_source_regression_v3.csv")))
cat(sprintf("  %s\n", pdf_path))
cat(sprintf("  %s\n", png_path))
cat("============================\n")
