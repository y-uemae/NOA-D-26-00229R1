# =============================================================================
# step23_standardized_effect.R
# GBM/Glioma TP53xLAG3 解析 - Step 23: 標準化効果量比較（Glass's delta / Cohen's d）
#
# 目的:
#   source間でLAG3のSD/分散が異なる中で、標準化効果量を用いて
#   「TP53効果が両sourceで一貫している」を示す。
#   Glass's delta（WT群SDで標準化）がsource間比較に最も適切。
#
# 解析内容:
#   Glass's delta = (mean_mut - mean_wt) / sd_wt
#   Cohen's d     = (mean_mut - mean_wt) / pooled_SD
#   Cliff's delta = step09_lag3_summary.csvから読み込み（再計算なし）
#   95%CI         = bootstrap（B=2000）
#
# 出力:
#   figS23a_effect_forest.pdf/png  : Glass's d / Cohen's d source別Forest風
#   figS23b_effect_bar.pdf/png     : 3指標横並び棒グラフ
#   figS23_combined.pdf/png        : A+B上下配置
#   step23_standardized_effects.csv: 全指標まとめ
#   step23_log.txt
#
# 出力先: 23_standardized_effect/
#
# 作成日: 2026-02-24
# 修正日: 2026-02-26  Panel A テキストラベル表示修正
#   - scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) 追加
#   - coord_cartesian(clip = "off") に変更
#   - plot.margin 右余白を縮小（scale_x で制御）
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
OUT_DIR    <- file.path(RESULT_DIR, "23_standardized_effect")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# bootstrap設定
set.seed(42)
B_BOOT <- 2000

# 色設定（引継書10章に準拠）
COL_TCGA  <- "#3C5488"
COL_CPTAC <- "#E07B54"
COL_ALL   <- "#8491B4"

# 入力ファイル
GDC_PATH    <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
STEP09_PATH <- file.path(RESULT_DIR, "09_statistics/step09_lag3_summary.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step23_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 23: 標準化効果量比較 開始 ===")
log_msg(sprintf("Bootstrap: B=%d, seed=42", B_BOOT))

# =============================================================================
# 2. データ読み込み
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

step09 <- read_csv(STEP09_PATH, show_col_types = FALSE)

log_msg(sprintf("GDC: %d行", nrow(gdc_base)))
log_msg(sprintf("Step09 summary: %d行", nrow(step09)))

# Grade4解析データ
gdc_g4 <- gdc_base %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    source_label = factor(source,
                          levels = c("TCGA", "CPTAC_HCMI"),
                          labels = c("TCGA", "CPTAC/HCMI"))
  )

log_msg(sprintf("GDC Grade4: n=%d (TCGA=%d, CPTAC/HCMI=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$source == "TCGA"),
                sum(gdc_g4$source == "CPTAC_HCMI")))

# =============================================================================
# 3. Bootstrap CI計算ヘルパー
# =============================================================================

boot_ci <- function(vec_mut, vec_wt, stat_fn, B = B_BOOT) {
  n_mut <- length(vec_mut)
  n_wt  <- length(vec_wt)
  boot_vals <- replicate(B, {
    s_mut <- sample(vec_mut, n_mut, replace = TRUE)
    s_wt  <- sample(vec_wt,  n_wt,  replace = TRUE)
    stat_fn(s_mut, s_wt)
  })
  quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE)
}

# Glass's delta関数
glass_delta_fn <- function(mut, wt) {
  (mean(mut, na.rm = TRUE) - mean(wt, na.rm = TRUE)) / sd(wt, na.rm = TRUE)
}

# Cohen's d関数（pooled SD）
cohen_d_fn <- function(mut, wt) {
  n1 <- length(mut); n2 <- length(wt)
  pooled_sd <- sqrt(((n1 - 1) * var(mut, na.rm = TRUE) +
                       (n2 - 1) * var(wt,  na.rm = TRUE)) / (n1 + n2 - 2))
  (mean(mut, na.rm = TRUE) - mean(wt, na.rm = TRUE)) / pooled_sd
}

# =============================================================================
# 4. 標準化効果量の計算（All / TCGA / CPTAC/HCMI）
# =============================================================================

log_msg("--- 標準化効果量計算 ---")

calc_effects <- function(df, subset_name) {
  vec_mut <- df %>% filter(tp53_status == "mutant")   %>% pull(LAG3_log2tpm)
  vec_wt  <- df %>% filter(tp53_status == "wildtype") %>% pull(LAG3_log2tpm)
  
  # 点推定
  gd <- glass_delta_fn(vec_mut, vec_wt)
  cd <- cohen_d_fn(vec_mut, vec_wt)
  
  # Bootstrap CI
  log_msg(sprintf("  [%s] bootstrap開始 (n_mut=%d, n_wt=%d)...",
                  subset_name, length(vec_mut), length(vec_wt)))
  ci_gd <- boot_ci(vec_mut, vec_wt, glass_delta_fn)
  ci_cd <- boot_ci(vec_mut, vec_wt, cohen_d_fn)
  
  log_msg(sprintf("  [%s] Glass's delta=%.4f [%.4f, %.4f]",
                  subset_name, gd, ci_gd[1], ci_gd[2]))
  log_msg(sprintf("  [%s] Cohen's d     =%.4f [%.4f, %.4f]",
                  subset_name, cd, ci_cd[1], ci_cd[2]))
  
  tibble(
    subset      = subset_name,
    n_mut       = length(vec_mut),
    n_wt        = length(vec_wt),
    glass_delta = round(gd,       4),
    glass_ci_lo = round(ci_gd[1], 4),
    glass_ci_hi = round(ci_gd[2], 4),
    cohen_d     = round(cd,       4),
    cohen_ci_lo = round(ci_cd[1], 4),
    cohen_ci_hi = round(ci_cd[2], 4)
  )
}

results_all   <- calc_effects(gdc_g4,                                    "All")
results_tcga  <- calc_effects(gdc_g4 %>% filter(source == "TCGA"),       "TCGA")
results_cptac <- calc_effects(gdc_g4 %>% filter(source == "CPTAC_HCMI"), "CPTAC/HCMI")

effect_results <- do.call(rbind, list(results_all, results_tcga, results_cptac))

# Step09からCliff's deltaを結合
cliffs_map <- step09 %>%
  filter(gene == "LAG3",
         subset %in% c("All", "TCGA", "CPTAC_HCMI")) %>%
  mutate(
    subset_key = dplyr::recode(subset,
                               "All"        = "All",
                               "TCGA"       = "TCGA",
                               "CPTAC_HCMI" = "CPTAC/HCMI")
  ) %>%
  select(subset_key, cliffs_delta)

effect_results <- effect_results %>%
  left_join(cliffs_map, by = c("subset" = "subset_key"))

log_msg("--- 全指標まとめ ---")
for (i in seq_len(nrow(effect_results))) {
  r <- effect_results[i, ]
  log_msg(sprintf("  [%s] Glass_d=%.4f | Cohen_d=%.4f | Cliff_d=%.4f",
                  r$subset, r$glass_delta, r$cohen_d, r$cliffs_delta))
}

write_csv(effect_results, file.path(OUT_DIR, "step23_standardized_effects.csv"))
log_msg("保存: step23_standardized_effects.csv")

# =============================================================================
# 5. Panel A: Glass's delta + Cohen's d Forest風
# =============================================================================

log_msg("--- Panel A: Forest風 作成 ---")

# long形式（Glass's delta / Cohen's d）
forest_data <- effect_results %>%
  mutate(
    subset = factor(subset, levels = c("CPTAC/HCMI", "TCGA", "All")),
    source_color = dplyr::recode(as.character(subset),
                                 "All"        = COL_ALL,
                                 "TCGA"       = COL_TCGA,
                                 "CPTAC/HCMI" = COL_CPTAC)
  ) %>%
  pivot_longer(
    cols      = c(glass_delta, cohen_d),
    names_to  = "metric",
    values_to = "estimate"
  ) %>%
  mutate(
    ci_lo = ifelse(metric == "glass_delta", glass_ci_lo, cohen_ci_lo),
    ci_hi = ifelse(metric == "glass_delta", glass_ci_hi, cohen_ci_hi),
    metric_label = factor(
      dplyr::recode(metric,
                    glass_delta = "Glass's delta\n(mean diff / SD_WT)",
                    cohen_d     = "Cohen's d\n(mean diff / pooled SD)"),
      levels = c("Glass's delta\n(mean diff / SD_WT)",
                 "Cohen's d\n(mean diff / pooled SD)")
    )
  )

# 色ベクター
col_vals <- c("All" = COL_ALL, "TCGA" = COL_TCGA, "CPTAC/HCMI" = COL_CPTAC)

fig_a <- ggplot(forest_data,
                aes(x = estimate, y = subset,
                    color = subset, fill = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                orientation = "y", width = 0.2, linewidth = 0.8) +
  geom_point(shape = 21, size = 4.5, stroke = 1.1) +
  geom_text(aes(x = ci_hi,
                label = sprintf("%.3f [%.3f, %.3f]", estimate, ci_lo, ci_hi)),
            hjust = -0.12, size = 2.9, color = "grey30") +
  facet_wrap(~ metric_label, nrow = 1, scales = "free_x") +
  scale_color_manual(values = col_vals, guide = "none") +
  scale_fill_manual(values  = col_vals, guide = "none") +
  # ★修正1: x軸右側に55%の余白を確保してテキストラベル用の描画領域を作る
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.55))) +
  labs(
    x        = "Effect size (95% bootstrap CI)",
    y        = NULL,
    title    = "A  Standardized effect sizes by source",
    subtitle = paste0(
      "GDC Grade4. Bootstrap CI (B=", B_BOOT, ", seed=42). ",
      "Both sources show consistent positive direction."
    )
  ) +
  # ★修正2: clip = "off" でパネル境界外のテキストも描画
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "#F0F0F0", color = NA),
    strip.text         = element_text(size = 9, face = "bold"),
    axis.text.y        = element_text(size = 10, face = "bold"),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40"),
    # ★修正3: 右余白はscale_xで制御するためmarginは縮小
    plot.margin        = margin(10, 20, 10, 10)
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS23a_effect_forest.pdf"), width = 11, height = 4.5)
print(fig_a)
dev.off()
agg_png(file.path(OUT_DIR, "figS23a_effect_forest_450dpi.png"),
        width = 11, height = 4.5, units = "in", res = 450)
print(fig_a)
dev.off()
log_msg("保存: figS23a_effect_forest.pdf/png")

# =============================================================================
# 6. Panel B: 3指標横並び棒グラフ
# =============================================================================

log_msg("--- Panel B: 3指標棒グラフ 作成 ---")

bar_data <- effect_results %>%
  mutate(subset = factor(subset, levels = c("All", "TCGA", "CPTAC/HCMI"))) %>%
  select(subset, glass_delta, cohen_d, cliffs_delta) %>%
  pivot_longer(
    cols      = c(glass_delta, cohen_d, cliffs_delta),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = factor(
      dplyr::recode(metric,
                    glass_delta  = "Glass's delta",
                    cohen_d      = "Cohen's d",
                    cliffs_delta = "Cliff's delta"),
      levels = c("Glass's delta", "Cohen's d", "Cliff's delta")
    )
  )

fig_b <- ggplot(bar_data,
                aes(x = subset, y = value, fill = subset)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", value)),
            vjust = -0.4, size = 3.0, color = "grey30") +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.4) +
  facet_wrap(~ metric_label, nrow = 1) +
  scale_fill_manual(values = col_vals, guide = "none") +
  scale_y_continuous(
    name   = "Effect size",
    limits = c(0, max(bar_data$value, na.rm = TRUE) * 1.25),
    breaks = seq(0, 1.0, by = 0.1)
  ) +
  labs(
    x        = NULL,
    title    = "B  Three effect size metrics across sources",
    subtitle = paste0(
      "All metrics positive across sources. ",
      "Glass's delta accounts for source-specific SD differences."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text       = element_text(size = 9.5, face = "bold"),
    axis.text.x      = element_text(size = 9, face = "bold"),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8.5, color = "grey40")
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS23b_effect_bar.pdf"), width = 9, height = 4.5)
print(fig_b)
dev.off()
agg_png(file.path(OUT_DIR, "figS23b_effect_bar_450dpi.png"),
        width = 9, height = 4.5, units = "in", res = 450)
print(fig_b)
dev.off()
log_msg("保存: figS23b_effect_bar.pdf/png")

# =============================================================================
# 7. Combined: A + B 上下配置
# =============================================================================

log_msg("--- Combined 上下配置 作成 ---")

fig_combined <- wrap_plots(fig_a) / wrap_plots(fig_b) +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title   = "TP53 effect on LAG3 is consistent across sources after standardization",
    caption = paste0(
      "GDC Grade4 (All n=442, TCGA n=245, CPTAC/HCMI n=197). ",
      "Glass's delta = (mean_Mut - mean_WT) / SD_WT; ",
      "Cohen's d = (mean_Mut - mean_WT) / pooled SD; ",
      "Cliff's delta from Step09 (Mann-Whitney). ",
      "Bootstrap 95% CI (B=", B_BOOT, ") for Glass's delta and Cohen's d."
    ),
    theme = theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

pdf(file.path(OUT_DIR, "figS23_combined.pdf"), width = 11, height = 10)
print(fig_combined)
dev.off()
agg_png(file.path(OUT_DIR, "figS23_combined_450dpi.png"),
        width = 11, height = 10, units = "in", res = 450)
print(fig_combined)
dev.off()
log_msg("保存: figS23_combined.pdf/png")

# =============================================================================
# 8. 完了
# =============================================================================

log_msg("=== Step 23: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 23 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/figS23a_effect_forest.pdf/png\n",   OUT_DIR))
cat(sprintf("  %s/figS23b_effect_bar.pdf/png\n",      OUT_DIR))
cat(sprintf("  %s/figS23_combined.pdf/png\n",         OUT_DIR))
cat(sprintf("  %s/step23_standardized_effects.csv\n", OUT_DIR))
cat(sprintf("  %s/step23_log.txt\n",                  OUT_DIR))
cat("============================\n")
