# ============================================================
# step10_glass.R  GLASS Grade4 検証解析（構文修正版）
# ============================================================

library(tidyverse)

BASE_DIR  <- here::here("results", "TP53", "20260221")
GLASS_DIR <- file.path(BASE_DIR, "05c_glass")
OUT_DIR   <- file.path(BASE_DIR, "10_glass_validation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 共通関数
# ============================================================
compare_two_groups <- function(df, expr_col,
                               group_col = "tp53_status",
                               mut_label = "Mut",
                               wt_label  = "WT") {
  df2 <- df %>%
    filter(.data[[group_col]] %in% c(mut_label, wt_label)) %>%
    select(expr = all_of(expr_col), group = all_of(group_col)) %>%
    drop_na()
  
  mut_vals <- df2$expr[df2$group == mut_label]
  wt_vals  <- df2$expr[df2$group == wt_label]
  n_mut <- length(mut_vals); n_wt <- length(wt_vals)
  
  if (n_mut < 3 || n_wt < 3) {
    return(tibble(
      n_mut=n_mut, n_wt=n_wt,
      median_mut=NA_real_, median_wt=NA_real_, median_diff=NA_real_,
      mean_mut=NA_real_, mean_wt=NA_real_,
      sd_mut=NA_real_, sd_wt=NA_real_,
      HL_estimate=NA_real_, p_wilcox=NA_real_,
      cliffs_delta=NA_real_, effect_size_label=NA_character_
    ))
  }
  
  wtest <- wilcox.test(mut_vals, wt_vals, exact = FALSE, conf.int = TRUE)
  dominance <- outer(mut_vals, wt_vals, FUN = function(a, b) sign(a - b))
  delta <- sum(dominance) / (n_mut * n_wt)
  abs_d <- abs(delta)
  label <- case_when(
    abs_d < 0.147 ~ "negligible",
    abs_d < 0.330 ~ "small",
    abs_d < 0.474 ~ "medium",
    TRUE          ~ "large"
  )
  
  tibble(
    n_mut=n_mut, n_wt=n_wt,
    median_mut=median(mut_vals), median_wt=median(wt_vals),
    median_diff=median(mut_vals) - median(wt_vals),
    mean_mut=mean(mut_vals), mean_wt=mean(wt_vals),
    sd_mut=sd(mut_vals), sd_wt=sd(wt_vals),
    HL_estimate=as.numeric(wtest$estimate),
    p_wilcox=wtest$p.value,
    cliffs_delta=delta,
    effect_size_label=label
  )
}

# ============================================================
# データ読込
# ============================================================
glass_main <- read_csv(
  file.path(GLASS_DIR, "glass_final_cohort_wxs_notcga.csv"),
  show_col_types = FALSE
)
glass_sens <- read_csv(
  file.path(GLASS_DIR, "glass_final_cohort_wxs_wgs_notcga.csv"),
  show_col_types = FALSE
)

cat("主解析 n=", nrow(glass_main), "| 感度 n=", nrow(glass_sens), "\n")
cat("tp53_status（主解析）:\n")
print(table(glass_main$tp53_status, useNA = "always"))

# ============================================================
# 1. 主解析: WXS-only, TCGA除外（n=79）
# ============================================================
cat("\n=== 主解析: WXS-only TCGA除外 ===\n")
res_main <- compare_two_groups(glass_main, "LAG3_log2tpm") %>%
  mutate(analysis="GLASS_main", cohort="GLASS",
         subset="WXS_noTCGA", gene="LAG3", p_BH=p_wilcox)
print(res_main %>% select(n_mut, n_wt, median_mut, median_wt,
                          median_diff, p_wilcox, cliffs_delta, effect_size_label))

# ============================================================
# 2. 感度解析: WXS+WGS, TCGA除外（n=95）
# ============================================================
cat("\n=== 感度解析: WXS+WGS TCGA除外 ===\n")
res_sens <- compare_two_groups(glass_sens, "LAG3_log2tpm") %>%
  mutate(analysis="GLASS_sensitivity", cohort="GLASS",
         subset="WXS_WGS_noTCGA", gene="LAG3", p_BH=p_wilcox)
print(res_sens %>% select(n_mut, n_wt, median_mut, median_wt,
                          median_diff, p_wilcox, cliffs_delta, effect_size_label))

# ============================================================
# 3. seq_type別（感度コホート内）
# ============================================================
cat("\n=== seq_type別（感度コホート内）===\n")
if ("seq_type" %in% names(glass_sens)) {
  seq_types <- unique(na.omit(glass_sens$seq_type))
  res_seqtype <- map_dfr(seq_types, function(st) {
    compare_two_groups(
      glass_sens %>% filter(seq_type == st),
      "LAG3_log2tpm"
    ) %>%
      mutate(analysis = "GLASS_by_seqtype", subset = st, gene = "LAG3")
  }) %>%
    mutate(p_BH = p.adjust(p_wilcox, method = "BH"))
  print(res_seqtype %>% select(subset, n_mut, n_wt, median_diff,
                               p_wilcox, p_BH, cliffs_delta, effect_size_label))
} else {
  cat("seq_type列なし → スキップ\n")
  res_seqtype <- tibble()
}

# ============================================================
# 4. 簡易回帰（主解析コホート: TP53 only）
# ============================================================
cat("\n=== GLASS 簡易回帰（主解析コホート）===\n")
glass_reg_main <- glass_main %>%
  mutate(tp53_bin = if_else(tp53_status == "Mut", 1L, 0L)) %>%
  filter(!is.na(tp53_bin), !is.na(LAG3_log2tpm))

cat("回帰用 n=", nrow(glass_reg_main),
    "| Mut:", sum(glass_reg_main$tp53_bin),
    "| WT:", sum(1L - glass_reg_main$tp53_bin), "\n")

fit_main <- lm(LAG3_log2tpm ~ tp53_bin, data = glass_reg_main)
sm <- summary(fit_main)
cat("全係数:\n"); print(sm$coefficients)
cat("95%CI:\n"); print(confint(fit_main))
cat(sprintf("R²=%.3f\n", sm$r.squared))

# ============================================================
# 5. 簡易回帰（感度コホート: TP53 + seq_type調整）
# ============================================================
cat("\n=== GLASS 簡易回帰（感度コホート: TP53 + seq_type）===\n")
glass_reg_sens <- glass_sens %>%
  mutate(
    tp53_bin    = if_else(tp53_status == "Mut", 1L, 0L),
    seqtype_bin = if_else(seq_type == "WXS", 1L, 0L)
  ) %>%
  filter(!is.na(tp53_bin), !is.na(seqtype_bin), !is.na(LAG3_log2tpm))

cat("回帰用 n=", nrow(glass_reg_sens),
    "| Mut:", sum(glass_reg_sens$tp53_bin),
    "| WT:", sum(1L - glass_reg_sens$tp53_bin), "\n")

fit_sens <- lm(LAG3_log2tpm ~ tp53_bin + seqtype_bin, data = glass_reg_sens)
ss <- summary(fit_sens)
cat("全係数:\n"); print(ss$coefficients)
cat("95%CI:\n"); print(confint(fit_sens))
cat(sprintf("R²=%.3f\n", ss$r.squared))

# ============================================================
# 6. GDC vs GLASS 効果量比較表
# ============================================================
cat("\n=== GDC vs GLASS 効果量比較 ===\n")

forest_df <- bind_rows(
  tibble(cohort="GDC", subset="Grade4_all",
         n_mut=147L, n_wt=295L, median_diff=0.288, cliffs_delta=0.289, p_wilcox=7.13e-7),
  tibble(cohort="GDC", subset="Grade4_TCGA",
         n_mut=79L,  n_wt=166L, median_diff=0.338, cliffs_delta=0.326, p_wilcox=3.78e-5),
  tibble(cohort="GDC", subset="Grade4_CPTAC_HCMI",
         n_mut=68L,  n_wt=129L, median_diff=0.193, cliffs_delta=0.234, p_wilcox=7.08e-3),
  tibble(cohort="GDC", subset="Grade4_IDH_WT",
         n_mut=120L, n_wt=290L, median_diff=0.254, cliffs_delta=0.251, p_wilcox=6.49e-5),
  tibble(cohort="GLASS", subset="WXS_noTCGA",
         n_mut=res_main$n_mut, n_wt=res_main$n_wt,
         median_diff=res_main$median_diff,
         cliffs_delta=res_main$cliffs_delta,
         p_wilcox=res_main$p_wilcox),
  tibble(cohort="GLASS", subset="WXS_WGS_noTCGA",
         n_mut=res_sens$n_mut, n_wt=res_sens$n_wt,
         median_diff=res_sens$median_diff,
         cliffs_delta=res_sens$cliffs_delta,
         p_wilcox=res_sens$p_wilcox)
)

print(forest_df)

# ============================================================
# 7. 保存
# ============================================================
res_all <- bind_rows(res_main, res_sens) %>%
  select(analysis, cohort, subset, gene, n_mut, n_wt,
         median_mut, median_wt, median_diff,
         HL_estimate, p_wilcox, p_BH, cliffs_delta, effect_size_label)

write_csv(res_all,   file.path(OUT_DIR, "step10_glass_lag3_results.csv"))
write_csv(forest_df, file.path(OUT_DIR, "step10_forest_summary.csv"))
cat("\n✅ 保存:\n")
cat(" ", file.path(OUT_DIR, "step10_glass_lag3_results.csv"), "\n")
cat(" ", file.path(OUT_DIR, "step10_forest_summary.csv"), "\n")

# ============================================================
# 8. 最終サマリ
# ============================================================
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("【Step10 完了サマリ】\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("\n■ GLASS 主解析（WXS-only, TCGA除外）\n")
with(res_main, cat(sprintf(
  "  Mut n=%d, WT n=%d | median_diff=%.3f | p=%.4e | δ=%.3f (%s)\n",
  n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label
)))

cat("\n■ GLASS 感度解析（WXS+WGS, TCGA除外）\n")
with(res_sens, cat(sprintf(
  "  Mut n=%d, WT n=%d | median_diff=%.3f | p=%.4e | δ=%.3f (%s)\n",
  n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label
)))

cat("\n■ δ一覧（GDC→GLASS 方向確認）\n")
print(forest_df %>% select(cohort, subset, n_mut, n_wt,
                           median_diff, cliffs_delta, p_wilcox))
