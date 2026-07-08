# =============================================================================
# 20260224_step12c_subgroup_analysis.R
# Step 12c: TP53 4群サブグループ解析
#
# 解析セット（事前固定）:
#   A) GDC Grade4 全体（n=442）
#   B) GDC Grade4 TCGA のみ（n=245）
#   C) GDC Grade4 CPTAC_HCMI のみ（n=197）
#   D) GDC Grade4 IDH-WT限定（n≈410）
#
# 解析内容:
#   1. Kruskal-Wallis（4群全体差）
#   2. 事前比較2本: Hotspot vs WT、Truncating vs WT（BH補正）
#      参考: Other_missense vs WT（同BH補正）
#   3. 効果量: Cliff's δ + median差（各MutサブグループvsWT）
#   4. 回帰: LAG3_log2tpm ~ tp53_class4 + source + idh_bin（全体のみ）
#
# 入力:
#   12_subgroup/step12b_subgroup_classified.csv
# 出力:
#   12_subgroup/step12c_kw_results.csv
#   12_subgroup/step12c_pairwise_results.csv
#   12_subgroup/step12c_effectsize.csv
#   12_subgroup/step12c_regression.csv
#   12_subgroup/step12c_summary.txt
# =============================================================================

library(tidyverse)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR <- here::here("results", "TP53", "20260221")
IN_CSV   <- file.path(BASE_DIR, "12_subgroup/step12b_subgroup_classified.csv")
OUT_DIR  <- file.path(BASE_DIR, "12_subgroup")

# ── 1. データ読み込み・整形 ───────────────────────────────────────────────────
df <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(
    tp53_class4 = factor(tp53_class4,
                         levels = c("WT", "Hotspot", "Truncating", "Other_missense")),
    idh_bin     = as.integer(idh_status == "mutant"),
    source_bin  = as.integer(source == "CPTAC_HCMI")
  )

cat("=== 解析データ確認 ===\n")
df %>% count(tp53_class4, source) %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric)))) %>%
  print()

# 解析サブセット定義
subsets <- list(
  "Grade4_all"      = df,
  "Grade4_TCGA"     = filter(df, source == "TCGA"),
  "Grade4_CPTAC"    = filter(df, source == "CPTAC_HCMI"),
  "Grade4_IDH_WT"   = filter(df, idh_status == "wildtype")
)

# ── 2. Cliff's δ 関数 ────────────────────────────────────────────────────────
cliffs_delta <- function(x, y) {
  # x: Mut群, y: WT群
  nx <- length(x); ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  sum(outer(x, y, function(a, b) sign(a - b))) / (nx * ny)
}

# ── 3. 解析関数（1サブセット分）────────────────────────────────────────────
run_subgroup_analysis <- function(data, subset_label) {
  
  wt_vals <- data$LAG3_log2tpm[data$tp53_class4 == "WT"]
  
  # --- 3a. Kruskal-Wallis ---
  kw <- kruskal.test(LAG3_log2tpm ~ tp53_class4, data = data)
  
  kw_result <- tibble(
    subset   = subset_label,
    n_total  = nrow(data),
    n_WT     = sum(data$tp53_class4 == "WT"),
    n_Hot    = sum(data$tp53_class4 == "Hotspot"),
    n_Trunc  = sum(data$tp53_class4 == "Truncating"),
    n_Other  = sum(data$tp53_class4 == "Other_missense"),
    KW_stat  = kw$statistic,
    KW_df    = kw$parameter,
    KW_p     = kw$p.value
  )
  
  # --- 3b. 事前比較（Wilcoxon 3本）+ BH補正 ---
  compare_groups <- c("Hotspot", "Truncating", "Other_missense")
  
  pairwise_raw <- map_dfr(compare_groups, function(grp) {
    mut_vals <- data$LAG3_log2tpm[data$tp53_class4 == grp]
    if (length(mut_vals) == 0) return(NULL)
    
    wt  <- wilcox.test(mut_vals, wt_vals, exact = FALSE)
    med_diff <- median(mut_vals, na.rm = TRUE) - median(wt_vals, na.rm = TRUE)
    delta    <- cliffs_delta(mut_vals, wt_vals)
    
    tibble(
      subset      = subset_label,
      comparison  = paste0(grp, "_vs_WT"),
      n_mut       = length(mut_vals),
      n_wt        = length(wt_vals),
      median_mut  = median(mut_vals, na.rm = TRUE),
      median_wt   = median(wt_vals,  na.rm = TRUE),
      median_diff = med_diff,
      cliffs_delta = delta,
      p_wilcox    = wt$p.value
    )
  })
  
  # BH補正（3本まとめて）
  pairwise_raw <- pairwise_raw %>%
    mutate(
      p_BH = p.adjust(p_wilcox, method = "BH"),
      prespecified = comparison %in% c("Hotspot_vs_WT", "Truncating_vs_WT"),
      effect_label = case_when(
        abs(cliffs_delta) >= 0.474 ~ "large",
        abs(cliffs_delta) >= 0.330 ~ "medium",
        abs(cliffs_delta) >= 0.147 ~ "small",
        TRUE ~ "negligible"
      )
    )
  
  list(kw = kw_result, pairwise = pairwise_raw)
}

# ── 4. 全サブセットで実行 ─────────────────────────────────────────────────────
kw_all       <- map_dfr(names(subsets), ~run_subgroup_analysis(subsets[[.]], .)$kw)
pairwise_all <- map_dfr(names(subsets), ~run_subgroup_analysis(subsets[[.]], .)$pairwise)

cat("\n=== Kruskal-Wallis 結果 ===\n")
kw_all %>% mutate(across(where(is.numeric), ~round(., 4))) %>% print()

cat("\n=== 事前比較結果（主要2本 + 参考1本）===\n")
pairwise_all %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print(n = Inf)

# ── 5. 回帰（GDC Grade4全体・共変量調整）────────────────────────────────────
cat("\n=== 回帰: LAG3 ~ tp53_class4 + source + IDH ===\n")

fit_reg <- lm(LAG3_log2tpm ~ tp53_class4 + source_bin + idh_bin, data = df)
reg_sum <- summary(fit_reg)
cat("R² =", round(reg_sum$r.squared, 4), "\n\n")
print(round(reg_sum$coefficients, 4))

reg_result <- as_tibble(reg_sum$coefficients, rownames = "term") %>%
  rename(beta = Estimate, se = `Std. Error`,
         t_val = `t value`, p_val = `Pr(>|t|)`) %>%
  mutate(
    ci_lo  = beta - 1.96 * se,
    ci_hi  = beta + 1.96 * se,
    subset = "Grade4_all",
    r_squared = reg_sum$r.squared
  )

# ── 6. 各群のmedian・IQR サマリー ────────────────────────────────────────────
cat("\n=== 群別 LAG3 記述統計 ===\n")
desc_stats <- df %>%
  group_by(tp53_class4) %>%
  summarise(
    n      = n(),
    median = round(median(LAG3_log2tpm, na.rm = TRUE), 3),
    q1     = round(quantile(LAG3_log2tpm, 0.25, na.rm = TRUE), 3),
    q3     = round(quantile(LAG3_log2tpm, 0.75, na.rm = TRUE), 3),
    mean   = round(mean(LAG3_log2tpm, na.rm = TRUE), 3),
    sd     = round(sd(LAG3_log2tpm, na.rm = TRUE), 3),
    .groups = "drop"
  )
print(desc_stats)

# ── 7. CSV出力 ────────────────────────────────────────────────────────────────
write_csv(kw_all,       file.path(OUT_DIR, "step12c_kw_results.csv"))
write_csv(pairwise_all, file.path(OUT_DIR, "step12c_pairwise_results.csv"))
write_csv(reg_result,   file.path(OUT_DIR, "step12c_regression.csv"))
write_csv(desc_stats,   file.path(OUT_DIR, "step12c_desc_stats.csv"))

# ── 8. サマリーテキスト ───────────────────────────────────────────────────────
sink(file.path(OUT_DIR, "step12c_summary.txt"))
cat("=== Step 12c Summary ===\n")
cat("実行日時:", format(Sys.time()), "\n\n")
cat("【KW結果】\n"); print(kw_all %>% mutate(across(where(is.numeric), ~round(., 4))))
cat("\n【事前比較（BH補正）】\n"); print(pairwise_all %>% mutate(across(where(is.numeric), ~round(., 4))))
cat("\n【回帰係数】\n"); print(round(reg_sum$coefficients, 4))
cat("\n【記述統計】\n"); print(desc_stats)
sink()

cat("\n=== Step 12c 完了 ===\n")
cat("  step12c_kw_results.csv\n")
cat("  step12c_pairwise_results.csv\n")
cat("  step12c_regression.csv\n")
cat("  step12c_desc_stats.csv\n")
cat("  step12c_summary.txt\n")
