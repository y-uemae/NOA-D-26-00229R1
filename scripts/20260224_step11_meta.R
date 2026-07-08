# =============================================================================
# Step 11: メタ解析（β統合版・修正版）
# Phase 1: GDC source別サブグループ回帰 → β・SE取得
# Phase 2: GLASS WXS_noTCGA 回帰 → β・SE取得
# Phase 3: 3コホート fixed + random effectメタ解析（metafor）
# 入力:
#   08_final_cohort/final_cohort.csv
#   05c_glass/glass_final_cohort_wxs_notcga.csv
# 出力:
#   11_meta_analysis/step11_subgroup_regression.csv
#   11_meta_analysis/step11_meta_results.csv
#   11_meta_analysis/step11_meta_forest_summary.csv
#   11_meta_analysis/step11_forest_check.png
# =============================================================================

library(tidyverse)
library(metafor)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR    <- here::here("results", "TP53", "20260221")
GDC_CSV     <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
GLASS_CSV   <- file.path(BASE_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")
OUT_DIR     <- file.path(BASE_DIR, "11_meta_analysis")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. GDCデータ読み込み・整形 ────────────────────────────────────────────────
gdc_raw <- read_csv(GDC_CSV, show_col_types = FALSE)

gdc_g4 <- gdc_raw %>%
  filter(include_flag == TRUE, grade == "Grade4") %>%
  mutate(
    tp53_bin = as.integer(tp53_status == "mutant"),
    idh_bin  = as.integer(idh_status  == "mutant")
  )

cat("=== GDC Grade4 サンプル数 ===\n")
gdc_g4 %>% count(source, tp53_status) %>% print()

# ── 2. GLASS データ読み込み・確認 ─────────────────────────────────────────────
glass_raw <- read_csv(GLASS_CSV, show_col_types = FALSE)

cat("\n=== GLASS WXS_noTCGA 列名確認 ===\n")
cat("列名:", paste(names(glass_raw), collapse = ", "), "\n")
cat("行数:", nrow(glass_raw), "\n")

cat("\n=== GLASS tp53_status の値 ===\n")
glass_raw %>% count(tp53_status) %>% print()

# tp53_bin: 引継書より GLASS は "Mut"/"WT"
glass_df <- glass_raw %>%
  mutate(
    tp53_bin = as.integer(tp53_status == "Mut"),
    # idh_status列があれば共変量に使う（なければ省略）
    idh_bin  = if ("idh_status" %in% names(.))
      as.integer(idh_status == "mutant" | idh_status == "Mut")
    else NA_integer_
  )

cat("\n=== GLASS tp53_bin 確認 ===\n")
glass_df %>% count(tp53_status, tp53_bin) %>% print()

# LAG3列の確認
if (!"LAG3_log2tpm" %in% names(glass_df)) {
  cat("⚠️  LAG3_log2tpm列なし。実際のLAG3列名:\n")
  print(names(glass_df)[str_detect(names(glass_df), "LAG3|lag3")])
  stop("LAG3列名を確認してください")
}

# ── 3. 回帰関数（共通）────────────────────────────────────────────────────────
run_lm <- function(df, study_label, use_idh = TRUE) {
  # idh_binが使えるか確認
  has_idh <- use_idh && "idh_bin" %in% names(df) && !all(is.na(df$idh_bin))
  
  formula_str <- if (has_idh) "LAG3_log2tpm ~ tp53_bin + idh_bin" else "LAG3_log2tpm ~ tp53_bin"
  
  fit <- lm(as.formula(formula_str), data = df)
  cf  <- summary(fit)$coefficients
  
  tibble(
    study      = study_label,
    model_used = formula_str,
    n_total    = nrow(df),
    n_mut      = sum(df$tp53_bin == 1, na.rm = TRUE),
    n_wt       = sum(df$tp53_bin == 0, na.rm = TRUE),
    beta_TP53  = cf["tp53_bin", "Estimate"],
    se_TP53    = cf["tp53_bin", "Std. Error"],
    t_val      = cf["tp53_bin", "t value"],
    p_val      = cf["tp53_bin", "Pr(>|t|)"],
    ci_lo      = beta_TP53 - 1.96 * se_TP53,
    ci_hi      = beta_TP53 + 1.96 * se_TP53,
    r_squared  = summary(fit)$r.squared
  )
}

# ── 4. 各コホートの回帰実行 ───────────────────────────────────────────────────
res_tcga   <- run_lm(filter(gdc_g4, source == "TCGA"),       "GDC_TCGA",       use_idh = TRUE)
res_cptac  <- run_lm(filter(gdc_g4, source == "CPTAC_HCMI"), "GDC_CPTAC_HCMI", use_idh = TRUE)
res_glass  <- run_lm(glass_df,                                "GLASS_WXS",      use_idh = FALSE)
# GLASS: IDH情報の有無に応じて自動切替（use_idh=FALSEで TP53 only）
# TP53 only にする理由: GLASS主解析との一貫性、IDH情報の信頼性

cat("\n=== 全コホート回帰結果 ===\n")
bind_rows(res_tcga, res_cptac, res_glass) %>%
  select(study, model_used, n_mut, n_wt, beta_TP53, se_TP53, ci_lo, ci_hi, p_val, r_squared) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

# ── 5. メタ解析投入データ ─────────────────────────────────────────────────────
meta_df <- bind_rows(res_tcga, res_cptac, res_glass)

# SE欠損チェック
if (any(is.na(meta_df$se_TP53))) {
  cat("⚠️  SE にNAがあります:\n")
  print(filter(meta_df, is.na(se_TP53)))
  stop("SE取得失敗。回帰結果を確認してください")
}

cat("\n=== メタ解析投入データ（3コホート）===\n")
meta_df %>%
  select(study, n_mut, n_wt, beta_TP53, se_TP53, ci_lo, ci_hi, p_val) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

# ── 6. メタ解析実行 ───────────────────────────────────────────────────────────
fit_fe <- rma(yi = beta_TP53, sei = se_TP53, data = meta_df,
              method = "FE",   slab = study)
fit_re <- rma(yi = beta_TP53, sei = se_TP53, data = meta_df,
              method = "REML", slab = study)

cat("\n============================\n")
cat("=== Fixed Effect Model  ===\n")
cat("============================\n")
print(summary(fit_fe))

cat("\n============================\n")
cat("=== Random Effect Model ===\n")
cat("============================\n")
print(summary(fit_re))

# ── 7. 結果整理 ───────────────────────────────────────────────────────────────
extract_meta <- function(fit, model_type) {
  tibble(
    model       = model_type,
    k           = fit$k,
    beta_pooled = as.numeric(fit$b),
    ci_lo       = fit$ci.lb,
    ci_hi       = fit$ci.ub,
    se_pooled   = fit$se,
    z_val       = fit$zval,
    p_val       = fit$pval,
    Q           = fit$QE,
    Q_df        = fit$QEdf,
    Q_p         = fit$QEp,
    I2          = fit$I2,
    tau2        = fit$tau2
  )
}

results_tbl <- bind_rows(
  extract_meta(fit_fe, "Fixed Effect"),
  extract_meta(fit_re, "Random Effect")
)

cat("\n=== 統合結果サマリー ===\n")
results_tbl %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>%
  print()

# 参照値比較
cat("\n=== 参照値比較 ===\n")
cat("  旧コホート Fixed  β=0.183, p=0.016\n")
cat("  旧コホート Random β=0.270, p=0.172, I²=79%\n")
cat(sprintf("  新コホート Fixed  β=%.3f (95%%CI %.3f–%.3f), p=%.4f\n",
            results_tbl$beta_pooled[1], results_tbl$ci_lo[1],
            results_tbl$ci_hi[1], results_tbl$p_val[1]))
cat(sprintf("  新コホート Random β=%.3f (95%%CI %.3f–%.3f), p=%.4f, I²=%.1f%%\n",
            results_tbl$beta_pooled[2], results_tbl$ci_lo[2],
            results_tbl$ci_hi[2], results_tbl$p_val[2], results_tbl$I2[2]))

# ── 8. Forest summary CSV（Step12可視化用）────────────────────────────────────
individual_rows <- meta_df %>%
  transmute(
    study, type = "individual",
    n_mut, n_wt,
    beta  = beta_TP53,
    se    = se_TP53,
    ci_lo, ci_hi, p_val,
    I2 = NA_real_, tau2 = NA_real_
  )

pooled_rows <- results_tbl %>%
  transmute(
    study = paste0("Pooled_", sub(" ", "_", model)),
    type  = "pooled",
    n_mut = NA_real_, n_wt = NA_real_,
    beta  = beta_pooled,
    se    = se_pooled,
    ci_lo, ci_hi, p_val, I2, tau2
  )

forest_summary <- bind_rows(individual_rows, pooled_rows)

# ── 9. CSV出力 ────────────────────────────────────────────────────────────────
write_csv(bind_rows(res_tcga, res_cptac, res_glass),
          file.path(OUT_DIR, "step11_subgroup_regression.csv"))
write_csv(results_tbl,
          file.path(OUT_DIR, "step11_meta_results.csv"))
write_csv(forest_summary,
          file.path(OUT_DIR, "step11_meta_forest_summary.csv"))

# ── 10. 確認用 forest plot ────────────────────────────────────────────────────
png(file.path(OUT_DIR, "step11_forest_check.png"),
    width = 1000, height = 450, res = 130)
forest(fit_re,
       xlab    = "β (TP53 Mutant vs Wildtype, LAG3 log2TPM)",
       main    = "Meta-analysis: TP53 mutation and LAG3 expression in Grade4 Glioma",
       header  = c("Study", "β [95% CI]"),
       refline = 0)
dev.off()

cat("\n=== 出力完了 ===\n")
cat("  step11_subgroup_regression.csv\n")
cat("  step11_meta_results.csv\n")
cat("  step11_meta_forest_summary.csv\n")
cat("  step11_forest_check.png\n")
cat("metafor version:", as.character(packageVersion("metafor")), "\n")
