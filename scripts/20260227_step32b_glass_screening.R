# ============================================================
# ファイル名: step32b_glass_screening.R
# 目的: GLASSで最小4モデルのΔβ計算
#   M1: Base（TP53 + Tcell + APM + IFNg）
#   M2: Base + TumorPurity
#   M3: Base + TumorPurity + G2M
#   M4: Base + TumorPurity + G2M + MYC + E2F（top3同時投入）
# 出力:
#   - 32_glass_validation_screening/glass_stage2_results.rds/.csv
# ============================================================

library(tidyverse)

base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "32_glass_validation_screening")

glass <- readRDS(file.path(out_dir, "glass_analysis_dataset.rds"))

# ============================================================
# 共通関数：モデルを回してβ_TP53と統計量を返す
# ============================================================
fit_model <- function(data, formula_str, model_label) {
  fit <- lm(as.formula(formula_str), data = data)
  cf  <- summary(fit)$coefficients
  
  # tp53_statusMut 行を取得
  tp53_row <- cf[grep("tp53_statusMut", rownames(cf)), , drop = FALSE]
  
  tibble(
    model       = model_label,
    formula     = formula_str,
    beta_TP53   = tp53_row[1, "Estimate"],
    se_TP53     = tp53_row[1, "Std. Error"],
    t_TP53      = tp53_row[1, "t value"],
    p_TP53      = tp53_row[1, "Pr(>|t|)"],
    adj_r2      = summary(fit)$adj.r.squared,
    n           = nobs(fit)
  )
}

# ============================================================
# 4モデルを順番に実行
# ============================================================
results <- bind_rows(
  fit_model(glass,
            "LAG3 ~ tp53_status + Tcell + APM + IFNg",
            "M1_Base"),
  fit_model(glass,
            "LAG3 ~ tp53_status + Tcell + APM + IFNg + TumorPurity",
            "M2_Base_Purity"),
  fit_model(glass,
            "LAG3 ~ tp53_status + Tcell + APM + IFNg + TumorPurity + G2M",
            "M3_Purity_G2M"),
  fit_model(glass,
            "LAG3 ~ tp53_status + Tcell + APM + IFNg + TumorPurity + G2M + MYC + E2F",
            "M4_Purity_top3")
)

# ============================================================
# Δβ% 計算（M1_Baseを基準）
# ============================================================
beta_base <- results$beta_TP53[results$model == "M1_Base"]

results <- results %>%
  mutate(
    delta_beta_pct = round((beta_TP53 - beta_base) / abs(beta_base) * 100, 1),
    beta_TP53      = round(beta_TP53, 4),
    se_TP53        = round(se_TP53,   4),
    p_TP53         = signif(p_TP53,   3),
    adj_r2         = round(adj_r2,    3)
  )

# ============================================================
# 結果表示
# ============================================================
cat("=== GLASS 最小4モデル Δβ結果 ===\n\n")
results %>%
  select(model, beta_TP53, se_TP53, p_TP53, delta_beta_pct, adj_r2, n) %>%
  print(n = Inf)

cat("\n=== GDCとの対比（参考） ===\n")
gdc_ref <- tibble(
  model          = c("M1_Base", "M2_Base_Purity", "M3_Purity_G2M", "M4_Purity_top3"),
  GDC_beta_TP53  = c(0.2665,    0.2406,           0.2007,          0.1973),
  GDC_delta_pct  = c(0,        -9.7,             -16.6,           -18.0)
)
print(gdc_ref)

# ============================================================
# 保存
# ============================================================
saveRDS(results, file.path(out_dir, "glass_stage2_results.rds"))
write_excel_csv(results, file.path(out_dir, "glass_stage2_results.csv"))

cat("\n✅ Step 32b 完了\n")
