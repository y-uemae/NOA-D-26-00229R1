# ============================================================
# ファイル名: step32c_make_suppl_table_s5.R
# 目的: GDC vs GLASS の Δβ比較表（Supplementary Table S5）作成
# 出力:
#   - 32_glass_validation_screening/suppltable_s5_purity_proliferation.csv
# ============================================================

library(tidyverse)

base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "32_glass_validation_screening")

# GLASS結果（Step32b）
glass_res <- readRDS(file.path(out_dir, "glass_stage2_results.rds"))

# GDC結果（引継書の確定値を直接入力）
gdc_res <- tibble(
  model          = c("M1_Base", "M2_Base_Purity",
                     "M3_Purity_G2M", "M4_Purity_top3"),
  GDC_formula    = c(
    "LAG3 ~ TP53 + Tcell + APM + IFNg",
    "LAG3 ~ TP53 + Tcell + APM + IFNg + TumorPurity",
    "LAG3 ~ TP53 + Tcell + APM + IFNg + TumorPurity + G2M",
    "LAG3 ~ TP53 + Tcell + APM + IFNg + TumorPurity + G2M + MYC + E2F"
  ),
  GDC_n          = 442,
  GDC_beta_TP53  = c(0.2665, 0.2406, 0.2007, 0.1973),
  GDC_p_TP53     = c(9.30e-7, 3.24e-6, 1.11e-4, 1.44e-4),
  GDC_delta_pct  = c(0, -9.7, -16.6, -18.0),
  GDC_adj_r2     = c(0.303, 0.371, 0.389, 0.390)
)

# GLASS結果を整形
glass_sel <- glass_res %>%
  select(
    model,
    GLASS_n          = n,
    GLASS_beta_TP53  = beta_TP53,
    GLASS_p_TP53     = p_TP53,
    GLASS_delta_pct  = delta_beta_pct,
    GLASS_adj_r2     = adj_r2
  )

# 結合
suppl_s5 <- gdc_res %>%
  left_join(glass_sel, by = "model") %>%
  mutate(
    # モデルの説明ラベル
    Model_description = c(
      "Base (TP53 + immune scores)",
      "Base + Tumor purity",
      "Base + Tumor purity + G2M checkpoint",
      "Base + Tumor purity + G2M + MYC targets v1 + E2F targets"
    ),
    # p値を有効数字3桁のsignif形式に
    GDC_p_TP53   = signif(GDC_p_TP53,   3),
    GLASS_p_TP53 = signif(GLASS_p_TP53, 3),
    # Δβを文字列化（基準モデルはハイフン）
    GDC_delta_pct_str   = ifelse(GDC_delta_pct   == 0, "—",
                                 paste0(GDC_delta_pct,   "%")),
    GLASS_delta_pct_str = ifelse(GLASS_delta_pct == 0, "—",
                                 paste0(GLASS_delta_pct, "%"))
  ) %>%
  select(
    Model = Model_description,
    `GDC n`            = GDC_n,
    `GDC β(TP53)`      = GDC_beta_TP53,
    `GDC p`            = GDC_p_TP53,
    `GDC Δβ%`          = GDC_delta_pct_str,
    `GDC adj.R²`       = GDC_adj_r2,
    `GLASS n`          = GLASS_n,
    `GLASS β(TP53)`    = GLASS_beta_TP53,
    `GLASS p`          = GLASS_p_TP53,
    `GLASS Δβ%`        = GLASS_delta_pct_str,
    `GLASS adj.R²`     = GLASS_adj_r2
  )

cat("=== Supplementary Table S5 プレビュー ===\n")
print(suppl_s5, width = Inf)

# 保存
write_excel_csv(suppl_s5,
                file.path(out_dir, "suppltable_s5_purity_proliferation.csv"))

cat("\n✅ Supplementary Table S5 完了\n")
