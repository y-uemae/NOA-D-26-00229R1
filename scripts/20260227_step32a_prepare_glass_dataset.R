# ============================================================
# ファイル名: step32a_prepare_glass_dataset_v3.R
# 目的: GLASS検証用データセット作成
#       glass_final_cohort_wxs_notcga.csv を起点に
#       TumorPurity・Hallmark3スコアを結合するだけ
# 出力:
#   - 32_glass_validation_screening/glass_analysis_dataset.rds/.csv
# ============================================================

library(tidyverse)

base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "32_glass_validation_screening")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. GLASSコホートファイルを読み込み（79行・全列揃い）
# ============================================================
glass_cohort <- read_csv(
  file.path(base_dir, "05c_glass/glass_final_cohort_wxs_notcga.csv"),
  show_col_types = FALSE
)

# ============================================================
# 2. 共変量スコアを行平均で計算（_log2tpm列が既存）
#    HLA-A等はHLA.A形式で保存されているためそのまま使用
# ============================================================
apm_cols   <- c("B2M_log2tpm", "TAP1_log2tpm", "TAP2_log2tpm", "TAPBP_log2tpm",
                "HLA.A_log2tpm", "HLA.B_log2tpm", "HLA.C_log2tpm", "NLRC5_log2tpm")
ifng_cols  <- c("STAT1_log2tpm", "IRF1_log2tpm", "IRF9_log2tpm",
                "CXCL9_log2tpm", "CXCL10_log2tpm", "CXCL11_log2tpm",
                "GBP1_log2tpm", "GBP2_log2tpm", "GBP4_log2tpm",
                "GBP5_log2tpm", "IDO1_log2tpm")
tcell_cols <- c("CD3D_log2tpm", "CD3E_log2tpm", "CD3G_log2tpm",
                "CD8A_log2tpm", "CD8B_log2tpm",
                "GZMA_log2tpm", "GZMB_log2tpm", "PRF1_log2tpm")

glass_base <- glass_cohort %>%
  mutate(
    tp53_status = factor(tp53_status, levels = c("WT", "Mut")),
    LAG3        = LAG3_log2tpm,
    APM         = rowMeans(select(., all_of(apm_cols)),   na.rm = TRUE),
    IFNg        = rowMeans(select(., all_of(ifng_cols)),  na.rm = TRUE),
    Tcell       = rowMeans(select(., all_of(tcell_cols)), na.rm = TRUE)
  ) %>%
  select(pair_id, tp53_status, LAG3, APM, IFNg, Tcell)

cat("=== glass_base ===\n")
print(head(glass_base, 5))
cat("行数:", nrow(glass_base), "\n")

# ============================================================
# 3. TumorPurity を結合
# ============================================================
est_glass <- readRDS(
  file.path(base_dir, "28_estimate_tumor_purity/estimate_scores_glass.rds")
)

# pair_idのドット→ハイフン変換（念のため）
est_glass <- est_glass %>%
  mutate(pair_id = gsub("\\.", "-", pair_id)) %>%
  select(pair_id, TumorPurity)

cat("\n=== TumorPurity 結合前マッチ確認 ===\n")
cat("est_glass件数:", nrow(est_glass), "\n")
cat("マッチ件数:", sum(glass_base$pair_id %in% est_glass$pair_id), "\n")

# ============================================================
# 4. Hallmark3スコアを結合（G2M / MYC_V1 / E2F）
# ============================================================
ssgsea_glass <- readRDS(
  file.path(base_dir, "29_ssgsea_hallmark/ssgsea_hallmark_glass.rds")
)

ssgsea_glass <- ssgsea_glass %>%
  mutate(pair_id = gsub("\\.", "-", pair_id)) %>%
  select(
    pair_id,
    G2M = HALLMARK_G2M_CHECKPOINT,
    MYC = HALLMARK_MYC_TARGETS_V1,
    E2F = HALLMARK_E2F_TARGETS
  )

cat("\n=== ssGSEA 結合前マッチ確認 ===\n")
cat("ssgsea_glass件数:", nrow(ssgsea_glass), "\n")
cat("マッチ件数:", sum(glass_base$pair_id %in% ssgsea_glass$pair_id), "\n")

# ============================================================
# 5. 全結合・最終確認
# ============================================================
glass_final <- glass_base %>%
  left_join(est_glass,    by = "pair_id") %>%
  left_join(ssgsea_glass, by = "pair_id")

cat("\n=== 最終データセット ===\n")
cat("行数:", nrow(glass_final), "  列数:", ncol(glass_final), "\n")
cat("列名:", colnames(glass_final), "\n")
cat("\n=== 欠損値確認 ===\n")
print(colSums(is.na(glass_final)))
cat("\n=== 先頭5行 ===\n")
print(head(glass_final, 5))
cat("\n=== tp53_status件数 ===\n")
print(table(glass_final$tp53_status))

# ============================================================
# 6. 保存
# ============================================================
saveRDS(glass_final, file.path(out_dir, "glass_analysis_dataset.rds"))
write_excel_csv(glass_final, file.path(out_dir, "glass_analysis_dataset.csv"))

cat("\n✅ Step 32a 完了\n")
cat("出力: 32_glass_validation_screening/glass_analysis_dataset.rds/.csv\n")
