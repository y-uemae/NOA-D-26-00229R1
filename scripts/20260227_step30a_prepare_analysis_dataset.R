# =============================================================================
# step30a_prepare_analysis_dataset.R 【v1.0】
# 目的: Step 28/29の結果をGrade4コホート（n=442）に結合し
#       Stage1・Stage2回帰の入力データセットを作成する
# 入力: 08_final_cohort/final_cohort.csv
#       28_estimate_tumor_purity/estimate_scores_gdc.rds
#       29_ssgsea_hallmark/ssgsea_hallmark_gdc.rds
# 出力先: 30_screening/
#   analysis_dataset_grade4.rds / .csv（n=442）
# 列構成:
#   pair_id, case_barcode, source, grade, tp53_status, idh_status,
#   LAG3（log2TPM）, Tcell, APM, IFNg,（既存共変量）
#   TumorPurity（Step28）, HALLMARK_xxx×50（Step29）
# 作成日: 2026-02-27
# バージョン履歴:
#   v1.0 - 初版
# =============================================================================

library(tidyverse)

# =============================================================================
# 0. パス設定
# =============================================================================
base_dir   <- here::here("results", "TP53", "20260221")
out_dir    <- file.path(base_dir, "30_screening")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("フォルダ作成:", out_dir, "\n")
} else {
  cat("出力フォルダ確認済み:", out_dir, "\n")
}

# --- 入力ファイル ---
input_cohort   <- file.path(base_dir, "08_final_cohort",        "final_cohort.csv")
input_estimate <- file.path(base_dir, "28_estimate_tumor_purity","estimate_scores_gdc.rds")
input_ssgsea   <- file.path(base_dir, "29_ssgsea_hallmark",      "ssgsea_hallmark_gdc.rds")

# --- 出力ファイル ---
output_rds <- file.path(out_dir, "analysis_dataset_grade4.rds")
output_csv <- file.path(out_dir, "analysis_dataset_grade4.csv")

# =============================================================================
# 1. コホートデータ読み込み
# =============================================================================
cat("コホートデータ読み込み中...\n")
cohort <- read_csv(input_cohort, show_col_types = FALSE)
cat("全サンプル数:", nrow(cohort), "\n")
cat("列名確認（最初の20列）:\n")
print(names(cohort)[1:20])

# =============================================================================
# 2. Grade4に絞る
# =============================================================================
cat("\nGrade4フィルタ中...\n")
cohort_g4 <- cohort %>% filter(grade == "Grade4")
cat("Grade4サンプル数:", nrow(cohort_g4), "\n")
cat("source別:\n")
print(table(cohort_g4$source))
cat("tp53_status別:\n")
print(table(cohort_g4$tp53_status))

# =============================================================================
# 3. ESTIMATE（TumorPurity）の結合
# =============================================================================
cat("\nESTIMATE結果結合中...\n")
estimate <- readRDS(input_estimate) %>%
  select(pair_id, TumorPurity, ImmuneScore, StromalScore, ESTIMATEScore)

cohort_g4 <- cohort_g4 %>%
  left_join(estimate, by = "pair_id")

cat("TumorPurity欠損:", sum(is.na(cohort_g4$TumorPurity)), "\n")
cat("TumorPurity summary:\n")
print(summary(cohort_g4$TumorPurity))

# =============================================================================
# 4. ssGSEA（Hallmark50）の結合
# =============================================================================
cat("\nssGSEA結果結合中...\n")
ssgsea <- readRDS(input_ssgsea)
hallmark_cols <- names(ssgsea)[grepl("^HALLMARK_", names(ssgsea))]
cat("Hallmarkセット数:", length(hallmark_cols), "\n")

cohort_g4 <- cohort_g4 %>%
  left_join(ssgsea %>% select(pair_id, all_of(hallmark_cols)), by = "pair_id")

cat("ssGSEA欠損チェック（最初の3セット）:\n")
print(colSums(is.na(cohort_g4[, hallmark_cols[1:3]])))

# =============================================================================
# 5. 必要列の確認
# =============================================================================
cat("\n最終データセット確認:\n")
cat("サンプル数:", nrow(cohort_g4), "\n")
cat("列数:", ncol(cohort_g4), "\n")

# 回帰に使う主要列の存在確認
required_cols <- c("LAG3", "tp53_status", "source", "idh_status",
                   "Tcell", "APM", "IFNg", "TumorPurity")
missing_cols  <- setdiff(required_cols, names(cohort_g4))

if (length(missing_cols) > 0) {
  cat("⚠ 以下の列が見つかりません（final_cohort.csvの列名を確認）:\n")
  print(missing_cols)
  cat("\nfinal_cohort.csvの全列名:\n")
  print(names(cohort_g4))
} else {
  cat("✔ 必要列すべて存在\n")
  cat("\n主要列のsummary:\n")
  print(summary(cohort_g4[, required_cols]))
}

# =============================================================================
# 6. 上書き保存
# =============================================================================
saveRDS(cohort_g4, output_rds)
write_csv(cohort_g4, output_csv)
cat("\n保存完了（RDS）:", output_rds, "\n")
cat("保存完了（CSV）:", output_csv, "\n")
