# =============================================================================
# step06b1_glass_wxs_checkpoint.R
# GBM/Glioma TP53×LAG3 解析 - Step 06b.1: GLASS主解析ファイルへのチェックポイント列追加
#
# 目的:
#   Step06bでは glass_final_cohort.csv（101例）と glass_final_cohort_all101.csv に
#   チェックポイント6遺伝子を追加したが、GLASS主解析個票である
#   glass_final_cohort_wxs_notcga.csv（WXS-only・TCGA除外・n=79）には未追加。
#   本スクリプトはその整合を取り、主解析個票を完成させる。
#
# 処理:
#   glass_final_cohort_wxs_notcga.csv の pair_id を使って、
#   更新済み glass_final_cohort.csv（101例）から対応する発現列を左結合する。
#   → 新規抽出は不要。すでにStep06bで抽出済みのデータを流用。
#
# 追加遺伝子（6遺伝子・Step06bと同一）:
#   PDCD1, CTLA4, TIGIT, HAVCR2, CD274, PDCD1LG2
#
# 出力:
#   05c_glass/glass_final_cohort_wxs_notcga.csv  ★上書き（6遺伝子×2列追加）
#
# 前提:
#   Step06b実行済み（glass_final_cohort.csv に6遺伝子列が追加されていること）
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)

# =============================================================================
# 0. 設定
# =============================================================================

RESULT_DIR    <- here::here("results", "TP53", "20260221")
GLASS_OUT_DIR <- file.path(RESULT_DIR, "05c_glass")

CHECKPOINT_GENES <- c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2")

# 追加する列名（_tpm と _log2tpm）
CP_COLS <- as.vector(rbind(
  paste0(CHECKPOINT_GENES, "_tpm"),
  paste0(CHECKPOINT_GENES, "_log2tpm")
))

# ファイルパス
WXS_PATH  <- file.path(GLASS_OUT_DIR, "glass_final_cohort_wxs_notcga.csv")
FULL_PATH  <- file.path(GLASS_OUT_DIR, "glass_final_cohort.csv")

# =============================================================================
# 1. ファイル存在確認
# =============================================================================

cat("=== Step 06b.1: GLASS主解析ファイルへのチェックポイント列追加 ===\n\n")

if (!file.exists(WXS_PATH)) {
  stop("ERROR: 主解析ファイルが存在しません: ", WXS_PATH)
}
if (!file.exists(FULL_PATH)) {
  stop("ERROR: 101例ファイルが存在しません（Step06b未実行？）: ", FULL_PATH)
}

# =============================================================================
# 2. 読み込み
# =============================================================================

wxs  <- read_csv(WXS_PATH, show_col_types = FALSE)
full <- read_csv(FULL_PATH, show_col_types = FALSE)

cat(sprintf("主解析ファイル(wxs_notcga): %d行 × %d列\n", nrow(wxs), ncol(wxs)))
cat(sprintf("101例ファイル(full):        %d行 × %d列\n", nrow(full), ncol(full)))

# =============================================================================
# 3. 事前チェック
# =============================================================================

# Step06bで追加されているはずの列が存在するか確認
missing_in_full <- setdiff(CP_COLS, names(full))
if (length(missing_in_full) > 0) {
  stop(
    "ERROR: glass_final_cohort.csv にチェックポイント列がありません。",
    "Step06bを先に実行してください。\n",
    "不足列: ", paste(missing_in_full, collapse = ", ")
  )
}
cat("✅ Step06b実行確認: 101例ファイルにチェックポイント列あり\n")

# wxs_notcga のn確認
cat(sprintf("\n主解析ファイル行数確認: %d行 （参照値: 79）\n", nrow(wxs)))
if (nrow(wxs) != 79) {
  cat("  WARNING: 行数が79ではありません。内容を確認してください。\n")
}

# pair_idの結合可能性を確認
id_col_wxs  <- if ("pair_id" %in% names(wxs))  "pair_id"  else "sample_id"
id_col_full <- if ("pair_id" %in% names(full)) "pair_id"  else "sample_id"

matched <- sum(wxs[[id_col_wxs]] %in% full[[id_col_full]])
cat(sprintf("ID結合確認: %d/%d件が101例ファイルと一致\n", matched, nrow(wxs)))
if (matched < nrow(wxs)) {
  unmatched <- wxs[[id_col_wxs]][!wxs[[id_col_wxs]] %in% full[[id_col_full]]]
  cat("  WARNING: 不一致ID:\n")
  cat(paste0("    ", head(unmatched, 10), collapse = "\n"), "\n")
}

# =============================================================================
# 4. 既存チェックポイント列の削除（再実行安全）
# =============================================================================

already_cols <- intersect(names(wxs), CP_COLS)
if (length(already_cols) > 0) {
  cat(sprintf("\nINFO: 既存のチェックポイント列を削除して再追加します（%d列）\n",
              length(already_cols)))
  wxs <- wxs %>% select(-all_of(already_cols))
}

# =============================================================================
# 5. 101例ファイルからチェックポイント列を抽出して結合
# =============================================================================

cp_from_full <- full %>%
  select(all_of(c(id_col_full, CP_COLS))) %>%
  rename(!!id_col_wxs := !!id_col_full)

wxs_updated <- wxs %>%
  left_join(cp_from_full, by = id_col_wxs)

# 列順：既存列の後ろに新列を追加
existing_cols  <- names(wxs)
add_cols_present <- intersect(CP_COLS, names(wxs_updated))
wxs_updated <- wxs_updated %>%
  select(all_of(c(existing_cols, add_cols_present)))

cat(sprintf("\n更新後: %d行 × %d列（追加列: %d）\n",
            nrow(wxs_updated), ncol(wxs_updated), length(add_cols_present)))

# =============================================================================
# 6. 追加列の品質確認
# =============================================================================

cat("\n--- チェックポイント遺伝子 発現統計（主解析79例） ---\n")
for (gene in CHECKPOINT_GENES) {
  col <- paste0(gene, "_log2tpm")
  if (col %in% names(wxs_updated)) {
    vals    <- wxs_updated[[col]]
    na_rate <- mean(is.na(vals)) * 100
    cat(sprintf("  %-12s median=%.3f, mean=%.3f, NA率=%.1f%%\n",
                gene,
                median(vals, na.rm = TRUE),
                mean(vals,   na.rm = TRUE),
                na_rate))
    if (na_rate > 0) cat(sprintf("    WARNING: NA あり（%d件）\n", sum(is.na(vals))))
  }
}

# LAG3との比較（参照）
if ("LAG3_log2tpm" %in% names(wxs_updated)) {
  lag3 <- wxs_updated$LAG3_log2tpm
  cat(sprintf("  %-12s median=%.3f（参照・Step10と一致確認用）\n",
              "LAG3", median(lag3, na.rm = TRUE)))
}

# =============================================================================
# 7. 保存
# =============================================================================

write_csv(wxs_updated, WXS_PATH)
cat(sprintf("\n✅ 上書き保存: glass_final_cohort_wxs_notcga.csv\n"))
cat(sprintf("   %d行 × %d列\n", nrow(wxs_updated), ncol(wxs_updated)))

# =============================================================================
# 8. 完了
# =============================================================================

cat("\n=== Step 06b.1 完了 ===\n")
cat("GLASS主解析個票（n=79）にチェックポイント列の追加が完了しました。\n")
cat("追加列:\n")
for (col in add_cols_present) cat(sprintf("  %s\n", col))
cat("\n次のステップ: Step16a（チェックポイント特異性解析）\n")
cat(sprintf("入力ファイル（確定）: %s\n", WXS_PATH))
