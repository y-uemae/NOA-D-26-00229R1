# =============================================================================
# step28_estimate_tumor_purity.R 【v2.5 warning修正版】
# 修正:
#   - pair_idをcharacterで統一し、join直前にのみinteger化（warning解消）
#   - select(-Description)をselect(-matches("^Description"))に変更
# 作成日: 2026-02-27
# バージョン履歴:
#   v1.0 - 初版
#   v1.1 - 転置時colnames修正
#   v2.0 - 全遺伝子GCT入力に変更
#   v2.1 - colnames明示設定
#   v2.2 - GCT→TSV形式に変更、pair_id↔rna_sample_id対応付けを追加
#   v2.3 - Pプレフィックス付与、列フィルタ厳密化、stopifnot追加
#   v2.4 - GLSS列名ドット変換対応、pair_id復元追加
#   v2.5 - pair_idをcharacter統一→join直前にinteger化（warning解消）
#         - Description列除去をmatches()で確実化
# =============================================================================

library(tidyverse)
library(data.table)
library(estimate)

# =============================================================================
# 0. パス設定
# =============================================================================
base_dir   <- here::here("results", "TP53", "20260221")
expr_dir   <- file.path(base_dir, "27_expression_matrix")
cohort_dir <- file.path(base_dir, "08_final_cohort")
out_dir    <- file.path(base_dir, "28_estimate_tumor_purity")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("フォルダ作成:", out_dir, "\n")
} else {
  cat("出力フォルダ確認済み:", out_dir, "\n")
}

input_gdc    <- file.path(expr_dir,   "expression_full_log2tpm_wide.csv")
input_glass  <- file.path(expr_dir,   "glass_expression_full_log2tpm_wide.csv")
input_cohort <- file.path(cohort_dir, "final_cohort.csv")

output_gdc_rds   <- file.path(out_dir, "estimate_scores_gdc.rds")
output_gdc_csv   <- file.path(out_dir, "estimate_scores_gdc.csv")
output_glass_rds <- file.path(out_dir, "estimate_scores_glass.rds")
output_glass_csv <- file.path(out_dir, "estimate_scores_glass.csv")

# =============================================================================
# 1. pair_id対応表の読み込み（GDC用）
# =============================================================================
cat("pair_id対応表読み込み中...\n")
cohort <- read_csv(input_cohort, show_col_types = FALSE) %>%
  select(pair_id, rna_sample_id, case_barcode, source, grade,
         tp53_status, idh_status)
cat("cohortサンプル数:", nrow(cohort), "\n")
cat("pair_id範囲:", min(cohort$pair_id), "〜", max(cohort$pair_id), "\n")

# =============================================================================
# 2. 関数定義: run_estimate_full()
# =============================================================================
run_estimate_full <- function(expr_wide_path, cohort_name, out_dir,
                              out_rds, out_csv,
                              id_map = NULL) {
  cat("\n========================================\n")
  cat(cohort_name, "処理開始\n")
  cat("========================================\n")
  
  # --- 読み込み ---
  cat("発現行列読み込み中...\n")
  expr_wide <- read_csv(expr_wide_path, show_col_types = FALSE)
  cat("読み込み完了:", nrow(expr_wide), "サンプル×",
      ncol(expr_wide) - 1, "遺伝子\n")
  
  # --- 転置: 行=遺伝子, 列=サンプル ---
  cat("転置中（全遺伝子）...\n")
  sample_ids <- as.character(expr_wide$pair_id)
  
  mat <- expr_wide %>%
    select(-pair_id) %>%
    as.matrix()
  rownames(mat) <- sample_ids
  mat <- t(mat)
  
  cat("転置完了:", nrow(mat), "遺伝子×", ncol(mat), "サンプル\n")
  cat("colnames先頭5件（変換前）:", paste(head(colnames(mat), 5), collapse = ", "), "\n")
  
  # --- GDC: 数値IDにPプレフィックス付与（X1変換回避）---
  is_numeric_ids <- all(grepl("^[0-9]+$", colnames(mat)))
  if (is_numeric_ids) {
    colnames(mat) <- paste0("P", colnames(mat))
    cat("NOTE: 数値IDにPプレフィックスを付与\n")
    cat("colnames先頭5件（変換後）:", paste(head(colnames(mat), 5), collapse = ", "), "\n")
  }
  
  # --- TSV書き出し ---
  tmp_tsv_raw  <- file.path(out_dir, paste0("tmp_", cohort_name, "_full.tsv"))
  tmp_tsv_filt <- file.path(out_dir, paste0("tmp_", cohort_name, "_filtered.gct"))
  tmp_gct_out  <- file.path(out_dir, paste0("tmp_", cohort_name, "_scored.gct"))
  
  cat("TSV書き出し中（全", nrow(mat), "遺伝子、数分かかります）...\n")
  
  con <- file(tmp_tsv_raw, open = "wt")
  writeLines(paste(c("NAME", "Description", colnames(mat)), collapse = "\t"), con)
  close(con)
  
  dt_gct <- as.data.table(mat, keep.rownames = "NAME")
  dt_gct[, Description := NAME]
  setcolorder(dt_gct, c("NAME", "Description",
                        setdiff(names(dt_gct), c("NAME", "Description"))))
  
  fwrite(dt_gct, file = tmp_tsv_raw, append = TRUE,
         sep = "\t", col.names = FALSE, quote = FALSE)
  
  cat("TSV書き出し完了\n")
  cat("TSV先頭3行確認:\n")
  for (l in readLines(tmp_tsv_raw, n = 3)) cat(" ", substr(l, 1, 120), "\n")
  
  # --- filterCommonGenes ---
  cat("filterCommonGenes実行中...\n")
  filterCommonGenes(
    input.f  = tmp_tsv_raw,
    output.f = tmp_tsv_filt,
    id       = "GeneSymbol"
  )
  
  # --- estimateScore ---
  cat("estimateScore実行中...\n")
  estimateScore(
    input.ds  = tmp_tsv_filt,
    output.ds = tmp_gct_out,
    platform  = "illumina"
  )
  
  # --- 結果読み込み ---
  cat("スコア整形中...\n")
  scored <- read_tsv(tmp_gct_out, skip = 2, show_col_types = FALSE) %>%
    rename(ScoreType = NAME) %>%
    select(-matches("^Description"))  # Description.1等を確実に除去
  
  cat("scored列数:", ncol(scored), "\n")
  cat("scored列名先頭10件:", paste(head(names(scored), 10), collapse = ", "), "\n")
  
  # --- サンプル列をパターンで厳密フィルタ ---
  if (cohort_name == "GDC") {
    sample_cols <- names(scored)[grepl("^P[0-9]+$", names(scored))]
  } else {
    sample_cols <- names(scored)[grepl("^GLSS[-.]", names(scored))]
  }
  
  cat("検出されたサンプル列数:", length(sample_cols), "\n")
  cat("サンプル列先頭5件:", paste(head(sample_cols, 5), collapse = ", "), "\n")
  
  if (length(sample_cols) == 0) {
    stop("サンプル列が検出されませんでした。scored GCTのヘッダーを確認してください。")
  }
  
  # --- ★pair_idはcharacterで統一（case_whenでのinteger変換を廃止）---
  score_df <- scored %>%
    select(ScoreType, all_of(sample_cols)) %>%
    pivot_longer(-ScoreType, names_to = "pair_id", values_to = "value") %>%
    pivot_wider(names_from = ScoreType, values_from = value) %>%
    mutate(
      pair_id = if (cohort_name == "GDC") sub("^P", "", pair_id)   # "P1"→"1"（文字列）
      else gsub("\\.", "-", pair_id),                      # "GLSS.19."→"GLSS-19-"
      TumorPurity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)
    )
  
  cat("サンプル数:", nrow(score_df), "\n")
  cat("pair_id先頭5件:", paste(head(score_df$pair_id, 5), collapse = ", "), "\n")
  
  # --- GDCの場合: integer化→安全装置→id_map結合 ---
  if (!is.null(id_map)) {
    # join直前にのみinteger化（warningが出る変換をここだけに限定）
    score_df <- score_df %>% mutate(pair_id = as.integer(pair_id))
    
    stopifnot("pair_idにNAがあります"            = !anyNA(score_df$pair_id))
    stopifnot("pair_idがid_mapにない行があります" =
                all(score_df$pair_id %in% id_map$pair_id))
    
    score_df <- score_df %>%
      left_join(id_map, by = "pair_id") %>%
      select(pair_id, rna_sample_id, case_barcode, source, grade,
             tp53_status, idh_status,
             StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity)
    
    cat("rna_sample_id結合完了\n")
    cat("結合確認（先頭3件）:\n")
    print(head(score_df %>% select(pair_id, rna_sample_id, TumorPurity), 3))
  }
  
  cat("ESTIMATEScore summary:\n")
  print(summary(score_df$ESTIMATEScore))
  cat("TumorPurity summary:\n")
  print(summary(score_df$TumorPurity))
  
  # --- 上書き保存 ---
  saveRDS(score_df, out_rds)
  write_csv(score_df, out_csv)
  cat("保存完了（RDS）:", out_rds, "\n")
  cat("保存完了（CSV）:", out_csv, "\n")
  
  # --- 一時ファイル削除 ---
  file.remove(c(tmp_tsv_raw, tmp_tsv_filt, tmp_gct_out))
  cat("一時ファイル削除完了\n")
  
  return(score_df)
}

# =============================================================================
# 3. GDC実行
# =============================================================================
scores_gdc <- run_estimate_full(
  expr_wide_path = input_gdc,
  cohort_name    = "GDC",
  out_dir        = out_dir,
  out_rds        = output_gdc_rds,
  out_csv        = output_gdc_csv,
  id_map         = cohort
)

# =============================================================================
# 4. GLASS実行
# =============================================================================
scores_glass <- run_estimate_full(
  expr_wide_path = input_glass,
  cohort_name    = "GLASS",
  out_dir        = out_dir,
  out_rds        = output_glass_rds,
  out_csv        = output_glass_csv,
  id_map         = NULL
)

# =============================================================================
# 5. 完了サマリー
# =============================================================================
cat("\n========================================\n")
cat("Step 28 完了\n")
cat("========================================\n")
cat("GDC:   ", nrow(scores_gdc),   "サンプル →", output_gdc_csv,   "\n")
cat("GLASS: ", nrow(scores_glass), "サンプル →", output_glass_csv, "\n")
cat("GDC列:   pair_id, rna_sample_id, case_barcode, source, grade,\n")
cat("         tp53_status, idh_status,\n")
cat("         StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity\n")
cat("GLASS列: pair_id, StromalScore, ImmuneScore, ESTIMATEScore, TumorPurity\n")
