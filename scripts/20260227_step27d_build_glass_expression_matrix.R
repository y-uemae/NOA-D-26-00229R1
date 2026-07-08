# =============================================================================
# step27d_build_glass_expression_matrix.R
# GBM/Glioma TP53xLAG3解析 - Step 27d: GLASS全遺伝子発現行列の整形
#
# 目的:
#   data_mrna_seq_tpm.txt（35438遺伝子×355サンプル）から
#   glass_final_cohort_wxs_notcga.csv（n=79）の
#   79サンプル分を抽出し、GDCと同形式のwide行列を構築する
#
# 確認済み仕様:
#   - 遺伝子ID列名: Hugo_Symbol（gene symbol形式）
#   - サンプルIDキー: pair_id（GLSS-xx-xxxx-TP形式）で79/79完全一致
#   - 区切り文字: タブ
#   - TPM値（非log）が格納されている
#
# 処理内容:
#   1. data_mrna_seq_tpm.txt を読み込み
#   2. pair_idで79サンプルを抽出
#   3. 重複Hugo_Symbol → 最大TPMを採用（GDCと同ルール）
#   4. log2(TPM+1)変換
#   5. wide形式（行=pair_id、列=Hugo_Symbol）でCSV保存
#   6. LAG3値の照合確認
#
# 出力先: results/TP53/20260221/27_expression_matrix/
#   glass_expression_full_log2tpm_wide.csv  # 行=pair_id, 列=遺伝子
#   glass_expression_full_tpm_wide.csv      # TPM生値
#   step27d_gene_list.csv                   # 抽出遺伝子名一覧
#   step27d_log.txt
#
# 作成日: 2026-02-27
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

# =============================================================================
# 0. パス設定
# =============================================================================

BASE_DIR    <- here::here()
RESULT_DIR  <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR     <- file.path(RESULT_DIR, "27_expression_matrix")

GLASS_RNA   <- file.path(BASE_DIR,
                         "data/raw/external_validation/difg_glass/data_mrna_seq_tpm.txt")
GLASS_COHORT <- file.path(RESULT_DIR,
                          "05c_glass/glass_final_cohort_wxs_notcga.csv")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step27d_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("=== Step 27d: GLASS全遺伝子発現行列の整形 開始 ===")

# =============================================================================
# 2. GLASSコホート読み込み（pair_idリスト取得）
# =============================================================================

log_msg("--- 2. コホート読み込み ---")

glass_cohort <- read_csv(GLASS_COHORT, show_col_types = FALSE)
log_msg(sprintf("コホートサンプル数: %d", nrow(glass_cohort)))

target_pair_ids <- glass_cohort$pair_id
log_msg(sprintf("抽出対象pair_id数: %d", length(target_pair_ids)))
log_msg(sprintf("pair_id例（先頭3件）: %s",
                paste(head(target_pair_ids, 3), collapse = ", ")))

# =============================================================================
# 3. 発現ファイルの読み込み
# =============================================================================

log_msg("--- 3. 発現ファイル読み込み（数分かかります）---")

# 全列を文字列で読み込む（数値精度担保のため）
rna_raw <- read.table(GLASS_RNA,
                      sep              = "\t",
                      header           = TRUE,
                      check.names      = FALSE,  # "-"を"."に変換しない
                      stringsAsFactors = FALSE)

log_msg(sprintf("読み込み完了: %d行 × %d列", nrow(rna_raw), ncol(rna_raw)))

# 遺伝子ID列名の確認
gene_col <- names(rna_raw)[1]
log_msg(sprintf("遺伝子ID列名: '%s'", gene_col))

# サンプル列名（1列目以外）
all_sample_cols <- names(rna_raw)[-1]
log_msg(sprintf("総サンプル列数: %d", length(all_sample_cols)))

# =============================================================================
# 4. 79サンプルを抽出
# =============================================================================

log_msg("--- 4. 79サンプルの抽出 ---")

# pair_idと列名のマッチング確認
n_match <- sum(target_pair_ids %in% all_sample_cols)
log_msg(sprintf("pair_idと列名の一致数: %d / %d", n_match, length(target_pair_ids)))

if (n_match == 0) {
  stop("pair_idが列名と一致しません。IDの形式を確認してください。")
}

# 一致しないpair_idを報告
unmatched <- target_pair_ids[!target_pair_ids %in% all_sample_cols]
if (length(unmatched) > 0) {
  log_msg(sprintf("[WARNING] 一致しないpair_id（%d件）: %s",
                  length(unmatched),
                  paste(head(unmatched, 5), collapse = ", ")))
}

# 79サンプル列を抽出
matched_cols <- target_pair_ids[target_pair_ids %in% all_sample_cols]
rna_79 <- rna_raw[, c(gene_col, matched_cols), drop = FALSE]

log_msg(sprintf("抽出後: %d行 × %d列（遺伝子ID列含む）",
                nrow(rna_79), ncol(rna_79)))

# =============================================================================
# 5. 重複Hugo_Symbolの処理（最大TPMを採用）
# =============================================================================

log_msg("--- 5. 重複Hugo_Symbolの処理 ---")

# 全TPM列を数値変換したlong形式で処理
rna_long <- rna_79 %>%
  # Hugo_Symbolが空やNAの行を除外
  filter(!is.na(.data[[gene_col]]),
         .data[[gene_col]] != "",
         .data[[gene_col]] != "Hugo_Symbol") %>%
  # 各列を数値に変換
  mutate(across(-all_of(gene_col),
                ~ suppressWarnings(as.numeric(.))))

# 重複確認
n_before   <- nrow(rna_long)
n_unique   <- n_distinct(rna_long[[gene_col]])
n_dup      <- n_before - n_unique
log_msg(sprintf("処理前: %d行, ユニーク遺伝子数: %d, 重複行数: %d",
                n_before, n_unique, n_dup))

if (n_dup > 0) {
  # 各遺伝子で全サンプルのTPM合計が最大の行を採用
  # （step06と同様に最大TPMを採用; 合計で代表行を選ぶ）
  rna_long <- rna_long %>%
    mutate(row_sum = rowSums(across(-all_of(gene_col)), na.rm = TRUE)) %>%
    group_by(.data[[gene_col]]) %>%
    slice_max(row_sum, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-row_sum)
  
  log_msg(sprintf("重複処理後: %d行（最大TPM行を採用）", nrow(rna_long)))
} else {
  log_msg("重複なし")
}

# =============================================================================
# 6. wide行列の構築（行=pair_id、列=Hugo_Symbol）
# =============================================================================

log_msg("--- 6. wide行列の構築（転置）---")

# 行=pair_id、列=Hugo_Symbolに転置
# 現在: 行=遺伝子、列=サンプル → 転置して: 行=サンプル、列=遺伝子

gene_names <- rna_long[[gene_col]]

# 数値行列として取得（遺伝子名列を除く）
tpm_matrix_t <- rna_long %>%
  select(-all_of(gene_col)) %>%
  as.matrix() %>%
  t()  # 転置: 行=サンプル、列=遺伝子

colnames(tpm_matrix_t) <- gene_names
rownames(tpm_matrix_t) <- matched_cols

log_msg(sprintf("TPM行列: %d サンプル × %d 遺伝子",
                nrow(tpm_matrix_t), ncol(tpm_matrix_t)))

# log2(TPM+1)変換
log2tpm_matrix_t <- log2(tpm_matrix_t + 1)

# =============================================================================
# 7. LAG3値の照合確認
# =============================================================================

log_msg("--- 7. LAG3照合確認 ---")

if ("LAG3" %in% colnames(log2tpm_matrix_t)) {
  lag3_from_matrix <- log2tpm_matrix_t[, "LAG3"]
  
  # コホートのLAG3_log2tpmと照合
  lag3_check <- glass_cohort %>%
    filter(pair_id %in% matched_cols) %>%
    select(pair_id, LAG3_log2tpm) %>%
    mutate(
      LAG3_matrix = lag3_from_matrix[pair_id],
      diff        = abs(LAG3_log2tpm - LAG3_matrix)
    )
  
  n_matched  <- sum(!is.na(lag3_check$LAG3_matrix))
  max_diff   <- max(lag3_check$diff, na.rm = TRUE)
  mean_diff  <- mean(lag3_check$diff, na.rm = TRUE)
  
  log_msg(sprintf("LAG3照合: マッチ数=%d / %d", n_matched, nrow(glass_cohort)))
  log_msg(sprintf("  最大誤差: %.6f", max_diff))
  log_msg(sprintf("  平均誤差: %.6f", mean_diff))
  
  if (max_diff < 0.01) {
    log_msg("  [OK] LAG3値が既存コホートと一致しています")
  } else if (max_diff < 0.1) {
    log_msg("  [注意] 軽微な誤差あり（正規化方法の差の可能性）")
    # 誤差が大きい例を表示
    large_diff <- lag3_check %>% arrange(desc(diff)) %>% head(3)
    for (i in seq_len(nrow(large_diff))) {
      log_msg(sprintf("    %s: cohort=%.4f, matrix=%.4f, diff=%.4f",
                      large_diff$pair_id[i],
                      large_diff$LAG3_log2tpm[i],
                      large_diff$LAG3_matrix[i],
                      large_diff$diff[i]))
    }
  } else {
    log_msg("  [WARNING] 誤差が大きいです。TPM値の確認が必要です")
  }
} else {
  log_msg("[WARNING] LAG3列が行列に存在しません")
}

# =============================================================================
# 8. GDCとの共通遺伝子数確認
# =============================================================================

log_msg("--- 8. GDC行列との共通遺伝子数確認 ---")

gdc_gene_list_path <- file.path(OUT_DIR, "step27b_gene_list.csv")
if (file.exists(gdc_gene_list_path)) {
  gdc_genes  <- read_csv(gdc_gene_list_path, show_col_types = FALSE)$gene_name
  glass_genes <- colnames(log2tpm_matrix_t)
  n_common   <- sum(glass_genes %in% gdc_genes)
  log_msg(sprintf("GDC遺伝子数: %d", length(gdc_genes)))
  log_msg(sprintf("GLASS遺伝子数: %d", length(glass_genes)))
  log_msg(sprintf("共通遺伝子数: %d", n_common))
  log_msg(sprintf("（ssGSEAはこの共通セットで実施予定）"))
} else {
  log_msg("[INFO] step27b_gene_list.csv が見つかりません（GDC未処理）")
}

# =============================================================================
# 9. CSV保存
# =============================================================================

log_msg("--- 9. CSV保存 ---")

# log2TPM wide
log2tpm_df <- as.data.frame(log2tpm_matrix_t) %>%
  rownames_to_column("pair_id")

write_csv(log2tpm_df,
          file.path(OUT_DIR, "glass_expression_full_log2tpm_wide.csv"))
log_msg(sprintf("保存: glass_expression_full_log2tpm_wide.csv  (%d行 × %d列)",
                nrow(log2tpm_df), ncol(log2tpm_df)))

# TPM生値 wide
tpm_df <- as.data.frame(tpm_matrix_t) %>%
  rownames_to_column("pair_id")

write_csv(tpm_df,
          file.path(OUT_DIR, "glass_expression_full_tpm_wide.csv"))
log_msg(sprintf("保存: glass_expression_full_tpm_wide.csv  (%d行 × %d列)",
                nrow(tpm_df), ncol(tpm_df)))

# 遺伝子リスト
gene_list_df <- tibble(gene_name = colnames(log2tpm_matrix_t))
write_csv(gene_list_df,
          file.path(OUT_DIR, "step27d_gene_list.csv"))
log_msg(sprintf("保存: step27d_gene_list.csv  (%d遺伝子)", nrow(gene_list_df)))

# =============================================================================
# 10. 完了サマリー
# =============================================================================

log_msg("=== Step 27d 完了 ===")

cat("\n====================================================\n")
cat("Step 27d 完了サマリー\n")
cat("====================================================\n")
cat(sprintf("GLASSサンプル数:  %d\n", nrow(log2tpm_matrix_t)))
cat(sprintf("GLASS遺伝子数:    %d\n", ncol(log2tpm_matrix_t)))
if (file.exists(gdc_gene_list_path)) {
  cat(sprintf("GDCとの共通遺伝子数: %d\n", n_common))
}
cat(sprintf("LAG3照合: 最大誤差=%.6f, 平均誤差=%.6f\n",
            max_diff, mean_diff))

cat("\n以下をコピーして返してください:\n")
cat(sprintf("  GLASSサンプル数: %d\n", nrow(log2tpm_matrix_t)))
cat(sprintf("  GLASS遺伝子数: %d\n", ncol(log2tpm_matrix_t)))
cat(sprintf("  GDCとの共通遺伝子数: %d\n", n_common))
cat(sprintf("  LAG3照合: マッチ数=%d, 最大誤差=%.6f\n",
            n_matched, max_diff))
cat("====================================================\n")
