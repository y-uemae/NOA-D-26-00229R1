# =============================================================================
# step27b_build_full_expression_matrix.R
# GBM/Glioma TP53xLAG3解析 - Step 27b: 全遺伝子発現行列の構築
#
# 目的:
#   GDC augmented_star_gene_counts.tsv（951サンプル分）から
#   全遺伝子のTPMを抽出し、サンプル×遺伝子のwide行列を構築する
#   （Step 28以降のssGSEA・ESTIMATE解析に使用）
#
# step06との違い:
#   - 抽出対象: 28遺伝子 → 全遺伝子（gene_nameが存在する全行）
#   - 出力形式: CSVのみ（RDSなし）
#   - ペアリング元: final_cohort.csv（step06はpairs_final.csv）
#
# 処理内容:
#   1. final_cohort.csvからrna_file_idを取得（951サンプル）
#   2. 各tsvからgene_name + tpm_unstrandedを全行抽出
#   3. 重複gene_name → 最大TPMを採用（step06と同ルール）
#   4. log2(TPM+1)変換
#   5. wide形式（行=pair_id、列=gene_name）でCSV保存
#
# 出力先: results/TP53/20260221/27_expression_matrix/
#   expression_full_log2tpm_wide.csv   # 行=pair_id, 列=遺伝子（log2TPM）
#   expression_full_tpm_wide.csv       # 行=pair_id, 列=遺伝子（TPM生値）
#   step27b_gene_list.csv              # 抽出できた遺伝子名一覧
#   step27b_sample_qc.csv              # サンプル別QC（欠損遺伝子数等）
#   step27b_log.txt
#
# 実行時間の目安: 20〜40分（951ファイル × 全遺伝子）
# 作成日: 2026-02-27
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)

# =============================================================================
# 0. パス設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
RNA_DIR    <- file.path(BASE_DIR, "data/raw/GDC/glioma/RNAseq/downloads")
OUT_DIR    <- file.path(RESULT_DIR, "27_expression_matrix")

GDC_COHORT <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TPM_COL      <- "tpm_unstranded"
GENE_NAME_COL <- "gene_name"

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step27b_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("=== Step 27b: 全遺伝子発現行列の構築 開始 ===")

# =============================================================================
# 2. コホート読み込み・rna_file_id取得
# =============================================================================

log_msg("--- 2. コホート読み込み ---")

cohort <- read_csv(GDC_COHORT, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

log_msg(sprintf("コホートサンプル数: %d", nrow(cohort)))

if (!"rna_file_id" %in% names(cohort)) {
  stop("rna_file_id 列が見つかりません")
}

# pair_id と rna_file_id の対応表
id_map <- cohort %>%
  select(pair_id, rna_file_id, source, grade) %>%
  filter(!is.na(rna_file_id)) %>%
  distinct(rna_file_id, .keep_all = TRUE)  # 同一rna_file_idの重複除去

log_msg(sprintf("抽出対象 rna_file_id（重複除去後）: %d件", nrow(id_map)))

# =============================================================================
# 3. tsvファイルの構造確認（1ファイルだけ先読み）
# =============================================================================

log_msg("--- 3. tsvファイル構造の確認 ---")

sample_fid  <- id_map$rna_file_id[1]
sample_path <- file.path(RNA_DIR, sample_fid,
                         list.files(file.path(RNA_DIR, sample_fid),
                                    pattern = "augmented_star_gene_counts\\.tsv$")[1])

sample_tsv <- read_tsv(sample_path, comment = "#",
                       show_col_types = FALSE, n_max = 10)

log_msg(sprintf("列名: %s", paste(names(sample_tsv), collapse = ", ")))

# 必要列の確認
if (!GENE_NAME_COL %in% names(sample_tsv)) {
  stop("gene_name 列が見つかりません: ", paste(names(sample_tsv), collapse = ", "))
}
if (!TPM_COL %in% names(sample_tsv)) {
  stop("tpm_unstranded 列が見つかりません: ", paste(names(sample_tsv), collapse = ", "))
}

log_msg("gene_name 列: OK")
log_msg("tpm_unstranded 列: OK")

# 全遺伝子数の確認（コメント除いた全行）
n_total_genes <- nrow(read_tsv(sample_path, comment = "#",
                               show_col_types = FALSE,
                               col_select = all_of(GENE_NAME_COL))) - 4
# augmented_starの先頭4行はN_xxx サマリー行
log_msg(sprintf("1ファイルあたりの遺伝子行数（概算）: %d", n_total_genes))

# =============================================================================
# 4. 全サンプルから全遺伝子TPMを抽出
# =============================================================================

log_msg("--- 4. 全サンプル抽出開始 ---")
log_msg(sprintf("対象: %d サンプル × 全遺伝子（20分〜40分かかります）", nrow(id_map)))

n_samples   <- nrow(id_map)
tpm_list    <- vector("list", n_samples)
sample_qc   <- vector("list", n_samples)

for (i in seq_len(n_samples)) {
  
  fid     <- id_map$rna_file_id[i]
  pair_id <- id_map$pair_id[i]
  
  if (i %% 50 == 0 || i == 1 || i == n_samples) {
    log_msg(sprintf("  処理中: %d / %d  (pair_id: %s)", i, n_samples, pair_id))
  }
  
  # ファイルパスの特定
  fdir <- file.path(RNA_DIR, fid)
  if (!dir.exists(fdir)) {
    sample_qc[[i]] <- tibble(pair_id = pair_id, rna_file_id = fid,
                             status = "dir_missing", n_genes = NA_integer_,
                             n_duplicated = NA_integer_)
    next
  }
  
  tsv_files <- list.files(fdir,
                          pattern = "augmented_star_gene_counts\\.tsv$",
                          full.names = TRUE)
  if (length(tsv_files) == 0) {
    sample_qc[[i]] <- tibble(pair_id = pair_id, rna_file_id = fid,
                             status = "file_missing", n_genes = NA_integer_,
                             n_duplicated = NA_integer_)
    next
  }
  
  # 読み込み（gene_name + tpm_unstranded の2列のみ取得→高速化）
  df <- tryCatch(
    read_tsv(tsv_files[1],
             comment        = "#",
             col_types       = cols(.default = "c"),
             show_col_types  = FALSE,
             col_select      = all_of(c(GENE_NAME_COL, TPM_COL))),
    error = function(e) {
      log_msg(sprintf("  [ERROR] %s: %s", pair_id, conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(df)) {
    sample_qc[[i]] <- tibble(pair_id = pair_id, rna_file_id = fid,
                             status = "parse_error", n_genes = NA_integer_,
                             n_duplicated = NA_integer_)
    next
  }
  
  # N_xxx サマリー行を除外（augmented_starの先頭4行）
  df <- df %>% filter(!grepl("^N_", .[[GENE_NAME_COL]]))
  
  # gene_nameがNAの行を除外
  df <- df %>% filter(!is.na(.data[[GENE_NAME_COL]]),
                      .data[[GENE_NAME_COL]] != "")
  
  # TPMを数値変換
  df <- df %>%
    mutate(tpm_val = suppressWarnings(as.numeric(.data[[TPM_COL]])))
  
  # 重複gene_name → 最大TPMを採用（step06と同ルール）
  n_before <- nrow(df)
  df <- df %>%
    group_by(.data[[GENE_NAME_COL]]) %>%
    slice_max(tpm_val, n = 1, with_ties = FALSE) %>%
    ungroup()
  n_dup <- n_before - nrow(df)
  
  # named vectorとして保存（名前=gene_name、値=TPM）
  tpm_vec        <- df$tpm_val
  names(tpm_vec) <- df[[GENE_NAME_COL]]
  tpm_list[[i]]  <- tpm_vec
  names(tpm_list)[i] <- pair_id
  
  sample_qc[[i]] <- tibble(
    pair_id      = pair_id,
    rna_file_id  = fid,
    status       = "ok",
    n_genes      = length(tpm_vec),
    n_duplicated = n_dup
  )
}

# NULLを除去
ok_idx    <- !sapply(tpm_list, is.null)
tpm_list  <- tpm_list[ok_idx]
n_ok      <- length(tpm_list)

log_msg(sprintf("読み込み完了: %d / %d サンプル成功", n_ok, n_samples))

# =============================================================================
# 5. 全サンプル共通の遺伝子セットでwide行列を構築
# =============================================================================

log_msg("--- 5. wide行列の構築 ---")

# 全サンプル共通の遺伝子（intersect）で行列化
# ※ 共通でない遺伝子はNA扱いにするunionも選択肢だが、
#    ssGSEA/ESTIMATEの品質担保のためintersectを使用
all_genes <- Reduce(intersect, lapply(tpm_list, names))
log_msg(sprintf("全サンプル共通の遺伝子数: %d", length(all_genes)))

# 行列化（行=サンプル、列=遺伝子）
log_msg("  行列化中（数分かかります）...")
tpm_mat <- do.call(rbind,
                   lapply(tpm_list, function(x) x[all_genes]))
rownames(tpm_mat) <- names(tpm_list)

log_msg(sprintf("TPM行列サイズ: %d サンプル × %d 遺伝子",
                nrow(tpm_mat), ncol(tpm_mat)))

# log2(TPM+1)変換
log2tpm_mat <- log2(tpm_mat + 1)

# =============================================================================
# 6. LAG3照合（既存コホートとの整合性確認）
# =============================================================================

log_msg("--- 6. LAG3照合（品質確認）---")

if ("LAG3" %in% colnames(log2tpm_mat)) {
  lag3_check <- cohort %>%
    filter(pair_id %in% rownames(log2tpm_mat)) %>%
    select(pair_id, LAG3_log2tpm) %>%
    mutate(
      LAG3_from_matrix = log2tpm_mat[pair_id, "LAG3"],
      diff = abs(LAG3_log2tpm - LAG3_from_matrix)
    )
  
  max_diff  <- max(lag3_check$diff, na.rm = TRUE)
  mean_diff <- mean(lag3_check$diff, na.rm = TRUE)
  log_msg(sprintf("LAG3照合: n=%d, 最大誤差=%.6f, 平均誤差=%.6f",
                  nrow(lag3_check), max_diff, mean_diff))
  
  if (max_diff < 0.01) {
    log_msg("[OK] LAG3値が既存コホートとほぼ一致しています")
  } else {
    log_msg("[WARNING] LAG3値に乖離があります。遺伝子ID/列の確認を推奨")
  }
} else {
  log_msg("[WARNING] LAG3列が行列に存在しません")
}

# =============================================================================
# 7. データフレーム化してCSV保存
# =============================================================================

log_msg("--- 7. CSV保存 ---")
log_msg("  （大きなファイルのため数分かかります）")

# log2TPM wide CSV（行=pair_id、列=遺伝子）
log2tpm_df <- as.data.frame(log2tpm_mat) %>%
  tibble::rownames_to_column("pair_id")

write_csv(log2tpm_df,
          file.path(OUT_DIR, "expression_full_log2tpm_wide.csv"))
log_msg(sprintf("保存: expression_full_log2tpm_wide.csv  (%d行 × %d列)",
                nrow(log2tpm_df), ncol(log2tpm_df)))

# TPM wide CSV
tpm_df <- as.data.frame(tpm_mat) %>%
  tibble::rownames_to_column("pair_id")

write_csv(tpm_df,
          file.path(OUT_DIR, "expression_full_tpm_wide.csv"))
log_msg(sprintf("保存: expression_full_tpm_wide.csv  (%d行 × %d列)",
                nrow(tpm_df), ncol(tpm_df)))

# 遺伝子リスト
gene_list_df <- tibble(gene_name = all_genes)
write_csv(gene_list_df,
          file.path(OUT_DIR, "step27b_gene_list.csv"))
log_msg(sprintf("保存: step27b_gene_list.csv  (%d遺伝子)", nrow(gene_list_df)))

# サンプルQC
sample_qc_df <- bind_rows(sample_qc)
write_csv(sample_qc_df,
          file.path(OUT_DIR, "step27b_sample_qc.csv"))
log_msg(sprintf("保存: step27b_sample_qc.csv  (%d行)", nrow(sample_qc_df)))

# =============================================================================
# 8. 完了サマリー
# =============================================================================

log_msg("=== Step 27b 完了 ===")

cat("\n====================================================\n")
cat("Step 27b 完了サマリー\n")
cat("====================================================\n")
cat(sprintf("サンプル数:    %d\n", nrow(log2tpm_mat)))
cat(sprintf("遺伝子数:      %d\n", ncol(log2tpm_mat)))
cat(sprintf("遺伝子ID形式（先頭5件）: %s\n",
            paste(head(colnames(log2tpm_mat), 5), collapse = ", ")))
cat(sprintf("コホートとの一致数: %d / %d\n",
            sum(cohort$pair_id %in% rownames(log2tpm_mat)), nrow(cohort)))

qc_summary <- sample_qc_df %>% count(status)
cat("サンプルQCサマリー:\n")
for (i in seq_len(nrow(qc_summary))) {
  cat(sprintf("  %-15s: %d件\n", qc_summary$status[i], qc_summary$n[i]))
}

cat("\n以下をコピーして返してください:\n")
cat(sprintf("  サンプル数: %d\n", nrow(log2tpm_mat)))
cat(sprintf("  遺伝子数: %d\n", ncol(log2tpm_mat)))
cat(sprintf("  遺伝子ID形式（先頭5件）: %s\n",
            paste(head(colnames(log2tpm_mat), 5), collapse = ", ")))
cat(sprintf("  コホートとの一致数: %d / %d\n",
            sum(cohort$pair_id %in% rownames(log2tpm_mat)), nrow(cohort)))
cat(sprintf("  LAG3照合結果: \n"))
cat("  QCサマリー（ok/error件数）:\n")
for (i in seq_len(nrow(qc_summary))) {
  cat(sprintf("    %-15s: %d件\n", qc_summary$status[i], qc_summary$n[i]))
}
cat("====================================================\n")
