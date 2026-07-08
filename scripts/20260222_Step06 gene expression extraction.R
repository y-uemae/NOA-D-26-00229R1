# =============================================================================
# step06_gene_expression_extraction.R
# GBM/Glioma TP53×LAG3 再解析 - Step 06: 解析対象遺伝子セット 発現量一括抽出
#
# 目的:
#   WXS-RNA ペアリング済みの RNA file_id に対応する
#   augmented_star_gene_counts.tsv から、解析対象遺伝子セット（28遺伝子）の
#   TPM を一括抽出し log2(TPM+1) に変換する。
#
# -----------------------------------------------------------------------
# 解析対象遺伝子セット（28遺伝子・固定）
#   MHC-I 提示経路   (8): B2M, TAP1, TAP2, TAPBP, HLA-A, HLA-B, HLA-C, NLRC5
#   IFN-γ シグナル   (3): STAT1, IRF1, IRF9
#   ケモカイン        (3): CXCL9, CXCL10, CXCL11
#   GBP ファミリー   (4): GBP1, GBP2, GBP4, GBP5
#   免疫抑制         (1): IDO1
#   T 細胞マーカー   (5): CD3D, CD3E, CD3G, CD8A, CD8B
#   細胞傷害分子     (3): GZMA, GZMB, PRF1
#   免疫チェックポイント (1): LAG3
#   合計: 8+3+3+4+1+5+3+1 = 28遺伝子
#
# -----------------------------------------------------------------------
# 抽出仕様（固定）
#   マッチ列（優先）: gene_name
#   フォールバック  : gene_name 列が存在しない/全 NA の場合 → gene_id（ENSG…）で照合
#                    ※ターゲットは遺伝子記号のため通常は gene_name で完結
#   match_column_used: 使用したマッチ列名を gene_extraction_log_file.csv に記録
#   値列            : tpm_unstranded
#   変換            : log2(tpm_unstranded + 1)
#   コメント行      : # で始まる行をスキップ（先頭: # gene-model: GENCODE v36）
#
# 重複行処理ルール（固定）:
#   同一 gene_name が複数行 → 最大 TPM の行を採用
#   status = "duplicated", selected_value_rule = "max_tpm_selected"
#   n_hits（マッチ行数）と gene_id（全候補）を note に記録
#
# -----------------------------------------------------------------------
# 出力
#   06_gene_expression/
#     gene_expression_wide.csv       : 解析用（1行=ペア、列=遺伝子）
#     gene_extraction_log_long.csv   : 監査用（1行=file_id×gene）
#     gene_extraction_log_file.csv   : ファイル単位の読み込みステータス
#     step06_log.txt
#
# 作成日: 2026-02-21
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
RNA_DL_DIR <- file.path(BASE_DIR, "data/raw/GDC/glioma/RNAseq/downloads")
OUT_DIR    <- file.path(RESULT_DIR, "06_gene_expression")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 実行識別子（監査用・wide テーブルに付加）
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")

# -----------------------------------------------------------------------
# 解析対象遺伝子セット（28遺伝子・固定）
# カテゴリコメントはカウント根拠として残す
# -----------------------------------------------------------------------
TARGET_GENES <- c(
  # MHC-I 提示経路 (8)
  "B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5",
  # IFN-γ シグナル (3)
  "STAT1", "IRF1", "IRF9",
  # ケモカイン (3)
  "CXCL9", "CXCL10", "CXCL11",
  # GBP ファミリー (4)
  "GBP1", "GBP2", "GBP4", "GBP5",
  # 免疫抑制 (1)
  "IDO1",
  # T 細胞マーカー (5)
  "CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
  # 細胞傷害分子 (3)
  "GZMA", "GZMB", "PRF1",
  # 免疫チェックポイント (1)
  "LAG3"
)
TARGET_GENES <- unique(TARGET_GENES)   # 重複除去（念のため）
stopifnot("TARGET_GENES は 28 遺伝子である必要があります" = length(TARGET_GENES) == 28)

TPM_COL        <- "tpm_unstranded"
GENE_NAME_COL  <- "gene_name"
GENE_ID_COL    <- "gene_id"   # フォールバック用（ENSG番号）

# ターゲット遺伝子の ENSG ID マッピング（gene_name フォールバック用）
# ※ GENCODE v36 に基づく主要 ENSG を記載。フォールバック発動時のみ使用。
TARGET_GENE_ENSG <- c(
  B2M    = "ENSG00000166710", TAP1   = "ENSG00000168394",
  TAP2   = "ENSG00000141505", TAPBP  = "ENSG00000231925",
  `HLA-A`= "ENSG00000206503", `HLA-B`= "ENSG00000234745",
  `HLA-C`= "ENSG00000204525", NLRC5  = "ENSG00000140853",
  STAT1  = "ENSG00000115415", IRF1   = "ENSG00000125347",
  IRF9   = "ENSG00000213928", CXCL9  = "ENSG00000138755",
  CXCL10 = "ENSG00000169245", CXCL11 = "ENSG00000169248",
  GBP1   = "ENSG00000117228", GBP2   = "ENSG00000163703",
  GBP4   = "ENSG00000162654", GBP5   = "ENSG00000154451",
  IDO1   = "ENSG00000131203", CD3D   = "ENSG00000167286",
  CD3E   = "ENSG00000198851", CD3G   = "ENSG00000160654",
  CD8A   = "ENSG00000153563", CD8B   = "ENSG00000172116",
  GZMA   = "ENSG00000145649", GZMB   = "ENSG00000100453",
  PRF1   = "ENSG00000180644", LAG3   = "ENSG00000089692"
)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step06_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 06: 解析対象遺伝子セット 発現量一括抽出 開始 ===")
log_msg(sprintf("run_id: %s", RUN_ID))
log_msg(sprintf("対象遺伝子数: %d（カウント確認済み）", length(TARGET_GENES)))
log_msg(sprintf("遺伝子リスト: %s", paste(TARGET_GENES, collapse = ", ")))
log_msg(sprintf("抽出列: %s  変換: log2(TPM+1)", TPM_COL))
log_msg("重複行ルール: 同一gene_nameが複数行 → 最大TPMを採用（max_tpm_selected）")
log_msg("マッチ列: gene_name 優先 → 欠損/全NA時は gene_id（ENSG）でフォールバック")

# =============================================================================
# 2. ペアリング済みリストの読み込み
# =============================================================================

log_msg("--- ペアリング済みリスト読み込み ---")

read_pairs <- function(fpath, label) {
  if (!file.exists(fpath)) {
    log_msg(sprintf("  WARNING [%s]: ファイルが存在しません: %s", label, fpath))
    return(NULL)
  }
  df <- read_csv(fpath, show_col_types = FALSE)
  df$source <- label
  log_msg(sprintf("  [%s]: %d件", label, nrow(df)))
  df
}

pairs_tcga  <- read_pairs(
  file.path(RESULT_DIR, "05a_wxs_rna_match_tcga/05a_tcga_pairs_final.csv"), "TCGA")
pairs_cptac <- read_pairs(
  file.path(RESULT_DIR, "05b_wxs_rna_match_cptac_hcmi/05b_cptac_hcmi_pairs_final.csv"),
  "CPTAC_HCMI")

all_pairs <- bind_rows(pairs_tcga, pairs_cptac) %>%
  mutate(pair_id = row_number())

log_msg(sprintf("全ペア数: %d件", nrow(all_pairs)))

rna_file_ids <- unique(all_pairs$rna_file_id[!is.na(all_pairs$rna_file_id)])
log_msg(sprintf("抽出対象 RNA file_id（重複除去）: %d件", length(rna_file_ids)))

# =============================================================================
# 3. 1ファイルから対象遺伝子セットを抽出する関数
# =============================================================================

extract_genes_from_file <- function(file_id) {
  
  dir_path <- file.path(RNA_DL_DIR, file_id)
  
  # ---- ディレクトリ・ファイル存在確認 ----
  if (!dir.exists(dir_path)) {
    return(list(
      file_status      = "file_missing",
      file_note        = sprintf("Directory not found: %s", dir_path),
      match_column_used = NA_character_,
      long_rows        = make_empty_long(file_id, "file_missing")
    ))
  }
  
  tsv_files <- list.files(dir_path,
                          pattern = "augmented_star_gene_counts\\.tsv$",
                          full.names = TRUE, recursive = FALSE)
  if (length(tsv_files) == 0) {
    return(list(
      file_status      = "file_missing",
      file_note        = sprintf("No TSV found in: %s", dir_path),
      match_column_used = NA_character_,
      long_rows        = make_empty_long(file_id, "file_missing")
    ))
  }
  
  tsv_path <- tsv_files[1]
  
  # ---- ファイル読み込み ----
  df <- tryCatch(
    read_tsv(tsv_path,
             comment       = "#",
             col_types     = cols(.default = "c"),
             show_col_types = FALSE),
    error = function(e) e
  )
  
  if (inherits(df, "error")) {
    return(list(
      file_status      = "parse_error",
      file_note        = conditionMessage(df),
      match_column_used = NA_character_,
      long_rows        = make_empty_long(file_id, "parse_error")
    ))
  }
  
  # ---- マッチ列の決定（gene_name 優先 → gene_id フォールバック） ----
  use_gene_name <- (GENE_NAME_COL %in% names(df)) &&
    (sum(!is.na(df[[GENE_NAME_COL]])) > 0)
  
  if (use_gene_name) {
    match_col  <- GENE_NAME_COL
    match_vals <- TARGET_GENES          # 記号でマッチ
  } else {
    # フォールバック: gene_id（ENSG）で照合
    match_col  <- GENE_ID_COL
    match_vals <- TARGET_GENE_ENSG      # ENSG番号ベクター
    log_msg(sprintf("  INFO [%s]: gene_name 欠損 → gene_id でフォールバック", file_id))
  }
  
  # TPM 列確認
  if (!TPM_COL %in% names(df)) {
    msg <- sprintf("Missing column '%s'. Available: %s",
                   TPM_COL, paste(names(df), collapse = ", "))
    return(list(
      file_status      = "parse_error",
      file_note        = msg,
      match_column_used = match_col,
      long_rows        = make_empty_long(file_id, "parse_error")
    ))
  }
  
  if (!match_col %in% names(df)) {
    msg <- sprintf("Match column '%s' not found. Available: %s",
                   match_col, paste(names(df), collapse = ", "))
    return(list(
      file_status      = "parse_error",
      file_note        = msg,
      match_column_used = match_col,
      long_rows        = make_empty_long(file_id, "parse_error")
    ))
  }
  
  # ---- 対象遺伝子を一括フィルタ ----
  df_target <- df %>% filter(.data[[match_col]] %in% match_vals)
  
  # ---- 遺伝子ごとに 1 値を決定 ----
  long_rows <- lapply(TARGET_GENES, function(gene) {
    
    # フォールバック時は ENSG で検索
    search_val <- if (use_gene_name) gene else TARGET_GENE_ENSG[gene]
    
    rows <- df_target %>%
      filter(!is.na(search_val) & .data[[match_col]] == search_val)
    n <- nrow(rows)
    
    if (n == 0) {
      return(tibble(
        file_id             = file_id,
        gene                = gene,
        tpm                 = NA_real_,
        log2tpm             = NA_real_,
        status              = "not_found",
        n_hits              = 0L,
        selected_value_rule = NA_character_,
        note                = NA_character_,
        match_column_used   = match_col
      ))
    }
    
    tpm_vals <- suppressWarnings(as.numeric(rows[[TPM_COL]]))
    valid    <- !is.na(tpm_vals)
    
    if (!any(valid)) {
      return(tibble(
        file_id             = file_id,
        gene                = gene,
        tpm                 = NA_real_,
        log2tpm             = NA_real_,
        status              = "parse_error",
        n_hits              = n,
        selected_value_rule = NA_character_,
        note                = sprintf("TPM not numeric: '%s'", rows[[TPM_COL]][1]),
        match_column_used   = match_col
      ))
    }
    
    tpm_use <- max(tpm_vals[valid])
    
    if (n == 1) {
      status_val <- "ok"
      rule_val   <- "single_row"
      note_val   <- NA_character_
    } else {
      # 重複: 最大 TPM を採用（max_tpm_selected）
      status_val <- "duplicated"
      rule_val   <- "max_tpm_selected"
      gene_ids   <- if (GENE_ID_COL %in% names(rows)) {
        paste(rows[[GENE_ID_COL]], collapse = ",")
      } else NA_character_
      note_val   <- sprintf("n_hits=%d; gene_ids=[%s]; max_tpm=%.4f",
                            n, gene_ids, tpm_use)
    }
    
    tibble(
      file_id             = file_id,
      gene                = gene,
      tpm                 = tpm_use,
      log2tpm             = log2(tpm_use + 1),
      status              = status_val,
      n_hits              = n,
      selected_value_rule = rule_val,
      note                = note_val,
      match_column_used   = match_col
    )
  })
  
  list(
    file_status      = "ok",
    file_note        = NA_character_,
    match_column_used = match_col,
    long_rows        = bind_rows(long_rows)
  )
}

# 空の long_rows を作るヘルパー
make_empty_long <- function(file_id, status_val) {
  tibble(
    file_id             = file_id,
    gene                = TARGET_GENES,
    tpm                 = NA_real_,
    log2tpm             = NA_real_,
    status              = status_val,
    n_hits              = NA_integer_,
    selected_value_rule = NA_character_,
    note                = NA_character_,
    match_column_used   = NA_character_
  )
}

# =============================================================================
# 4. 全 file_id に対して抽出実行
# =============================================================================

log_msg("--- 遺伝子発現量抽出開始 ---")
n_total <- length(rna_file_ids)

file_status_rows <- vector("list", n_total)
all_long_rows    <- vector("list", n_total)

for (i in seq_along(rna_file_ids)) {
  fid    <- rna_file_ids[i]
  result <- extract_genes_from_file(fid)
  
  file_status_rows[[i]] <- tibble(
    file_id           = fid,
    file_status       = result$file_status,
    match_column_used = result$match_column_used,
    file_note         = result$file_note
  )
  all_long_rows[[i]] <- result$long_rows
  
  if (i %% 100 == 0 || i == n_total) {
    n_ok <- sum(sapply(file_status_rows[seq_len(i)],
                       function(r) !is.null(r) && r$file_status == "ok"))
    log_msg(sprintf("  進捗: %d/%d 完了 (file_ok=%d, 問題=%d)",
                    i, n_total, n_ok, i - n_ok))
  }
}

# =============================================================================
# 5. 結果整形
# =============================================================================

log_msg("--- 結果整形 ---")

file_log <- bind_rows(file_status_rows)
long_log <- bind_rows(all_long_rows)

# フォールバック使用状況をログ
fallback_n <- sum(file_log$match_column_used == GENE_ID_COL, na.rm = TRUE)
if (fallback_n > 0) {
  log_msg(sprintf("  INFO: gene_id フォールバック使用ファイル: %d件", fallback_n))
  fb_files <- file_log %>% filter(match_column_used == GENE_ID_COL) %>% pull(file_id)
  log_msg(sprintf("    %s", paste(head(fb_files, 10), collapse = ", ")))
} else {
  log_msg("  INFO: 全ファイルで gene_name マッチ（フォールバックなし）")
}

# ステータス集計
log_msg("ファイル読み込みステータス:")
file_log %>% count(file_status) %>% arrange(desc(n)) %>%
  { for (i in seq_len(nrow(.))) log_msg(sprintf("  %-15s: %d件", .[[1]][i], .[[2]][i])) }

log_msg("遺伝子抽出ステータス（全 file_id × gene）:")
long_log %>% count(status) %>% arrange(desc(n)) %>%
  { for (i in seq_len(nrow(.))) log_msg(sprintf("  %-15s: %d件", .[[1]][i], .[[2]][i])) }

# 問題ファイル詳細（最大20件）
prob_files <- file_log %>% filter(file_status != "ok")
if (nrow(prob_files) > 0) {
  log_msg(sprintf("  問題ファイル（%d件）:", nrow(prob_files)))
  for (i in seq_len(min(nrow(prob_files), 20))) {
    log_msg(sprintf("    [%s] %s: %s",
                    prob_files$file_status[i],
                    prob_files$file_id[i],
                    prob_files$file_note[i]))
  }
}

# =============================================================================
# 6. wide 形式への変換
# =============================================================================

log_msg("--- wide 形式への変換 ---")

wide_tpm <- long_log %>%
  select(file_id, gene, tpm) %>%
  pivot_wider(names_from = gene, values_from = tpm,   names_glue = "{gene}_tpm")

wide_log2 <- long_log %>%
  select(file_id, gene, log2tpm) %>%
  pivot_wider(names_from = gene, values_from = log2tpm, names_glue = "{gene}_log2tpm")

gene_status_summary <- long_log %>%
  group_by(file_id) %>%
  summarise(
    n_ok          = sum(status == "ok"),
    n_not_found   = sum(status == "not_found"),
    n_duplicated  = sum(status == "duplicated"),
    n_parse_error = sum(status == "parse_error"),
    .groups = "drop"
  )

# ペアテーブルに結合 + run_id を付加
wide_expr <- all_pairs %>%
  left_join(wide_tpm,            by = c("rna_file_id" = "file_id")) %>%
  left_join(wide_log2,           by = c("rna_file_id" = "file_id")) %>%
  left_join(gene_status_summary, by = c("rna_file_id" = "file_id")) %>%
  mutate(run_id = RUN_ID)

log_msg(sprintf("wide テーブル: %d行 × %d列 (run_id=%s)",
                nrow(wide_expr), ncol(wide_expr), RUN_ID))

# long にサンプルメタ情報を付加
sample_key <- all_pairs %>%
  select(rna_file_id, grade, source,
         any_of(c("wxs_sample_id", "rna_sample_id",
                  "case_barcode", "case_id", "sample_id"))) %>%
  distinct(rna_file_id, .keep_all = TRUE)

long_log_meta <- long_log %>%
  left_join(sample_key, by = c("file_id" = "rna_file_id"))

# =============================================================================
# 7. 出力
# =============================================================================

log_msg("--- ファイル出力 ---")

write_csv(wide_expr,     file.path(OUT_DIR, "gene_expression_wide.csv"))
log_msg(sprintf("保存: gene_expression_wide.csv (%d行 × %d列)",
                nrow(wide_expr), ncol(wide_expr)))

write_csv(long_log_meta, file.path(OUT_DIR, "gene_extraction_log_long.csv"))
log_msg(sprintf("保存: gene_extraction_log_long.csv (%d行)",
                nrow(long_log_meta)))

write_csv(file_log,      file.path(OUT_DIR, "gene_extraction_log_file.csv"))
log_msg(sprintf("保存: gene_extraction_log_file.csv (%d行)", nrow(file_log)))

# =============================================================================
# 8. サマリー
# =============================================================================

log_msg("=== Step 06 サマリー ===")

# LAG3 を代表として Grade × Source 別統計
lag3_ok <- long_log_meta %>% filter(gene == "LAG3", status == "ok")
if (nrow(lag3_ok) > 0) {
  grade_stats <- lag3_ok %>%
    group_by(grade, source) %>%
    summarise(n      = n(),
              median = round(median(log2tpm, na.rm = TRUE), 3),
              mean   = round(mean(log2tpm,   na.rm = TRUE), 3),
              min    = round(min(log2tpm,    na.rm = TRUE), 3),
              max    = round(max(log2tpm,    na.rm = TRUE), 3),
              .groups = "drop")
  log_msg("LAG3 log2(TPM+1) 基本統計（Grade × Source）:")
  for (i in seq_len(nrow(grade_stats))) {
    r <- grade_stats[i, ]
    log_msg(sprintf("  [%-7s][%-10s] n=%d, median=%.3f, mean=%.3f, range=[%.3f, %.3f]",
                    r$grade, r$source, r$n, r$median, r$mean, r$min, r$max))
  }
}

# 遺伝子別ステータス一覧
gene_nf <- long_log %>%
  group_by(gene) %>%
  summarise(n_ok        = sum(status == "ok"),
            n_not_found = sum(status == "not_found"),
            n_dup       = sum(status == "duplicated"),
            n_error     = sum(status %in% c("parse_error","file_missing")),
            .groups = "drop") %>%
  arrange(desc(n_not_found), desc(n_dup))

log_msg("遺伝子別 抽出ステータス（not_found 多い順）:")
for (i in seq_len(nrow(gene_nf))) {
  r    <- gene_nf[i, ]
  flag <- if (r$n_not_found > 0 || r$n_dup > 0) " ←要確認" else ""
  log_msg(sprintf("  %-10s ok=%d, not_found=%d, dup=%d, error=%d%s",
                  r$gene, r$n_ok, r$n_not_found, r$n_dup, r$n_error, flag))
}

log_msg("=== Step 06: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 06 完了\n")
cat(sprintf("run_id: %s\n", RUN_ID))
cat(sprintf("対象遺伝子: 28個\n"))
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    gene_expression_wide.csv        ← 解析用（1行=ペア、列=遺伝子、run_id付き）\n")
cat("    gene_extraction_log_long.csv    ← 監査用（1行=file_id×gene、selected_value_rule付き）\n")
cat("    gene_extraction_log_file.csv    ← ファイル単位ステータス（match_column_used付き）\n")
cat("    step06_log.txt\n")
cat("============================\n")
