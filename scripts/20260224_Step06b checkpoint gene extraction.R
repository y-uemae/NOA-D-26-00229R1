# =============================================================================
# step06b_checkpoint_gene_extraction.R
# GBM/Glioma TP53×LAG3 解析 - Step 06b: チェックポイント遺伝子 追加抽出
#
# 目的:
#   チェックポイント特異性解析（Step16a予定）のために、
#   Step06の28遺伝子セットに含まれていない免疫チェックポイント分子6遺伝子を
#   GDC（生TSVファイル）とGLASS（発現マトリックス）から追加抽出する。
#
# 追加対象遺伝子（6遺伝子・事前固定）:
#   PDCD1   (PD-1)
#   CTLA4
#   TIGIT
#   HAVCR2  (TIM-3)
#   CD274   (PD-L1)
#   PDCD1LG2(PD-L2)
#
# 処理内容:
#   [GDC]   Step06と同じロジックで augmented_star_gene_counts.tsv から抽出
#            → gene_expression_wide.csv に6列を追加して上書き保存
#   [GLASS] data_mrna_seq_tpm.txt（横持ちマトリックス）から列抽出
#            → glass_final_cohort.csv に6列×2形式（_tpm, _log2tpm）を追加して上書き保存
#            → glass_final_cohort_all101.csv も同様に更新
#
# 出力:
#   06_gene_expression/
#     gene_expression_wide.csv          ★上書き（6列追加）
#     step06b_gdc_extraction_log.csv    追加抽出監査ログ（GDC用）
#     step06b_log.txt
#   05c_glass/
#     glass_final_cohort.csv            ★上書き（6列×2追加）
#     glass_final_cohort_all101.csv     ★上書き（6列×2追加）
#
# 注意:
#   Step06の gene_expression_wide.csv が存在することが前提。
#   GLASSの glass_final_cohort.csv / glass_final_cohort_all101.csv が存在することが前提。
#   既存列（_tpm, _log2tpm）が存在する場合は上書き（再実行安全）。
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR    <- here::here()
RESULT_DIR  <- file.path(BASE_DIR, "results/TP53/20260221")
RNA_DL_DIR  <- file.path(BASE_DIR, "data/raw/GDC/glioma/RNAseq/downloads")
GLASS_DIR   <- file.path(BASE_DIR, "data/raw/external_validation/difg_glass")
GDC_OUT_DIR <- file.path(RESULT_DIR, "06_gene_expression")
GLASS_OUT_DIR <- file.path(RESULT_DIR, "05c_glass")

# 追加対象遺伝子（6遺伝子・事前固定）
CHECKPOINT_GENES <- c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2")

# GENCODE v36 ENSG（gene_nameフォールバック用）
CHECKPOINT_ENSG <- c(
  PDCD1    = "ENSG00000188389",
  CTLA4    = "ENSG00000163599",
  TIGIT    = "ENSG00000181847",
  HAVCR2   = "ENSG00000135077",
  CD274    = "ENSG00000120217",
  PDCD1LG2 = "ENSG00000197646"
)

TPM_COL       <- "tpm_unstranded"
GENE_NAME_COL <- "gene_name"
GENE_ID_COL   <- "gene_id"

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(GDC_OUT_DIR, "step06b_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 06b: チェックポイント遺伝子 追加抽出 開始 ===")
log_msg(sprintf("追加遺伝子: %s", paste(CHECKPOINT_GENES, collapse = ", ")))

# =============================================================================
# 2. GDC 側：gene_expression_wide.csv からfile_idリストを取得
# =============================================================================

log_msg("--- [GDC] gene_expression_wide.csv 読み込み ---")

wide_path <- file.path(GDC_OUT_DIR, "gene_expression_wide.csv")
if (!file.exists(wide_path)) {
  log_msg("ERROR: gene_expression_wide.csv が見つかりません。Step06を先に実行してください。")
  close(log_con)
  stop("gene_expression_wide.csv が存在しません: ", wide_path)
}

wide_existing <- read_csv(wide_path, show_col_types = FALSE)
log_msg(sprintf("既存wide: %d行 × %d列", nrow(wide_existing), ncol(wide_existing)))

# すでに追加済みの列がある場合は削除して再抽出（再実行安全）
already_cols <- intersect(
  names(wide_existing),
  c(paste0(CHECKPOINT_GENES, "_tpm"), paste0(CHECKPOINT_GENES, "_log2tpm"))
)
if (length(already_cols) > 0) {
  log_msg(sprintf("INFO: 既存列を削除して再抽出します: %s", paste(already_cols, collapse = ", ")))
  wide_existing <- wide_existing %>% select(-all_of(already_cols))
}

rna_file_ids <- unique(wide_existing$rna_file_id[!is.na(wide_existing$rna_file_id)])
log_msg(sprintf("抽出対象 rna_file_id: %d件", length(rna_file_ids)))

# =============================================================================
# 3. GDC 側：1ファイルからチェックポイント遺伝子を抽出する関数
# =============================================================================

make_empty_long_cp <- function(file_id, status_val) {
  tibble(
    file_id             = file_id,
    gene                = CHECKPOINT_GENES,
    tpm                 = NA_real_,
    log2tpm             = NA_real_,
    status              = status_val,
    n_hits              = NA_integer_,
    selected_value_rule = NA_character_,
    note                = NA_character_,
    match_column_used   = NA_character_
  )
}

extract_checkpoint_from_file <- function(file_id) {
  dir_path  <- file.path(RNA_DL_DIR, file_id)
  
  if (!dir.exists(dir_path)) {
    return(list(file_status = "file_missing",
                long_rows   = make_empty_long_cp(file_id, "file_missing")))
  }
  
  tsv_files <- list.files(dir_path,
                          pattern    = "augmented_star_gene_counts\\.tsv$",
                          full.names = TRUE, recursive = FALSE)
  if (length(tsv_files) == 0) {
    return(list(file_status = "file_missing",
                long_rows   = make_empty_long_cp(file_id, "file_missing")))
  }
  
  df <- tryCatch(
    read_tsv(tsv_files[1], comment = "#",
             col_types = cols(.default = "c"), show_col_types = FALSE),
    error = function(e) e
  )
  if (inherits(df, "error")) {
    return(list(file_status = "parse_error",
                long_rows   = make_empty_long_cp(file_id, "parse_error")))
  }
  
  # マッチ列決定（gene_name優先・ENSGフォールバック）
  use_gene_name <- (GENE_NAME_COL %in% names(df)) &&
    (sum(!is.na(df[[GENE_NAME_COL]])) > 0)
  match_col  <- if (use_gene_name) GENE_NAME_COL else GENE_ID_COL
  match_vals <- if (use_gene_name) CHECKPOINT_GENES else CHECKPOINT_ENSG
  
  if (!TPM_COL %in% names(df) || !match_col %in% names(df)) {
    return(list(file_status = "parse_error",
                long_rows   = make_empty_long_cp(file_id, "parse_error")))
  }
  
  df_target <- df %>% filter(.data[[match_col]] %in% match_vals)
  
  long_rows <- lapply(CHECKPOINT_GENES, function(gene) {
    search_val <- if (use_gene_name) gene else CHECKPOINT_ENSG[gene]
    rows <- df_target %>% filter(.data[[match_col]] == search_val)
    n    <- nrow(rows)
    
    if (n == 0) {
      return(tibble(file_id = file_id, gene = gene,
                    tpm = NA_real_, log2tpm = NA_real_,
                    status = "not_found", n_hits = 0L,
                    selected_value_rule = NA_character_, note = NA_character_,
                    match_column_used = match_col))
    }
    
    tpm_vals <- suppressWarnings(as.numeric(rows[[TPM_COL]]))
    valid    <- !is.na(tpm_vals)
    if (!any(valid)) {
      return(tibble(file_id = file_id, gene = gene,
                    tpm = NA_real_, log2tpm = NA_real_,
                    status = "parse_error", n_hits = n,
                    selected_value_rule = NA_character_,
                    note = sprintf("TPM not numeric: '%s'", rows[[TPM_COL]][1]),
                    match_column_used = match_col))
    }
    
    tpm_use    <- max(tpm_vals[valid])
    status_val <- if (n == 1) "ok" else "duplicated"
    rule_val   <- if (n == 1) "single_row" else "max_tpm_selected"
    note_val   <- if (n > 1) {
      gene_ids <- if (GENE_ID_COL %in% names(rows)) paste(rows[[GENE_ID_COL]], collapse = ",") else NA
      sprintf("n_hits=%d; gene_ids=[%s]; max_tpm=%.4f", n, gene_ids, tpm_use)
    } else NA_character_
    
    tibble(file_id = file_id, gene = gene,
           tpm = tpm_use, log2tpm = log2(tpm_use + 1),
           status = status_val, n_hits = n,
           selected_value_rule = rule_val, note = note_val,
           match_column_used = match_col)
  })
  
  list(file_status = "ok", long_rows = bind_rows(long_rows))
}

# =============================================================================
# 4. GDC 側：全ファイルに対して抽出実行
# =============================================================================

log_msg("--- [GDC] 全ファイル抽出開始 ---")
n_total       <- length(rna_file_ids)
all_long_rows <- vector("list", n_total)
file_statuses <- character(n_total)

for (i in seq_along(rna_file_ids)) {
  fid              <- rna_file_ids[i]
  result           <- extract_checkpoint_from_file(fid)
  all_long_rows[[i]] <- result$long_rows
  file_statuses[i]   <- result$file_status
  
  if (i %% 100 == 0 || i == n_total) {
    n_ok <- sum(file_statuses[seq_len(i)] == "ok")
    log_msg(sprintf("  進捗: %d/%d (ok=%d, 問題=%d)", i, n_total, n_ok, i - n_ok))
  }
}

long_all <- bind_rows(all_long_rows)

# ステータス集計
log_msg("遺伝子別抽出ステータス（not_found多い順）:")
gene_summary <- long_all %>%
  group_by(gene) %>%
  summarise(n_ok        = sum(status == "ok"),
            n_not_found = sum(status == "not_found"),
            n_dup       = sum(status == "duplicated"),
            n_error     = sum(status %in% c("parse_error","file_missing")),
            .groups = "drop") %>%
  arrange(desc(n_not_found), desc(n_dup))

for (i in seq_len(nrow(gene_summary))) {
  r <- gene_summary[i,]
  log_msg(sprintf("  %-12s ok=%d, not_found=%d, dup=%d, error=%d",
                  r$gene, r$n_ok, r$n_not_found, r$n_dup, r$n_error))
}

# 監査ログ保存
write_csv(long_all, file.path(GDC_OUT_DIR, "step06b_gdc_extraction_log.csv"))
log_msg("保存: step06b_gdc_extraction_log.csv")

# =============================================================================
# 5. GDC 側：wide形式に変換してgene_expression_wide.csvに列追加
# =============================================================================

log_msg("--- [GDC] wide形式変換・列追加 ---")

wide_tpm <- long_all %>%
  select(file_id, gene, tpm) %>%
  pivot_wider(names_from = gene, values_from = tpm,
              names_glue = "{gene}_tpm")

wide_log2 <- long_all %>%
  select(file_id, gene, log2tpm) %>%
  pivot_wider(names_from = gene, values_from = log2tpm,
              names_glue = "{gene}_log2tpm")

# 列の順序：_tpm と _log2tpm を遺伝子ごとにペアで並べる
new_cols_ordered <- as.vector(rbind(
  paste0(CHECKPOINT_GENES, "_tpm"),
  paste0(CHECKPOINT_GENES, "_log2tpm")
))

wide_new <- wide_tpm %>%
  left_join(wide_log2, by = c("file_id" = "file_id"))

# gene_expression_wide.csv に結合
wide_updated <- wide_existing %>%
  left_join(wide_new, by = c("rna_file_id" = "file_id"))

# 列順整理：既存列の後ろに新列を追加
existing_cols <- names(wide_existing)
add_cols      <- intersect(new_cols_ordered, names(wide_updated))
wide_updated  <- wide_updated %>% select(all_of(c(existing_cols, add_cols)))

log_msg(sprintf("wide更新後: %d行 × %d列（追加列: %d）",
                nrow(wide_updated), ncol(wide_updated), length(add_cols)))

# 追加列のNA率確認
for (col in add_cols) {
  na_rate <- mean(is.na(wide_updated[[col]])) * 100
  if (na_rate > 5) log_msg(sprintf("  WARNING: %s のNA率 %.1f%%", col, na_rate))
}

write_csv(wide_updated, wide_path)
log_msg(sprintf("上書き保存: gene_expression_wide.csv"))

# =============================================================================
# 6. GLASS 側：data_mrna_seq_tpm.txt から追加抽出
# =============================================================================

log_msg("--- [GLASS] data_mrna_seq_tpm.txt からチェックポイント遺伝子抽出 ---")

glass_rna_path <- file.path(GLASS_DIR, "data_mrna_seq_tpm.txt")
if (!file.exists(glass_rna_path)) {
  log_msg(sprintf("ERROR: GLASSのRNAファイルが見つかりません: %s", glass_rna_path))
  close(log_con)
  stop("GLASSのRNAファイルが存在しません")
}

rna_glass <- read.table(glass_rna_path, sep = "\t", header = TRUE,
                        comment.char = "#", stringsAsFactors = FALSE, quote = "")
log_msg(sprintf("GLASS RNAマトリックス: %d行 × %d列", nrow(rna_glass), ncol(rna_glass)))

# チェックポイント遺伝子を抽出
rna_cp <- rna_glass %>%
  filter(Hugo_Symbol %in% CHECKPOINT_GENES) %>%
  group_by(Hugo_Symbol) %>%
  slice(1) %>%   # 重複は先頭採用（Step05c_02と同じルール）
  ungroup()

found_genes <- rna_cp$Hugo_Symbol
missing_genes <- setdiff(CHECKPOINT_GENES, found_genes)
log_msg(sprintf("GLASS: 発見された遺伝子: %s", paste(found_genes, collapse = ", ")))
if (length(missing_genes) > 0) {
  log_msg(sprintf("GLASS: WARNING 見つからない遺伝子: %s", paste(missing_genes, collapse = ", ")))
}

# =============================================================================
# 7. GLASS 側：glass_final_cohort.csv と _all101.csv に列追加
# =============================================================================

update_glass_cohort <- function(cohort_path, rna_cp_df, label) {
  if (!file.exists(cohort_path)) {
    log_msg(sprintf("WARNING: %s が見つかりません（スキップ）: %s", label, cohort_path))
    return(invisible(NULL))
  }
  
  cohort <- read_csv(cohort_path, show_col_types = FALSE)
  log_msg(sprintf("[GLASS %s] 読み込み: %d行 × %d列", label, nrow(cohort), ncol(cohort)))
  
  # mappingからrna_col_dotを取得（sample_idとrna_col_dotの対応）
  mapping_path <- file.path(GLASS_OUT_DIR, "step05c_01_glass_mapping_included.csv")
  if (!file.exists(mapping_path)) {
    log_msg(sprintf("WARNING: GLASSのmappingファイルが見つかりません: %s", mapping_path))
    log_msg("  → pair_idをそのままdotスタイルに変換して試みます")
    # pair_id（sample_id）からrna_col_dotへの変換を試みる
    sample_ids   <- cohort$pair_id
    rna_col_dots <- str_replace_all(sample_ids, "-", ".")
  } else {
    mapping_df   <- read_csv(mapping_path, show_col_types = FALSE)
    # pair_id = sample_id のケース
    id_col       <- if ("pair_id" %in% names(cohort)) "pair_id" else "sample_id"
    sample_ids   <- cohort[[id_col]]
    rna_col_dots <- mapping_df$rna_col_dot[match(sample_ids, mapping_df$sample_id)]
  }
  
  # 既存のチェックポイント列を削除（再実行安全）
  drop_cols <- intersect(names(cohort),
                         c(paste0(CHECKPOINT_GENES, "_tpm"),
                           paste0(CHECKPOINT_GENES, "_log2tpm")))
  if (length(drop_cols) > 0) {
    cohort <- cohort %>% select(-all_of(drop_cols))
    log_msg(sprintf("  既存列を削除して再抽出: %s", paste(drop_cols, collapse = ", ")))
  }
  
  # サンプルごとにTPMを取得してlong形式に
  n_samples <- nrow(cohort)
  expr_rows <- lapply(seq_len(n_samples), function(i) {
    rna_col <- rna_col_dots[i]
    sid     <- sample_ids[i]
    
    if (is.na(rna_col) || !rna_col %in% names(rna_cp_df)) {
      # カラムが存在しない場合はNA
      base_row <- tibble(pair_id = sid)
      for (g in CHECKPOINT_GENES) {
        base_row[[paste0(g, "_tpm")]]     <- NA_real_
        base_row[[paste0(g, "_log2tpm")]] <- NA_real_
      }
      return(base_row)
    }
    
    base_row <- tibble(pair_id = sid)
    for (g in found_genes) {
      tpm_val <- suppressWarnings(
        as.numeric(rna_cp_df[rna_cp_df$Hugo_Symbol == g, rna_col, drop = TRUE])
      )
      tpm_val <- if (length(tpm_val) == 0 || is.na(tpm_val)) NA_real_ else tpm_val[1]
      base_row[[paste0(g, "_tpm")]]     <- tpm_val
      base_row[[paste0(g, "_log2tpm")]] <- if (!is.na(tpm_val)) log2(tpm_val + 1) else NA_real_
    }
    # 見つからなかった遺伝子はNA
    for (g in missing_genes) {
      base_row[[paste0(g, "_tpm")]]     <- NA_real_
      base_row[[paste0(g, "_log2tpm")]] <- NA_real_
    }
    base_row
  })
  
  expr_wide <- bind_rows(expr_rows)
  
  # pair_idで結合
  id_col_cohort <- if ("pair_id" %in% names(cohort)) "pair_id" else "sample_id"
  cohort_updated <- cohort %>%
    left_join(expr_wide %>% rename(!!id_col_cohort := pair_id),
              by = id_col_cohort)
  
  # 列順：既存列の後ろに新列
  new_cols <- as.vector(rbind(paste0(CHECKPOINT_GENES, "_tpm"),
                              paste0(CHECKPOINT_GENES, "_log2tpm")))
  add_cols_present <- intersect(new_cols, names(cohort_updated))
  cohort_updated <- cohort_updated %>%
    select(all_of(c(names(cohort), add_cols_present)))
  
  write_csv(cohort_updated, cohort_path)
  log_msg(sprintf("  上書き保存: %s (%d行 × %d列)",
                  basename(cohort_path), nrow(cohort_updated), ncol(cohort_updated)))
  
  # 追加列のNA率確認
  for (col in add_cols_present[grepl("_log2tpm$", add_cols_present)]) {
    na_rate <- mean(is.na(cohort_updated[[col]])) * 100
    gene_name <- str_remove(col, "_log2tpm$")
    log_msg(sprintf("  %s: median=%.3f, NA率=%.1f%%",
                    gene_name,
                    median(cohort_updated[[col]], na.rm = TRUE),
                    na_rate))
  }
  
  invisible(cohort_updated)
}

# glass_final_cohort.csv（主解析セット：coverage あり）
update_glass_cohort(
  file.path(GLASS_OUT_DIR, "glass_final_cohort.csv"),
  rna_cp, "main"
)

# glass_final_cohort_all101.csv（全101例・監査用）
update_glass_cohort(
  file.path(GLASS_OUT_DIR, "glass_final_cohort_all101.csv"),
  rna_cp, "all101"
)

# =============================================================================
# 8. 完了サマリー
# =============================================================================

log_msg("=== Step 06b サマリー ===")
log_msg(sprintf("追加遺伝子: %s", paste(CHECKPOINT_GENES, collapse = ", ")))
log_msg("[GDC]   gene_expression_wide.csv に6遺伝子×2形式（_tpm, _log2tpm）を追加")
log_msg("[GLASS] glass_final_cohort.csv に6遺伝子×2形式を追加")
log_msg("[GLASS] glass_final_cohort_all101.csv に6遺伝子×2形式を追加")
log_msg("次のステップ: Step16a（チェックポイント特異性解析）")
log_msg("=== Step 06b: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 06b 完了\n")
cat(sprintf("追加遺伝子: %s\n", paste(CHECKPOINT_GENES, collapse = ", ")))
cat("更新されたファイル:\n")
cat(sprintf("  [GDC]   %s/gene_expression_wide.csv\n", GDC_OUT_DIR))
cat(sprintf("  [GLASS] %s/glass_final_cohort.csv\n", GLASS_OUT_DIR))
cat(sprintf("  [GLASS] %s/glass_final_cohort_all101.csv\n", GLASS_OUT_DIR))
cat("新規ファイル:\n")
cat(sprintf("  [GDC]   %s/step06b_gdc_extraction_log.csv\n", GDC_OUT_DIR))
cat(sprintf("  [GDC]   %s/step06b_log.txt\n", GDC_OUT_DIR))
cat("============================\n")
