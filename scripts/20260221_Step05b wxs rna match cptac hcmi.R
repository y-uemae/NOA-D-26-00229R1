# =============================================================================
# step05b_wxs_rna_match_cptac_hcmi.R
# GBM/Glioma TP53×LAG3 再解析 - Step 05b: WXS-RNA 対応付け（CPTAC/HCMI）改訂版
#
# 診断結果（step05b_diagnose.R）による確定事項:
#   - CPTAC/HCMI の "Case ID" 列は sample sheet の展開形式で値が重複しており
#     GDC UUID ではない（例: "C1230738, C1230738"）
#   - WXS[Case ID] × RNA[Case ID] の共通 = 0件（フォーマット不一致）
#   - WXS[Sample ID] × RNA[Sample ID] の共通 = 197件 ← 正しい結合キー
#   → CPTAC/HCMI の結合キーは Sample ID を使用する
#
# -----------------------------------------------------------------------
# マッチングキー（確定）
#   主キー: Sample ID（GDC sample sheet の "Sample ID" 列）
#   補助記録: Case ID（ログ・監査用に保持するが結合には使わない）
#
# マッチング優先順位
#   Level 1 [sample_1to1] : Sample ID で WXS/RNA が 1:1 一致 → 採用（主解析）
#   除外 [no_wxs_after_rep]    : WXS に対応 Sample なし
#   除外 [no_rna_after_rep]    : RNA に対応 Sample なし
#   除外 [multiple_candidates] : Step04 後も複数候補残存（保険）
#
# -----------------------------------------------------------------------
# 出力
#   05b_pairing_table.csv              : 全 Sample の採否決定表（監査用）
#   05b_cptac_hcmi_pairs_final.csv     : 解析に入る最終 1:1 ペア
#   step05b_log.txt
#
# 作成日: 2026-02-21
# =============================================================================

library(dplyr)
library(readr)
library(stringr)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
IN_DIR     <- file.path(RESULT_DIR, "04_representative")
OUT_DIR    <- file.path(RESULT_DIR, "05b_wxs_rna_match_cptac_hcmi")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step05b_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 05b: WXS-RNA 対応付け（CPTAC/HCMI）改訂版 開始 ===")
log_msg("結合キー確定: Sample ID（診断スクリプトで WXS×RNA 共通197件を確認）")
log_msg("Case ID 列はフォーマット不一致のため結合には使用しない（監査記録のみ）")

# =============================================================================
# 2. ヘルパー関数
# =============================================================================

find_col <- function(df, candidates, label = "", required = TRUE) {
  col <- candidates[candidates %in% names(df)]
  if (length(col) == 0) {
    if (required) {
      log_msg(sprintf("  ERROR [%s]: 必須列が見つかりません。候補: [%s]",
                      label, paste(candidates, collapse = ", ")))
      log_msg(sprintf("  利用可能な列名: %s", paste(names(df), collapse = ", ")))
      stop("Required column not found.")
    }
    return(NA_character_)
  }
  col[1]
}

is_non_tcga <- function(v) !str_starts(v, "TCGA")

# =============================================================================
# 3. 入力ファイル読み込み & CPTAC/HCMI 絞り込み
# =============================================================================

read_and_prep <- function(fname, modality_label) {
  fpath <- file.path(IN_DIR, fname)
  if (!file.exists(fpath)) stop(sprintf("File not found: %s", fpath))
  df <- read_csv(fpath, show_col_types = FALSE)
  log_msg(sprintf("読み込み: %s (%d行 × %d列)", fname, nrow(df), ncol(df)))
  
  proj_col   <- find_col(df, c("Project ID", "project_id"),   label = modality_label)
  fileid_col <- find_col(df, c("File ID", "file_id"),         label = modality_label)
  sample_col <- find_col(df, c("Sample ID", "sample_id",
                               "sample_submitter_id",
                               "Sample Submitter ID"),        label = modality_label)
  case_col   <- find_col(df, c("Case ID", "case_id",
                               "Case Submitter ID",
                               "case_submitter_id"),          label = modality_label,
                         required = FALSE)
  
  df <- df %>%
    filter(is_non_tcga(.[[proj_col]])) %>%
    rename(
      project_id = all_of(proj_col),
      file_id    = all_of(fileid_col),
      sample_id  = all_of(sample_col)
    ) %>%
    mutate(modality = modality_label)
  
  # case_id 列は存在すれば保持（結合には使わないが監査用に残す）
  if (!is.na(case_col) && case_col %in% names(df)) {
    df <- df %>% rename(case_id_raw = all_of(case_col))
  } else {
    df$case_id_raw <- NA_character_
  }
  
  log_msg(sprintf("  [%s] CPTAC/HCMI 行数: %d件", modality_label, nrow(df)))
  
  # プロジェクト内訳
  proj_tbl <- df %>% count(project_id) %>% arrange(desc(n))
  for (i in seq_len(nrow(proj_tbl))) {
    log_msg(sprintf("    %s: %d件", proj_tbl$project_id[i], proj_tbl$n[i]))
  }
  
  df
}

log_msg("--- 入力ファイル読み込み ---")
wxs_df <- read_and_prep("wxs_representative.csv",        "WXS")
rna_g2 <- read_and_prep("rna_grade2_representative.csv", "RNA_grade2")
rna_g3 <- read_and_prep("rna_grade3_representative.csv", "RNA_grade3")
rna_g4 <- read_and_prep("rna_grade4_representative.csv", "RNA_grade4")

# =============================================================================
# 4. 診断ブロック: Sample ID の重複・欠損チェック
# =============================================================================

diagnose_sample_ids <- function(df, label) {
  n_rows    <- nrow(df)
  n_na      <- sum(is.na(df$sample_id))
  n_valid   <- n_rows - n_na
  n_distinct <- n_distinct(df$sample_id[!is.na(df$sample_id)])
  n_dup     <- n_valid - n_distinct
  
  log_msg(sprintf("  診断 [%s]: n_rows=%d, n_distinct(sample_id)=%d, n_missing(sample_id)=%d",
                  label, n_rows, n_distinct, n_na))
  
  if (n_dup > 0) {
    dup_samples <- df %>%
      filter(!is.na(sample_id)) %>%
      count(sample_id) %>%
      filter(n > 1) %>%
      arrange(desc(n))
    log_msg(sprintf("  WARNING [%s]: sample_id 重複 %d件（%d超過行）",
                    label, nrow(dup_samples), n_dup))
    for (i in seq_len(min(nrow(dup_samples), 5))) {
      log_msg(sprintf("    sample_id=%s: %d行",
                      dup_samples$sample_id[i], dup_samples$n[i]))
    }
  } else if (n_na > 0) {
    log_msg(sprintf("  WARNING [%s]: sample_id NA %d件 → ペアリング対象外",
                    label, n_na))
  } else {
    log_msg(sprintf("  OK [%s]: 全行で sample_id が一意（重複・欠損なし）", label))
  }
}

log_msg("--- 診断: Sample ID の重複・欠損チェック ---")
diagnose_sample_ids(wxs_df, "WXS")
diagnose_sample_ids(rna_g2, "RNA_grade2")
diagnose_sample_ids(rna_g3, "RNA_grade3")
diagnose_sample_ids(rna_g4, "RNA_grade4")

# =============================================================================
# 5. WXS-RNA ペアリング関数（Sample ID ベース）
# =============================================================================

pair_wxs_rna_cptac <- function(wxs_df, rna_df, grade_label) {
  
  log_msg(sprintf("--- ペアリング: %s (WXS=%d, RNA=%d) ---",
                  grade_label, nrow(wxs_df), nrow(rna_df)))
  
  wxs_valid <- wxs_df %>% filter(!is.na(sample_id))
  rna_valid <- rna_df %>% filter(!is.na(sample_id))
  
  wxs_samples <- unique(wxs_valid$sample_id)
  rna_samples <- unique(rna_valid$sample_id)
  all_samples  <- union(wxs_samples, rna_samples)
  
  log_msg(sprintf("  Sample数: WXS=%d / RNA=%d / 共通=%d / WXSのみ=%d / RNAのみ=%d",
                  length(wxs_samples), length(rna_samples),
                  length(intersect(wxs_samples, rna_samples)),
                  length(setdiff(wxs_samples, rna_samples)),
                  length(setdiff(rna_samples, wxs_samples))))
  
  pairing_rows <- vector("list", length(all_samples))
  final_rows   <- vector("list", length(all_samples))
  
  for (idx in seq_along(all_samples)) {
    sid <- all_samples[idx]
    
    w <- wxs_valid %>% filter(sample_id == sid)
    r <- rna_valid  %>% filter(sample_id == sid)
    
    has_wxs <- nrow(w) > 0
    has_rna <- nrow(r) > 0
    
    # ---- WXS のみ ----
    if (!has_wxs) {
      pairing_rows[[idx]] <- tibble(
        sample_id           = sid,
        grade               = grade_label,
        project_id          = r$project_id[1],
        case_id_raw         = r$case_id_raw[1],
        wxs_sample_id       = NA_character_,
        rna_sample_id       = r$sample_id[1],
        wxs_file_id         = NA_character_,
        rna_file_id         = r$file_id[1],
        n_wxs_candidates    = 0L,
        n_rna_candidates    = nrow(r),
        match_type          = "RNA_only",
        include_flag        = FALSE,
        exclude_reason      = "no_wxs_after_rep"
      )
      next
    }
    
    # ---- RNA のみ ----
    if (!has_rna) {
      pairing_rows[[idx]] <- tibble(
        sample_id           = sid,
        grade               = grade_label,
        project_id          = w$project_id[1],
        case_id_raw         = w$case_id_raw[1],
        wxs_sample_id       = w$sample_id[1],
        rna_sample_id       = NA_character_,
        wxs_file_id         = w$file_id[1],
        rna_file_id         = NA_character_,
        n_wxs_candidates    = nrow(w),
        n_rna_candidates    = 0L,
        match_type          = "WXS_only",
        include_flag        = FALSE,
        exclude_reason      = "no_rna_after_rep"
      )
      next
    }
    
    # ---- 複数候補残存（保険） ----
    if (nrow(w) > 1 || nrow(r) > 1) {
      log_msg(sprintf("  WARNING [%s] sample_id=%s: 複数候補残存 (WXS=%d, RNA=%d) → 除外",
                      grade_label, sid, nrow(w), nrow(r)))
      pairing_rows[[idx]] <- tibble(
        sample_id           = sid,
        grade               = grade_label,
        project_id          = w$project_id[1],
        case_id_raw         = w$case_id_raw[1],
        wxs_sample_id       = paste(w$sample_id, collapse = "|"),
        rna_sample_id       = paste(r$sample_id, collapse = "|"),
        wxs_file_id         = paste(w$file_id, collapse = "|"),
        rna_file_id         = paste(r$file_id, collapse = "|"),
        n_wxs_candidates    = nrow(w),
        n_rna_candidates    = nrow(r),
        match_type          = "multiple_candidates",
        include_flag        = FALSE,
        exclude_reason      = "multiple_candidates"
      )
      next
    }
    
    # ---- 正常: 各 1 本 → sample_1to1 で採用 ----
    w1 <- w[1, ]; r1 <- r[1, ]
    
    pairing_rows[[idx]] <- tibble(
      sample_id           = sid,
      grade               = grade_label,
      project_id          = w1$project_id,
      case_id_raw         = w1$case_id_raw,
      wxs_sample_id       = w1$sample_id,
      rna_sample_id       = r1$sample_id,
      wxs_file_id         = w1$file_id,
      rna_file_id         = r1$file_id,
      n_wxs_candidates    = 1L,
      n_rna_candidates    = 1L,
      match_type          = "sample_1to1",
      include_flag        = TRUE,
      exclude_reason      = NA_character_
    )
    
    final_rows[[idx]] <- tibble(
      sample_id     = sid,
      grade         = grade_label,
      project_id    = w1$project_id,
      case_id_raw   = w1$case_id_raw,
      wxs_sample_id = w1$sample_id,
      rna_sample_id = r1$sample_id,
      wxs_file_id   = w1$file_id,
      rna_file_id   = r1$file_id,
      match_type    = "sample_1to1"
    )
  }
  
  pairing_table <- bind_rows(pairing_rows)
  final_df      <- bind_rows(final_rows)
  
  tbl <- pairing_table %>%
    count(match_type, include_flag) %>%
    arrange(desc(include_flag), match_type)
  
  log_msg(sprintf("  [%s] 結果内訳:", grade_label))
  for (i in seq_len(nrow(tbl))) {
    ro <- tbl[i, ]
    log_msg(sprintf("    %-25s %s: %d件",
                    ro$match_type, if (ro$include_flag) "採用" else "除外", ro$n))
  }
  log_msg(sprintf("  [%s] 最終採用: %d件", grade_label, nrow(final_df)))
  
  list(pairing_table = pairing_table, final = final_df)
}

# =============================================================================
# 6. Grade 別実行
# =============================================================================

log_msg("--- Grade 別ペアリング実行 ---")
res_g2 <- pair_wxs_rna_cptac(wxs_df, rna_g2, "Grade2")
res_g3 <- pair_wxs_rna_cptac(wxs_df, rna_g3, "Grade3")
res_g4 <- pair_wxs_rna_cptac(wxs_df, rna_g4, "Grade4")

# =============================================================================
# 7. 出力
# =============================================================================

log_msg("--- ファイル出力 ---")

pairing_table_all <- bind_rows(
  res_g2$pairing_table,
  res_g3$pairing_table,
  res_g4$pairing_table
)
write_csv(pairing_table_all, file.path(OUT_DIR, "05b_pairing_table.csv"))
log_msg(sprintf("保存: 05b_pairing_table.csv (%d行)", nrow(pairing_table_all)))

final_all <- bind_rows(res_g2$final, res_g3$final, res_g4$final)
write_csv(final_all, file.path(OUT_DIR, "05b_cptac_hcmi_pairs_final.csv"))
log_msg(sprintf("保存: 05b_cptac_hcmi_pairs_final.csv (%d行)", nrow(final_all)))

# =============================================================================
# 8. サマリー
# =============================================================================

log_msg("=== Step 05b サマリー ===")

smry <- pairing_table_all %>%
  count(grade, match_type, include_flag) %>%
  arrange(grade, desc(include_flag), match_type)

log_msg("全 Grade の match_type 集計:")
for (i in seq_len(nrow(smry))) {
  ro <- smry[i, ]
  log_msg(sprintf("  [%-7s] %-25s %s: %d件",
                  ro$grade, ro$match_type,
                  if (ro$include_flag) "採用" else "除外", ro$n))
}

log_msg(sprintf("最終採用ペア総計: %d件", nrow(final_all)))
log_msg(sprintf("  Grade2: %d件 / Grade3: %d件 / Grade4: %d件",
                nrow(res_g2$final), nrow(res_g3$final), nrow(res_g4$final)))

# プロジェクト別内訳
if (nrow(final_all) > 0) {
  proj_smry <- final_all %>% count(project_id, grade) %>% arrange(grade, desc(n))
  log_msg("採用ペアのプロジェクト × Grade 内訳:")
  for (i in seq_len(nrow(proj_smry))) {
    log_msg(sprintf("  [%s] %s: %d件",
                    proj_smry$grade[i], proj_smry$project_id[i], proj_smry$n[i]))
  }
}

# 除外理由内訳
excl_smry <- pairing_table_all %>%
  filter(!include_flag, !is.na(exclude_reason)) %>%
  count(exclude_reason) %>%
  arrange(desc(n))
if (nrow(excl_smry) > 0) {
  log_msg("除外理由内訳:")
  for (i in seq_len(nrow(excl_smry))) {
    log_msg(sprintf("  %-35s: %d件",
                    excl_smry$exclude_reason[i], excl_smry$n[i]))
  }
}

# 前回（診断時）との比較メモ
log_msg("【診断結果との照合】")
log_msg("  診断スクリプトで確認した共通 Sample ID 数: 197件")
log_msg(sprintf("  今回の sample_1to1 採用数 (Grade4): %d件", nrow(res_g4$final)))
log_msg("  ※ 差分がある場合は重複 sample_id 等の影響を確認してください")

log_msg("=== Step 05b: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 05b（改訂版）完了\n")
cat("結合キー: Sample ID（Case IDから変更）\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    05b_pairing_table.csv              ← 全 Sample 採否決定表（監査用）\n")
cat("    05b_cptac_hcmi_pairs_final.csv     ← 解析採用ペア\n")
cat("    step05b_log.txt\n")
cat("============================\n")
