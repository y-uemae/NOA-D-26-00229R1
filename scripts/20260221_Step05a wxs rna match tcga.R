# =============================================================================
# step05a_wxs_rna_match_tcga.R
# GBM/Glioma TP53×LAG3 再解析 - Step 05a: WXS-RNA 対応付け（TCGA）
#
# 目的:
#   TCGA サンプルについて、代表選択済みの WXS ファイルと RNA ファイルを
#   1:1 対応付けする。
#
# -----------------------------------------------------------------------
# マッチングキー（固定）
#   Case キー  : case_id（GDC UUID）
#   Sample キー: sample_submitter_id（例: TCGA-XX-YYYY-01A/01B）
#
# マッチング優先順位（事前固定・恣意性排除）
#   Level 1 [A-A]              : -01A 同士で 1:1 → 採用（主解析）
#   Level 2 [B-B]              : -01A 不可、-01B 同士で 1:1 → 採用（主解析）
#   Level 3 [case_only]        : A/B 問わず Case 単位で WXS/RNA が各 1 本
#                                かつ 1:1 が明確 → 採用（主解析、フラグ付き）
#   Level 4 [cross_vial_rescue]: A-B / B-A 跨ぎ、上記全て不可の場合のみ
#                                → 感度解析枠（is_sensitivity_only=TRUE）
#   除外 [no_wxs_after_rep]      : WXS に対応 Case なし
#   除外 [no_rna_after_rep]      : RNA に対応 Case なし
#   除外 [multiple_candidates]   : Step04 後も複数候補残存（保険）
#   除外 [no_01A_or_01B_match]   : A/B どちらでも整合不能かつ case_only 非該当
#
# -----------------------------------------------------------------------
# 出力
#   05a_pairing_table.csv    : 全 Case の候補数・match_type・採否・除外理由
#                              （監査用完全記録）
#   05a_tcga_pairs_final.csv : 解析に入る最終 1:1 ペアのみ
#                              （感度解析枠は is_sensitivity_only=TRUE で区別）
#   step05a_log.txt
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
OUT_DIR    <- file.path(RESULT_DIR, "05a_wxs_rna_match_tcga")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step05a_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 05a: WXS-RNA 対応付け（TCGA） 開始 ===")
log_msg("マッチング優先順位: A-A > B-B > case_only > cross_vial_rescue（感度解析枠）")
log_msg("A-B / B-A 跨ぎは cross_vial_rescue として感度解析枠に限定")

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

extract_case_barcode <- function(x) {
  str_extract(x, "^TCGA-[A-Z0-9]+-[A-Z0-9]+")
}

extract_vial <- function(x) {
  # TCGA-XX-YYYY-01A-... の "01A" 部分を抽出
  str_extract(x, "(?<=-)(\\d{2}[A-Z])(?=-|$)")
}

# =============================================================================
# 3. 入力ファイル読み込み & TCGA 絞り込み
# =============================================================================

read_and_prep <- function(fname, modality_label) {
  fpath <- file.path(IN_DIR, fname)
  if (!file.exists(fpath)) stop(sprintf("File not found: %s", fpath))
  df <- read_csv(fpath, show_col_types = FALSE)
  log_msg(sprintf("読み込み: %s (%d行 × %d列)", fname, nrow(df), ncol(df)))
  
  proj_col   <- find_col(df, c("Project ID", "project_id"),        label = modality_label)
  fileid_col <- find_col(df, c("File ID", "file_id"),              label = modality_label)
  sample_col <- find_col(df, c("Sample ID", "sample_id",
                               "sample_submitter_id",
                               "Sample Submitter ID"),            label = modality_label)
  case_col   <- find_col(df, c("Case ID", "case_id",
                               "Case Submitter ID",
                               "case_submitter_id"),              label = modality_label)
  
  df <- df %>%
    rename(
      project_id = all_of(proj_col),
      file_id    = all_of(fileid_col),
      sample_id  = all_of(sample_col),
      case_id    = all_of(case_col)
    ) %>%
    filter(str_starts(project_id, "TCGA")) %>%
    mutate(
      case_barcode = extract_case_barcode(sample_id),
      vial         = extract_vial(sample_id),
      is_01A       = (!is.na(vial) & vial == "01A"),
      is_01B       = (!is.na(vial) & vial == "01B"),
      modality     = modality_label
    )
  
  na_cb <- sum(is.na(df$case_barcode))
  if (na_cb > 0) {
    log_msg(sprintf("  WARNING [%s]: case_barcode が NA の行: %d件",
                    modality_label, na_cb))
  }
  log_msg(sprintf("  [%s] TCGA 行数: %d件", modality_label, nrow(df)))
  df
}

log_msg("--- 入力ファイル読み込み ---")
wxs_tcga <- read_and_prep("wxs_representative.csv",        "WXS")
g2_tcga  <- read_and_prep("rna_grade2_representative.csv", "RNA_grade2")
g3_tcga  <- read_and_prep("rna_grade3_representative.csv", "RNA_grade3")
g4_tcga  <- read_and_prep("rna_grade4_representative.csv", "RNA_grade4")

# =============================================================================
# 4. WXS-RNA ペアリング関数
# =============================================================================

pair_wxs_rna <- function(wxs_df, rna_df, grade_label) {
  
  log_msg(sprintf("--- ペアリング: %s (WXS=%d, RNA=%d) ---",
                  grade_label, nrow(wxs_df), nrow(rna_df)))
  
  wxs_cases <- unique(wxs_df$case_barcode[!is.na(wxs_df$case_barcode)])
  rna_cases <- unique(rna_df$case_barcode[!is.na(rna_df$case_barcode)])
  all_cases  <- union(wxs_cases, rna_cases)
  
  log_msg(sprintf("  Case数: WXS=%d / RNA=%d / 共通=%d / WXSのみ=%d / RNAのみ=%d",
                  length(wxs_cases), length(rna_cases),
                  length(intersect(wxs_cases, rna_cases)),
                  length(setdiff(wxs_cases, rna_cases)),
                  length(setdiff(rna_cases, wxs_cases))))
  
  pairing_rows <- vector("list", length(all_cases))
  final_rows   <- vector("list", length(all_cases))
  
  for (idx in seq_along(all_cases)) {
    cb <- all_cases[idx]
    
    w <- wxs_df %>% filter(case_barcode == cb)
    r <- rna_df  %>% filter(case_barcode == cb)
    
    has_wxs <- nrow(w) > 0
    has_rna <- nrow(r) > 0
    
    # ---- WXS / RNA 片方のみ ----
    if (!has_wxs) {
      pairing_rows[[idx]] <- tibble(
        case_barcode        = cb,       grade = grade_label,
        wxs_sample_id       = NA_character_,
        rna_sample_id       = r$sample_id[1],
        wxs_file_id         = NA_character_,
        rna_file_id         = r$file_id[1],
        n_wxs_candidates    = 0L,
        n_rna_candidates    = nrow(r),
        match_type          = "RNA_only",
        include_flag        = FALSE,
        is_sensitivity_only = FALSE,
        exclude_reason      = "no_wxs_after_rep"
      )
      next
    }
    if (!has_rna) {
      pairing_rows[[idx]] <- tibble(
        case_barcode        = cb,       grade = grade_label,
        wxs_sample_id       = w$sample_id[1],
        rna_sample_id       = NA_character_,
        wxs_file_id         = w$file_id[1],
        rna_file_id         = NA_character_,
        n_wxs_candidates    = nrow(w),
        n_rna_candidates    = 0L,
        match_type          = "WXS_only",
        include_flag        = FALSE,
        is_sensitivity_only = FALSE,
        exclude_reason      = "no_rna_after_rep"
      )
      next
    }
    
    # ---- Step04 後も複数候補残存（保険） ----
    if (nrow(w) > 1 || nrow(r) > 1) {
      log_msg(sprintf("  WARNING [%s] %s: 複数候補残存 (WXS=%d, RNA=%d) → 除外",
                      grade_label, cb, nrow(w), nrow(r)))
      pairing_rows[[idx]] <- tibble(
        case_barcode        = cb,       grade = grade_label,
        wxs_sample_id       = w$sample_id[1],
        rna_sample_id       = r$sample_id[1],
        wxs_file_id         = w$file_id[1],
        rna_file_id         = r$file_id[1],
        n_wxs_candidates    = nrow(w),
        n_rna_candidates    = nrow(r),
        match_type          = "multiple_candidates",
        include_flag        = FALSE,
        is_sensitivity_only = FALSE,
        exclude_reason      = "multiple_candidates"
      )
      next
    }
    
    # ---- 以下は WXS/RNA それぞれ 1 行の Case ----
    w1 <- w[1, ]; r1 <- r[1, ]
    
    if (w1$is_01A && r1$is_01A) {
      # Level 1: A-A（主解析）
      match_type   <- "A-A"
      include_flag <- TRUE
      is_sens      <- FALSE
      excl         <- NA_character_
      
    } else if (w1$is_01B && r1$is_01B) {
      # Level 2: B-B（主解析）
      match_type   <- "B-B"
      include_flag <- TRUE
      is_sens      <- FALSE
      excl         <- NA_character_
      
    } else {
      # A/B 跨ぎ（A-B または B-A）
      # Case 単位で各 1 本しかない = case_only として主解析採用
      # ※ nrow(w)==1 && nrow(r)==1 は上記で保証済み
      is_cross_vial <- (w1$is_01A != r1$is_01A)  # 異なる vial
      
      if (is_cross_vial) {
        # 跨ぎだが各 1 本しかないので case_only（主解析）
        match_type   <- "case_only"
        include_flag <- TRUE
        is_sens      <- FALSE
        excl         <- NA_character_
        log_msg(sprintf("  INFO [%s] %s: case_only採用 (WXS=%s, RNA=%s)",
                        grade_label, cb, w1$vial, r1$vial))
      } else {
        # 両方 vial 不明など整合不能
        match_type   <- "no_01A_or_01B_match"
        include_flag <- FALSE
        is_sens      <- FALSE
        excl         <- "no_01A_or_01B_match"
      }
    }
    
    pairing_rows[[idx]] <- tibble(
      case_barcode        = cb,       grade = grade_label,
      wxs_sample_id       = w1$sample_id,
      rna_sample_id       = r1$sample_id,
      wxs_file_id         = w1$file_id,
      rna_file_id         = r1$file_id,
      n_wxs_candidates    = 1L,
      n_rna_candidates    = 1L,
      match_type          = match_type,
      include_flag        = include_flag,
      is_sensitivity_only = is_sens,
      exclude_reason      = excl
    )
    
    if (include_flag) {
      final_rows[[idx]] <- tibble(
        case_barcode        = cb,
        grade               = grade_label,
        project_id          = w1$project_id,
        wxs_sample_id       = w1$sample_id,
        rna_sample_id       = r1$sample_id,
        wxs_file_id         = w1$file_id,
        rna_file_id         = r1$file_id,
        match_type          = match_type,
        is_sensitivity_only = is_sens
      )
    }
  }
  
  pairing_table <- bind_rows(pairing_rows)
  final_df      <- bind_rows(final_rows)
  
  # ---- ログ集計 ----
  tbl <- pairing_table %>%
    count(match_type, include_flag, is_sensitivity_only) %>%
    arrange(desc(include_flag), match_type)
  
  log_msg(sprintf("  [%s] 結果内訳:", grade_label))
  for (i in seq_len(nrow(tbl))) {
    ro <- tbl[i, ]
    sens_note <- if (isTRUE(ro$is_sensitivity_only)) " [感度解析枠]" else ""
    incl_note <- if (ro$include_flag) "採用" else "除外"
    log_msg(sprintf("    %-25s %s%s: %d件",
                    ro$match_type, incl_note, sens_note, ro$n))
  }
  n_sens <- sum(isTRUE(final_df$is_sensitivity_only))
  log_msg(sprintf("  [%s] 最終採用: %d件（うち感度解析枠: %d件）",
                  grade_label, nrow(final_df), n_sens))
  
  list(pairing_table = pairing_table, final = final_df)
}

# =============================================================================
# 5. Grade 別実行
# =============================================================================

# =============================================================================
# 4b. 診断ブロック: 行数 vs 重複 case_barcode の不一致を事前に検出・ログ出力
#     （例: Grade4の284行 → 283 distinct caseの差分を説明可能にする）
# =============================================================================

diagnose_duplicates <- function(df, label) {
  n_rows    <- nrow(df)
  n_na      <- sum(is.na(df$case_barcode))
  n_valid   <- n_rows - n_na
  n_distinct <- n_distinct(df$case_barcode[!is.na(df$case_barcode)])
  n_dup_cases <- n_valid - n_distinct  # 同一case_barcodeに2行以上存在する行の超過数
  
  log_msg(sprintf("  診断 [%s]: n_rows=%d, n_distinct(case_barcode)=%d, n_missing(case_barcode)=%d",
                  label, n_rows, n_distinct, n_na))
  
  if (n_dup_cases > 0) {
    # 重複しているcase_barcodeを特定してWARNING出力
    dup_cases <- df %>%
      filter(!is.na(case_barcode)) %>%
      count(case_barcode) %>%
      filter(n > 1) %>%
      arrange(desc(n))
    log_msg(sprintf("  WARNING [%s]: case_barcode重複 %d件（%d超過行）。詳細:",
                    label, nrow(dup_cases), n_dup_cases))
    for (i in seq_len(min(nrow(dup_cases), 10))) {
      log_msg(sprintf("    %s: %d行", dup_cases$case_barcode[i], dup_cases$n[i]))
    }
    log_msg(sprintf("  → ペアリング関数内でfirst-row選択により1:1に解決されます（保険コード適用）"))
  } else if (n_na > 0) {
    log_msg(sprintf("  WARNING [%s]: case_barcode NA %d件 → ペアリング対象外（集計から除外）",
                    label, n_na))
  } else {
    log_msg(sprintf("  OK [%s]: 全行でcase_barcodeが一意（重複・欠損なし）", label))
  }
}

log_msg("--- 診断: 行数 vs distinct case_barcode ---")
diagnose_duplicates(wxs_tcga, "WXS")
diagnose_duplicates(g2_tcga,  "RNA_grade2")
diagnose_duplicates(g3_tcga,  "RNA_grade3")
diagnose_duplicates(g4_tcga,  "RNA_grade4")

# =============================================================================
# 5. Grade 別ペアリング実行
# =============================================================================

log_msg("--- Grade 別ペアリング実行 ---")
res_g2 <- pair_wxs_rna(wxs_tcga, g2_tcga, "Grade2")
res_g3 <- pair_wxs_rna(wxs_tcga, g3_tcga, "Grade3")
res_g4 <- pair_wxs_rna(wxs_tcga, g4_tcga, "Grade4")

# =============================================================================
# 6. 出力
# =============================================================================

log_msg("--- ファイル出力 ---")

pairing_table_all <- bind_rows(
  res_g2$pairing_table,
  res_g3$pairing_table,
  res_g4$pairing_table
)
write_csv(pairing_table_all, file.path(OUT_DIR, "05a_pairing_table.csv"))
log_msg(sprintf("保存: 05a_pairing_table.csv (%d行)", nrow(pairing_table_all)))

final_all <- bind_rows(res_g2$final, res_g3$final, res_g4$final)
write_csv(final_all, file.path(OUT_DIR, "05a_tcga_pairs_final.csv"))
log_msg(sprintf("保存: 05a_tcga_pairs_final.csv (%d行)", nrow(final_all)))

# =============================================================================
# 7. サマリー
# =============================================================================

log_msg("=== Step 05a サマリー ===")

smry <- pairing_table_all %>%
  count(grade, match_type, include_flag, is_sensitivity_only) %>%
  arrange(grade, desc(include_flag), match_type)

log_msg("全 Grade の match_type 集計:")
for (i in seq_len(nrow(smry))) {
  ro <- smry[i, ]
  sens <- if (isTRUE(ro$is_sensitivity_only)) " [感度解析枠]" else ""
  log_msg(sprintf("  [%-7s] %-25s %s%s: %d件",
                  ro$grade, ro$match_type,
                  if (ro$include_flag) "採用" else "除外",
                  sens, ro$n))
}

log_msg(sprintf("最終採用ペア総計: %d件", nrow(final_all)))
log_msg(sprintf("  Grade2: %d件 / Grade3: %d件 / Grade4: %d件",
                nrow(res_g2$final), nrow(res_g3$final), nrow(res_g4$final)))

total_sens <- sum(isTRUE(final_all$is_sensitivity_only))
log_msg(sprintf("  うち感度解析枠 (cross_vial_rescue): %d件", total_sens))

# case_only の件数と内訳
case_only_df <- final_all %>% filter(match_type == "case_only")
log_msg(sprintf("  うち case_only 採用: %d件", nrow(case_only_df)))
if (nrow(case_only_df) > 0) {
  for (i in seq_len(nrow(case_only_df))) {
    log_msg(sprintf("    %s (WXS=%s, RNA=%s, grade=%s)",
                    case_only_df$case_barcode[i],
                    case_only_df$wxs_sample_id[i],
                    case_only_df$rna_sample_id[i],
                    case_only_df$grade[i]))
  }
}

# 除外理由トップ
excl_summary <- pairing_table_all %>%
  filter(!include_flag, !is.na(exclude_reason)) %>%
  count(exclude_reason) %>%
  arrange(desc(n))
log_msg("除外理由内訳:")
for (i in seq_len(nrow(excl_summary))) {
  log_msg(sprintf("  %-35s: %d件", excl_summary$exclude_reason[i], excl_summary$n[i]))
}

# multiple_candidates が存在する場合はcase_barcodeを明示
mc_cases <- pairing_table_all %>%
  filter(match_type == "multiple_candidates") %>%
  select(case_barcode, grade, n_wxs_candidates, n_rna_candidates)
if (nrow(mc_cases) > 0) {
  log_msg(sprintf("  multiple_candidates 除外ケース（%d件）:", nrow(mc_cases)))
  for (i in seq_len(nrow(mc_cases))) {
    log_msg(sprintf("    %s [%s] WXS=%d, RNA=%d",
                    mc_cases$case_barcode[i], mc_cases$grade[i],
                    mc_cases$n_wxs_candidates[i], mc_cases$n_rna_candidates[i]))
  }
}

log_msg("=== Step 05a: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 05a 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    05a_pairing_table.csv       ← 全 Case 採否決定表（監査用）\n")
cat("    05a_tcga_pairs_final.csv    ← 解析採用ペア（感度解析枠はフラグ付き）\n")
cat("    step05a_log.txt\n")
cat("============================\n")
