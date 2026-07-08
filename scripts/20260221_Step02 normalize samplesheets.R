# ============================================================
# Step 02: sample sheet 正規化スクリプト
# ファイル名: step02_normalize_samplesheets.R
#
# 目的:
#   GDC sample sheet（WXS・RNA Grade別）を正規化し、
#   Step 03以降で使える「1行1ファイル・1行1サンプル」形式に変換する。
#
# Step 01で判明した重要事項（設計に反映済み）:
#   1. aligned fields（Sample ID/Tissue Type/Tumor Descriptor/Specimen Type）は
#      全sheetで整合100% → 安全に展開可能
#   2. Grade別sheetにWXS（MAF）行が混入している
#      → Data Category AND Data Type の二重フィルタでRNAのみ抽出
#        RNA : Data Category == "Transcriptome Profiling"
#              AND Data Type == "Gene Expression Quantification"
#        WXS : Data Category == "Simple Nucleotide Variation"
#              AND Data Type == "Masked Somatic Mutation"
#              AND File Name が ".maf.gz" で終わる
#   3. 複数Gradeに出現するCase IDが3例（TCGA-CS-6670, TCGA-DH-5141, TCGA-HT-A619）
#      → 主解析から除外（ambiguous grade）しログに記録
#
# ---------------------------------------------------------------
# [Methods固定文 / 査読対応コメント]
# multi-grade 3例の除外根拠（論文 Methods または補足に記載）:
#
# "A small number of TCGA-LGG cases appeared in multiple
#  grade-filtered downloads. As grade could not be unambiguously
#  reconciled from the available GDC sample sheets and API
#  clinical fields, these cases were excluded from
#  grade-stratified analyses to avoid misclassification."
#
# 背景:
#   GDC APIのclinical fieldsにもgrade情報が存在しなかったため、
#   どちらのgradeが正しいかを一次情報で確定できない。
#   誤分類リスクを避けるため、保守的に除外する（査読に強い選択）。
# ---------------------------------------------------------------
#
# 処理フロー:
#   A. aligned fields展開（カンマ区切り → 行展開）
#   B. Data Category + Data Type 二重フィルタ
#      （Grade別RNAシートからWXS行を除去）
#   C. 複数Grade出現Case IDの除外（ambiguous grade）
#   D. 正規化後のsheetをCSVとして保存
#   E. 各sheetの基本統計をログ出力
#
# 入力:
#   WXS : GDC/glioma/WXS/gdc_sample_sheet.2025-09-29.tsv
#   RNA : GDC/glioma/gdc_sample_sheet.2025-09-30_grade{2,3,4}.tsv
#   除外: 01_inspect/multi_grade_cases.csv
#
# 出力先: D:/Projects/GBM_Analysis/results/TP53/20260221/02_normalize/
# 作成日: 2026-02-21
# ============================================================

library(tidyverse)

# ============================================================
# 0. 設定
# ============================================================

BASE_GDC    <- here::here("data", "raw", "GDC", "glioma")
INSPECT_DIR <- here::here("results", "TP53", "20260221", "01_inspect")
OUT_DIR     <- here::here("results", "TP53", "20260221", "02_normalize")

WXS_SHEET_PATH <- file.path(BASE_GDC, "WXS/gdc_sample_sheet.2025-09-29.tsv")

RNA_GRADE_PATHS <- list(
  Grade2 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade2.tsv"),
  Grade3 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade3.tsv"),
  Grade4 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade4.tsv")
)

MULTI_GRADE_PATH <- file.path(INSPECT_DIR, "multi_grade_cases.csv")

# aligned fields: この4列がカンマ区切りで連動している
ALIGNED_COLS <- c("Sample ID", "Tissue Type", "Tumor Descriptor", "Specimen Type")

# ---------------------------------------------------------------
# フィルタ条件（表記ゆれ・混入に強い二重フィルタ）
# ---------------------------------------------------------------
RNA_DATA_CATEGORY <- "Transcriptome Profiling"
RNA_DATA_TYPE     <- "Gene Expression Quantification"
RNA_FILE_EXT      <- "\\.tsv$"          # augmented_star_gene_counts.tsv

WXS_DATA_CATEGORY <- "Simple Nucleotide Variation"
WXS_DATA_TYPE     <- "Masked Somatic Mutation"
WXS_FILE_EXT      <- "\\.maf\\.gz$"    # aliquot_ensemble_masked.maf.gz

# ============================================================
# 1. 出力先・ログ開始
# ============================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "step02_log.txt")
log_con  <- file(LOG_FILE, open = "wt")

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., sep = ""))
  message(msg)
  writeLines(msg, log_con)
}

log_msg("=== Step 02: sample sheet 正規化 開始 ===")
log_msg("出力先: ", OUT_DIR)
log_msg("R version: ", R.version$version.string)
log_msg("実行日時: ", format(Sys.time()))
log_msg("")
log_msg("フィルタ条件（RNA）: Data Category='", RNA_DATA_CATEGORY,
        "' AND Data Type='", RNA_DATA_TYPE, "'")
log_msg("フィルタ条件（WXS）: Data Category='", WXS_DATA_CATEGORY,
        "' AND Data Type='", WXS_DATA_TYPE, "' AND File Name=~'", WXS_FILE_EXT, "'")

# ============================================================
# ヘルパー関数
# ============================================================

# aligned fields をカンマ区切りで展開する
# 4列（Sample ID / Tissue Type / Tumor Descriptor / Specimen Type）を
# 同じ順番で同時に分割し、1要素ごとに1行へ展開する
expand_aligned_fields <- function(df, label) {
  log_msg("--- aligned fields 展開: ", label, " ---")
  log_msg("  展開前行数: ", nrow(df))
  
  present <- intersect(ALIGNED_COLS, colnames(df))
  missing <- setdiff(ALIGNED_COLS, colnames(df))
  if (length(missing) > 0)
    log_msg("  [警告] aligned cols が存在しない: ",
            paste(missing, collapse = ", "))
  
  non_aligned <- setdiff(colnames(df), ALIGNED_COLS)
  result_list <- vector("list", nrow(df))
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    
    splits <- lapply(present, function(col) {
      trimws(strsplit(as.character(row[[col]]), ",", fixed = TRUE)[[1]])
    })
    names(splits) <- present
    
    # 分割数の整合確認
    n_splits <- sapply(splits, length)
    if (length(unique(n_splits)) > 1) {
      log_msg("  [警告] 行", i, ": 分割数不一致 ",
              paste(names(n_splits), "=", n_splits, collapse = ", "),
              " → 最小数で切り詰め")
      min_n  <- min(n_splits)
      splits <- lapply(splits, function(x) x[seq_len(min_n)])
    }
    
    n <- length(splits[[present[1]]])
    base_row   <- row[rep(1, n), non_aligned, drop = FALSE]
    rownames(base_row) <- NULL
    aligned_df <- as.data.frame(splits, check.names = FALSE,
                                stringsAsFactors = FALSE)
    
    result_list[[i]] <- bind_cols(base_row, aligned_df)
  }
  
  expanded <- bind_rows(result_list)
  log_msg("  展開後行数: ", nrow(expanded))
  
  # 展開後の各aligned colのユニーク値と件数
  for (col in present) {
    tbl <- sort(table(expanded[[col]]), decreasing = TRUE)
    log_msg("  [", col, "] (",  length(tbl), "種): ",
            paste(paste0(names(tbl), "=", as.integer(tbl)), collapse = " | "))
  }
  
  return(expanded)
}

# 基本統計の出力
report_stats <- function(df, label) {
  log_msg("--- 基本統計: ", label, " ---")
  log_msg("  行数=", nrow(df), "  列数=", ncol(df))
  
  for (col in c("Project ID", "Data Category", "Data Type",
                "Tissue Type", "Tumor Descriptor", "Specimen Type")) {
    if (!col %in% colnames(df)) next
    tbl <- sort(table(df[[col]]), decreasing = TRUE)
    log_msg("  [", col, "]: ",
            paste(paste0(names(tbl), "=", as.integer(tbl)), collapse = " | "))
  }
  
  # Project × Tissue Type クロス集計
  if (all(c("Project ID", "Tissue Type") %in% colnames(df))) {
    cross <- df %>%
      count(`Project ID`, `Tissue Type`) %>%
      arrange(`Project ID`, `Tissue Type`)
    log_msg("  [Project × Tissue Type]:")
    for (i in seq_len(nrow(cross)))
      log_msg("    ", cross$`Project ID`[i], " / ",
              cross$`Tissue Type`[i], " : ", cross$n[i], "件")
  }
}

# ============================================================
# 2. 複数Grade出現Case IDの読み込み
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("2. 複数Grade出現 Case ID の読み込み")
log_msg("========================================")
log_msg("  除外方針: ambiguous grade のため主解析から除外（保守的・査読に強い）")
log_msg("  Methods固定文（コメント参照）を根拠として除外ログに記録する")

multi_grade_cases <- character(0)
if (file.exists(MULTI_GRADE_PATH)) {
  mg <- read_csv(MULTI_GRADE_PATH, show_col_types = FALSE)
  multi_grade_cases <- mg$Case_ID
  log_msg("  除外対象 Case ID（", length(multi_grade_cases), "件）:")
  for (cid in multi_grade_cases) {
    grades <- mg$Grades[mg$Case_ID == cid]
    log_msg("    ", cid, "  出現Grade=（", grades, "）")
  }
} else {
  log_msg("  [警告] ", MULTI_GRADE_PATH, " が見つかりません")
  log_msg("  複数Grade除外をスキップします")
}

# ============================================================
# 3. WXS sample sheet の正規化
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("3. WXS sample sheet 正規化")
log_msg("========================================")

wxs_raw <- read_tsv(WXS_SHEET_PATH, show_col_types = FALSE)
log_msg("  読み込み: ", nrow(wxs_raw), "行")

# Data Category + Data Type + File Name 二重フィルタ（WXS確認）
n_before <- nrow(wxs_raw)
wxs_filtered <- wxs_raw %>%
  filter(
    `Data Category` == WXS_DATA_CATEGORY,
    `Data Type`     == WXS_DATA_TYPE,
    grepl(WXS_FILE_EXT, `File Name`, ignore.case = TRUE)
  )
n_removed <- n_before - nrow(wxs_filtered)

log_msg("  WXS フィルタ後: ", nrow(wxs_filtered), "行",
        " / 除去=", n_removed, "行")
if (n_removed > 0) {
  log_msg("  [注意] WXS sheetに想定外の行が含まれていました")
  unexpected <- wxs_raw %>%
    filter(!(
      `Data Category` == WXS_DATA_CATEGORY &
        `Data Type`     == WXS_DATA_TYPE &
        grepl(WXS_FILE_EXT, `File Name`, ignore.case = TRUE)
    ))
  log_msg("  除去された行の Data Category: ",
          paste(unique(unexpected$`Data Category`), collapse = " | "))
  log_msg("  除去された行の Data Type: ",
          paste(unique(unexpected$`Data Type`),     collapse = " | "))
}

# aligned fields展開
wxs_expanded <- expand_aligned_fields(wxs_filtered, "WXS")

# 統計
report_stats(wxs_expanded, "WXS（展開後）")

write_csv(wxs_expanded, file.path(OUT_DIR, "wxs_normalized.csv"))
log_msg("  → wxs_normalized.csv に出力（", nrow(wxs_expanded), "行）")

# ============================================================
# 4. RNA Grade別 sample sheet の正規化
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("4. RNA Grade別 sample sheet 正規化")
log_msg("========================================")

rna_normalized   <- list()
exclusion_records <- list()

for (gr in names(RNA_GRADE_PATHS)) {
  log_msg("")
  log_msg("--- ", gr, " ---")
  path <- RNA_GRADE_PATHS[[gr]]
  
  if (!file.exists(path)) {
    log_msg("  [NOT FOUND] ", path); next
  }
  
  df_raw <- read_tsv(path, show_col_types = FALSE)
  log_msg("  読み込み: ", nrow(df_raw), "行")
  
  # --- Step B: Data Category + Data Type + File Name 二重フィルタ ---
  # Grade別sheetにWXS行が混入しているため必須
  n_before_dt <- nrow(df_raw)
  
  df_rna <- df_raw %>%
    filter(
      `Data Category` == RNA_DATA_CATEGORY,
      `Data Type`     == RNA_DATA_TYPE,
      grepl(RNA_FILE_EXT, `File Name`, ignore.case = TRUE)
    )
  
  # 除去された行（WXS混入行）
  df_removed_dt <- df_raw %>%
    filter(!(
      `Data Category` == RNA_DATA_CATEGORY &
        `Data Type`     == RNA_DATA_TYPE &
        grepl(RNA_FILE_EXT, `File Name`, ignore.case = TRUE)
    ))
  
  n_removed_dt <- nrow(df_removed_dt)
  log_msg("  [Data Type/Category フィルタ] RNA残存=", nrow(df_rna),
          " / 除去=", n_removed_dt)
  
  if (n_removed_dt > 0) {
    # 除去行の内訳をログ
    removed_types <- df_removed_dt %>%
      count(`Data Category`, `Data Type`) %>%
      arrange(desc(n))
    for (i in seq_len(nrow(removed_types)))
      log_msg("    除去内訳: Category='", removed_types$`Data Category`[i],
              "' Type='", removed_types$`Data Type`[i],
              "' =", removed_types$n[i], "行")
    
    exclusion_records[[paste0(gr, "_dtype")]] <- df_removed_dt %>%
      mutate(exclusion_step   = "DataType_filter",
             exclusion_reason = paste0(
               "Non-RNA row in grade sheet (Category=",
               `Data Category`, ", Type=", `Data Type`, ")"),
             Grade = gr)
  }
  
  # --- Step A: aligned fields展開 ---
  df_expanded <- expand_aligned_fields(df_rna, paste0("RNA_", gr))
  
  # --- Step C: 複数Grade出現 Case ID の除外 ---
  n_before_mg <- nrow(df_expanded)
  
  # Case IDの先頭要素（展開済みなのでカンマなしのはずだが念のため）
  df_expanded <- df_expanded %>%
    mutate(.case_first = trimws(sapply(strsplit(`Case ID`, ","),
                                       function(x) x[1])))
  
  df_excl_mg <- df_expanded %>% filter(.case_first %in% multi_grade_cases)
  df_expanded <- df_expanded %>%
    filter(!.case_first %in% multi_grade_cases) %>%
    select(-.case_first)
  
  n_removed_mg <- n_before_mg - nrow(df_expanded)
  log_msg("  [multi-grade 除外] 除去=", n_removed_mg,
          " / 残存=", nrow(df_expanded))
  if (n_removed_mg > 0) {
    removed_cids <- unique(df_excl_mg$.case_first)
    log_msg("  除外 Case ID: ", paste(removed_cids, collapse = ", "))
    
    exclusion_records[[paste0(gr, "_multigrade")]] <- df_excl_mg %>%
      select(-.case_first) %>%
      mutate(
        exclusion_step   = "multi_grade_exclusion",
        exclusion_reason = paste0(
          "Ambiguous grade: case appears in multiple grade-filtered downloads.",
          " Grade could not be reconciled from GDC sample sheets or API.",
          " Excluded to avoid misclassification."),
        Grade = gr
      )
  }
  
  # Grade列を付与
  df_expanded <- df_expanded %>% mutate(Grade = gr)
  
  # 統計
  report_stats(df_expanded, paste0("RNA_", gr, "（展開後・除外後）"))
  
  rna_normalized[[gr]] <- df_expanded
  
  fname <- paste0("rna_", tolower(gr), "_normalized.csv")
  write_csv(df_expanded, file.path(OUT_DIR, fname))
  log_msg("  → ", fname, " に出力（", nrow(df_expanded), "行）")
}

# ============================================================
# 5. RNA 全Grade統合
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("5. RNA 全Grade統合（正規化後）")
log_msg("========================================")

rna_all_normalized <- bind_rows(rna_normalized)
log_msg("  統合後行数: ", nrow(rna_all_normalized))

# Project × Grade クロス集計
log_msg("  [Project × Grade]:")
rna_all_normalized %>%
  count(`Project ID`, Grade) %>%
  arrange(`Project ID`, Grade) %>%
  { for (i in seq_len(nrow(.)))
    log_msg("    ", .$`Project ID`[i], " / ", .$Grade[i],
            " : ", .$n[i], "件")
    invisible(.) }

# Specimen Type × Grade（cell-ish / Solid Tissue / Unknown の確認）
# これがStep 03フィルタ設計に直結する
log_msg("  [Specimen Type × Grade]:")
rna_all_normalized %>%
  count(`Specimen Type`, Grade) %>%
  arrange(Grade, desc(n)) %>%
  { for (i in seq_len(nrow(.)))
    log_msg("    Grade=", .$Grade[i], " / ",
            .$`Specimen Type`[i], " : ", .$n[i], "件")
    invisible(.) }

write_csv(rna_all_normalized, file.path(OUT_DIR, "rna_all_normalized.csv"))
log_msg("  → rna_all_normalized.csv に出力（", nrow(rna_all_normalized), "行）")

# ============================================================
# 6. 除外記録の出力
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("6. 除外記録の出力")
log_msg("========================================")

if (length(exclusion_records) > 0) {
  keep_cols <- c("File ID", "File Name", "Data Category", "Data Type",
                 "Project ID", "Case ID", "Sample ID",
                 "Tissue Type", "Tumor Descriptor", "Specimen Type",
                 "exclusion_step", "exclusion_reason", "Grade")
  
  excl_all <- bind_rows(lapply(exclusion_records, function(df) {
    df[, intersect(keep_cols, colnames(df))]
  }))
  
  write_csv(excl_all, file.path(OUT_DIR, "step02_excluded_records.csv"))
  log_msg("  除外記録: ", nrow(excl_all), "行 → step02_excluded_records.csv")
  
  excl_all %>%
    count(exclusion_step, exclusion_reason) %>%
    { for (i in seq_len(nrow(.)))
      log_msg("  [", .$exclusion_step[i], "] ",
              .$exclusion_reason[i], " : ", .$n[i], "件")
      invisible(.) }
} else {
  log_msg("  除外記録なし")
}

# ============================================================
# 7. 全体サマリー
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("7. 全体サマリー")
log_msg("========================================")

raw_rows <- c(
  WXS      = nrow(wxs_raw),
  Grade2   = nrow(read_tsv(RNA_GRADE_PATHS$Grade2, show_col_types = FALSE)),
  Grade3   = nrow(read_tsv(RNA_GRADE_PATHS$Grade3, show_col_types = FALSE)),
  Grade4   = nrow(read_tsv(RNA_GRADE_PATHS$Grade4, show_col_types = FALSE))
)

norm_rows <- c(
  WXS      = nrow(wxs_expanded),
  Grade2   = if (!is.null(rna_normalized$Grade2)) nrow(rna_normalized$Grade2) else NA_integer_,
  Grade3   = if (!is.null(rna_normalized$Grade3)) nrow(rna_normalized$Grade3) else NA_integer_,
  Grade4   = if (!is.null(rna_normalized$Grade4)) nrow(rna_normalized$Grade4) else NA_integer_
)

summary_df <- tibble(
  sheet            = names(raw_rows),
  rows_raw         = as.integer(raw_rows),
  rows_normalized  = as.integer(norm_rows),
  rows_diff        = as.integer(raw_rows) - as.integer(norm_rows)
)

write_csv(summary_df, file.path(OUT_DIR, "step02_summary.csv"))
log_msg("  サマリーを step02_summary.csv に出力")
log_msg(sprintf("  %-10s %8s %12s %10s", "sheet", "raw", "normalized", "diff"))
for (i in seq_len(nrow(summary_df)))
  log_msg(sprintf("  %-10s %8d %12d %10d",
                  summary_df$sheet[i],
                  summary_df$rows_raw[i],
                  summary_df$rows_normalized[i],
                  summary_df$rows_diff[i]))

# ============================================================
# 終了
# ============================================================

log_msg("")
log_msg("=== Step 02 完了 ===")
log_msg("Step 03 へ進む前に確認すべき点:")
log_msg("  1. raw→normalized の行数差（サマリー参照）")
log_msg("     Grade2/3は WXS混入行が約半数あったはずなので大幅減が正常")
log_msg("  2. Specimen Type × Grade の集計（Step 03の cell-ish 除外設計に使用）")
log_msg("  3. step02_excluded_records.csv で除外内容を確認")
log_msg("出力ファイル一覧:")
for (f in list.files(OUT_DIR, full.names = FALSE)) log_msg("  ", f)

close(log_con)
message("Step 02 完了。ログ: ", LOG_FILE)
