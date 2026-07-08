# ============================================================
# Step 01: データ構造確認スクリプト（修正版）
# ファイル名: step01_inspect_data_structure.R
#
# 判明したファイル構造:
#   WXS sample sheet : GDC/glioma/WXS/gdc_sample_sheet.2025-09-29.tsv
#   RNA sample sheet : GDC/glioma/RNAseq/gdc_sample_sheet.2025-09-29.tsv
#                      （全Grade混在、WXSとのSample IDマッチ確認用）
#   RNA Grade別sheet : GDC/glioma/gdc_sample_sheet.2025-09-30_grade2/3/4.tsv
#                      （3ファイル、Grade別正規化の主軸）
#
# 重要設計方針:
#   ・RNAseq は Grade別に正規化するため Grade別sheetを主軸とする
#   ・Grade別3ファイルを統合後、複数Gradeに出現するCase IDを検出→除外予定
#   ・WXS-RNA対応はSample ID（TCGA-XX-XXXX-01A）で行う（TCGA）
#   ・代表選択（updated_datetime）はStep 04で実施。ここでは分布のみ。
#
# 確認項目:
#    1. WXS sample sheet の構造・aligned fields・重複分布
#    2. RNA sample sheet（全Grade）の構造・aligned fields・重複分布
#    3. RNA Grade別sheetの構造・統合後の複数Grade出現Case検出
#    4. WXS/RNA 列名比較
#    5. aligned fields 分割数整合チェック（全sheet）
#    6. RNA ダウンロードファイル構造確認
#    7. WXS MAF ファイル構造確認
#    8. GLASS データ確認
#    9. GLASS フォルダ内ファイル一覧
#   10. 全体サマリー
#
# 出力先: D:/Projects/GBM_Analysis/results/TP53/20260221/01_inspect/
# 作成日: 2026-02-21
# ============================================================

library(tidyverse)
library(tools)  # file_ext

# ============================================================
# 0. 設定
# ============================================================

OUT_DIR <- here::here("results", "TP53", "20260221", "01_inspect")

BASE_GDC  <- here::here("data", "raw", "GDC", "glioma")

# WXS sample sheet
WXS_SHEET_PATH <- file.path(BASE_GDC, "WXS/gdc_sample_sheet.2025-09-29.tsv")

# RNA sample sheet（全Grade混在）
RNA_SHEET_PATH <- file.path(BASE_GDC, "RNAseq/gdc_sample_sheet.2025-09-29.tsv")

# RNA Grade別 sample sheet（3ファイル）
RNA_GRADE_PATHS <- list(
  Grade2 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade2.tsv"),
  Grade3 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade3.tsv"),
  Grade4 = file.path(BASE_GDC, "gdc_sample_sheet.2025-09-30_grade4.tsv")
)

# ダウンロードフォルダ
RNA_DOWNLOADS_ROOT <- file.path(BASE_GDC, "RNAseq/downloads")
WXS_DOWNLOADS_ROOT <- file.path(BASE_GDC, "WXS/downloads")

# GLASS
GLASS_DATA_PATH <- here::here("data", "raw", "external_validation", "difg_glass", "data_mrna_seq_tpm.txt")

# ============================================================
# 1. 出力先・ログ開始
# ============================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "step01_log.txt")
log_con  <- file(LOG_FILE, open = "wt")

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., sep = ""))
  message(msg)
  writeLines(msg, log_con)
}

log_msg("=== Step 01: データ構造確認 開始 ===")
log_msg("出力先: ", OUT_DIR)
log_msg("R version: ", R.version$version.string)
log_msg("実行日時: ", format(Sys.time()))

# ============================================================
# ヘルパー関数
# ============================================================

check_file <- function(path, label = NULL) {
  lbl <- if (!is.null(label)) label else basename(path)
  if (file.exists(path)) {
    sz <- file.info(path)$size
    log_msg("  [FOUND] ", lbl, "  (", format(sz, big.mark = ","), " bytes)")
    return(TRUE)
  } else {
    log_msg("  [NOT FOUND] ", lbl, " -- ", path)
    return(FALSE)
  }
}

read_sheet <- function(path, label) {
  log_msg("--- ", label, " 読み込み ---")
  if (!file.exists(path)) {
    log_msg("  ファイルが見つかりません: ", path)
    return(NULL)
  }
  raw_lines <- readLines(path, n = 3)
  log_msg("  先頭3行: ")
  for (ln in raw_lines) log_msg("    ", ln)
  
  tryCatch({
    df <- read_tsv(path, show_col_types = FALSE)
    log_msg("  行数=", nrow(df), "  列数=", ncol(df))
    log_msg("  列名: ", paste(colnames(df), collapse = " | "))
    return(df)
  }, error = function(e) {
    log_msg("  エラー: ", conditionMessage(e)); return(NULL)
  })
}

report_col <- function(df, col, indent = "  ") {
  if (!col %in% colnames(df)) return(invisible(NULL))
  vals  <- df[[col]]
  n_na  <- sum(is.na(vals))
  n_mc  <- sum(grepl(",", vals, fixed = TRUE), na.rm = TRUE)
  uvals <- sort(unique(vals[!is.na(vals)]))
  log_msg(indent, "[", col, "] 総行=", length(vals),
          " カンマ含む=", n_mc, " NA=", n_na,
          " ユニーク=", length(uvals))
  if (length(uvals) <= 30) {
    log_msg(indent, "  ユニーク値: ", paste(uvals, collapse = " | "))
  } else {
    log_msg(indent, "  ユニーク値（先頭30）: ",
            paste(head(uvals, 30), collapse = " | "), " ...")
  }
}

# aligned fields 分割数整合チェック
check_aligned <- function(df, label) {
  log_msg("--- aligned fields チェック: ", label, " ---")
  aligned_cols <- c("Sample ID", "Tissue Type", "Tumor Descriptor", "Specimen Type")
  present <- intersect(aligned_cols, colnames(df))
  log_msg("  対象列: ", paste(present, collapse = ", "))
  if (length(present) < 2) {
    log_msg("  ※ 列が1つ以下のためスキップ"); return(invisible(NULL))
  }
  
  sc <- lapply(present, function(col) {
    vals <- df[[col]]
    ifelse(is.na(vals), 0L,
           vapply(strsplit(as.character(vals), ",", fixed = TRUE),
                  length, integer(1L)))
  })
  names(sc) <- present
  sc_df <- as.data.frame(sc)
  
  ok  <- apply(sc_df, 1, function(r) length(unique(r)) == 1)
  log_msg("  整合=", sum(ok), " / 不整合=", sum(!ok), " / 全", nrow(df), "行")
  
  for (col in present) {
    tbl <- sort(table(sc_df[[col]]), decreasing = TRUE)
    log_msg("  [", col, "] 分割数分布: ",
            paste(paste0(names(tbl), "個=", as.integer(tbl), "行"),
                  collapse = " | "))
  }
  
  if (sum(!ok) > 0) {
    bad <- which(!ok)
    log_msg("  [警告] 不整合行（先頭10件）:")
    for (i in head(bad, 10)) {
      rv <- paste(present, "=",
                  sapply(present, function(c) df[[c]][i]), collapse = " | ")
      log_msg("    行", i, ": ", rv)
    }
    key <- intersect(c("File ID", "Case ID", "Project ID"), colnames(df))
    write_csv(df[bad, c(key, present)],
              file.path(OUT_DIR,
                        paste0("aligned_inconsistent_",
                               gsub("[ /]", "_", label), ".csv")))
    log_msg("  不整合行CSV出力済み")
  }
}

# 重複分布の集計（代表選択はStep 04）
dup_dist <- function(df, label) {
  if (is.null(df)) return(invisible(NULL))
  log_msg("--- 重複分布: ", label, " ---")
  log_msg("  ※ 代表選択はStep 04。ここでは分布のみ出力。")
  
  sid_col  <- "Sample ID"
  case_col <- "Case ID"
  if (!sid_col %in% colnames(df)) {
    log_msg("  Sample ID列なし"); return(invisible(NULL))
  }
  
  ids_first <- sapply(strsplit(as.character(df[[sid_col]]), ","),
                      function(x) trimws(x[1]))
  
  # TCGA: -01A 単位
  is_tcga <- grepl("^TCGA", ids_first)
  if (any(is_tcga)) {
    ids_01a <- ids_first[is_tcga & grepl("-01A$", ids_first)]
    cnt     <- table(ids_01a)
    dist    <- as.data.frame(table(cnt))
    colnames(dist) <- c("files_per_sampleid", "n_sampleids")
    log_msg("  TCGA -01A Sample ID ファイル数分布:")
    for (i in seq_len(nrow(dist))) {
      log_msg("    ", dist$files_per_sampleid[i],
              "本: ", dist$n_sampleids[i], "件")
    }
    dup_ex <- names(cnt[cnt >= 2])
    if (length(dup_ex) > 0)
      log_msg("  重複例（先頭5）: ", paste(head(dup_ex, 5), collapse = ", "))
    
    fname <- paste0("dup_dist_tcga01a_", gsub("[ /]", "_", label), ".csv")
    write_csv(dist, file.path(OUT_DIR, fname))
    log_msg("  → ", fname)
  }
  
  # CPTAC/HCMI: Case ID 単位
  is_ch <- grepl("^(C3[LN]|HCM)", ids_first)
  if (any(is_ch) && case_col %in% colnames(df)) {
    case_ids <- sapply(strsplit(as.character(df[[case_col]][is_ch]), ","),
                       function(x) trimws(x[1]))
    cnt2  <- table(case_ids)
    dist2 <- as.data.frame(table(cnt2))
    colnames(dist2) <- c("files_per_caseid", "n_caseids")
    log_msg("  CPTAC/HCMI Case ID ファイル数分布:")
    for (i in seq_len(nrow(dist2))) {
      log_msg("    ", dist2$files_per_caseid[i],
              "本: ", dist2$n_caseids[i], "件")
    }
    fname2 <- paste0("dup_dist_cptac_hcmi_", gsub("[ /]", "_", label), ".csv")
    write_csv(dist2, file.path(OUT_DIR, fname2))
    log_msg("  → ", fname2)
  }
}

# ============================================================
# 2. WXS sample sheet
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("2. WXS sample sheet 確認")
log_msg("========================================")
check_file(WXS_SHEET_PATH, "WXS sample sheet")
wxs <- read_sheet(WXS_SHEET_PATH, "WXS sample sheet")
if (!is.null(wxs)) {
  for (col in c("Data Type", "Experimental Strategy", "Project ID",
                "Case ID", "Sample ID", "Tissue Type",
                "Tumor Descriptor", "Specimen Type",
                "File ID", "File Name")) {
    report_col(wxs, col)
  }
  check_aligned(wxs, "WXS sample sheet")
  dup_dist(wxs, "WXS")
  write_csv(head(wxs, 5), file.path(OUT_DIR, "wxs_samplesheet_head5.csv"))
  log_msg("  先頭5行 → wxs_samplesheet_head5.csv")
}

# ============================================================
# 3. RNA sample sheet（全Grade混在）
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("3. RNA sample sheet（全Grade混在）確認")
log_msg("========================================")
check_file(RNA_SHEET_PATH, "RNA sample sheet（全Grade）")
rna_all <- read_sheet(RNA_SHEET_PATH, "RNA sample sheet（全Grade）")
if (!is.null(rna_all)) {
  for (col in c("Data Type", "Experimental Strategy", "Project ID",
                "Case ID", "Sample ID", "Tissue Type",
                "Tumor Descriptor", "Specimen Type",
                "File ID", "File Name")) {
    report_col(rna_all, col)
  }
  check_aligned(rna_all, "RNA_all")
  dup_dist(rna_all, "RNA_all")
  write_csv(head(rna_all, 5), file.path(OUT_DIR, "rna_all_samplesheet_head5.csv"))
  log_msg("  先頭5行 → rna_all_samplesheet_head5.csv")
}

# ============================================================
# 4. RNA Grade別 sample sheet（3ファイル）
#    重要: 複数Gradeに出現するCase IDを検出する
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("4. RNA Grade別 sample sheet 確認")
log_msg("   + 複数Grade出現Case ID検出（除外候補）")
log_msg("========================================")

rna_grade <- list()
for (gr in names(RNA_GRADE_PATHS)) {
  path <- RNA_GRADE_PATHS[[gr]]
  log_msg("")
  log_msg("--- ", gr, " ---")
  check_file(path, paste0("RNA ", gr))
  df_g <- read_sheet(path, paste0("RNA ", gr))
  rna_grade[[gr]] <- df_g
  
  if (!is.null(df_g)) {
    for (col in c("Project ID", "Case ID", "Sample ID",
                  "Tissue Type", "Tumor Descriptor", "Specimen Type")) {
      report_col(df_g, col, "  ")
    }
    check_aligned(df_g, paste0("RNA_", gr))
    dup_dist(df_g, paste0("RNA_", gr))
    write_csv(head(df_g, 3),
              file.path(OUT_DIR,
                        paste0("rna_", tolower(gr), "_samplesheet_head3.csv")))
    log_msg("  先頭3行 → rna_", tolower(gr), "_samplesheet_head3.csv")
  }
}

# 複数Gradeに出現するCase IDの検出
log_msg("")
log_msg("--- 複数Grade出現 Case ID の検出 ---")
case_col <- "Case ID"

grade_case_lists <- lapply(names(rna_grade), function(gr) {
  df_g <- rna_grade[[gr]]
  if (is.null(df_g) || !case_col %in% colnames(df_g)) return(NULL)
  # aligned fields 展開前のため先頭要素を取り出す
  ids <- sapply(strsplit(as.character(df_g[[case_col]]), ","),
                function(x) trimws(x[1]))
  data.frame(Grade = gr, Case_ID = unique(ids), stringsAsFactors = FALSE)
})
grade_case_df <- bind_rows(grade_case_lists)

if (nrow(grade_case_df) > 0) {
  multi_grade <- grade_case_df %>%
    group_by(Case_ID) %>%
    summarise(Grades = paste(sort(Grade), collapse = ","),
              n_grades = n(), .groups = "drop") %>%
    filter(n_grades >= 2) %>%
    arrange(desc(n_grades), Case_ID)
  
  log_msg("  複数Grade出現Case ID数: ", nrow(multi_grade))
  if (nrow(multi_grade) > 0) {
    log_msg("  [除外候補] 以下のCase IDは複数Gradeに出現:")
    for (i in seq_len(nrow(multi_grade))) {
      log_msg("    ", multi_grade$Case_ID[i],
              "  Grades=", multi_grade$Grades[i])
    }
    write_csv(multi_grade,
              file.path(OUT_DIR, "multi_grade_cases.csv"))
    log_msg("  → multi_grade_cases.csv に出力（Step 02で除外処理に使用）")
  } else {
    log_msg("  複数Grade出現なし（全Case IDはいずれか1つのGradeにのみ存在）")
  }
}

# ============================================================
# 5. WXS / RNA 列名比較
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("5. WXS / RNA（全Grade）列名比較")
log_msg("========================================")
if (!is.null(wxs) && !is.null(rna_all)) {
  both     <- intersect(colnames(wxs), colnames(rna_all))
  wxs_only <- setdiff(colnames(wxs),    colnames(rna_all))
  rna_only <- setdiff(colnames(rna_all), colnames(wxs))
  log_msg("  共通列（", length(both),     "）: ", paste(both,     collapse = " | "))
  log_msg("  WXSのみ（", length(wxs_only), "）: ", paste(wxs_only, collapse = " | "))
  log_msg("  RNAのみ（", length(rna_only), "）: ", paste(rna_only, collapse = " | "))
}

# ============================================================
# 6. RNA ダウンロードファイル構造確認
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("6. RNA ダウンロードファイル構造確認")
log_msg("========================================")

if (dir.exists(RNA_DOWNLOADS_ROOT)) {
  subdirs <- list.dirs(RNA_DOWNLOADS_ROOT, full.names = TRUE, recursive = FALSE)
  log_msg("  サブフォルダ数（file_id数）: ", length(subdirs))
  checked <- 0
  for (sd in subdirs) {
    if (checked >= 3) break
    for (f in list.files(sd, full.names = TRUE)) {
      if (!tolower(file_ext(f)) %in% c("tsv", "txt")) next
      sz <- file.info(f)$size
      if (sz == 0) next
      log_msg("  [", basename(sd), "] ", basename(f),
              " (", format(sz, big.mark = ","), " bytes)")
      tryCatch({
        lns <- readLines(f, n = 12)
        log_msg("  先頭12行:")
        for (ln in lns) log_msg("    ", ln)
        data_lns <- lns[!grepl("^#", lns)]
        if (length(data_lns) >= 1) {
          hdr <- strsplit(data_lns[1], "\t")[[1]]
          log_msg("  ヘッダ列数=", length(hdr))
          log_msg("  ヘッダ（先頭10）: ", paste(head(hdr, 10), collapse = " | "))
          log_msg("  tpm_unstranded列: ",
                  if ("tpm_unstranded" %in% hdr) "存在" else "不在")
          log_msg("  gene_name列: ",
                  if ("gene_name" %in% hdr) "存在" else "不在")
        }
        # コメント行数の確認
        n_comment <- sum(grepl("^#", lns))
        log_msg("  先頭12行中のコメント行（#）数: ", n_comment)
      }, error = function(e) log_msg("  読み込みエラー: ", conditionMessage(e)))
      checked <- checked + 1
      break
    }
  }
} else {
  log_msg("  RNA_DOWNLOADS_ROOT が見つかりません: ", RNA_DOWNLOADS_ROOT)
}

# ============================================================
# 7. WXS MAF ファイル構造確認
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("7. WXS MAF ファイル構造確認")
log_msg("========================================")

if (dir.exists(WXS_DOWNLOADS_ROOT)) {
  subdirs_wxs <- list.dirs(WXS_DOWNLOADS_ROOT, full.names = TRUE, recursive = FALSE)
  log_msg("  サブフォルダ数（file_id数）: ", length(subdirs_wxs))
  checked_maf <- 0
  for (sd in subdirs_wxs) {
    if (checked_maf >= 2) break
    for (f in list.files(sd, full.names = TRUE)) {
      if (!grepl("\\.maf(\\.gz)?$", f, ignore.case = TRUE)) next
      sz <- file.info(f)$size
      log_msg("  [", basename(sd), "] ", basename(f),
              " (", format(sz, big.mark = ","), " bytes)")
      tryCatch({
        lns <- if (grepl("\\.gz$", f)) {
          con <- gzcon(file(f, "rb")); x <- readLines(con, n = 15); close(con); x
        } else readLines(f, n = 15)
        
        log_msg("  先頭15行:")
        for (ln in lns) log_msg("    ", ln)
        data_lns <- lns[!grepl("^#", lns)]
        if (length(data_lns) >= 1) {
          hdr <- strsplit(data_lns[1], "\t")[[1]]
          log_msg("  MAFヘッダ列数=", length(hdr))
          log_msg("  MAFヘッダ（先頭15列）: ",
                  paste(head(hdr, 15), collapse = " | "))
          key_cols <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type",
                        "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                        "GDC_FILTER", "HGVSp_Short", "FILTER")
          for (kc in key_cols) {
            log_msg("  重要列[", kc, "]: ",
                    if (kc %in% hdr) "存在" else "不在")
          }
        }
      }, error = function(e) log_msg("  MAFエラー: ", conditionMessage(e)))
      checked_maf <- checked_maf + 1
      break
    }
  }
} else {
  log_msg("  WXS_DOWNLOADS_ROOT が見つかりません: ", WXS_DOWNLOADS_ROOT)
}

# ============================================================
# 8. GLASS データ確認
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("8. GLASS データ確認")
log_msg("========================================")
check_file(GLASS_DATA_PATH, "GLASS RNA-seq TPM")
if (file.exists(GLASS_DATA_PATH)) {
  tryCatch({
    raw_lines <- readLines(GLASS_DATA_PATH, n = 5)
    log_msg("  先頭5行（生テキスト）:")
    for (ln in raw_lines) log_msg("    ", ln)
    glass_raw <- read_tsv(GLASS_DATA_PATH, show_col_types = FALSE)
    log_msg("  行数=", nrow(glass_raw), "  列数=", ncol(glass_raw))
    log_msg("  先頭5列名: ", paste(head(colnames(glass_raw), 5), collapse = " | "))
    log_msg("  末尾3列名: ", paste(tail(colnames(glass_raw), 3), collapse = " | "))
    log_msg("  推定サンプル数（列数-1）: ", ncol(glass_raw) - 1)
    gene_col <- colnames(glass_raw)[1]
    log_msg("  遺伝子名列（1列目）: ", gene_col)
    lag3 <- glass_raw %>% filter(.data[[gene_col]] == "LAG3")
    log_msg("  LAG3行数: ", nrow(lag3))
    if (nrow(lag3) > 0) {
      v <- as.numeric(lag3[1, -1])
      log_msg("  LAG3 TPM: min=", round(min(v, na.rm=T), 3),
              " median=", round(median(v, na.rm=T), 3),
              " max=", round(max(v, na.rm=T), 3),
              " NA=", sum(is.na(v)))
    }
    write_csv(glass_raw[1:5, 1:min(6, ncol(glass_raw))],
              file.path(OUT_DIR, "glass_data_head5x6.csv"))
    log_msg("  先頭5行×6列 → glass_data_head5x6.csv")
  }, error = function(e) log_msg("  エラー: ", conditionMessage(e)))
}

# ============================================================
# 9. GLASS フォルダ内ファイル一覧
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("9. GLASS フォルダ内ファイル一覧")
log_msg("========================================")
glass_dir <- dirname(GLASS_DATA_PATH)
if (dir.exists(glass_dir)) {
  gf_all <- list.files(glass_dir, full.names = TRUE,
                       pattern = ".*\\.(txt|csv|tsv|maf|gz|xlsx)$")
  log_msg("  ファイル数: ", length(gf_all))
  for (gf in gf_all) {
    sz <- file.info(gf)$size
    log_msg("    ", basename(gf), "  (", format(sz, big.mark = ","), " bytes)")
  }
}

# ============================================================
# 10. 全体サマリー
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("10. サマリー")
log_msg("========================================")

summary_df <- tibble(
  item = c(
    "wxs_sheet_found",
    "rna_sheet_all_found",
    "rna_sheet_grade2_found",
    "rna_sheet_grade3_found",
    "rna_sheet_grade4_found",
    "wxs_rows", "rna_all_rows",
    "rna_grade2_rows", "rna_grade3_rows", "rna_grade4_rows",
    "multi_grade_cases_detected",
    "rna_downloads_exist", "wxs_downloads_exist",
    "glass_found"
  ),
  value = as.character(c(
    file.exists(WXS_SHEET_PATH),
    file.exists(RNA_SHEET_PATH),
    file.exists(RNA_GRADE_PATHS$Grade2),
    file.exists(RNA_GRADE_PATHS$Grade3),
    file.exists(RNA_GRADE_PATHS$Grade4),
    if (!is.null(wxs))      nrow(wxs)      else NA,
    if (!is.null(rna_all))  nrow(rna_all)  else NA,
    if (!is.null(rna_grade$Grade2)) nrow(rna_grade$Grade2) else NA,
    if (!is.null(rna_grade$Grade3)) nrow(rna_grade$Grade3) else NA,
    if (!is.null(rna_grade$Grade4)) nrow(rna_grade$Grade4) else NA,
    if (exists("multi_grade") && nrow(multi_grade) > 0)
      paste(multi_grade$Case_ID, collapse = ",") else "none",
    dir.exists(RNA_DOWNLOADS_ROOT),
    dir.exists(WXS_DOWNLOADS_ROOT),
    file.exists(GLASS_DATA_PATH)
  ))
)

write_csv(summary_df, file.path(OUT_DIR, "step01_summary.csv"))
for (i in seq_len(nrow(summary_df)))
  log_msg("  ", summary_df$item[i], ": ", summary_df$value[i])

# ============================================================
# 終了
# ============================================================

log_msg("")
log_msg("=== Step 01 完了 ===")
log_msg("Step 02 へ進む前に確認すべき点:")
log_msg("  1. aligned fields の列名と分割数（WXS/RNA 各sheet）")
log_msg("  2. 複数Grade出現Case ID（multi_grade_cases.csv）→ Step 02で除外")
log_msg("  3. RNA ダウンロードファイルのコメント行数・tpm_unstranded列の位置")
log_msg("  4. MAFの重要列（GDC_FILTER等）の存在確認")
log_msg("  5. TCGA/CPTAC/HCMI の重複分布（Step 04の代表選択の規模感）")
log_msg("出力ファイル一覧:")
for (f in list.files(OUT_DIR, full.names = FALSE)) log_msg("  ", f)

close(log_con)
message("Step 01 完了。ログ: ", LOG_FILE)
