# ============================================================
# Step 03: 腫瘍Primary サンプル確定フィルタ
# ファイル名: step03_select_tumor_primary.R
#
# 目的:
#   Step 02で正規化したsample sheetから、
#   「腫瘍Primary」サンプルのみを確定する。
#   引継書 Step B の実装。
#
# Step 02で判明した重要事項（設計に反映済み）:
#   ・Grade2/3: Solid Tissue(正常) + Unknown(要検討) のみ
#   ・Grade4: Solid Tissue + Unknown(340件=50%) + cell-ish(35件) + 1 ER
#   ・EXCEPTIONAL_RESPONDERS-ER が Grade4 に1件混入 → 対象外プロジェクトとして除外
#   ・Unknown は即除外しない（引継書方針）→ Solid TissueとともにQC候補として保持
#
# フィルタ規則（事前固定・恣意性排除）:
#
#   [共通 Step B-1] 対象プロジェクトのみ残す
#     TCGA-GBM, TCGA-LGG, CPTAC-3, HCMI-CMDC
#     ※ EXCEPTIONAL_RESPONDERS-ER は除外（対象外コホート）
#
#   [RNA Step B-2] Tissue Type / Tumor Descriptor フィルタ
#     Tissue Type     == "Tumor"
#     AND Tumor Descriptor == "Primary"
#     ※ RNA Grade別sheetは全行Tumor/Primaryなので確認のみ
#
#   [WXS Step B-2] Tissue Type / Tumor Descriptor フィルタ
#     aligned fields展開後、腫瘍Primary行を抽出
#     Tissue Type     == "Tumor"
#     AND Tumor Descriptor == "Primary"
#     → TCGAなら Sample ID が TCGA-..-01A 形式であることを確認
#
#   [WXS/RNA Step B-3] Specimen Type フィルタ
#     除外（cell-ish）:
#       Adherent Cell Line, 3D Neurosphere,
#       2D Modified Conditionally Reprogrammed Cells,
#       Mixed Adherent Suspension
#     保持（Solid Tissue / Unknown）:
#       Solid Tissue → 採用
#       Unknown      → 採用（「Unknown群」として保持、QCでレスキュー候補）
#                      ※ 引継書: "Unknown＝除外ではなく、Tumor Primary担保＋QC通過で採用"
#     その他（血液系: Peripheral Blood*, Buffy Coat, Buccal Cells）:
#       Tissue Type == "Normal" のものは既にB-2で除外済みのため
#       ここでは Tumor 側に残る血液系 Specimen Type を確認のみ
#       （WXS Normal行はB-2除外済み）
#
#   [TCGA Step B-4] Sample ID形式の確認
#     腫瘍Primary: Sample ID が TCGA-..-XXXX-01[A-Z] 形式
#     ※ 01B等も腫瘍Primaryなので 01[A-Z] で一括確認
#     ※ WXS-RNA対応の主キーは TCGA-..-01A（Step 05aで使用）
#
# 出力:
#   WXS腫瘍Primary確定リスト  : wxs_tumor_primary.csv
#   RNA Grade別腫瘍Primary確定 : rna_grade{2,3,4}_tumor_primary.csv
#   RNA全Grade統合             : rna_all_tumor_primary.csv
#   除外記録                   : step03_excluded_records.csv
#   サマリー                   : step03_summary.csv
#
# 出力先: D:/Projects/GBM_Analysis/results/TP53/20260221/03_tumor_primary/
# 作成日: 2026-02-21
# ============================================================

library(tidyverse)

# ============================================================
# 0. 設定
# ============================================================

NORMALIZE_DIR <- here::here("results", "TP53", "20260221", "02_normalize")
OUT_DIR       <- here::here("results", "TP53", "20260221", "03_tumor_primary")

# 対象プロジェクト（事前固定）
TARGET_PROJECTS <- c("TCGA-GBM", "TCGA-LGG", "CPTAC-3", "HCMI-CMDC")

# cell-ish Specimen Type（除外対象、事前固定）
CELLINE_SPECIMEN_TYPES <- c(
  "Adherent Cell Line",
  "3D Neurosphere",
  "2D Modified Conditionally Reprogrammed Cells",
  "Mixed Adherent Suspension"
)

# 採用 Specimen Type
SOLID_SPECIMEN_TYPE   <- "Solid Tissue"
UNKNOWN_SPECIMEN_TYPE <- "Unknown"

# ============================================================
# 1. 出力先・ログ開始
# ============================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "step03_log.txt")
log_con  <- file(LOG_FILE, open = "wt")

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., sep = ""))
  message(msg)
  writeLines(msg, log_con)
}

log_msg("=== Step 03: 腫瘍Primary サンプル確定 開始 ===")
log_msg("出力先: ", OUT_DIR)
log_msg("R version: ", R.version$version.string)
log_msg("実行日時: ", format(Sys.time()))
log_msg("対象プロジェクト: ", paste(TARGET_PROJECTS, collapse = ", "))
log_msg("cell-ish除外対象: ", paste(CELLINE_SPECIMEN_TYPES, collapse = " | "))

# ============================================================
# ヘルパー関数
# ============================================================

# フィルタを段階的に適用し、各ステップの除去件数をログに残す
apply_filters_with_log <- function(df, sheet_label) {
  log_msg("--- フィルタ適用: ", sheet_label, " ---")
  log_msg("  入力行数: ", nrow(df))
  
  exclusions <- list()
  
  # Step B-1: 対象プロジェクトフィルタ
  excl_proj <- df %>% filter(!`Project ID` %in% TARGET_PROJECTS)
  df        <- df %>% filter(`Project ID`  %in% TARGET_PROJECTS)
  if (nrow(excl_proj) > 0) {
    log_msg("  [B-1 プロジェクト除外] 除去=", nrow(excl_proj),
            " 除外プロジェクト: ",
            paste(unique(excl_proj$`Project ID`), collapse = ", "))
    exclusions[["B1_project"]] <- excl_proj %>%
      mutate(exclusion_step   = "B1_project_filter",
             exclusion_reason = paste0("Not in target projects: ", `Project ID`))
  }
  log_msg("  B-1後: ", nrow(df), "行")
  
  # Step B-2: Tissue Type / Tumor Descriptor フィルタ
  excl_tt <- df %>%
    filter(!(`Tissue Type` == "Tumor" & `Tumor Descriptor` == "Primary"))
  df      <- df %>%
    filter(`Tissue Type` == "Tumor" & `Tumor Descriptor` == "Primary")
  
  if (nrow(excl_tt) > 0) {
    tt_summary <- excl_tt %>%
      count(`Tissue Type`, `Tumor Descriptor`) %>%
      arrange(desc(n))
    for (i in seq_len(nrow(tt_summary)))
      log_msg("  [B-2 除外内訳] Tissue=", tt_summary$`Tissue Type`[i],
              " Descriptor=", tt_summary$`Tumor Descriptor`[i],
              " : ", tt_summary$n[i], "件")
    exclusions[["B2_tissue"]] <- excl_tt %>%
      mutate(exclusion_step   = "B2_tissue_descriptor_filter",
             exclusion_reason = paste0("Not Tumor/Primary: Tissue=",
                                       `Tissue Type`, ", Descriptor=",
                                       `Tumor Descriptor`))
  }
  log_msg("  B-2後（Tumor/Primary）: ", nrow(df), "行")
  
  # Step B-3: Specimen Type フィルタ
  # cell-ish を除外、Solid Tissue / Unknown を保持
  # その他（血液系等）はB-2でTumor/Primaryとして残っているものだけ確認
  specimen_tbl <- sort(table(df$`Specimen Type`), decreasing = TRUE)
  log_msg("  [B-3前 Specimen Type分布]: ",
          paste(paste0(names(specimen_tbl), "=", as.integer(specimen_tbl)),
                collapse = " | "))
  
  # cell-ish 除外
  excl_cell <- df %>% filter(`Specimen Type` %in% CELLINE_SPECIMEN_TYPES)
  df        <- df %>% filter(!`Specimen Type` %in% CELLINE_SPECIMEN_TYPES)
  
  if (nrow(excl_cell) > 0) {
    cell_summary <- excl_cell %>%
      count(`Specimen Type`, `Project ID`) %>%
      arrange(desc(n))
    for (i in seq_len(nrow(cell_summary)))
      log_msg("  [B-3 cell-ish除外] ",
              cell_summary$`Specimen Type`[i], " / ",
              cell_summary$`Project ID`[i], " : ",
              cell_summary$n[i], "件")
    exclusions[["B3_celline"]] <- excl_cell %>%
      mutate(exclusion_step   = "B3_cellline_filter",
             exclusion_reason = paste0("Cell-ish Specimen Type: ",
                                       `Specimen Type`))
  }
  log_msg("  B-3後（cell-ish除外）: ", nrow(df), "行")
  
  # Unknown 群の確認（除外しない・ログのみ）
  n_solid   <- sum(df$`Specimen Type` == SOLID_SPECIMEN_TYPE,   na.rm = TRUE)
  n_unknown <- sum(df$`Specimen Type` == UNKNOWN_SPECIMEN_TYPE, na.rm = TRUE)
  n_other   <- nrow(df) - n_solid - n_unknown
  log_msg("  [B-3後 内訳] Solid Tissue=", n_solid,
          " / Unknown=", n_unknown,
          " / その他=", n_other)
  
  if (n_unknown > 0)
    log_msg("  [Unknown群] ", n_unknown,
            "件を保持（QC通過で採用候補）。",
            "B-3時点では除外しない（引継書方針）。")
  
  if (n_other > 0) {
    other_tbl <- df %>%
      filter(!`Specimen Type` %in% c(SOLID_SPECIMEN_TYPE, UNKNOWN_SPECIMEN_TYPE)) %>%
      count(`Specimen Type`, `Project ID`) %>%
      arrange(desc(n))
    log_msg("  [その他 Specimen Type（確認）]:")
    for (i in seq_len(nrow(other_tbl)))
      log_msg("    ", other_tbl$`Specimen Type`[i], " / ",
              other_tbl$`Project ID`[i], " : ", other_tbl$n[i], "件")
  }
  
  # Specimen_Type_Group 列を付与（Step 04以降で使用）
  df <- df %>%
    mutate(Specimen_Type_Group = case_when(
      `Specimen Type` == SOLID_SPECIMEN_TYPE   ~ "Solid",
      `Specimen Type` == UNKNOWN_SPECIMEN_TYPE ~ "Unknown",
      TRUE                                      ~ "Other"
    ))
  
  log_msg("  最終残存: ", nrow(df), "行")
  
  return(list(df = df, exclusions = exclusions))
}

# TCGA Sample ID の形式チェック（-01[A-Z] 腫瘍Primary確認）
check_tcga_sample_id <- function(df, sheet_label) {
  log_msg("--- TCGA Sample ID 形式確認: ", sheet_label, " ---")
  
  tcga_rows <- df %>% filter(grepl("^TCGA", `Project ID`))
  if (nrow(tcga_rows) == 0) {
    log_msg("  TCGA行なし"); return(invisible(NULL))
  }
  
  # -01[A-Z] パターン（腫瘍Primary）
  is_01x <- grepl("-01[A-Z]$", tcga_rows$`Sample ID`)
  # -01A のみ（WXS-RNA対応の主キー）
  is_01a <- grepl("-01A$",     tcga_rows$`Sample ID`)
  # それ以外のSample Type
  is_other_tumor <- !is_01x
  
  log_msg("  TCGA行数=", nrow(tcga_rows))
  log_msg("  -01A（推奨主キー）=", sum(is_01a),
          " / -01B以降=", sum(is_01x) - sum(is_01a),
          " / それ以外=", sum(is_other_tumor))
  
  # -01A 以外の Sample ID を確認
  not_01a <- tcga_rows %>%
    filter(!grepl("-01A$", `Sample ID`)) %>%
    count(`Sample ID`) %>%
    arrange(`Sample ID`)
  if (nrow(not_01a) > 0) {
    log_msg("  -01A以外のSample ID（先頭10件）:")
    for (i in seq_len(min(10, nrow(not_01a))))
      log_msg("    ", not_01a$`Sample ID`[i], " : ", not_01a$n[i], "件")
  }
}

# ============================================================
# 2. WXS 腫瘍Primary確定
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("2. WXS 腫瘍Primary確定")
log_msg("========================================")

wxs_norm <- read_csv(file.path(NORMALIZE_DIR, "wxs_normalized.csv"),
                     show_col_types = FALSE)
log_msg("  WXS正規化済み読み込み: ", nrow(wxs_norm), "行")

wxs_result   <- apply_filters_with_log(wxs_norm, "WXS")
wxs_filtered <- wxs_result$df

check_tcga_sample_id(wxs_filtered, "WXS")

write_csv(wxs_filtered,
          file.path(OUT_DIR, "wxs_tumor_primary.csv"))
log_msg("  → wxs_tumor_primary.csv に出力（", nrow(wxs_filtered), "行）")

# WXS除外記録
wxs_excl_all <- bind_rows(wxs_result$exclusions)

# ============================================================
# 3. RNA Grade別 腫瘍Primary確定
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("3. RNA Grade別 腫瘍Primary確定")
log_msg("========================================")

rna_filtered <- list()
rna_excl_all <- list()

for (gr in c("Grade2", "Grade3", "Grade4")) {
  log_msg("")
  log_msg("--- ", gr, " ---")
  
  fname_in <- paste0("rna_", tolower(gr), "_normalized.csv")
  path_in  <- file.path(NORMALIZE_DIR, fname_in)
  
  if (!file.exists(path_in)) {
    log_msg("  [NOT FOUND] ", path_in); next
  }
  
  df_norm <- read_csv(path_in, show_col_types = FALSE)
  log_msg("  正規化済み読み込み: ", nrow(df_norm), "行")
  
  result    <- apply_filters_with_log(df_norm, paste0("RNA_", gr))
  df_out    <- result$df
  
  check_tcga_sample_id(df_out, paste0("RNA_", gr))
  
  rna_filtered[[gr]]   <- df_out
  rna_excl_all[[gr]]   <- bind_rows(result$exclusions)
  
  fname_out <- paste0("rna_", tolower(gr), "_tumor_primary.csv")
  write_csv(df_out, file.path(OUT_DIR, fname_out))
  log_msg("  → ", fname_out, " に出力（", nrow(df_out), "行）")
}

# ============================================================
# 4. RNA 全Grade統合（腫瘍Primary確定後）
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("4. RNA 全Grade統合（腫瘍Primary確定後）")
log_msg("========================================")

rna_all_tp <- bind_rows(rna_filtered)
log_msg("  統合後行数: ", nrow(rna_all_tp))

# Project × Grade × Specimen_Type_Group クロス集計
log_msg("  [Project × Grade × Specimen_Type_Group]:")
rna_all_tp %>%
  count(`Project ID`, Grade, Specimen_Type_Group) %>%
  arrange(`Project ID`, Grade, Specimen_Type_Group) %>%
  { for (i in seq_len(nrow(.)))
    log_msg("    ", .$`Project ID`[i], " / ", .$Grade[i],
            " / ", .$Specimen_Type_Group[i], " : ", .$n[i], "件")
    invisible(.) }

# Unknown群のSample ID一覧（Step 04 QC判定に備えて）
unknown_samples <- rna_all_tp %>%
  filter(Specimen_Type_Group == "Unknown") %>%
  select(`Project ID`, Grade, `Case ID`, `Sample ID`, `File ID`, `File Name`)
log_msg("  Unknown群: ", nrow(unknown_samples), "件")
write_csv(unknown_samples,
          file.path(OUT_DIR, "rna_unknown_group_list.csv"))
log_msg("  → rna_unknown_group_list.csv に出力（Step 04 QC参照用）")

write_csv(rna_all_tp, file.path(OUT_DIR, "rna_all_tumor_primary.csv"))
log_msg("  → rna_all_tumor_primary.csv に出力（", nrow(rna_all_tp), "行）")

# ============================================================
# 5. 除外記録の統合出力
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("5. 除外記録の統合出力")
log_msg("========================================")

keep_cols <- c("File ID", "File Name", "Data Category", "Data Type",
               "Project ID", "Case ID", "Sample ID",
               "Tissue Type", "Tumor Descriptor", "Specimen Type",
               "Grade", "exclusion_step", "exclusion_reason")

all_excl <- bind_rows(
  wxs_excl_all,
  bind_rows(rna_excl_all)
) %>%
  { .[, intersect(keep_cols, colnames(.))] }

write_csv(all_excl, file.path(OUT_DIR, "step03_excluded_records.csv"))
log_msg("  除外記録: ", nrow(all_excl), "行 → step03_excluded_records.csv")

all_excl %>%
  count(exclusion_step, exclusion_reason) %>%
  arrange(exclusion_step) %>%
  { for (i in seq_len(nrow(.)))
    log_msg("  [", .$exclusion_step[i], "] ",
            .$exclusion_reason[i], " : ", .$n[i], "件")
    invisible(.) }

# ============================================================
# 6. 全体サマリー
# ============================================================

log_msg("")
log_msg("========================================")
log_msg("6. 全体サマリー")
log_msg("========================================")

summary_rows <- list(
  list(sheet="WXS",    n_in=nrow(wxs_norm),   n_out=nrow(wxs_filtered)),
  list(sheet="RNA_G2", n_in=nrow(read_csv(file.path(NORMALIZE_DIR,"rna_grade2_normalized.csv"), show_col_types=FALSE)),
       n_out=if(!is.null(rna_filtered$Grade2)) nrow(rna_filtered$Grade2) else NA),
  list(sheet="RNA_G3", n_in=nrow(read_csv(file.path(NORMALIZE_DIR,"rna_grade3_normalized.csv"), show_col_types=FALSE)),
       n_out=if(!is.null(rna_filtered$Grade3)) nrow(rna_filtered$Grade3) else NA),
  list(sheet="RNA_G4", n_in=nrow(read_csv(file.path(NORMALIZE_DIR,"rna_grade4_normalized.csv"), show_col_types=FALSE)),
       n_out=if(!is.null(rna_filtered$Grade4)) nrow(rna_filtered$Grade4) else NA),
  list(sheet="RNA_all",n_in=NA, n_out=nrow(rna_all_tp))
)

summary_df <- bind_rows(lapply(summary_rows, as_tibble))
write_csv(summary_df, file.path(OUT_DIR, "step03_summary.csv"))
log_msg("  サマリーを step03_summary.csv に出力")
log_msg(sprintf("  %-10s %8s %8s", "sheet", "n_in", "n_out"))
for (i in seq_len(nrow(summary_df)))
  log_msg(sprintf("  %-10s %8s %8s",
                  summary_df$sheet[i],
                  ifelse(is.na(summary_df$n_in[i]),  "-", summary_df$n_in[i]),
                  ifelse(is.na(summary_df$n_out[i]), "-", summary_df$n_out[i])))

# ============================================================
# 終了
# ============================================================

log_msg("")
log_msg("=== Step 03 完了 ===")
log_msg("Step 04 へ進む前に確認すべき点:")
log_msg("  1. WXS腫瘍Primary行数（Normal行が除去されているか）")
log_msg("  2. RNA Grade4: Unknown群（340件前後）が保持されているか")
log_msg("  3. cell-ish除外数（Grade4: 35件前後が除去されているか）")
log_msg("  4. TCGA以外（CPTAC/HCMI）に -01A 以外のSample IDがあるか")
log_msg("  5. rna_unknown_group_list.csv の内容（Project・Gradeの内訳）")
log_msg("     → Step 04のGDC API代表選択で、Unknownのupdated_datetimeも取得対象")
log_msg("出力ファイル一覧:")
for (f in list.files(OUT_DIR, full.names = FALSE)) log_msg("  ", f)

close(log_con)
message("Step 03 完了。ログ: ", LOG_FILE)
