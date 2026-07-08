# =============================================================================
# step04_representative_selection.R
# GBM/Glioma TP53×LAG3 再解析 - Step 04: GDC API 代表選択
#
# 目的:
#   同一Sample ID（TCGA）または同一Case ID（CPTAC/HCMI）内で
#   複数ファイルが存在する場合、GDC API の analysis.updated_datetime
#   が最新のファイルを代表として選択する。
#
# 代表選択ルール（事前固定・恣意性排除）:
#   1. analysis.updated_datetime が最大（最新）のファイルを採用
#   2. 同点の場合: file_id 昇順（辞書順最小）でtiebreak
#   3. analysis.updated_datetime が欠損の場合: 最後尾扱い（選ばない）
#   4. グループ内全欠損の場合: file_id 昇順で代表を決め、
#      reason_code = "updated_datetime_missing_all" を記録
#
# グループキー定義:
#   TCGA   : sample_submitter_id（例: TCGA-XX-YYYY-01A）
#   CPTAC/HCMI: case_id（GDC UUID）
#
# 入力:
#   results/TP53/20260221/03_tumor_primary/wxs_tumor_primary.csv
#   results/TP53/20260221/03_tumor_primary/rna_grade2_tumor_primary.csv
#   results/TP53/20260221/03_tumor_primary/rna_grade3_tumor_primary.csv
#   results/TP53/20260221/03_tumor_primary/rna_grade4_tumor_primary.csv
#
# 出力:
#   results/TP53/20260221/04_representative/
#     wxs_representative.csv
#     rna_grade2_representative.csv
#     rna_grade3_representative.csv
#     rna_grade4_representative.csv
#     gdc_api_metadata.json          （APIから取得した生JSONの監査ログ）
#     step04_selection_log.csv       （各ファイルの選択根拠）
#     step04_log.txt
#
# 作成日: 2026-02-21
# =============================================================================

library(httr)
library(jsonlite)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR    <- here::here()
RESULT_DIR  <- file.path(BASE_DIR, "results/TP53/20260221")
IN_DIR      <- file.path(RESULT_DIR, "03_tumor_primary")
OUT_DIR     <- file.path(RESULT_DIR, "04_representative")

GDC_API_URL <- "https://api.gdc.cancer.gov/files"

# APIリクエスト設定
CHUNK_SIZE   <- 200   # 1リクエストあたりのfile_id数
MAX_RETRY    <- 3     # 最大リトライ回数
RETRY_WAIT   <- c(5, 15, 30)  # リトライ間隔（秒）：指数バックオフ

# GDC APIで取得するフィールド
GDC_FIELDS <- paste(c(
  "file_id",
  "file_name",
  "file_size",
  "md5sum",
  "data_type",
  "data_category",
  "experimental_strategy",
  "platform",
  "access",
  "analysis.workflow_type",
  "analysis.workflow_version",
  "analysis.updated_datetime",
  "cases.case_id",
  "cases.submitter_id",
  "cases.samples.portions.analytes.aliquots.aliquot_id",
  "cases.samples.portions.analytes.aliquots.submitter_id"
), collapse = ",")

# 出力ディレクトリ作成
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step04_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 04: GDC API 代表選択 開始 ===")
log_msg(sprintf("入力ディレクトリ: %s", IN_DIR))
log_msg(sprintf("出力ディレクトリ: %s", OUT_DIR))

# =============================================================================
# 2. 入力ファイル読み込み
# =============================================================================

read_input <- function(fname) {
  fpath <- file.path(IN_DIR, fname)
  if (!file.exists(fpath)) {
    log_msg(sprintf("ERROR: ファイルが見つかりません: %s", fpath))
    stop(sprintf("File not found: %s", fpath))
  }
  df <- read_csv(fpath, show_col_types = FALSE)
  log_msg(sprintf("読み込み完了: %s (%d行 × %d列)", fname, nrow(df), ncol(df)))
  return(df)
}

wxs_df    <- read_input("wxs_tumor_primary.csv")
rna_g2_df <- read_input("rna_grade2_tumor_primary.csv")
rna_g3_df <- read_input("rna_grade3_tumor_primary.csv")
rna_g4_df <- read_input("rna_grade4_tumor_primary.csv")

# =============================================================================
# 3. ユーティリティ関数
# =============================================================================

# --- file_id 列を特定するヘルパー ---
get_file_id_col <- function(df) {
  candidates <- c("File ID", "file_id", "FileID")
  col <- candidates[candidates %in% names(df)]
  if (length(col) == 0) stop("file_id列が見つかりません。列名を確認してください。")
  col[1]
}

# --- sample/case ID列を特定するヘルパー ---
get_sample_id_col <- function(df) {
  # GDC sample sheetの列名候補
  candidates <- c("Sample ID", "sample_id", "sample_submitter_id",
                  "Sample Submitter ID", "Aliquot ID")
  col <- candidates[candidates %in% names(df)]
  if (length(col) == 0) {
    # 全列名をログに出力してデバッグ支援
    log_msg(sprintf("  利用可能な列名: %s", paste(names(df), collapse = ", ")))
    stop("Sample ID列が見つかりません。")
  }
  col[1]
}

get_case_id_col <- function(df) {
  candidates <- c("Case ID", "case_id", "Case Submitter ID",
                  "case_submitter_id", "cases.case_id")
  col <- candidates[candidates %in% names(df)]
  if (length(col) == 0) {
    log_msg(sprintf("  利用可能な列名: %s", paste(names(df), collapse = ", ")))
    stop("Case ID列が見つかりません。")
  }
  col[1]
}

# --- プロジェクト判定（TCGAか否か） ---
is_tcga_project <- function(project_vec) {
  str_detect(project_vec, "^TCGA")
}

# --- グループキーの決定 ---
# TCGA: Sample ID（-01A, -01B などを含む）
# CPTAC/HCMI: Case ID
assign_group_key <- function(df) {
  file_id_col  <- get_file_id_col(df)
  sample_id_col <- tryCatch(get_sample_id_col(df), error = function(e) NA)
  case_id_col   <- tryCatch(get_case_id_col(df),   error = function(e) NA)
  
  # Project列の特定
  proj_col_candidates <- c("Project ID", "project_id", "Project")
  proj_col <- proj_col_candidates[proj_col_candidates %in% names(df)]
  if (length(proj_col) == 0) {
    log_msg("  WARNING: Project ID列が見つかりません。全行をCase IDベースで処理します。")
    df$group_key      <- df[[case_id_col]]
    df$group_key_type <- "case_id"
    return(df)
  }
  proj_col <- proj_col[1]
  
  # TCGA行とそれ以外を振り分け
  tcga_flag <- is_tcga_project(df[[proj_col]])
  
  df$group_key      <- NA_character_
  df$group_key_type <- NA_character_
  
  # TCGA: sample_submitter_id を使用
  if (!is.na(sample_id_col) && any(tcga_flag)) {
    df$group_key[tcga_flag]      <- df[[sample_id_col]][tcga_flag]
    df$group_key_type[tcga_flag] <- "sample_id"
  }
  
  # 非TCGA: case_id を使用
  if (!is.na(case_id_col) && any(!tcga_flag)) {
    df$group_key[!tcga_flag]      <- df[[case_id_col]][!tcga_flag]
    df$group_key_type[!tcga_flag] <- "case_id"
  }
  
  # group_key が NA になった行はケース IDで救済
  na_flag <- is.na(df$group_key)
  if (any(na_flag) && !is.na(case_id_col)) {
    df$group_key[na_flag]      <- df[[case_id_col]][na_flag]
    df$group_key_type[na_flag] <- "case_id_fallback"
    log_msg(sprintf("  WARNING: %d行でgroup_keyがNAになり、case_idで補完しました。",
                    sum(na_flag)))
  }
  
  return(df)
}

# =============================================================================
# 4. GDC API 呼び出し関数
# =============================================================================

#' GDC API に POST リクエストを送り、ファイルメタデータを取得する
#'
#' @param file_ids character vector of GDC file UUIDs
#' @return list (raw API response hits) or NULL on failure
gdc_api_post <- function(file_ids) {
  
  body <- list(
    filters  = list(
      op      = "in",
      content = list(
        field = "file_id",
        value = as.list(file_ids)
      )
    ),
    fields = GDC_FIELDS,
    format = "JSON",
    size   = length(file_ids)
  )
  
  for (attempt in seq_len(MAX_RETRY)) {
    tryCatch({
      resp <- POST(
        url         = GDC_API_URL,
        body        = toJSON(body, auto_unbox = TRUE),
        content_type_json(),
        timeout(120)
      )
      
      if (http_error(resp)) {
        err_msg <- sprintf("HTTP %d: %s", status_code(resp),
                           rawToChar(resp$content[1:min(200, length(resp$content))]))
        log_msg(sprintf("  API エラー (attempt %d/%d): %s", attempt, MAX_RETRY, err_msg))
        if (attempt < MAX_RETRY) {
          log_msg(sprintf("  %d秒後にリトライします...", RETRY_WAIT[attempt]))
          Sys.sleep(RETRY_WAIT[attempt])
          next
        } else {
          log_msg("  最大リトライ回数に達しました。このチャンクをスキップします。")
          return(NULL)
        }
      }
      
      parsed <- content(resp, as = "parsed", type = "application/json")
      return(parsed$data$hits)
      
    }, error = function(e) {
      log_msg(sprintf("  例外発生 (attempt %d/%d): %s", attempt, MAX_RETRY,
                      conditionMessage(e)))
      if (attempt < MAX_RETRY) {
        log_msg(sprintf("  %d秒後にリトライします...", RETRY_WAIT[min(attempt, length(RETRY_WAIT))]))
        Sys.sleep(RETRY_WAIT[min(attempt, length(RETRY_WAIT))])
      } else {
        log_msg("  最大リトライ回数に達しました。このチャンクをスキップします。")
        return(NULL)
      }
    })
  }
  return(NULL)
}

#' file_id リスト全体に対してAPIを叩き、フラットなメタデータ data.frame を返す
#'
#' @param file_ids character vector of all file UUIDs to query
#' @param label    ラベル文字列（ログ用）
#' @return data.frame of metadata (one row per file)
fetch_gdc_metadata <- function(file_ids, label = "") {
  
  file_ids <- unique(file_ids)
  n_total  <- length(file_ids)
  chunks   <- split(file_ids, ceiling(seq_along(file_ids) / CHUNK_SIZE))
  n_chunks <- length(chunks)
  
  log_msg(sprintf("  [%s] 合計 %d file_id を %d チャンクに分割して取得開始",
                  label, n_total, n_chunks))
  
  all_hits     <- list()
  failed_ids   <- character(0)
  
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    log_msg(sprintf("  チャンク %d/%d (%d件) を取得中...", i, n_chunks, length(chunk)))
    
    hits <- gdc_api_post(chunk)
    
    if (is.null(hits)) {
      log_msg(sprintf("  WARNING: チャンク %d の取得に失敗。%d件をfailed_idsへ",
                      i, length(chunk)))
      failed_ids <- c(failed_ids, chunk)
      next
    }
    
    all_hits <- c(all_hits, hits)
    log_msg(sprintf("  チャンク %d: %d件取得完了", i, length(hits)))
    
    # APIへの礼儀：短いインターバル
    if (i < n_chunks) Sys.sleep(0.5)
  }
  
  log_msg(sprintf("  [%s] 取得完了: %d件成功 / %d件失敗",
                  label, length(all_hits), length(failed_ids)))
  
  if (length(failed_ids) > 0) {
    log_msg(sprintf("  失敗したfile_id: %s", paste(failed_ids, collapse = ", ")))
  }
  
  # hits をフラットな data.frame に変換
  meta_df <- parse_hits_to_df(all_hits, failed_ids)
  
  return(meta_df)
}

#' API hits（list）をフラットな data.frame に変換する
parse_hits_to_df <- function(hits, failed_ids = character(0)) {
  
  if (length(hits) == 0) {
    log_msg("  WARNING: hits が空です。空のdata.frameを返します。")
    return(data.frame())
  }
  
  rows <- lapply(hits, function(h) {
    # analysis フィールドの安全な取得
    analysis <- h[["analysis"]]
    workflow_type    <- if (!is.null(analysis[["workflow_type"]]))    analysis[["workflow_type"]]    else NA_character_
    workflow_version <- if (!is.null(analysis[["workflow_version"]])) analysis[["workflow_version"]] else NA_character_
    updated_datetime <- if (!is.null(analysis[["updated_datetime"]])) analysis[["updated_datetime"]] else NA_character_
    
    # cases → 最初の case のみ取得（通常1件）
    cases <- h[["cases"]]
    case_id       <- NA_character_
    case_submitter <- NA_character_
    aliquot_submitter <- NA_character_
    
    if (!is.null(cases) && length(cases) > 0) {
      case_id        <- cases[[1]][["case_id"]]        %||% NA_character_
      case_submitter <- cases[[1]][["submitter_id"]]   %||% NA_character_
      
      # aliquot submitter_id を深掘り
      samps <- cases[[1]][["samples"]]
      if (!is.null(samps) && length(samps) > 0) {
        portions <- samps[[1]][["portions"]]
        if (!is.null(portions) && length(portions) > 0) {
          analytes <- portions[[1]][["analytes"]]
          if (!is.null(analytes) && length(analytes) > 0) {
            aliquots <- analytes[[1]][["aliquots"]]
            if (!is.null(aliquots) && length(aliquots) > 0) {
              aliquot_submitter <- aliquots[[1]][["submitter_id"]] %||% NA_character_
            }
          }
        }
      }
    }
    
    data.frame(
      file_id              = h[["file_id"]]           %||% NA_character_,
      file_name            = h[["file_name"]]          %||% NA_character_,
      file_size            = h[["file_size"]]          %||% NA_real_,
      md5sum               = h[["md5sum"]]             %||% NA_character_,
      data_type            = h[["data_type"]]          %||% NA_character_,
      data_category        = h[["data_category"]]      %||% NA_character_,
      experimental_strategy = h[["experimental_strategy"]] %||% NA_character_,
      platform             = h[["platform"]]           %||% NA_character_,
      access               = h[["access"]]             %||% NA_character_,
      workflow_type        = workflow_type,
      workflow_version     = workflow_version,
      analysis_updated_datetime = updated_datetime,
      api_case_id          = case_id,
      api_case_submitter_id = case_submitter,
      api_aliquot_submitter_id = aliquot_submitter,
      api_fetch_status     = "ok",
      stringsAsFactors     = FALSE
    )
  })
  
  meta_df <- bind_rows(rows)
  
  # 失敗したfile_idをNAで追記
  if (length(failed_ids) > 0) {
    failed_df <- data.frame(
      file_id              = failed_ids,
      file_name            = NA_character_,
      file_size            = NA_real_,
      md5sum               = NA_character_,
      data_type            = NA_character_,
      data_category        = NA_character_,
      experimental_strategy = NA_character_,
      platform             = NA_character_,
      access               = NA_character_,
      workflow_type        = NA_character_,
      workflow_version     = NA_character_,
      analysis_updated_datetime = NA_character_,
      api_case_id          = NA_character_,
      api_case_submitter_id = NA_character_,
      api_aliquot_submitter_id = NA_character_,
      api_fetch_status     = "api_fetch_failed",
      stringsAsFactors     = FALSE
    )
    meta_df <- bind_rows(meta_df, failed_df)
  }
  
  return(meta_df)
}

# NULL合体演算子（Rの古いバージョン対応）
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# =============================================================================
# 5. 代表選択関数
# =============================================================================

#' グループキー単位で代表ファイルを選択する
#'
#' @param df       sample sheet data.frame（group_key, group_key_type 列付き）
#' @param meta_df  GDC APIから取得したメタデータ data.frame
#' @param file_id_col  df内のfile_id列名
#' @return list(representative = data.frame, selection_log = data.frame)
select_representative <- function(df, meta_df, file_id_col) {
  
  # メタデータを結合
  joined <- left_join(df, meta_df, by = setNames("file_id", file_id_col))
  
  # workflow が同一グループ内で混在するか確認
  workflow_check <- joined %>%
    group_by(group_key) %>%
    summarise(
      n_workflows = n_distinct(workflow_type, na.rm = TRUE),
      workflow_mixed_flag = n_distinct(workflow_type, na.rm = TRUE) > 1,
      .groups = "drop"
    )
  
  joined <- left_join(joined, workflow_check, by = "group_key")
  
  # updated_datetime を POSIXct に変換（欠損はNAのまま）
  joined <- joined %>%
    mutate(
      dt_parsed = suppressWarnings(
        ymd_hms(analysis_updated_datetime, quiet = TRUE)
      )
    )
  
  # 選択ロジック
  selection_log <- joined %>%
    group_by(group_key, group_key_type) %>%
    mutate(
      n_candidates = n(),
      n_with_datetime = sum(!is.na(dt_parsed))
    ) %>%
    arrange(group_key,
            # updated_datetime が欠損のものを最後尾に
            is.na(dt_parsed),
            # 最新を先頭に
            desc(dt_parsed),
            # tiebreak: file_id 昇順
            .data[[file_id_col]]) %>%
    mutate(
      rank_in_group = row_number(),
      is_selected   = (rank_in_group == 1),
      tie_break_applied = (n_with_datetime > 0 &
                             n_distinct(dt_parsed[!is.na(dt_parsed)]) < n_with_datetime),
      reason_code = case_when(
        n_candidates == 1                         ~ "single_candidate",
        n_with_datetime == 0 & is_selected        ~ "updated_datetime_missing_all",
        is.na(dt_parsed) & is_selected            ~ "updated_datetime_missing_selected",
        tie_break_applied & is_selected           ~ "tie_break_by_file_id",
        is_selected                               ~ "ok",
        TRUE                                      ~ "not_selected"
      )
    ) %>%
    ungroup()
  
  # 代表ファイルのみ抽出
  representative <- selection_log %>%
    filter(is_selected) %>%
    select(-rank_in_group, -is_selected, -dt_parsed,
           -n_candidates, -n_with_datetime, -n_workflows)
  
  # 選択ログ用に整形
  log_out <- selection_log %>%
    select(
      group_key, group_key_type,
      file_id        = all_of(file_id_col),
      n_candidates,
      n_with_datetime,
      is_selected,
      analysis_updated_datetime,
      workflow_type,
      workflow_version,
      workflow_mixed_flag,
      tie_break_applied,
      reason_code,
      api_fetch_status
    )
  
  return(list(representative = representative, selection_log = log_out))
}

# =============================================================================
# 6. 処理実行
# =============================================================================

# 入力データにグループキーを付与
log_msg("--- グループキー付与 ---")

wxs_df    <- assign_group_key(wxs_df)
rna_g2_df <- assign_group_key(rna_g2_df)
rna_g3_df <- assign_group_key(rna_g3_df)
rna_g4_df <- assign_group_key(rna_g4_df)

# 全file_idを収集
file_id_col_wxs <- get_file_id_col(wxs_df)
file_id_col_rna <- get_file_id_col(rna_g2_df)  # grade間で同じはず

all_file_ids <- unique(c(
  wxs_df[[file_id_col_wxs]],
  rna_g2_df[[file_id_col_rna]],
  rna_g3_df[[file_id_col_rna]],
  rna_g4_df[[file_id_col_rna]]
))

log_msg(sprintf("全file_id数（重複除去後）: %d", length(all_file_ids)))

# =============================================================================
# 7. GDC API からメタデータ取得
# =============================================================================

log_msg("--- GDC API メタデータ取得開始 ---")

meta_df <- fetch_gdc_metadata(all_file_ids, label = "all_files")

log_msg(sprintf("取得完了: %d件のメタデータ", nrow(meta_df)))

# 取得したJSONを含む生データをリストとして保存（監査ログ）
api_metadata_save <- list(
  retrieval_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  gdc_api_url         = GDC_API_URL,
  fields_requested    = GDC_FIELDS,
  chunk_size          = CHUNK_SIZE,
  total_file_ids      = length(all_file_ids),
  metadata            = as.list(meta_df)
)

json_path <- file.path(OUT_DIR, "gdc_api_metadata.json")
write_json(api_metadata_save, json_path, pretty = TRUE, auto_unbox = TRUE)
log_msg(sprintf("監査ログ保存: %s", json_path))

# =============================================================================
# 8. 代表選択の実行
# =============================================================================

log_msg("--- 代表選択開始 ---")

process_one <- function(df, label) {
  fid_col <- get_file_id_col(df)
  log_msg(sprintf("  処理中: %s (%d行)", label, nrow(df)))
  result <- select_representative(df, meta_df, fid_col)
  n_rep  <- nrow(result$representative)
  n_total <- nrow(df)
  n_dup  <- n_total - n_rep
  log_msg(sprintf("  %s: %d → %d件（%d件を代表外として除外）",
                  label, n_total, n_rep, n_dup))
  return(result)
}

wxs_result    <- process_one(wxs_df,    "WXS")
rna_g2_result <- process_one(rna_g2_df, "RNA Grade2")
rna_g3_result <- process_one(rna_g3_df, "RNA Grade3")
rna_g4_result <- process_one(rna_g4_df, "RNA Grade4")

# =============================================================================
# 9. 選択ログの結合・保存
# =============================================================================

combined_log <- bind_rows(
  rna_g2_result$selection_log %>% mutate(dataset = "RNA_grade2"),
  rna_g3_result$selection_log %>% mutate(dataset = "RNA_grade3"),
  rna_g4_result$selection_log %>% mutate(dataset = "RNA_grade4"),
  wxs_result$selection_log    %>% mutate(dataset = "WXS")
)

log_path <- file.path(OUT_DIR, "step04_selection_log.csv")
write_csv(combined_log, log_path)
log_msg(sprintf("選択ログ保存: %s (%d行)", log_path, nrow(combined_log)))

# =============================================================================
# 10. 代表ファイルリストの保存
# =============================================================================

save_representative <- function(result, fname) {
  out_path <- file.path(OUT_DIR, fname)
  # 追加した作業列を除去して保存
  out_df <- result$representative %>%
    select(-any_of(c("group_key", "group_key_type", "workflow_mixed_flag",
                     "tie_break_applied", "reason_code",
                     "analysis_updated_datetime", "workflow_type",
                     "workflow_version", "api_fetch_status",
                     "api_case_id", "api_case_submitter_id",
                     "api_aliquot_submitter_id", "file_name", "file_size",
                     "md5sum", "data_type", "data_category",
                     "experimental_strategy", "platform", "access")))
  write_csv(out_df, out_path)
  log_msg(sprintf("代表ファイルリスト保存: %s (%d行)", fname, nrow(out_df)))
  return(invisible(out_df))
}

save_representative(wxs_result,    "wxs_representative.csv")
save_representative(rna_g2_result, "rna_grade2_representative.csv")
save_representative(rna_g3_result, "rna_grade3_representative.csv")
save_representative(rna_g4_result, "rna_grade4_representative.csv")

# =============================================================================
# 11. 選択ログのサマリー出力
# =============================================================================

log_msg("=== Step 04 サマリー ===")

summary_df <- combined_log %>%
  filter(is_selected) %>%
  count(dataset, reason_code) %>%
  arrange(dataset, reason_code)

log_msg("代表選択の理由コード別集計:")
for (i in seq_len(nrow(summary_df))) {
  log_msg(sprintf("  [%s] %s: %d件",
                  summary_df$dataset[i],
                  summary_df$reason_code[i],
                  summary_df$n[i]))
}

# workflow混在フラグ確認
mixed_count <- combined_log %>%
  filter(is_selected, isTRUE(workflow_mixed_flag)) %>%
  nrow()
if (mixed_count > 0) {
  log_msg(sprintf("  WARNING: workflow_mixed_flag=TRUE の代表選択: %d件", mixed_count))
} else {
  log_msg("  INFO: workflow_mixed_flag=TRUE の案件なし（正常）")
}

# API取得失敗確認
failed_count <- combined_log %>%
  filter(api_fetch_status == "api_fetch_failed") %>%
  nrow()
if (failed_count > 0) {
  log_msg(sprintf("  WARNING: APIフェッチ失敗のファイル: %d件", failed_count))
  log_msg("  ※ 該当file_idを確認し、手動取得またはStep05以降で要確認")
} else {
  log_msg("  INFO: APIフェッチ失敗なし（正常）")
}

log_msg("=== Step 04: 完了 ===")

# ログファイルを閉じる
close(log_con)

# =============================================================================
# 12. 出力ファイル一覧の表示
# =============================================================================

cat("\n============================\n")
cat("Step 04 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    wxs_representative.csv\n")
cat("    rna_grade2_representative.csv\n")
cat("    rna_grade3_representative.csv\n")
cat("    rna_grade4_representative.csv\n")
cat("    gdc_api_metadata.json\n")
cat("    step04_selection_log.csv\n")
cat("    step04_log.txt\n")
cat("============================\n")
