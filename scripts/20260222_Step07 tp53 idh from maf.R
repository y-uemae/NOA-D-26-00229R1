# =============================================================================
# step07_tp53_idh_from_maf.R
# GBM/Glioma TP53×LAG3 再解析 - Step 07: TP53 / IDH1/2 変異抽出（MAF）
#
# 目的:
#   WXS ペアリング済みの MAF ファイルから TP53 および IDH1/2 の変異を
#   アミノ酸レベルで抽出し、サブグループ分類に使用するテーブルを作成する。
#
# -----------------------------------------------------------------------
# MAF ファイル仕様（確定済み）
#   形式   : .wxs.aliquot_ensemble_masked.maf.gz
#   コメント: 先頭複数行（#version, #annotation.spec, #contigs 等）
#   ヘッダ : # で始まらない最初の行（140列）
#   重要列 : Hugo_Symbol, Variant_Classification, Variant_Type,
#            Tumor_Sample_Barcode, HGVSp_Short,
#            GDC_FILTER, case_id, callers
#   ※ FILTER 列は存在しない（GDC_FILTER を使用）
#
# -----------------------------------------------------------------------
# フィルタ仕様（固定・補足I項）
#   GDC_FILTER 許容値: "", "PASS", "pass", "Pass"（正規化後に判定）
#   → それ以外は除外（reason: gdc_filter_fail）
#
# -----------------------------------------------------------------------
# 変異定義（固定・補足J項）
#
# TP53 非サイレント変異（全て対象）:
#   Missense_Mutation, Nonsense_Mutation,
#   Frame_Shift_Del, Frame_Shift_Ins,
#   Splice_Site,
#   In_Frame_Del, In_Frame_Ins,
#   Translation_Start_Site, Nonstop_Mutation
#
# IDH1/2 対象変異:
#   Hugo_Symbol %in% c("IDH1", "IDH2")
#   Variant_Classification == "Missense_Mutation"
#   主要ホットスポット: IDH1 p.R132H/C/S/G/L, IDH2 p.R172K/M/W/S/G
#   → ホットスポット以外の Missense も記録するが、IDH_mut 判定は全 Missense を含む
#
# -----------------------------------------------------------------------
# 出力
#   07_mutations/
#     tp53_mutation_table.csv     : 症例単位の TP53 変異サマリー
#     idh_mutation_table.csv      : 症例単位の IDH1/2 変異サマリー
#     mutation_combined.csv       : 両方を結合した解析用テーブル
#     mutation_raw_variants.csv   : 変異行レベルの生データ（監査用）
#     step07_log.txt
#
# 作成日: 2026-02-21
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
WXS_DL_DIR <- file.path(BASE_DIR, "data/raw/GDC/glioma/WXS/downloads")
OUT_DIR    <- file.path(RESULT_DIR, "07_mutations")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------
# GDC_FILTER 許容値（正規化後の値で比較）
# 診断結果（step07_diagnose.R）より:
#   このコホートの MAF は GDC_FILTER 列が全行 NA（値なし）
#   NA = フィルタ指定なし = PASS 相当として許容する
#   許容: NA, "", "pass"（tolower後）
#   除外: それ以外の文字列（"oxog", "strand_bias" 等が入る場合）
GDC_FILTER_PASS     <- c("", "pass")   # tolower() 後に比較
GDC_FILTER_ALLOW_NA <- TRUE            # NA を許容（診断結果に基づき TRUE）

# TP53 非サイレント Variant_Classification（固定）
TP53_NONSILENT <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

# IDH ホットスポット（記録用・判定は全 Missense を採用）
IDH_HOTSPOT <- c(
  "p.R132H", "p.R132C", "p.R132S", "p.R132G", "p.R132L",  # IDH1
  "p.R172K", "p.R172M", "p.R172W", "p.R172S", "p.R172G"   # IDH2
)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step07_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 07: TP53/IDH 変異抽出（MAF）開始 ===")
log_msg(sprintf("GDC_FILTER 許容値: %s%s",
                paste(ifelse(GDC_FILTER_PASS == "", '""', GDC_FILTER_PASS),
                      collapse = ", "),
                if (GDC_FILTER_ALLOW_NA) ", NA（診断結果に基づき許容）" else ""))
log_msg(sprintf("TP53 非サイレント分類: %s", paste(TP53_NONSILENT, collapse = ", ")))
log_msg("IDH 判定: IDH1/IDH2 の Missense_Mutation 全て（WT/Mut二値化）")

# =============================================================================
# 2. ペアリング済み WXS リストの読み込み
# =============================================================================

log_msg("--- WXS ペアリング済みリスト読み込み ---")

pairs_tcga  <- read_csv(
  file.path(RESULT_DIR, "05a_wxs_rna_match_tcga/05a_tcga_pairs_final.csv"),
  show_col_types = FALSE) %>% mutate(source = "TCGA")

pairs_cptac <- read_csv(
  file.path(RESULT_DIR, "05b_wxs_rna_match_cptac_hcmi/05b_cptac_hcmi_pairs_final.csv"),
  show_col_types = FALSE) %>% mutate(source = "CPTAC_HCMI")

all_pairs <- bind_rows(pairs_tcga, pairs_cptac)
log_msg(sprintf("全ペア数: %d件 (TCGA=%d / CPTAC_HCMI=%d)",
                nrow(all_pairs), nrow(pairs_tcga), nrow(pairs_cptac)))

# wxs_file_id の一覧（重複除去）
wxs_file_ids <- unique(all_pairs$wxs_file_id[!is.na(all_pairs$wxs_file_id)])
log_msg(sprintf("処理対象 WXS file_id（重複除去）: %d件", length(wxs_file_ids)))

# =============================================================================
# 3. MAF 読み込み関数
# =============================================================================

#' 1つの WXS file_id から MAF を読み込み、TP53/IDH 変異行を返す
#'
#' @return list(status, note, variants_df)
read_maf_variants <- function(file_id) {
  
  dir_path <- file.path(WXS_DL_DIR, file_id)
  
  empty_df <- data.frame(
    file_id = character(), Hugo_Symbol = character(),
    Variant_Classification = character(), Variant_Type = character(),
    Tumor_Sample_Barcode = character(), HGVSp_Short = character(),
    GDC_FILTER = character(), case_id = character(),
    callers = character(), stringsAsFactors = FALSE
  )
  
  if (!dir.exists(dir_path)) {
    return(list(status = "file_missing",
                note   = sprintf("Directory not found: %s", dir_path),
                gdc_filter_values = NULL,
                variants_df = empty_df))
  }
  
  maf_files <- list.files(dir_path, pattern = "\\.maf\\.gz$",
                          full.names = TRUE, recursive = FALSE)
  if (length(maf_files) == 0) {
    return(list(status = "file_missing",
                note   = sprintf("No .maf.gz in: %s", dir_path),
                gdc_filter_values = NULL,
                variants_df = empty_df))
  }
  
  maf_path <- maf_files[1]
  
  # MAF 読み込み: # コメント行をスキップ
  df <- tryCatch(
    read_tsv(maf_path,
             comment  = "#",
             col_types = cols(.default = "c"),
             show_col_types = FALSE),
    error = function(e) e
  )
  
  if (inherits(df, "error")) {
    return(list(status = "parse_error",
                note   = conditionMessage(df),
                gdc_filter_values = NULL,
                variants_df = empty_df))
  }
  
  # 必須列確認
  required_cols <- c("Hugo_Symbol", "Variant_Classification",
                     "Variant_Type", "Tumor_Sample_Barcode",
                     "HGVSp_Short", "GDC_FILTER")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    return(list(
      status = "parse_error",
      note   = sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")),
      gdc_filter_values = NULL,
      variants_df = empty_df
    ))
  }
  
  # -----------------------------------------------------------------------
  # GDC_FILTER 正規化・フィルタ（補足I項 + 診断結果に基づく修正）
  # 正規化: tolower(trimws()) + na_if() で文字列"NA"を本物NAに統一
  # 許容: 本物NA（GDC_FILTER_ALLOW_NA=TRUE時）, "", "pass"（正規化後）
  # 根拠: このデータセットでは GDC_FILTER が一貫して NA（診断で確認済み）
  #       aliquot_ensemble_masked MAF はキュレーション済み出力であり、
  #       NA はフィルタ指定なし＝pass-equivalent と解釈する
  # -----------------------------------------------------------------------
  df <- df %>%
    mutate(
      GDC_FILTER_norm = na_if(tolower(trimws(GDC_FILTER)), "na"),
      filter_pass     = GDC_FILTER_norm %in% GDC_FILTER_PASS |
        ((GDC_FILTER_ALLOW_NA) &
           (is.na(GDC_FILTER) | is.na(GDC_FILTER_norm)))
    )
  
  n_total  <- nrow(df)
  n_pass   <- sum(df$filter_pass, na.rm = TRUE)
  n_fail   <- n_total - n_pass
  
  # 非NA値を集積（全体分布ログ用）
  nonna_vals <- df$GDC_FILTER[!is.na(df$GDC_FILTER) &
                                trimws(df$GDC_FILTER) != "" &
                                tolower(trimws(df$GDC_FILTER)) != "na"]
  
  df_pass <- df %>% filter(filter_pass)
  
  # -----------------------------------------------------------------------
  # TP53 + IDH1/IDH2 の行のみ抽出
  # -----------------------------------------------------------------------
  target_rows <- df_pass %>%
    filter(Hugo_Symbol %in% c("TP53", "IDH1", "IDH2")) %>%
    select(
      Hugo_Symbol, Variant_Classification, Variant_Type,
      Tumor_Sample_Barcode, HGVSp_Short, GDC_FILTER,
      any_of(c("case_id", "callers"))
    ) %>%
    mutate(
      file_id        = file_id,
      n_total_maf    = n_total,
      n_filter_pass  = n_pass,
      n_filter_fail  = n_fail
    )
  
  list(
    status             = "ok",
    note               = sprintf("total=%d, pass=%d, fail=%d, target_rows=%d",
                                 n_total, n_pass, n_fail, nrow(target_rows)),
    gdc_filter_values  = nonna_vals,   # 全体分布ログ用（非NA・非空値）
    variants_df        = target_rows
  )
}

# =============================================================================
# 4. 全 file_id に対して MAF 読み込み実行
# =============================================================================

log_msg("--- MAF 読み込み開始 ---")
n_total_files <- length(wxs_file_ids)

file_status_rows <- vector("list", n_total_files)
all_variant_rows <- vector("list", n_total_files)

# GDC_FILTER の全体集計用（全ファイル合算）
gdc_filter_all_values <- character(0)

for (i in seq_along(wxs_file_ids)) {
  fid    <- wxs_file_ids[i]
  result <- read_maf_variants(fid)
  
  file_status_rows[[i]] <- data.frame(
    file_id     = fid,
    maf_status  = result$status,
    maf_note    = result$note,
    stringsAsFactors = FALSE
  )
  all_variant_rows[[i]] <- result$variants_df
  
  # GDC_FILTER 生値を集積（非NA値のみ）
  if (!is.null(result$gdc_filter_values)) {
    gdc_filter_all_values <- c(gdc_filter_all_values, result$gdc_filter_values)
  }
  
  if (i %% 100 == 0 || i == n_total_files) {
    n_ok <- sum(sapply(file_status_rows[seq_len(i)],
                       function(r) !is.null(r) && r$maf_status == "ok"))
    log_msg(sprintf("  進捗: %d/%d 完了 (ok=%d, 問題=%d)",
                    i, n_total_files, n_ok, i - n_ok))
  }
}

file_log     <- bind_rows(file_status_rows)
raw_variants <- bind_rows(all_variant_rows)

log_msg(sprintf("MAF 読み込み完了: ok=%d / file_missing=%d / parse_error=%d",
                sum(file_log$maf_status == "ok"),
                sum(file_log$maf_status == "file_missing"),
                sum(file_log$maf_status == "parse_error")))

# ------------------------------------------------------------------
# GDC_FILTER の全体分布ログ（根拠の明示・監査用）
# ------------------------------------------------------------------
log_msg("--- GDC_FILTER 全体分布（全ファイル合算・フィルタ前）---")
log_msg(sprintf("  非NA値の総数: %d件", length(gdc_filter_all_values)))
if (length(gdc_filter_all_values) == 0) {
  log_msg("  → 全ファイルで GDC_FILTER は一貫して NA（値なし）")
  log_msg("  → このデータセットでは GDC_FILTER=NA を pass-equivalent として扱う")
  log_msg("  （根拠: aliquot_ensemble_masked MAF はキュレーション済み出力のため、")
  log_msg("         NA はフィルタ指定なし＝採用相当と解釈。GDC_FILTER_ALLOW_NA=TRUE）")
} else {
  tbl <- sort(table(gdc_filter_all_values), decreasing = TRUE)
  log_msg(sprintf("  非NA値のユニーク値 (%d種):", length(tbl)))
  for (nm in names(head(tbl, 15))) {
    log_msg(sprintf("    %-40s: %d件", nm, tbl[nm]))
  }
  log_msg("  → 非NA値が存在するため GDC_FILTER_ALLOW_NA の影響を確認してください")
}

log_msg(sprintf("抽出変異行数（TP53/IDH1/IDH2）: %d行", nrow(raw_variants)))

# =============================================================================
# 5. TP53 変異サマリー（症例単位）
# =============================================================================

log_msg("--- TP53 変異サマリー作成 ---")

tp53_variants <- raw_variants %>%
  filter(Hugo_Symbol == "TP53",
         Variant_Classification %in% TP53_NONSILENT)

# Tumor_Sample_Barcode → file_id 経由でペアに結合するため、
# まず file_id 単位でサマリーを作る
tp53_by_file <- tp53_variants %>%
  group_by(file_id) %>%
  summarise(
    tp53_mut            = TRUE,
    tp53_n_variants     = n(),
    tp53_classifications = paste(sort(unique(Variant_Classification)),
                                 collapse = "|"),
    tp53_HGVSp          = paste(sort(unique(HGVSp_Short[!is.na(HGVSp_Short) &
                                                          HGVSp_Short != ""])),
                                collapse = "|"),
    tp53_Tumor_Sample_Barcode = first(Tumor_Sample_Barcode),
    .groups = "drop"
  )

log_msg(sprintf("TP53 変異あり（非サイレント）: %d / %d ファイル",
                nrow(tp53_by_file), sum(file_log$maf_status == "ok")))

# HGVSp 別頻度（上位20件）
tp53_hgvsp_freq <- tp53_variants %>%
  filter(!is.na(HGVSp_Short), HGVSp_Short != "") %>%
  count(HGVSp_Short) %>%
  arrange(desc(n))

log_msg("TP53 HGVSp 上位20件:")
for (i in seq_len(min(20, nrow(tp53_hgvsp_freq)))) {
  log_msg(sprintf("  %-20s: %d件", tp53_hgvsp_freq$HGVSp_Short[i],
                  tp53_hgvsp_freq$n[i]))
}

# =============================================================================
# 6. IDH1/IDH2 変異サマリー（症例単位）
# =============================================================================

log_msg("--- IDH1/IDH2 変異サマリー作成 ---")

idh_variants <- raw_variants %>%
  filter(Hugo_Symbol %in% c("IDH1", "IDH2"),
         Variant_Classification == "Missense_Mutation") %>%
  mutate(
    is_hotspot = HGVSp_Short %in% IDH_HOTSPOT
  )

idh_by_file <- idh_variants %>%
  group_by(file_id) %>%
  summarise(
    idh_mut         = TRUE,
    idh_gene        = paste(sort(unique(Hugo_Symbol)), collapse = "|"),
    idh_n_variants  = n(),
    idh_HGVSp       = paste(sort(unique(HGVSp_Short[!is.na(HGVSp_Short)])),
                            collapse = "|"),
    idh_is_hotspot  = any(is_hotspot, na.rm = TRUE),
    .groups = "drop"
  )

log_msg(sprintf("IDH 変異あり（Missense）: %d / %d ファイル",
                nrow(idh_by_file), sum(file_log$maf_status == "ok")))

# IDH HGVSp 内訳
idh_hgvsp_freq <- idh_variants %>%
  filter(!is.na(HGVSp_Short)) %>%
  count(Hugo_Symbol, HGVSp_Short) %>%
  arrange(Hugo_Symbol, desc(n))

log_msg("IDH HGVSp 内訳:")
for (i in seq_len(nrow(idh_hgvsp_freq))) {
  log_msg(sprintf("  [%-5s] %-20s: %d件",
                  idh_hgvsp_freq$Hugo_Symbol[i],
                  idh_hgvsp_freq$HGVSp_Short[i],
                  idh_hgvsp_freq$n[i]))
}

# =============================================================================
# 7. ペアリングテーブルへの結合・最終テーブル作成
# =============================================================================

log_msg("--- ペアリングテーブルへの結合 ---")

# all_pairs に file_id ステータスを付加
mutation_combined <- all_pairs %>%
  left_join(file_log %>% rename(wxs_file_id = file_id),
            by = "wxs_file_id") %>%
  # TP53
  left_join(tp53_by_file %>% rename(wxs_file_id = file_id),
            by = "wxs_file_id") %>%
  mutate(
    tp53_mut        = replace_na(tp53_mut, FALSE),
    tp53_status     = case_when(
      maf_status != "ok"  ~ "maf_error",
      tp53_mut            ~ "mutant",
      TRUE                ~ "wildtype"
    )
  ) %>%
  # IDH
  left_join(idh_by_file %>% rename(wxs_file_id = file_id),
            by = "wxs_file_id") %>%
  mutate(
    idh_mut     = replace_na(idh_mut, FALSE),
    idh_status  = case_when(
      maf_status != "ok"  ~ "maf_error",
      idh_mut             ~ "mutant",
      TRUE                ~ "wildtype"
    )
  )

# TP53 サブグループ分類（アミノ酸変異ベース）
# 主要ホットスポットを個別に、それ以外は分類でまとめる
mutation_combined <- mutation_combined %>%
  mutate(
    tp53_subgroup = case_when(
      tp53_status == "maf_error" ~ "maf_error",
      tp53_status == "wildtype"  ~ "WT",
      # 主要ホットスポット（頻度上位の変異を個別分類）
      str_detect(tp53_HGVSp, "p\\.R175H") ~ "p.R175H",
      str_detect(tp53_HGVSp, "p\\.R248W") ~ "p.R248W",
      str_detect(tp53_HGVSp, "p\\.R248Q") ~ "p.R248Q",
      str_detect(tp53_HGVSp, "p\\.R273H") ~ "p.R273H",
      str_detect(tp53_HGVSp, "p\\.R273C") ~ "p.R273C",
      str_detect(tp53_HGVSp, "p\\.R249S") ~ "p.R249S",
      str_detect(tp53_HGVSp, "p\\.G245S") ~ "p.G245S",
      str_detect(tp53_HGVSp, "p\\.R282W") ~ "p.R282W",
      str_detect(tp53_HGVSp, "p\\.H179R") ~ "p.H179R",
      # 分類別グループ
      str_detect(tp53_classifications, "Nonsense_Mutation|Frame_Shift")
      ~ "Truncating",
      str_detect(tp53_classifications, "Splice_Site")
      ~ "Splice_Site",
      str_detect(tp53_classifications, "Missense_Mutation")
      ~ "Missense_other",
      TRUE                               ~ "Other_nonsilent"
    ),
    # IDH サブグループ
    idh_subgroup = case_when(
      idh_status == "maf_error"  ~ "maf_error",
      idh_status == "wildtype"   ~ "WT",
      str_detect(idh_gene, "IDH1") & str_detect(idh_HGVSp, "p\\.R132H") ~ "IDH1_R132H",
      str_detect(idh_gene, "IDH1") ~ "IDH1_other",
      str_detect(idh_gene, "IDH2") & str_detect(idh_HGVSp, "p\\.R172K") ~ "IDH2_R172K",
      str_detect(idh_gene, "IDH2") ~ "IDH2_other",
      TRUE ~ "IDH_other"
    )
  )

# 集計ログ
log_msg("--- mutation_combined サマリー ---")
log_msg(sprintf("総ペア数: %d件", nrow(mutation_combined)))

tp53_tbl <- mutation_combined %>% count(grade, tp53_status) %>% arrange(grade, tp53_status)
log_msg("TP53 WT/Mut by Grade:")
for (i in seq_len(nrow(tp53_tbl))) {
  r <- tp53_tbl[i, ]
  log_msg(sprintf("  [%s] %-10s: %d件", r$grade, r$tp53_status, r$n))
}

idh_tbl <- mutation_combined %>% count(grade, idh_status) %>% arrange(grade, idh_status)
log_msg("IDH WT/Mut by Grade:")
for (i in seq_len(nrow(idh_tbl))) {
  r <- idh_tbl[i, ]
  log_msg(sprintf("  [%s] %-10s: %d件", r$grade, r$idh_status, r$n))
}

tp53_sub_tbl <- mutation_combined %>%
  count(grade, tp53_subgroup) %>% arrange(grade, desc(n))
log_msg("TP53 サブグループ by Grade:")
for (i in seq_len(nrow(tp53_sub_tbl))) {
  r <- tp53_sub_tbl[i, ]
  log_msg(sprintf("  [%s] %-20s: %d件", r$grade, r$tp53_subgroup, r$n))
}

# =============================================================================
# 8. 出力
# =============================================================================

log_msg("--- ファイル出力 ---")

# 症例単位の TP53 変異テーブル
write_csv(
  mutation_combined %>% select(
    any_of(c("case_barcode","case_id","sample_id","grade","source",
             "wxs_file_id","wxs_sample_id")),
    tp53_status, tp53_subgroup, tp53_n_variants,
    tp53_classifications, tp53_HGVSp
  ),
  file.path(OUT_DIR, "tp53_mutation_table.csv")
)
log_msg(sprintf("保存: tp53_mutation_table.csv (%d行)",  nrow(mutation_combined)))

# 症例単位の IDH 変異テーブル
write_csv(
  mutation_combined %>% select(
    any_of(c("case_barcode","case_id","sample_id","grade","source",
             "wxs_file_id","wxs_sample_id")),
    idh_status, idh_subgroup, idh_gene, idh_n_variants,
    idh_HGVSp, idh_is_hotspot
  ),
  file.path(OUT_DIR, "idh_mutation_table.csv")
)
log_msg(sprintf("保存: idh_mutation_table.csv (%d行)", nrow(mutation_combined)))

# 解析用結合テーブル
write_csv(mutation_combined, file.path(OUT_DIR, "mutation_combined.csv"))
log_msg(sprintf("保存: mutation_combined.csv (%d行 × %d列)",
                nrow(mutation_combined), ncol(mutation_combined)))

# 変異行レベルの生データ（監査用）
write_csv(raw_variants, file.path(OUT_DIR, "mutation_raw_variants.csv"))
log_msg(sprintf("保存: mutation_raw_variants.csv (%d行)", nrow(raw_variants)))

log_msg("=== Step 07: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 07 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    tp53_mutation_table.csv     ← TP53 WT/Mut・サブグループ\n")
cat("    idh_mutation_table.csv      ← IDH1/2 WT/Mut・サブグループ\n")
cat("    mutation_combined.csv       ← 解析用結合テーブル\n")
cat("    mutation_raw_variants.csv   ← 変異行レベル生データ（監査用）\n")
cat("    step07_log.txt\n")
cat("============================\n")
