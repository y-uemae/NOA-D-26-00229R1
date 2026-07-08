# =============================================================================
# step08_final_cohort.R
# GBM/Glioma TP53×LAG3 再解析 - Step 08: 最終コホート構築
#
# 目的:
#   Step 06（遺伝子発現）・Step 07（変異）の結果を
#   WXS-RNA ペアリングテーブルに結合し、解析に使用する最終コホートを確定する。
#
# -----------------------------------------------------------------------
# 入力
#   06_gene_expression/gene_expression_wide.csv       （951行）
#   07_mutations/mutation_combined.csv                （951行）
#   ※ ペア情報は gene_expression_wide に含まれている
#
# -----------------------------------------------------------------------
# 結合キー
#   gene_expression_wide と mutation_combined は共に pair_id を持つが、
#   生成順序の差異を避けるため wxs_file_id × rna_file_id × grade の
#   複合キーで結合する（より安全）。
#
# -----------------------------------------------------------------------
# 除外基準（最終コホートへの採用条件）
#   1. lag3_status == "ok"         : LAG3 発現が正常抽出されている
#   2. maf_status == "ok"          : MAF が正常読み込み済み
#   ※ 上記以外は exclude_reason に記録して selection_table に残す
#   ※ 発現ステータスは LAG3 のみ採用条件とする
#     （他の 27 遺伝子は not_found でも解析対象外にはしない）
#
# -----------------------------------------------------------------------
# 出力
#   08_final_cohort/
#     final_cohort.csv           : 解析採用サンプル（主解析用）
#     selection_table.csv        : 全候補の採否決定表（監査用）
#     cohort_summary.txt         : コホート記述統計（Methods 記載用）
#     step08_log.txt
#
# 作成日: 2026-02-22
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
OUT_DIR    <- file.path(RESULT_DIR, "08_final_cohort")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step08_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 08: 最終コホート構築 開始 ===")
log_msg("採用条件: lag3_status=ok AND maf_status=ok")

# =============================================================================
# 2. 入力読み込み
# =============================================================================

log_msg("--- 入力ファイル読み込み ---")

expr_wide <- read_csv(
  file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv"),
  show_col_types = FALSE)
log_msg(sprintf("  gene_expression_wide: %d行 × %d列", nrow(expr_wide), ncol(expr_wide)))

mut_combined <- read_csv(
  file.path(RESULT_DIR, "07_mutations/mutation_combined.csv"),
  show_col_types = FALSE)
log_msg(sprintf("  mutation_combined: %d行 × %d列", nrow(mut_combined), ncol(mut_combined)))

# =============================================================================
# 3. 結合
# =============================================================================

log_msg("--- テーブル結合 ---")

# Step 06 の long ログから rna_file_id × gene_name=LAG3 の status を取得
lag3_status_df <- read_csv(
  file.path(RESULT_DIR, "06_gene_expression/gene_extraction_log_long.csv"),
  show_col_types = FALSE) %>%
  filter(gene == "LAG3") %>%
  select(rna_file_id = file_id, lag3_status = status) %>%
  distinct(rna_file_id, .keep_all = TRUE)

log_msg(sprintf("  lag3_status 取得: %d件", nrow(lag3_status_df)))

# mutation_combined から発現テーブルに既にある列を除去（重複防止）
mut_cols_to_add <- mut_combined %>%
  select(
    wxs_file_id, rna_file_id, grade,
    maf_status, maf_note,
    tp53_status, tp53_subgroup, tp53_n_variants,
    tp53_classifications, tp53_HGVSp,
    idh_status, idh_subgroup, idh_gene,
    idh_n_variants, idh_HGVSp, idh_is_hotspot
  )

# 結合キー: wxs_file_id × rna_file_id × grade
combined <- expr_wide %>%
  left_join(lag3_status_df,  by = "rna_file_id") %>%
  left_join(mut_cols_to_add, by = c("wxs_file_id", "rna_file_id", "grade"))

log_msg(sprintf("  結合後: %d行 × %d列", nrow(combined), ncol(combined)))

# 結合漏れチェック
n_missing_maf <- sum(is.na(combined$maf_status))
if (n_missing_maf > 0) {
  log_msg(sprintf("  WARNING: maf_status が NA の行: %d件（結合キー不一致の可能性）",
                  n_missing_maf))
} else {
  log_msg("  OK: 全行で maf_status が結合済み")
}

# =============================================================================
# 4. 採否判定
# =============================================================================

log_msg("--- 採否判定 ---")

# 採否判定前に lag3_status / maf_status の分布を確認（採用根拠の明示）
log_msg("採用根拠の確認:")
lag3_tbl <- combined %>% count(lag3_status) %>% arrange(desc(n))
for (i in seq_len(nrow(lag3_tbl))) {
  log_msg(sprintf("  lag3_status %-15s: %d件", lag3_tbl$lag3_status[i], lag3_tbl$n[i]))
}
maf_tbl <- combined %>% count(maf_status) %>% arrange(desc(n))
for (i in seq_len(nrow(maf_tbl))) {
  log_msg(sprintf("  maf_status  %-15s: %d件", maf_tbl$maf_status[i], maf_tbl$n[i]))
}

selection_table <- combined %>%
  mutate(
    exclude_reason = case_when(
      lag3_status != "ok" & maf_status != "ok" ~ sprintf("lag3_%s_AND_maf_%s",
                                                         lag3_status, maf_status),
      lag3_status != "ok" ~ sprintf("lag3_%s", lag3_status),
      maf_status  != "ok" ~ sprintf("maf_%s",  maf_status),
      TRUE                ~ "included_ok"
    ),
    include_flag = (exclude_reason == "included_ok"),
    decision     = if_else(include_flag, "include", "exclude")
  )

# 採否内訳ログ
incl_tbl <- selection_table %>%
  count(decision, exclude_reason) %>%
  arrange(decision, exclude_reason)

log_msg("採否内訳:")
for (i in seq_len(nrow(incl_tbl))) {
  r <- incl_tbl[i, ]
  log_msg(sprintf("  %-7s %-40s: %d件",
                  r$decision, r$exclude_reason, r$n))
}

n_adopted <- sum(selection_table$include_flag)
log_msg(sprintf("最終採用: %d / %d件", n_adopted, nrow(selection_table)))

# =============================================================================
# 5. 最終コホートテーブルの整形
# =============================================================================

log_msg("--- 最終コホートテーブル整形 ---")

# 遺伝子列の順序を整理（メタ列 → TPM列 → log2TPM列 → ステータス列）
TARGET_GENES <- c(
  "B2M","TAP1","TAP2","TAPBP","HLA-A","HLA-B","HLA-C","NLRC5",
  "STAT1","IRF1","IRF9",
  "CXCL9","CXCL10","CXCL11",
  "GBP1","GBP2","GBP4","GBP5",
  "IDO1",
  "CD3D","CD3E","CD3G","CD8A","CD8B",
  "GZMA","GZMB","PRF1",
  "LAG3"
)

meta_cols <- c(
  "pair_id", "case_barcode", "grade", "source", "project_id",
  "wxs_sample_id", "rna_sample_id", "wxs_file_id", "rna_file_id",
  "match_type", "is_sensitivity_only",
  "sample_id", "case_id_raw", "run_id"
)
tpm_cols   <- paste0(TARGET_GENES, "_tpm")
log2_cols  <- paste0(TARGET_GENES, "_log2tpm")
status_cols <- c("n_ok","n_not_found","n_duplicated","n_parse_error",
                 "lag3_status",
                 "maf_status",
                 "tp53_status","tp53_subgroup","tp53_n_variants",
                 "tp53_classifications","tp53_HGVSp",
                 "idh_status","idh_subgroup","idh_gene",
                 "idh_n_variants","idh_HGVSp","idh_is_hotspot")

# 存在する列のみ選択
select_cols <- c(
  intersect(meta_cols,   names(selection_table)),
  intersect(tpm_cols,    names(selection_table)),
  intersect(log2_cols,   names(selection_table)),
  intersect(status_cols, names(selection_table)),
  "include_flag", "exclude_reason"
)

selection_table_out <- selection_table %>% select(all_of(select_cols))
final_cohort        <- selection_table_out %>% filter(include_flag)

log_msg(sprintf("final_cohort: %d行 × %d列", nrow(final_cohort), ncol(final_cohort)))

# =============================================================================
# 6. コホート記述統計
# =============================================================================

log_msg("--- コホート記述統計 ---")

summary_lines <- character(0)
add_line <- function(x) {
  summary_lines <<- c(summary_lines, x)
  log_msg(x)
}

add_line("=== コホート記述統計（Methods 記載用）===")
add_line(sprintf("解析採用ペア総数: %d件", nrow(final_cohort)))
add_line("")

# Grade × Source 内訳
grade_src <- final_cohort %>%
  count(grade, source) %>%
  arrange(grade, source)
add_line("Grade × Source 内訳:")
for (i in seq_len(nrow(grade_src))) {
  r <- grade_src[i, ]
  add_line(sprintf("  [%-7s][%-10s]: %d件", r$grade, r$source, r$n))
}
add_line("")

# Grade × TP53 WT/Mut
tp53_grade <- final_cohort %>%
  count(grade, tp53_status) %>%
  arrange(grade, tp53_status)
add_line("Grade × TP53 WT/Mut:")
for (i in seq_len(nrow(tp53_grade))) {
  r <- tp53_grade[i, ]
  add_line(sprintf("  [%-7s] %-10s: %d件", r$grade, r$tp53_status, r$n))
}
add_line("")

# Grade × IDH WT/Mut
idh_grade <- final_cohort %>%
  count(grade, idh_status) %>%
  arrange(grade, idh_status)
add_line("Grade × IDH WT/Mut:")
for (i in seq_len(nrow(idh_grade))) {
  r <- idh_grade[i, ]
  add_line(sprintf("  [%-7s] %-10s: %d件", r$grade, r$idh_status, r$n))
}
add_line("")

# TP53 × IDH の組み合わせ（Grade4 のみ）
add_line("TP53 × IDH 組み合わせ（Grade4）:")
cross_g4 <- final_cohort %>%
  filter(grade == "Grade4") %>%
  count(tp53_status, idh_status) %>%
  arrange(tp53_status, idh_status)
for (i in seq_len(nrow(cross_g4))) {
  r <- cross_g4[i, ]
  add_line(sprintf("  TP53=%-8s × IDH=%-8s: %d件",
                   r$tp53_status, r$idh_status, r$n))
}
add_line("")

# LAG3 log2(TPM+1) 記述統計
add_line("LAG3 log2(TPM+1) 記述統計:")
lag3_stats <- final_cohort %>%
  group_by(grade, source) %>%
  summarise(
    n      = n(),
    median = round(median(LAG3_log2tpm, na.rm = TRUE), 3),
    mean   = round(mean(LAG3_log2tpm,   na.rm = TRUE), 3),
    sd     = round(sd(LAG3_log2tpm,     na.rm = TRUE), 3),
    min    = round(min(LAG3_log2tpm,    na.rm = TRUE), 3),
    max    = round(max(LAG3_log2tpm,    na.rm = TRUE), 3),
    .groups = "drop"
  )
for (i in seq_len(nrow(lag3_stats))) {
  r <- lag3_stats[i, ]
  add_line(sprintf("  [%-7s][%-10s] n=%d, median=%.3f, mean±SD=%.3f±%.3f, range=[%.3f,%.3f]",
                   r$grade, r$source, r$n,
                   r$median, r$mean, r$sd, r$min, r$max))
}
add_line("")

# TP53 サブグループ全体
add_line("TP53 サブグループ内訳（全 Grade）:")
tp53_sub_all <- final_cohort %>%
  count(tp53_subgroup) %>%
  arrange(desc(n))
for (i in seq_len(nrow(tp53_sub_all))) {
  add_line(sprintf("  %-20s: %d件", tp53_sub_all$tp53_subgroup[i], tp53_sub_all$n[i]))
}
add_line("")

# IDH サブグループ全体
add_line("IDH サブグループ内訳（全 Grade）:")
idh_sub_all <- final_cohort %>%
  count(idh_subgroup) %>%
  arrange(desc(n))
for (i in seq_len(nrow(idh_sub_all))) {
  add_line(sprintf("  %-20s: %d件", idh_sub_all$idh_subgroup[i], idh_sub_all$n[i]))
}

# =============================================================================
# 7. 出力
# =============================================================================

log_msg("--- ファイル出力 ---")

write_csv(final_cohort, file.path(OUT_DIR, "final_cohort.csv"))
log_msg(sprintf("保存: final_cohort.csv (%d行 × %d列)",
                nrow(final_cohort), ncol(final_cohort)))

write_csv(selection_table_out, file.path(OUT_DIR, "selection_table.csv"))
log_msg(sprintf("保存: selection_table.csv (%d行 × %d列)",
                nrow(selection_table_out), ncol(selection_table_out)))

writeLines(summary_lines, file.path(OUT_DIR, "cohort_summary.txt"))
log_msg("保存: cohort_summary.txt")

log_msg("=== Step 08: 完了 ===")
close(log_con)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

cat("\n============================\n")
cat("Step 08 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/\n", OUT_DIR))
cat("    final_cohort.csv       ← 解析採用サンプル（主解析用）\n")
cat("    selection_table.csv    ← 全候補の採否決定表（監査用）\n")
cat("    cohort_summary.txt     ← 記述統計（Methods 記載用）\n")
cat("    step08_log.txt\n")
cat("============================\n")
