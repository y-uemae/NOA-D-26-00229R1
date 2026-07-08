# =============================================================================
# step16b_checkpoint_regression.R
# GBM/Glioma TP53×LAG3 解析 - Step 16b: チェックポイント特異性 回帰解析
#
# 目的:
#   Step16a（Wilcoxon/Cliff's δ）の補完として、共変量調整後も
#   LAG3のみがTP53変異と独立して関連するかを線形回帰で確認する。
#   "ノンパラ + 回帰の二枚腰"でチェックポイント特異性の査読耐性を高める。
#
# 解析対象遺伝子（7遺伝子・事前固定）:
#   LAG3, PDCD1, CTLA4, TIGIT, HAVCR2, CD274, PDCD1LG2
#
# モデル:
#   [GDC Grade4]  gene_log2tpm ~ tp53_bin + source_bin + idh_bin
#   [GLASS WXS]   gene_log2tpm ~ tp53_bin
#   ※ tp53_bin: Mut=1/WT=0（引継書 4章 参照カテゴリ定義に準拠）
#   ※ source_bin: CPTAC_HCMI=1/TCGA=0
#   ※ idh_bin: IDH Mut=1/IDH WT=0
#
# 多重補正: BH（コホート内・7遺伝子）
#
# 出力先: 16_checkpoint_specificity/
#   step16b_regression_results.csv   主結果
#   step16b_log.txt                  実行ログ
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(broom)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "16_checkpoint_specificity")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TARGET_GENES <- c("LAG3", "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2")

GDC_PATH   <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH  <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")
GLASS_PATH <- file.path(RESULT_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step16b_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 16b: チェックポイント特異性 回帰解析 開始 ===")
log_msg(sprintf("解析遺伝子: %s", paste(TARGET_GENES, collapse = ", ")))
log_msg("モデル [GDC]: gene ~ tp53_bin + source_bin + idh_bin")
log_msg("モデル [GLASS]: gene ~ tp53_bin")
log_msg("変数定義: tp53_bin Mut=1/WT=0, source_bin CPTAC_HCMI=1/TCGA=0, idh_bin IDH_Mut=1/IDH_WT=0")

# =============================================================================
# 2. データ読み込み・結合（Step16aと同じロジック）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

wide <- read_csv(WIDE_PATH, show_col_types = FALSE)

cp_log2_cols <- paste0(TARGET_GENES, "_log2tpm")
missing_cols <- setdiff(cp_log2_cols, names(wide))
if (length(missing_cols) > 0) {
  log_msg(sprintf("ERROR: wide に列が不足: %s", paste(missing_cols, collapse = ", ")))
  close(log_con); stop("Step06b未実行")
}

# GDC結合（TCGA: case_barcode、CPTAC: wxs_sample_id）
wide_sub <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id", "pair_id")),
         all_of(cp_log2_cols))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_sub %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".wide"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_sub %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".wide"))

gdc <- bind_rows(gdc_tcga, gdc_cptac)

# wide側の重複列をマージ
for (col in cp_log2_cols) {
  wide_col <- paste0(col, ".wide")
  if (wide_col %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wide_col]], gdc[[col]])
    gdc[[wide_col]] <- NULL
  }
}

log_msg(sprintf("GDC結合後: %d行（TCGA=%d, CPTAC=%d）",
                nrow(gdc), sum(gdc$source=="TCGA"), sum(gdc$source=="CPTAC_HCMI")))

# GDC Grade4 + 共変量ダミー変数作成
gdc_g4 <- gdc %>%
  filter(grade == "Grade4",
         tp53_status %in% c("mutant", "wildtype"),
         idh_status  %in% c("wildtype", "mutant")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),   # Mut=1 / WT=0
    source_bin = as.integer(source == "CPTAC_HCMI"),    # CPTAC_HCMI=1 / TCGA=0
    idh_bin    = as.integer(idh_status == "mutant")     # IDH Mut=1 / IDH WT=0
  )

log_msg(sprintf("GDC Grade4 解析対象: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$tp53_bin == 1),
                sum(gdc_g4$tp53_bin == 0)))
log_msg(sprintf("  source_bin: CPTAC=%d, TCGA=%d",
                sum(gdc_g4$source_bin == 1), sum(gdc_g4$source_bin == 0)))
log_msg(sprintf("  idh_bin: IDH_Mut=%d, IDH_WT=%d",
                sum(gdc_g4$idh_bin == 1), sum(gdc_g4$idh_bin == 0)))

# GLASS
glass <- read_csv(GLASS_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE,
         tp53_status %in% c("Mut", "WT")) %>%
  mutate(tp53_bin = as.integer(tp53_status == "Mut"))  # Mut=1 / WT=0

log_msg(sprintf("GLASS WXS: n=%d (Mut=%d, WT=%d)",
                nrow(glass),
                sum(glass$tp53_bin == 1),
                sum(glass$tp53_bin == 0)))

# =============================================================================
# 3. 回帰関数
# =============================================================================

run_regression <- function(df, gene, formula_str, cohort_name) {
  col <- paste0(gene, "_log2tpm")
  if (!col %in% names(df)) {
    return(tibble(cohort = cohort_name, gene = gene,
                  term = "tp53_bin", beta = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_,
                  se = NA_real_, p_value = NA_real_,
                  n = NA_integer_, r_squared = NA_real_,
                  note = "column_missing"))
  }
  
  df_use <- df %>% select(all_of(c(col, all.vars(as.formula(formula_str))[-1]))) %>%
    na.omit()
  n_use  <- nrow(df_use)
  
  if (n_use < 10) {
    return(tibble(cohort = cohort_name, gene = gene,
                  term = "tp53_bin", beta = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_,
                  se = NA_real_, p_value = NA_real_,
                  n = n_use, r_squared = NA_real_,
                  note = "insufficient_n"))
  }
  
  fml <- as.formula(paste0("`", col, "`", sub("^[^~]+", "", formula_str)))
  
  fit <- tryCatch(lm(fml, data = df_use), error = function(e) e)
  if (inherits(fit, "error")) {
    return(tibble(cohort = cohort_name, gene = gene,
                  term = "tp53_bin", beta = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_,
                  se = NA_real_, p_value = NA_real_,
                  n = n_use, r_squared = NA_real_,
                  note = conditionMessage(fit)))
  }
  
  ci  <- confint(fit, "tp53_bin", level = 0.95)
  tbl <- tidy(fit) %>% filter(term == "tp53_bin")
  r2  <- summary(fit)$r.squared
  
  tibble(
    cohort    = cohort_name,
    gene      = gene,
    term      = "tp53_bin",
    beta      = round(tbl$estimate, 4),
    ci_lower  = round(ci[1], 4),
    ci_upper  = round(ci[2], 4),
    se        = round(tbl$std.error, 4),
    p_value   = tbl$p.value,
    n         = n_use,
    r_squared = round(r2, 4),
    note      = "ok"
  )
}

# =============================================================================
# 4. 解析実行
# =============================================================================

log_msg("--- 回帰解析実行 ---")

# GDC Grade4：全体・TCGA・CPTAC・IDH_WT
gdc_configs <- list(
  list(data = gdc_g4,                               label = "GDC_Grade4_all",
       fml  = "y ~ tp53_bin + source_bin + idh_bin"),
  list(data = gdc_g4 %>% filter(source == "TCGA"),  label = "GDC_Grade4_TCGA",
       fml  = "y ~ tp53_bin + idh_bin"),
  list(data = gdc_g4 %>% filter(source == "CPTAC_HCMI"), label = "GDC_Grade4_CPTAC",
       fml  = "y ~ tp53_bin + idh_bin"),
  list(data = gdc_g4 %>% filter(idh_status == "wildtype"), label = "GDC_Grade4_IDH_WT",
       fml  = "y ~ tp53_bin + source_bin")
)

all_results <- list()

for (cfg in gdc_configs) {
  log_msg(sprintf("  [%s] モデル: %s, n=%d",
                  cfg$label, cfg$fml, nrow(cfg$data)))
  res_list <- lapply(TARGET_GENES, function(g) {
    run_regression(cfg$data, g, cfg$fml, cfg$label)
  })
  all_results <- c(all_results, res_list)
}

# GLASS WXS
log_msg(sprintf("  [GLASS_WXS] モデル: y ~ tp53_bin, n=%d", nrow(glass)))
glass_fml <- "y ~ tp53_bin"
res_glass <- lapply(TARGET_GENES, function(g) {
  run_regression(glass, g, glass_fml, "GLASS_WXS")
})
all_results <- c(all_results, res_glass)

results_raw <- bind_rows(all_results)
log_msg(sprintf("回帰完了: %d行（コホート×遺伝子）", nrow(results_raw)))

# =============================================================================
# 5. 多重補正（BH法：コホート内・7遺伝子）
# =============================================================================

log_msg("--- BH補正（コホート内・7遺伝子） ---")

results <- results_raw %>%
  group_by(cohort) %>%
  mutate(p_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(significant = p_BH < 0.05) %>%
  arrange(cohort, gene)

# =============================================================================
# 6. 保存
# =============================================================================

write_csv(results, file.path(OUT_DIR, "step16b_regression_results.csv"))
log_msg("保存: step16b_regression_results.csv")

# =============================================================================
# 7. サマリー出力
# =============================================================================

cohort_order <- c("GDC_Grade4_all", "GDC_Grade4_TCGA",
                  "GDC_Grade4_CPTAC", "GDC_Grade4_IDH_WT", "GLASS_WXS")

summary_lines <- c(
  "=== Step 16b: チェックポイント特異性 回帰解析 サマリー ===",
  sprintf("実行日時: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "変数定義: tp53_bin Mut=1/WT=0, source_bin CPTAC_HCMI=1/TCGA=0, idh_bin IDH_Mut=1/IDH_WT=0",
  ""
)

for (co in cohort_order) {
  res_co <- results %>% filter(cohort == co, note == "ok")
  if (nrow(res_co) == 0) next
  
  # モデル文字列
  model_str <- switch(co,
                      GDC_Grade4_all    = "gene ~ tp53_bin + source_bin + idh_bin",
                      GDC_Grade4_TCGA   = "gene ~ tp53_bin + idh_bin",
                      GDC_Grade4_CPTAC  = "gene ~ tp53_bin + idh_bin",
                      GDC_Grade4_IDH_WT = "gene ~ tp53_bin + source_bin",
                      GLASS_WXS         = "gene ~ tp53_bin",
                      "gene ~ tp53_bin + ..."
  )
  
  summary_lines <- c(summary_lines,
                     sprintf("【%s】 モデル: %s  n=%d", co, model_str, res_co$n[1])
  )
  
  for (i in seq_len(nrow(res_co))) {
    r <- res_co[i,]
    sig_mark <- if (!is.na(r$significant) && r$significant) "★" else "  "
    summary_lines <- c(summary_lines,
                       sprintf("  %s %-12s β=%+.4f [%+.4f, %+.4f]  p_BH=%.4f  R²=%.3f",
                               sig_mark, r$gene,
                               r$beta, r$ci_lower, r$ci_upper,
                               ifelse(is.na(r$p_BH), NA, r$p_BH),
                               ifelse(is.na(r$r_squared), NA, r$r_squared))
    )
  }
  summary_lines <- c(summary_lines, "")
}

# LAG3のみ有意かチェック（Step16aと同じ指標）
lag3_sig <- results %>%
  filter(note == "ok") %>%
  group_by(cohort) %>%
  summarise(
    lag3_sig    = any(gene == "LAG3" & significant, na.rm = TRUE),
    other_sig   = any(gene != "LAG3" & significant, na.rm = TRUE),
    n_sig_genes = sum(significant, na.rm = TRUE),
    sig_genes   = paste(gene[significant & !is.na(significant)], collapse = ", "),
    .groups = "drop"
  )

summary_lines <- c(summary_lines, "--- LAG3特異性チェック（回帰） ---")
for (i in seq_len(nrow(lag3_sig))) {
  r <- lag3_sig[i,]
  pattern <- dplyr::case_when(
    r$lag3_sig & !r$other_sig ~ "✅ LAG3のみ有意（特異性あり）",
    r$lag3_sig &  r$other_sig ~ "⚠️  LAG3と他も有意（特異性一部）",
    !r$lag3_sig               ~ "❌ LAG3非有意"
  )
  summary_lines <- c(summary_lines,
                     sprintf("  %-25s %s  有意: [%s]", r$cohort, pattern, r$sig_genes)
  )
}

# Step16aとの比較メモ
summary_lines <- c(summary_lines, "",
                   "--- Step16a（Wilcoxon）との整合確認 ---",
                   "  LAG3: Wilcoxonで有意だったコホートで回帰βも正方向・有意であれば一貫性OK",
                   "  他遺伝子: Wilcoxon非有意と回帰非有意が一致しているか確認"
)

writeLines(summary_lines, file.path(OUT_DIR, "step16b_summary.txt"))
log_msg("保存: step16b_summary.txt")

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n")

# =============================================================================
# 8. ログへの主要結果記録
# =============================================================================

log_msg("--- 主要結果（LAG3 特異性確認・回帰） ---")
for (i in seq_len(nrow(lag3_sig))) {
  r <- lag3_sig[i,]
  log_msg(sprintf("  %s: 有意遺伝子=%d個 [%s]", r$cohort, r$n_sig_genes, r$sig_genes))
}

# LAG3のβをログに記録
lag3_betas <- results %>%
  filter(gene == "LAG3", note == "ok") %>%
  select(cohort, beta, ci_lower, ci_upper, p_BH, significant)

log_msg("LAG3 β across cohorts:")
for (i in seq_len(nrow(lag3_betas))) {
  r <- lag3_betas[i,]
  sig <- if (!is.na(r$significant) && r$significant) "★" else "  "
  log_msg(sprintf("  %s %-25s β=%+.4f [%+.4f, %+.4f]  p_BH=%.4f",
                  sig, r$cohort, r$beta, r$ci_lower, r$ci_upper,
                  ifelse(is.na(r$p_BH), NA, r$p_BH)))
}

log_msg("=== Step 16b: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))
close(log_con)

cat("\n============================\n")
cat("Step 16b 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/step16b_regression_results.csv\n", OUT_DIR))
cat(sprintf("  %s/step16b_summary.txt\n", OUT_DIR))
cat(sprintf("  %s/step16b_log.txt\n", OUT_DIR))
cat("============================\n")
