# =============================================================================
# step16a_checkpoint_specificity_v2.R
# GBM/Glioma TP53×LAG3 解析 - Step 16a v2: チェックポイント特異性解析
#
# v2 変更点:
#   1. 出力ファイル名を step16a_v2_* に変更（既存ファイルを上書きしない）
#   2. p_wilcox（生p値）を結果CSVに明示的に保持
#   3. GLASS コホートで tied ranks 診断を追加出力
#   4. ログに全遺伝子の p_raw / p_BH を記録（LAG3以外のp値を可視化）
#   5. p_raw = 1.0 になった遺伝子に tied_warning フラグを付与
#   6. サマリーに tied 遺伝子数と注釈を追記
#
# 診断結果（v1の問題）:
#   GLASS_WXS コホートでは、低発現遺伝子（CTLA4, PDCD1LG2等）の
#   log2(TPM+1) が多数ゼロでほぼ完全な tied ranks になり、
#   wilcox.test(exact=FALSE, correct=TRUE) が p=1.0 を返す。
#   BH補正後も複数遺伝子が同じ p_BH=1.0 になる。
#   これはコードバグではなく、統計的・生物学的な現象（検出力不足）。
#   → v2 では tied 診断を追加し、結果の解釈に必要な情報を明示する。
#
# 作成日: 2026-02-24（v1）→ v2 修正
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "16_checkpoint_specificity")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TARGET_GENES <- c("LAG3", "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2")

GDC_PATH   <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
GLASS_PATH <- file.path(RESULT_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")
WIDE_PATH  <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step16a_v2_log.txt")   # ★ v2
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 16a v2: チェックポイント特異性解析 開始 ===")
log_msg(sprintf("解析遺伝子: %s", paste(TARGET_GENES, collapse = ", ")))
log_msg(sprintf("出力先: %s", OUT_DIR))
log_msg("v2変更点: tied診断追加・全遺伝子p_raw/p_BHをログ出力・ファイル名をv2に変更")

# =============================================================================
# 2. 補助関数
# =============================================================================

calc_cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  if (nx == 0 || ny == 0) return(NA_real_)
  mat <- outer(x, y, FUN = function(a, b) sign(a - b))
  sum(mat) / (nx * ny)
}

calc_hl <- function(x, y) {
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  diffs <- outer(x, y, FUN = "-")
  median(diffs, na.rm = TRUE)
}

# ★ v2追加: Tied ranks 診断
# ゼロ値（log2(TPM+1)≦ε）の割合と、全体のtied率を返す
calc_tied_diagnostics <- function(x_mut, x_wt, eps = 1e-6) {
  x_all    <- c(x_mut, x_wt)
  zero_rate_mut <- mean(x_mut <= eps)
  zero_rate_wt  <- mean(x_wt  <= eps)
  zero_rate_all <- mean(x_all <= eps)
  
  # Wilcoxon W統計量の期待値（帰無仮説下）
  n_mut <- length(x_mut); n_wt <- length(x_wt)
  W_expected <- n_mut * n_wt / 2
  
  # 実際のW統計量
  W_obs <- tryCatch(
    wilcox.test(x_mut, x_wt, exact = FALSE, correct = TRUE)$statistic,
    error = function(e) NA_real_
  )
  
  # W が期待値に近いほど tied の影響が大きい
  W_deviation <- if (!is.na(W_obs)) abs(W_obs - W_expected) / W_expected else NA_real_
  
  list(
    zero_rate_mut = round(zero_rate_mut, 4),
    zero_rate_wt  = round(zero_rate_wt,  4),
    zero_rate_all = round(zero_rate_all, 4),
    W_expected    = W_expected,
    W_obs         = W_obs,
    W_deviation   = round(W_deviation, 4)
  )
}

# =============================================================================
# 3. データ読み込み・前処理（v1と同一）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
log_msg(sprintf("GDC final_cohort: %d行", nrow(gdc_base)))

wide <- read_csv(WIDE_PATH, show_col_types = FALSE)
log_msg(sprintf("GDC wide: %d行 × %d列", nrow(wide), ncol(wide)))

cp_log2_cols <- paste0(TARGET_GENES, "_log2tpm")
missing_cols <- setdiff(cp_log2_cols, names(wide))
if (length(missing_cols) > 0) {
  log_msg(sprintf("ERROR: wide に列が不足: %s", paste(missing_cols, collapse = ", ")))
  close(log_con); stop("必要な発現量列が存在しません")
}
log_msg("✅ GDC: チェックポイント列確認OK")

id_cols_wide <- intersect(c("case_barcode", "wxs_sample_id", "pair_id"), names(wide))
log_msg(sprintf("wideのID列: %s", paste(id_cols_wide, collapse = ", ")))

expr_cols <- c(id_cols_wide, cp_log2_cols)
wide_sub  <- wide %>% select(all_of(intersect(expr_cols, names(wide))))

gdc_tcga  <- gdc_base %>%
  filter(source == "TCGA") %>%
  left_join(wide_sub %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".wide"))

gdc_cptac <- gdc_base %>%
  filter(source == "CPTAC_HCMI") %>%
  left_join(wide_sub %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".wide"))

gdc <- bind_rows(gdc_tcga, gdc_cptac)
log_msg(sprintf("GDC結合後: %d行（TCGA=%d, CPTAC=%d）",
                nrow(gdc), sum(gdc$source=="TCGA"), sum(gdc$source=="CPTAC_HCMI")))

for (col in cp_log2_cols) {
  wide_col <- paste0(col, ".wide")
  if (wide_col %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wide_col]], gdc[[col]])
    gdc[[wide_col]] <- NULL
  }
}

log_msg("GDC発現量 NA率確認:")
for (col in cp_log2_cols) {
  na_r <- mean(is.na(gdc[[col]])) * 100
  log_msg(sprintf("  %s NA率=%.1f%%", col, na_r))
}

glass <- read_csv(GLASS_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
log_msg(sprintf("GLASS WXS(n=79): %d行（include_flag==TRUE）", nrow(glass)))
log_msg(sprintf("GLASS TP53 Mut=%d, WT=%d",
                sum(glass$tp53_status=="Mut", na.rm=TRUE),
                sum(glass$tp53_status=="WT",  na.rm=TRUE)))

glass_missing <- setdiff(cp_log2_cols, names(glass))
if (length(glass_missing) > 0) {
  log_msg(sprintf("ERROR: GLASSファイルに列が不足: %s", paste(glass_missing, collapse=", ")))
  close(log_con); stop("Step06b.1を先に実行してください")
}
log_msg("✅ GLASS: チェックポイント列確認OK")

# =============================================================================
# 4. コホート定義
# =============================================================================

log_msg("--- コホート定義 ---")

cohorts <- list(
  GDC_Grade4_all = gdc %>%
    filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")),
  GDC_Grade4_TCGA = gdc %>%
    filter(grade == "Grade4", source == "TCGA",
           tp53_status %in% c("mutant", "wildtype")),
  GDC_Grade4_CPTAC = gdc %>%
    filter(grade == "Grade4", source == "CPTAC_HCMI",
           tp53_status %in% c("mutant", "wildtype")),
  GDC_Grade4_IDH_WT = gdc %>%
    filter(grade == "Grade4", idh_status == "wildtype",
           tp53_status %in% c("mutant", "wildtype")),
  GLASS_WXS = glass %>%
    filter(tp53_status %in% c("Mut", "WT"))
)

for (nm in names(cohorts)) {
  df     <- cohorts[[nm]]
  mut_val <- if (grepl("^GLASS", nm)) "Mut" else "mutant"
  wt_val  <- if (grepl("^GLASS", nm)) "WT"  else "wildtype"
  n_mut  <- sum(df$tp53_status == mut_val, na.rm = TRUE)
  n_wt   <- sum(df$tp53_status == wt_val,  na.rm = TRUE)
  log_msg(sprintf("  %-25s n=%d (Mut=%d, WT=%d)", nm, nrow(df), n_mut, n_wt))
}

# =============================================================================
# 5. 解析関数（★ v2: tied診断を追加）
# =============================================================================

run_checkpoint_analysis <- function(df, cohort_name) {
  
  is_glass <- grepl("^GLASS", cohort_name)
  mut_val  <- if (is_glass) "Mut" else "mutant"
  wt_val   <- if (is_glass) "WT"  else "wildtype"
  
  df_use <- df %>%
    filter(tp53_status %in% c(mut_val, wt_val)) %>%
    mutate(tp53_bin = ifelse(tp53_status == mut_val, "Mut", "WT"))
  
  n_mut <- sum(df_use$tp53_bin == "Mut")
  n_wt  <- sum(df_use$tp53_bin == "WT")
  
  results <- lapply(TARGET_GENES, function(gene) {
    col <- paste0(gene, "_log2tpm")
    
    if (!col %in% names(df_use)) {
      return(tibble(cohort = cohort_name, gene = gene,
                    n_mut = n_mut, n_wt = n_wt,
                    median_mut = NA_real_, median_wt = NA_real_,
                    median_diff = NA_real_, hl_estimate = NA_real_,
                    cliffs_delta = NA_real_,
                    p_wilcox = NA_real_,         # ★ v2: 生p値
                    zero_rate_mut = NA_real_,    # ★ v2
                    zero_rate_wt  = NA_real_,    # ★ v2
                    zero_rate_all = NA_real_,    # ★ v2
                    tied_warning  = FALSE,       # ★ v2
                    note = "column_missing"))
    }
    
    x_mut <- df_use %>% filter(tp53_bin == "Mut") %>% pull(!!col) %>%
      na.omit() %>% as.numeric()
    x_wt  <- df_use %>% filter(tp53_bin == "WT")  %>% pull(!!col) %>%
      na.omit() %>% as.numeric()
    
    n_mut_valid <- length(x_mut)
    n_wt_valid  <- length(x_wt)
    
    if (n_mut_valid < 3 || n_wt_valid < 3) {
      return(tibble(cohort = cohort_name, gene = gene,
                    n_mut = n_mut_valid, n_wt = n_wt_valid,
                    median_mut = median(x_mut), median_wt = median(x_wt),
                    median_diff = NA_real_, hl_estimate = NA_real_,
                    cliffs_delta = NA_real_,
                    p_wilcox = NA_real_,
                    zero_rate_mut = NA_real_, zero_rate_wt = NA_real_,
                    zero_rate_all = NA_real_, tied_warning = FALSE,
                    note = "insufficient_n"))
    }
    
    # Wilcoxon検定
    wt_result <- tryCatch(
      wilcox.test(x_mut, x_wt, exact = FALSE, correct = TRUE),
      error = function(e) list(p.value = NA_real_)
    )
    p_raw <- wt_result$p.value
    
    # ★ v2: Tied診断
    tied_diag <- calc_tied_diagnostics(x_mut, x_wt)
    
    # ★ v2: p=1.0（または極めて1.0に近い）かつゼロ率が高い場合に tied_warning=TRUE
    tied_warning <- (!is.na(p_raw) && p_raw > 0.999) ||
      (tied_diag$zero_rate_all > 0.70)
    
    tibble(
      cohort        = cohort_name,
      gene          = gene,
      n_mut         = n_mut_valid,
      n_wt          = n_wt_valid,
      median_mut    = round(median(x_mut), 4),
      median_wt     = round(median(x_wt),  4),
      median_diff   = round(median(x_mut) - median(x_wt), 4),
      hl_estimate   = round(calc_hl(x_mut, x_wt), 4),
      cliffs_delta  = round(calc_cliffs_delta(x_mut, x_wt), 4),
      p_wilcox      = p_raw,                          # ★ v2: 生p値を明示保持
      zero_rate_mut = tied_diag$zero_rate_mut,        # ★ v2
      zero_rate_wt  = tied_diag$zero_rate_wt,         # ★ v2
      zero_rate_all = tied_diag$zero_rate_all,        # ★ v2
      tied_warning  = tied_warning,                   # ★ v2
      note          = "ok"
    )
  })
  
  bind_rows(results)
}

# =============================================================================
# 6. 全コホートで解析実行
# =============================================================================

log_msg("--- 解析実行 ---")

all_results <- lapply(names(cohorts), function(nm) {
  log_msg(sprintf("  解析中: %s", nm))
  run_checkpoint_analysis(cohorts[[nm]], nm)
})

results_raw <- bind_rows(all_results)
log_msg(sprintf("解析完了: %d行（コホート×遺伝子）", nrow(results_raw)))

# =============================================================================
# 7. 多重補正（BH法：コホート内で7遺伝子）
# =============================================================================

log_msg("--- BH補正（コホート内・7遺伝子） ---")

results <- results_raw %>%
  group_by(cohort) %>%
  mutate(
    p_BH        = p.adjust(p_wilcox, method = "BH"),
    significant = p_BH < 0.05
  ) %>%
  ungroup() %>%
  mutate(
    effect_size_cat = case_when(
      abs(cliffs_delta) >= 0.474 ~ "large",
      abs(cliffs_delta) >= 0.330 ~ "medium",
      abs(cliffs_delta) >= 0.147 ~ "small",
      TRUE                        ~ "negligible"
    )
  ) %>%
  select(cohort, gene, n_mut, n_wt,
         median_mut, median_wt, median_diff, hl_estimate,
         cliffs_delta, effect_size_cat,
         p_wilcox, p_BH, significant,        # ★ v2: p_wilcox（生p値）を明示
         zero_rate_mut, zero_rate_wt,        # ★ v2
         zero_rate_all, tied_warning,        # ★ v2
         note) %>%
  arrange(cohort, gene)

# =============================================================================
# 8. 結果保存（★ v2: ファイル名を step16a_v2_* に変更）
# =============================================================================

write_csv(results, file.path(OUT_DIR, "step16a_v2_results.csv"))
log_msg("保存: step16a_v2_results.csv")

# =============================================================================
# 9. ★ v2: 全遺伝子の p_raw / p_BH をログに記録（診断用）
# =============================================================================

log_msg("--- 全遺伝子 p_raw / p_BH 詳細ログ（v2追加） ---")

for (nm in names(cohorts)) {
  res_co <- results %>% filter(cohort == nm, note == "ok")
  n_tied <- sum(res_co$tied_warning, na.rm = TRUE)
  log_msg(sprintf("  [%s]  (tied_warning遺伝子数: %d/7)", nm, n_tied))
  
  for (i in seq_len(nrow(res_co))) {
    r       <- res_co[i, ]
    sig_mark <- if (!is.na(r$significant) && r$significant) "★" else "  "
    tied_tag <- if (!is.na(r$tied_warning) && r$tied_warning) " [TIED]" else ""
    log_msg(sprintf("    %s %-12s p_raw=%s  p_BH=%s  δ=%+.3f  zero_all=%.1f%%%s",
                    sig_mark,
                    r$gene,
                    ifelse(is.na(r$p_wilcox), "   NA  ",
                           formatC(r$p_wilcox, format="f", digits=6)),
                    ifelse(is.na(r$p_BH), "   NA  ",
                           formatC(r$p_BH, format="f", digits=6)),
                    ifelse(is.na(r$cliffs_delta), 0, r$cliffs_delta),
                    ifelse(is.na(r$zero_rate_all), 0, r$zero_rate_all * 100),
                    tied_tag))
  }
}

# =============================================================================
# 10. ★ v2: GLASS専用 tied 診断サマリー
# =============================================================================

log_msg("--- GLASS_WXS: Tied診断サマリー（v2追加） ---")
glass_res <- results %>% filter(cohort == "GLASS_WXS", note == "ok")

n_tied <- sum(glass_res$tied_warning, na.rm = TRUE)
log_msg(sprintf("  tied_warning遺伝子数: %d / %d", n_tied, nrow(glass_res)))
log_msg("  （p_raw > 0.999 またはゼロ率 > 70%% の遺伝子）")

if (n_tied > 0) {
  tied_genes <- glass_res %>% filter(tied_warning) %>% pull(gene)
  log_msg(sprintf("  該当遺伝子: %s", paste(tied_genes, collapse=", ")))
  log_msg("  解釈: GLASSコホート（n=79）での低発現遺伝子は検出力が不足しており、")
  log_msg("        p値が1.0（=完全tied）になるのはバグではなく統計的限界です。")
  log_msg("        これらの遺伝子はTableS2でp_BH=1.00と記載し、")
  log_msg("        脚注に tied ranks の注釈を追加することを推奨します。")
}

# =============================================================================
# 11. サマリーテキスト出力（★ v2: tied情報を追記）
# =============================================================================

summary_lines <- c(
  "=== Step 16a v2: チェックポイント特異性解析 サマリー ===",
  sprintf("実行日時: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("解析遺伝子: %s", paste(TARGET_GENES, collapse=", ")),
  "v2変更点: tied診断追加・全遺伝子p_raw/p_BHを記録・ファイル名をv2に変更",
  ""
)

for (nm in names(cohorts)) {
  res_co <- results %>% filter(cohort == nm, note == "ok")
  n_tied <- sum(res_co$tied_warning, na.rm = TRUE)
  summary_lines <- c(summary_lines,
                     sprintf("【%s】（tied_warning: %d遺伝子）", nm, n_tied),
                     sprintf("  n_mut=%d, n_wt=%d", res_co$n_mut[1], res_co$n_wt[1])
  )
  
  for (i in seq_len(nrow(res_co))) {
    r <- res_co[i, ]
    sig_mark  <- if (!is.na(r$significant) && r$significant) "★" else "  "
    tied_tag  <- if (!is.na(r$tied_warning) && r$tied_warning) " [TIED]" else ""
    p_raw_str <- ifelse(is.na(r$p_wilcox), "NA",
                        formatC(r$p_wilcox, format = "e", digits = 3))
    p_bh_str  <- ifelse(is.na(r$p_BH), "NA",
                        formatC(r$p_BH, format = "f", digits = 4))
    summary_lines <- c(summary_lines,
                       sprintf("  %s %-12s p_raw=%-12s p_BH=%-8s δ=%+.3f (%s)  zero=%.0f%%%s",
                               sig_mark, r$gene,
                               p_raw_str, p_bh_str,
                               ifelse(is.na(r$cliffs_delta), 0, r$cliffs_delta),
                               r$effect_size_cat,
                               ifelse(is.na(r$zero_rate_all), 0,
                                      r$zero_rate_all * 100),
                               tied_tag)
    )
  }
  summary_lines <- c(summary_lines, "")
}

# LAG3特異性チェック
lag3_sig <- results %>%
  filter(note == "ok") %>%
  group_by(cohort) %>%
  summarise(
    lag3_sig    = any(gene == "LAG3" & significant, na.rm = TRUE),
    other_sig   = any(gene != "LAG3" & significant, na.rm = TRUE),
    n_sig_genes = sum(significant, na.rm = TRUE),
    sig_genes   = paste(gene[significant & !is.na(significant)], collapse=", "),
    n_tied      = sum(tied_warning, na.rm = TRUE),
    tied_genes  = paste(gene[tied_warning & !is.na(tied_warning)], collapse=", "),
    .groups = "drop"
  )

summary_lines <- c(summary_lines, "--- LAG3特異性チェック ---")
for (i in seq_len(nrow(lag3_sig))) {
  r <- lag3_sig[i, ]
  pattern <- dplyr::case_when(
    r$lag3_sig & !r$other_sig ~ "✅ LAG3のみ有意（特異性あり）",
    r$lag3_sig &  r$other_sig ~ "⚠️  LAG3と他も有意（特異性一部）",
    !r$lag3_sig               ~ "❌ LAG3非有意"
  )
  summary_lines <- c(summary_lines,
                     sprintf("  %-25s %s  有意遺伝子: [%s]", r$cohort, pattern, r$sig_genes))
  if (r$n_tied > 0) {
    summary_lines <- c(summary_lines,
                       sprintf("  %-25s ⚠️  Tied genes (%d個): [%s]",
                               "", r$n_tied, r$tied_genes),
                       sprintf("  %-25s    これらのp値はtied ranksにより1.0（検出力不足）",
                               ""))
  }
}

# ★ v2: TableS2 注釈推奨
summary_lines <- c(summary_lines,
                   "",
                   "--- TableS2 注釈推奨（v2追加） ---",
                   "  GLASS_WXS の tied 遺伝子については TableS2 の脚注に以下を追記してください:",
                   "  \"For the GLASS cohort, several checkpoint genes with low expression",
                   "   (high proportion of zero TPM values) produced p = 1.00 owing to",
                   "   complete tied ranks in the Wilcoxon test, indicating insufficient",
                   "   power to detect differences in these genes rather than a code error.\""
)

writeLines(summary_lines, file.path(OUT_DIR, "step16a_v2_summary.txt"))
log_msg("保存: step16a_v2_summary.txt")

cat("\n")
cat(paste(summary_lines, collapse="\n"))
cat("\n")

# =============================================================================
# 12. ログへの主要結果記録
# =============================================================================

log_msg("--- 主要結果（LAG3 特異性確認） ---")
for (i in seq_len(nrow(lag3_sig))) {
  r <- lag3_sig[i, ]
  log_msg(sprintf("  %s: 有意遺伝子=%d個 [%s]  tied=%d個 [%s]",
                  r$cohort, r$n_sig_genes, r$sig_genes, r$n_tied, r$tied_genes))
}

log_msg("=== Step 16a v2: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))
close(log_con)

cat("\n============================\n")
cat("Step 16a v2 完了\n")
cat("出力ファイル（v2・既存ファイルを上書きしません）:\n")
cat(sprintf("  %s/step16a_v2_results.csv\n",  OUT_DIR))
cat(sprintf("  %s/step16a_v2_summary.txt\n",  OUT_DIR))
cat(sprintf("  %s/step16a_v2_log.txt\n",      OUT_DIR))
cat("============================\n")
