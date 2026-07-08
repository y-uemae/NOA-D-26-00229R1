# =============================================================================
# step26_permutation.R
# GBM/Glioma TP53xLAG3 解析 - Step 26: 置換検定（TP53×source交互作用）
#
# 目的:
#   Step24の交互作用p値（M_int: p=0.059、M_int_immune: p=0.017）の安定性を
#   置換検定で確認する。
#   "少なくとも方向は安定"を示し、交互作用検定の頑健性を補強する。
#
# 手順:
#   1. 実データでM_int/M_int_immuneを実行 → 交互作用項のt統計量（観測値）
#   2. tp53_binをシャッフルしてN_PERM回再実行 → 置換分布
#   3. 置換p値 = |t_perm| >= |t_obs| の割合
#
# モデル:
#   M_int     : LAG3 ~ tp53_bin * source_bin + idh_bin
#   M_int_imm : LAG3 ~ tp53_bin * source_bin + idh_bin + tcell_score + apm_score
#
# 出力:
#   figS26_permutation.pdf/png   : 2モデルの置換分布（観測値赤線）
#   step26_permutation_results.csv
#   step26_log.txt
#
# 出力先: 26_permutation/
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "26_permutation")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

N_PERM <- 20000
set.seed(42)

# 免疫スコア定義（Step17/20/21/24と完全に同一）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

# 色設定
COL_MUT  <- "#E64B35"
COL_HIST <- "#AAAAAA"

# 入力ファイル
GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step26_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 26: 置換検定（TP53×source交互作用）開始 ===")
log_msg(sprintf("N_PERM=%d, seed=42", N_PERM))

# =============================================================================
# 2. データ読み込み・結合（Step24と同一ロジック）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
wide     <- read_csv(WIDE_PATH, show_col_types = FALSE)

log_msg(sprintf("GDC: %d行 / wide: %d行x%d列", nrow(gdc_base), nrow(wide), ncol(wide)))

# 列名正規化ヘルパー
get_log2_col <- function(gene, df) {
  c1 <- paste0(gene, "_log2tpm")
  c2 <- paste0(gsub("-", ".", gene), "_log2tpm")
  if (c1 %in% names(df)) return(c1)
  if (c2 %in% names(df)) return(c2)
  return(NA_character_)
}

calc_score <- function(df, genes, score_name) {
  cols  <- sapply(genes, get_log2_col, df = df)
  found <- cols[!is.na(cols)]
  missing <- genes[is.na(cols)]
  if (length(missing) > 0)
    log_msg(sprintf("  %s: 不足遺伝子 %s", score_name, paste(missing, collapse = ",")))
  log_msg(sprintf("  %s: %d/%d遺伝子で計算", score_name, length(found), length(genes)))
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

# GDC結合
all_score_genes  <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols_found <- intersect(
  c(paste0(all_score_genes, "_log2tpm"),
    paste0(gsub("-", ".", all_score_genes), "_log2tpm")),
  names(wide)
)

wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")), any_of(score_cols_found))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_score %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_score %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)

for (col in score_cols_found) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]])
    gdc[[wcol]] <- NULL
  }
}

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "IFN-gamma")
  )

# Grade4解析データ（2モデル用に別々に用意）
gdc_g4_base <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant")
  ) %>%
  select(LAG3_log2tpm, tp53_bin, source_bin, idh_bin,
         tcell_score, apm_score) %>%
  na.omit()

log_msg(sprintf("GDC Grade4: n=%d", nrow(gdc_g4_base)))

# =============================================================================
# 3. 置換検定ヘルパー
# =============================================================================

# 交互作用項のt統計量を取得
get_t_stat <- function(df, formula_str) {
  fit  <- lm(as.formula(formula_str), data = df)
  smr  <- coef(summary(fit))
  terms <- rownames(smr)
  int_term <- grep("tp53_bin:source_bin|source_bin:tp53_bin", terms, value = TRUE)
  if (length(int_term) == 0) return(NA_real_)
  smr[int_term[1], "t value"]
}

# 置換検定本体
run_permutation <- function(df, formula_str, model_name, n_perm) {
  log_msg(sprintf("  [%s] 観測t統計量を計算...", model_name))
  t_obs <- get_t_stat(df, formula_str)
  log_msg(sprintf("  [%s] t_obs=%.4f", model_name, t_obs))
  
  log_msg(sprintf("  [%s] 置換開始 (n=%d)...", model_name, n_perm))
  
  t_perm <- replicate(n_perm, {
    df_perm <- df
    df_perm$tp53_bin <- sample(df$tp53_bin)
    get_t_stat(df_perm, formula_str)
  })
  t_perm <- t_perm[is.finite(t_perm)]
  
  # 両側置換p値
  p_perm <- mean(abs(t_perm) >= abs(t_obs), na.rm = TRUE)
  
  log_msg(sprintf("  [%s] 置換p値=%.4f (n_perm=%d)",
                  model_name, p_perm, n_perm))
  log_msg(sprintf("  [%s] 置換分布: mean=%.4f, SD=%.4f, |t_obs| percentile=%.1f%%",
                  model_name,
                  mean(t_perm, na.rm = TRUE),
                  sd(t_perm, na.rm = TRUE),
                  mean(abs(t_perm) < abs(t_obs), na.rm = TRUE) * 100))
  
  list(
    model_name = model_name,
    formula    = formula_str,
    t_obs      = t_obs,
    p_perm     = p_perm,
    t_perm     = t_perm
  )
}

# =============================================================================
# 4. 置換検定実行
# =============================================================================

log_msg("--- 置換検定実行 ---")

# M_int
res_int <- run_permutation(
  df           = gdc_g4_base,
  formula_str  = "LAG3_log2tpm ~ tp53_bin * source_bin + idh_bin",
  model_name   = "M_int",
  n_perm       = N_PERM
)

# M_int_immune
res_int_imm <- run_permutation(
  df           = gdc_g4_base,
  formula_str  = "LAG3_log2tpm ~ tp53_bin * source_bin + idh_bin + tcell_score + apm_score",
  model_name   = "M_int_immune",
  n_perm       = N_PERM
)

# 結果をCSVに保存
perm_summary <- tibble(
  model       = c(res_int$model_name,     res_int_imm$model_name),
  formula     = c(res_int$formula,        res_int_imm$formula),
  t_obs       = round(c(res_int$t_obs,    res_int_imm$t_obs),   4),
  p_perm      = round(c(res_int$p_perm,   res_int_imm$p_perm),  4),
  n_perm      = N_PERM,
  sig_label   = case_when(
    p_perm < 0.001 ~ "***",
    p_perm < 0.01  ~ "**",
    p_perm < 0.05  ~ "*",
    TRUE           ~ "ns"
  )
)

write_csv(perm_summary, file.path(OUT_DIR, "step26_permutation_results.csv"))
log_msg("保存: step26_permutation_results.csv")

# =============================================================================
# 5. figS26: 置換分布ヒストグラム（2モデル並列）
# =============================================================================

log_msg("--- figS26 作成 ---")

make_perm_plot <- function(res, panel_label) {
  
  t_df   <- data.frame(t_val = res$t_perm)
  t_obs  <- res$t_obs
  p_perm <- res$p_perm
  sig    <- ifelse(p_perm < 0.001, "***",
                   ifelse(p_perm < 0.01,  "**",
                          ifelse(p_perm < 0.05,  "*", "ns")))
  
  # x軸範囲（対称）
  x_lim <- max(abs(c(t_df$t_val, t_obs)), na.rm = TRUE) * 1.15
  
  ggplot(t_df, aes(x = t_val)) +
    geom_histogram(bins = 60, fill = COL_HIST, color = "white",
                   linewidth = 0.2, alpha = 0.85) +
    # 観測t統計量（赤線）
    geom_vline(xintercept = t_obs,  color = COL_MUT, linewidth = 1.0,
               linetype = "solid") +
    geom_vline(xintercept = -t_obs, color = COL_MUT, linewidth = 0.6,
               linetype = "dashed") +
    # p値アノテーション
    annotate("text",
             x = t_obs * 0.85, y = Inf,
             label = sprintf("t_obs=%.3f\np_perm=%.4f (%s)", t_obs, p_perm, sig),
             hjust = ifelse(t_obs > 0, 1.05, -0.05),
             vjust = 1.4, size = 3.2, color = COL_MUT, fontface = "bold") +
    scale_x_continuous(
      name = "Permuted t-statistic (tp53:source interaction)"
    ) +
    coord_cartesian(xlim = c(-x_lim, x_lim)) +
    labs(
      y        = "Count",
      title    = sprintf("%s  %s", panel_label, res$model_name),
      subtitle = sprintf(
        "N_perm=%d. Red solid: observed t=%.3f. Red dashed: -|t_obs|.",
        N_PERM, t_obs
      )
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 8.5, color = "grey40")
    )
}

fig_top <- make_perm_plot(res_int,     "A")
fig_bot <- make_perm_plot(res_int_imm, "B")

fig_s26 <- wrap_plots(fig_top) / wrap_plots(fig_bot) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title   = "Permutation test: TP53 x source interaction",
    caption = paste0(
      "GDC Grade4 (n=", nrow(gdc_g4_base), "). ",
      "tp53_bin permuted ", N_PERM, " times (seed=42). ",
      "Two-sided permutation p: proportion of |t_perm| >= |t_obs|. ",
      "A: LAG3 ~ tp53*source + IDH. ",
      "B: LAG3 ~ tp53*source + IDH + T-cell + APM."
    ),
    theme = theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

pdf(file.path(OUT_DIR, "figS26_permutation.pdf"), width = 9, height = 8)
print(fig_s26)
dev.off()
log_msg("保存: figS26_permutation.pdf")

agg_png(file.path(OUT_DIR, "figS26_permutation_450dpi.png"),
        width = 9, height = 8, units = "in", res = 450)
print(fig_s26)
dev.off()
log_msg("保存: figS26_permutation_450dpi.png")

# =============================================================================
# 6. 完了
# =============================================================================

log_msg("=== Step 26: 完了 ===")
for (i in seq_len(nrow(perm_summary))) {
  r <- perm_summary[i, ]
  log_msg(sprintf("  [%s] t_obs=%.4f, p_perm=%.4f (%s)",
                  r$model, r$t_obs, r$p_perm, r$sig_label))
}
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 26 完了\n")
for (i in seq_len(nrow(perm_summary))) {
  r <- perm_summary[i, ]
  cat(sprintf("  [%s] p_perm=%.4f (%s)\n", r$model, r$p_perm, r$sig_label))
}
cat("出力ファイル:\n")
cat(sprintf("  %s/figS26_permutation.pdf\n",         OUT_DIR))
cat(sprintf("  %s/figS26_permutation_450dpi.png\n",  OUT_DIR))
cat(sprintf("  %s/step26_permutation_results.csv\n", OUT_DIR))
cat(sprintf("  %s/step26_log.txt\n",                 OUT_DIR))
cat("============================\n")
