# =============================================================================
# step20_residual_analysis.R
# GBM/Glioma TP53×LAG3 解析 - Step 20: 残差アプローチ（Checkpoint特異性）
#
# 目的:
#   7チェックポイント遺伝子それぞれを免疫スコアで残差化し、
#   「免疫量で説明されない成分」に対するTP53効果を示す。
#   LAG3だけが残差においてもTP53と関連することを可視化。
#
# 解析手順:
#   1. 各遺伝子を残差化
#      gene_resid = residuals(lm(gene ~ tcell_score + apm_score + ifng_score
#                                     + source_bin + idh_bin))
#   2. gene_resid ~ tp53_bin を7遺伝子で実行（BH補正）
#   3. Forest風プロット（β + 95%CI、LAG3だけ有意を一目で）
#
# スコア定義（Step17と同一）:
#   T-cell score : CD3D, CD3E, CD3G, CD8A, CD8B, GZMA, GZMB, PRF1  (8遺伝子)
#   APM score    : B2M, TAP1, TAP2, TAPBP, HLA-A, HLA-B, HLA-C, NLRC5  (8遺伝子)
#   IFNγ score   : STAT1, IRF1, IRF9, CXCL9, CXCL10, CXCL11,
#                  GBP1, GBP2, GBP4, GBP5, IDO1  (11遺伝子)
#
# 遺伝子順序（Step16/Fig2と同一・上から）:
#   LAG3 > PDCD1 > CD274 > PDCD1LG2 > CTLA4 > TIGIT > HAVCR2
#
# 出力先: 20_residual/
#   step20_residual_results.csv   : 7遺伝子の回帰結果（BH補正済み）
#   figS20_residual_forest.pdf/png: Forest風プロット
#   step20_log.txt
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "20_residual")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# スコア定義（Step17と完全に同一）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

# 解析対象遺伝子（Step16/Fig2と同一順序）
CP_GENES <- c("LAG3", "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "TIGIT", "HAVCR2")

# ggplot用factor順序（上→下 = LAG3が最上部）
CP_FACTOR_LEVELS <- rev(CP_GENES)

# 色設定（引継書10章に準拠）
COL_LAG3 <- "#E64B35"   # LAG3強調色
COL_OTHER <- "#AAAAAA"  # 他遺伝子

# 入力ファイル
GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定（Step17と同様にcat()方式で安全に）
# =============================================================================

log_file <- file.path(OUT_DIR, "step20_log.txt")
cat("", file = log_file, append = FALSE)  # 初期化

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 20: 残差アプローチ（Checkpoint特異性）開始 ===")
log_msg(sprintf("解析対象遺伝子: %s", paste(CP_GENES, collapse = ", ")))
log_msg(sprintf("T-cell genes (%d): %s", length(TCELL_GENES), paste(TCELL_GENES, collapse = ", ")))
log_msg(sprintf("APM genes    (%d): %s", length(APM_GENES),   paste(APM_GENES,   collapse = ", ")))
log_msg(sprintf("IFNg genes   (%d): %s", length(IFNG_GENES),  paste(IFNG_GENES,  collapse = ", ")))

# =============================================================================
# 2. データ読み込み・結合
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
wide     <- read_csv(WIDE_PATH, show_col_types = FALSE)

log_msg(sprintf("GDC: %d行 / wide: %d行×%d列", nrow(gdc_base), nrow(wide), ncol(wide)))

# =============================================================================
# 3. 列名正規化ヘルパー（HLA-A等のハイフン対応、Step17と同一）
# =============================================================================

get_log2_col <- function(gene, df) {
  c1 <- paste0(gene, "_log2tpm")
  c2 <- paste0(gsub("-", ".", gene), "_log2tpm")
  if (c1 %in% names(df)) return(c1)
  if (c2 %in% names(df)) return(c2)
  return(NA_character_)
}

calc_score <- function(df, genes, score_name, log_fn) {
  cols    <- sapply(genes, get_log2_col, df = df)
  found   <- cols[!is.na(cols)]
  missing <- genes[is.na(cols)]
  if (length(missing) > 0)
    log_fn(sprintf("  %s: 不足遺伝子 %s", score_name, paste(missing, collapse = ",")))
  log_fn(sprintf("  %s: %d/%d遺伝子で計算", score_name, length(found), length(genes)))
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

# =============================================================================
# 4. GDC結合・スコア計算（Step17と同一ロジック）
# =============================================================================

log_msg("--- GDC結合 ---")

all_score_genes  <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols_found <- intersect(
  c(paste0(all_score_genes, "_log2tpm"),
    paste0(gsub("-", ".", all_score_genes), "_log2tpm")),
  names(wide)
)

wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")),
         any_of(score_cols_found))

# チェックポイント遺伝子列も抽出
cp_cols_found <- intersect(
  c(paste0(CP_GENES, "_log2tpm"),
    paste0(gsub("-", ".", CP_GENES), "_log2tpm")),
  names(wide)
)
log_msg(sprintf("CP遺伝子列（wide内に存在）: %d / %d", length(cp_cols_found), length(CP_GENES)))

wide_all <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")),
         any_of(c(score_cols_found, cp_cols_found)))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_all %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_all %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)

# .w列マージ
all_extra_cols <- c(score_cols_found, cp_cols_found)
for (col in all_extra_cols) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]])
    gdc[[wcol]] <- NULL
  }
}
log_msg(sprintf("GDC結合後: %d行", nrow(gdc)))

# =============================================================================
# 5. 免疫スコア計算
# =============================================================================

log_msg("--- 免疫スコア計算 ---")

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell", log_msg),
    apm_score   = calc_score(., APM_GENES,   "APM",    log_msg),
    ifng_score  = calc_score(., IFNG_GENES,  "IFNg",   log_msg)
  )

# Grade4解析データ
gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant")
  )

log_msg(sprintf("GDC Grade4 解析対象: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$tp53_bin == 1),
                sum(gdc_g4$tp53_bin == 0)))

# スコアNA確認
for (sc in c("tcell_score", "apm_score", "ifng_score")) {
  log_msg(sprintf("  %s: NA=%d, median=%.3f",
                  sc, sum(is.na(gdc_g4[[sc]])),
                  median(gdc_g4[[sc]], na.rm = TRUE)))
}

# =============================================================================
# 6. 残差化 + TP53回帰（7遺伝子）
# =============================================================================

log_msg("--- 残差化 + TP53回帰 ---")

results_list <- lapply(CP_GENES, function(gene) {
  
  # log2tpm列名の解決
  gene_col <- get_log2_col(gene, gdc_g4)
  if (is.na(gene_col)) {
    log_msg(sprintf("  [%s] 列が見つかりません。スキップ。", gene))
    return(NULL)
  }
  
  # 残差化モデル: gene ~ tcell + apm + ifng + source + idh
  df_resid <- gdc_g4 %>%
    select(all_of(c(gene_col, "tcell_score", "apm_score", "ifng_score",
                    "source_bin", "idh_bin", "tp53_bin"))) %>%
    rename(gene_expr = all_of(gene_col)) %>%
    na.omit()
  
  # ステップ1: 残差化
  fit_resid <- lm(gene_expr ~ tcell_score + apm_score + ifng_score +
                    source_bin + idh_bin,
                  data = df_resid)
  df_resid$gene_resid <- residuals(fit_resid)
  
  # ステップ2: 残差 ~ tp53_bin
  fit_tp53 <- lm(gene_resid ~ tp53_bin, data = df_resid)
  ci  <- as.double(confint(fit_tp53, "tp53_bin"))
  smr <- summary(fit_tp53)
  coef_tbl <- coef(smr)
  
  beta  <- coef_tbl["tp53_bin", "Estimate"]
  se    <- coef_tbl["tp53_bin", "Std. Error"]
  pval  <- coef_tbl["tp53_bin", "Pr(>|t|)"]
  
  log_msg(sprintf("  [%s] n=%d, beta=%+.4f [%+.4f, %+.4f], p=%.4f",
                  gene, nrow(df_resid), beta, ci[1], ci[2], pval))
  
  tibble(
    gene     = gene,
    n        = nrow(df_resid),
    beta     = round(beta,  4),
    ci_lower = round(ci[1], 4),
    ci_upper = round(ci[2], 4),
    se       = round(se,    4),
    p_raw    = pval
  )
})

results <- do.call(rbind, Filter(Negate(is.null), results_list))

# BH補正
results <- results %>%
  mutate(
    p_BH      = p.adjust(p_raw, method = "BH"),
    sig_BH    = p_BH < 0.05,
    gene      = factor(gene, levels = CP_FACTOR_LEVELS),
    is_lag3   = (as.character(gene) == "LAG3")
  )

log_msg("--- BH補正後 結果 ---")
for (i in order(results$gene)) {
  r <- results[i, ]
  log_msg(sprintf("  [%s] beta=%+.4f [%+.4f, %+.4f] p_raw=%.4f p_BH=%.4f %s",
                  as.character(r$gene), r$beta, r$ci_lower, r$ci_upper,
                  r$p_raw, r$p_BH,
                  ifelse(r$sig_BH, "★", "")))
}

# 保存
write_csv(results %>% mutate(gene = as.character(gene)),
          file.path(OUT_DIR, "step20_residual_results.csv"))
log_msg("保存: step20_residual_results.csv")

# =============================================================================
# 7. figS20: Forest風プロット
# =============================================================================

log_msg("--- figS20 作成 ---")

# p値ラベル作成
results <- results %>%
  mutate(
    p_label = case_when(
      p_BH < 0.001 ~ sprintf("p[BH]=%.4f ***", p_BH),
      p_BH < 0.01  ~ sprintf("p[BH]=%.4f **",  p_BH),
      p_BH < 0.05  ~ sprintf("p[BH]=%.4f *",   p_BH),
      TRUE          ~ sprintf("p[BH]=%.3f",     p_BH)
    ),
    # 表示用遺伝子ラベル
    gene_label = case_when(
      as.character(gene) == "PDCD1"    ~ "PDCD1 (PD-1)",
      as.character(gene) == "CD274"    ~ "CD274 (PD-L1)",
      as.character(gene) == "PDCD1LG2" ~ "PDCD1LG2 (PD-L2)",
      as.character(gene) == "HAVCR2"   ~ "HAVCR2 (TIM-3)",
      TRUE ~ as.character(gene)
    ),
    gene_label = factor(gene_label,
                        levels = rev(c(
                          "LAG3",
                          "PDCD1 (PD-1)",
                          "CD274 (PD-L1)",
                          "PDCD1LG2 (PD-L2)",
                          "CTLA4",
                          "TIGIT",
                          "HAVCR2 (TIM-3)"
                        )))
  )

# x軸範囲
x_min <- min(results$ci_lower, na.rm = TRUE) * 1.2
x_max <- max(results$ci_upper, na.rm = TRUE) * 1.3

fig_s20 <- ggplot(results,
                  aes(x = beta, y = gene_label,
                      color = is_lag3, fill = is_lag3)) +
  # 参照線
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  # CI
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.3, linewidth = 0.8) +
  # 点
  geom_point(shape = 21, size = 4, stroke = 1.1) +
  # p値ラベル（CI右端の外側）
  geom_text(aes(x = ci_upper, label = p_label),
            hjust = -0.1, size = 3.2, color = "grey30") +
  # 色設定
  scale_color_manual(
    values = c("TRUE" = COL_LAG3, "FALSE" = "#555555"),
    guide  = "none"
  ) +
  scale_fill_manual(
    values = c("TRUE" = COL_LAG3, "FALSE" = "#BBBBBB"),
    guide  = "none"
  ) +
  # 軸
  scale_x_continuous(
    name   = expression(paste(beta[TP53], " in residual model (95% CI)")),
    limits = c(x_min, x_max),
    breaks = seq(-0.2, 0.5, by = 0.1)
  ) +
  labs(
    y        = NULL,
    title    = "TP53 effect on immune-adjusted residuals: checkpoint gene specificity",
    subtitle = paste0(
      "GDC Grade4 (n=442). Residuals from: gene ~ T-cell + APM + IFN\u03b3 + source + IDH.\n",
      "Then: residual ~ tp53_bin. BH correction across 7 genes."
    ),
    caption  = paste0(
      "Red: LAG3 (highlighted). \u2605 p[BH] < 0.05. ",
      "T-cell: CD3D/E/G, CD8A/B, GZMA/B, PRF1. ",
      "APM: B2M, TAP1/2, TAPBP, HLA-A/B/C, NLRC5. ",
      "IFN\u03b3: STAT1, IRF1/9, CXCL9/10/11, GBP1/2/4/5, IDO1."
    )
  ) +
  coord_cartesian(xlim = c(x_min, x_max)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    plot.caption       = element_text(size = 7.5, color = "grey50", hjust = 0),
    plot.margin        = margin(10, 80, 10, 10)
  )

# PDF出力
pdf_path <- file.path(OUT_DIR, "figS20_residual_forest.pdf")
pdf(pdf_path, width = 9, height = 5.5)
print(fig_s20)
dev.off()
log_msg("保存: figS20_residual_forest.pdf")

# PNG出力
png_path <- file.path(OUT_DIR, "figS20_residual_forest_450dpi.png")
agg_png(png_path, width = 9, height = 5.5, units = "in", res = 450)
print(fig_s20)
dev.off()
log_msg("保存: figS20_residual_forest_450dpi.png")

# =============================================================================
# 8. 完了
# =============================================================================

log_msg("=== Step 20: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 20 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/figS20_residual_forest.pdf\n",        OUT_DIR))
cat(sprintf("  %s/figS20_residual_forest_450dpi.png\n", OUT_DIR))
cat(sprintf("  %s/step20_residual_results.csv\n",       OUT_DIR))
cat(sprintf("  %s/step20_log.txt\n",                    OUT_DIR))
cat("============================\n")
