# =============================================================================
# step24_interaction_source.R
# GBM/Glioma TP53xLAG3 解析 - Step 24: TP53×source交互作用検定
#
# 目的:
#   source（TCGA vs CPTAC/HCMI）でTP53効果が統計的に異なるかを正式に検定する。
#   p非有意 → "差は統計的に明確でない"として査読対応
#   p有意   → "source差あり、ただし方向は一致"と記載
#
# モデル:
#   M_base    : LAG3 ~ tp53_bin + source_bin + idh_bin          （基準）
#   M_int     : LAG3 ~ tp53_bin * source_bin + idh_bin          （交互作用）
#   M_int_imm : LAG3 ~ tp53_bin * source_bin + idh_bin
#                    + tcell_score + apm_score                   （免疫調整後）
#
# 出力:
#   figS24a_interaction_forest.pdf/png  : source別β+CI（M_int）
#   figS24b_emm.pdf/png                 : EMM（source×TP53の組み合わせ）
#   figS24_combined.pdf/png             : A+B上下配置
#   step24_interaction_results.csv      : 交互作用項の検定結果
#   step24_emm_results.csv              : EMM結果
#   step24_log.txt
#
# 出力先: 24_interaction_source/
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
OUT_DIR    <- file.path(RESULT_DIR, "24_interaction_source")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 免疫スコア定義（Step17/20/21と完全に同一）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

# 色設定（引継書に準拠）
COL_MUT   <- "#E64B35"
COL_WT    <- "#AAAAAA"
COL_TCGA  <- "#3C5488"
COL_CPTAC <- "#E07B54"
COL_ALL   <- "#8491B4"

# 入力ファイル
GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step24_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 24: TP53 x source 交互作用検定 開始 ===")

# =============================================================================
# 2. データ読み込み・結合（Step20/21と同一ロジック）
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
log_msg(sprintf("GDC結合後: %d行", nrow(gdc)))

# =============================================================================
# 3. 免疫スコア計算
# =============================================================================

log_msg("--- 免疫スコア計算 ---")

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "IFN-gamma")
  )

# Grade4解析データ
gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = factor(
      ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT"),
      levels = c("TP53 WT", "TP53 Mut")
    ),
    source_label = factor(source,
                          levels = c("TCGA", "CPTAC_HCMI"),
                          labels = c("TCGA", "CPTAC/HCMI"))
  )

log_msg(sprintf("GDC Grade4: n=%d (Mut=%d, WT=%d | TCGA=%d, CPTAC/HCMI=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$tp53_bin == 1), sum(gdc_g4$tp53_bin == 0),
                sum(gdc_g4$source == "TCGA"), sum(gdc_g4$source == "CPTAC_HCMI")))

# =============================================================================
# 4. 交互作用モデル
# =============================================================================

log_msg("--- 交互作用モデル ---")

# M_base: 基準モデル（Step09b再現）
fit_base <- lm(LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin,
               data = gdc_g4)

# M_int: 交互作用項追加
fit_int  <- lm(LAG3_log2tpm ~ tp53_bin * source_bin + idh_bin,
               data = gdc_g4)

# M_int_imm: 免疫調整後の交互作用
fit_int_imm <- lm(LAG3_log2tpm ~ tp53_bin * source_bin + idh_bin
                  + tcell_score + apm_score,
                  data = gdc_g4 %>% filter(!is.na(tcell_score), !is.na(apm_score)))

# モデル比較（ANOVA）
anova_result <- anova(fit_base, fit_int)
log_msg("--- ANOVA: M_base vs M_int ---")
log_msg(sprintf("  F=%.4f, df=%d, p=%.4f",
                anova_result$F[2],
                anova_result$Df[2],
                anova_result$`Pr(>F)`[2]))

# 交互作用項の係数
extract_int_coef <- function(fit, model_name) {
  smr   <- coef(summary(fit))
  terms <- rownames(smr)
  
  # 交互作用項の行を取得
  int_term <- grep("tp53_bin:source_bin|source_bin:tp53_bin", terms, value = TRUE)
  
  results <- lapply(c("tp53_bin", "source_bin", int_term), function(term) {
    if (!term %in% terms) return(NULL)
    ci <- as.double(confint(fit, term))
    tibble(
      model    = model_name,
      term     = term,
      beta     = round(smr[term, "Estimate"],   4),
      ci_lower = round(ci[1],                   4),
      ci_upper = round(ci[2],                   4),
      se       = round(smr[term, "Std. Error"], 4),
      p_value  = smr[term, "Pr(>|t|)"]
    )
  })
  do.call(rbind, Filter(Negate(is.null), results))
}

coef_int     <- extract_int_coef(fit_int,     "M_int")
coef_int_imm <- extract_int_coef(fit_int_imm, "M_int_immune")

all_coefs <- do.call(rbind, list(coef_int, coef_int_imm)) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    r_squared = c(
      rep(round(summary(fit_int)$r.squared,     4), nrow(coef_int)),
      rep(round(summary(fit_int_imm)$r.squared, 4), nrow(coef_int_imm))
    )
  )

log_msg("--- 交互作用項の係数 ---")
for (i in seq_len(nrow(all_coefs))) {
  r <- all_coefs[i, ]
  log_msg(sprintf("  [%s] %s: beta=%+.4f [%+.4f, %+.4f] p=%.4f (%s)",
                  r$model, r$term, r$beta, r$ci_lower, r$ci_upper,
                  r$p_value, r$sig_label))
}

write_csv(all_coefs, file.path(OUT_DIR, "step24_interaction_results.csv"))
log_msg("保存: step24_interaction_results.csv")

# =============================================================================
# 5. EMM: source × TP53の組み合わせ別予測値
# =============================================================================

log_msg("--- EMM計算（手動）---")

# 各source内でのTP53効果（EMM差 = Mut - WT）
# M_int: LAG3 ~ tp53_bin * source_bin + idh_bin
# TCGA (source_bin=0): TP53効果 = beta_tp53
# CPTAC (source_bin=1): TP53効果 = beta_tp53 + beta_interaction

b  <- coef(fit_int)
ci_tp53  <- as.double(confint(fit_int, "tp53_bin"))
ci_inter <- as.double(confint(fit_int, "tp53_bin:source_bin"))

# TCGA TP53効果
eff_tcga  <- b["tp53_bin"]
# CPTAC TP53効果
eff_cptac <- b["tp53_bin"] + b["tp53_bin:source_bin"]

# SE（delta法でCIを近似）
vcov_mat  <- vcov(fit_int)
se_tcga   <- sqrt(vcov_mat["tp53_bin", "tp53_bin"])
se_cptac  <- sqrt(vcov_mat["tp53_bin",  "tp53_bin"] +
                    vcov_mat["tp53_bin:source_bin", "tp53_bin:source_bin"] +
                    2 * vcov_mat["tp53_bin", "tp53_bin:source_bin"])

emm_results <- tibble(
  source    = c("TCGA", "CPTAC/HCMI"),
  tp53_effect = round(c(eff_tcga, eff_cptac), 4),
  se          = round(c(se_tcga,  se_cptac),  4),
  ci_lower    = round(c(eff_tcga  - 1.96 * se_tcga,
                        eff_cptac - 1.96 * se_cptac), 4),
  ci_upper    = round(c(eff_tcga  + 1.96 * se_tcga,
                        eff_cptac + 1.96 * se_cptac), 4),
  n           = c(sum(gdc_g4$source == "TCGA"),
                  sum(gdc_g4$source == "CPTAC_HCMI")),
  model       = "M_int"
)

log_msg("--- EMM: source別TP53効果（M_int） ---")
for (i in seq_len(nrow(emm_results))) {
  r <- emm_results[i, ]
  log_msg(sprintf("  [%s] TP53 effect=%.4f [%.4f, %.4f] (n=%d)",
                  r$source, r$tp53_effect, r$ci_lower, r$ci_upper, r$n))
}

write_csv(emm_results, file.path(OUT_DIR, "step24_emm_results.csv"))
log_msg("保存: step24_emm_results.csv")

# =============================================================================
# 6. Panel A: source別TP53効果 Forest風
# =============================================================================

log_msg("--- Panel A: Forest風 作成 ---")

# 全体効果（M_base）も追加
beta_base <- coef(fit_base)["tp53_bin"]
ci_base   <- as.double(confint(fit_base, "tp53_bin"))
se_base   <- sqrt(vcov(fit_base)["tp53_bin", "tp53_bin"])

forest_data <- do.call(rbind, list(
  emm_results %>% select(source, tp53_effect, ci_lower, ci_upper, se, n),
  tibble(
    source      = "All (adjusted)",
    tp53_effect = round(beta_base, 4),
    ci_lower    = round(ci_base[1], 4),
    ci_upper    = round(ci_base[2], 4),
    se          = round(se_base, 4),
    n           = nrow(gdc_g4)
  )
)) %>%
  mutate(
    source = factor(source, levels = c("CPTAC/HCMI", "TCGA", "All (adjusted)")),
    col_group = dplyr::recode(as.character(source),
                              "All (adjusted)" = "all",
                              "TCGA"           = "tcga",
                              "CPTAC/HCMI"     = "cptac"),
    n_label = sprintf("n=%d", n)
  )

col_vals <- c("all" = COL_ALL, "tcga" = COL_TCGA, "cptac" = COL_CPTAC)

# 交互作用p値アノテーション用
int_p     <- all_coefs %>%
  filter(model == "M_int", grepl(":", term)) %>%
  pull(p_value)
int_p_imm <- all_coefs %>%
  filter(model == "M_int_immune", grepl(":", term)) %>%
  pull(p_value)

int_label <- sprintf(
  "Interaction p=%.3f (M_int)\nInteraction p=%.3f (M_int+immune)",
  int_p, int_p_imm
)

fig_a <- ggplot(forest_data,
                aes(x = tp53_effect, y = source,
                    color = col_group, fill = col_group)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.2, linewidth = 0.9) +
  geom_point(shape = 21, size = 5, stroke = 1.1) +
  geom_text(aes(x = ci_upper,
                label = sprintf("b=%+.3f [%.3f, %.3f], n=%d",
                                tp53_effect, ci_lower, ci_upper, n)),
            hjust = -0.08, size = 2.9, color = "grey30") +
  # 交互作用p値
  annotate("text",
           x = min(forest_data$ci_lower) - 0.02,
           y = 0.55,
           label = int_label,
           hjust = 0, vjust = 0, size = 3.0, color = "#333333",
           fontface = "italic") +
  scale_color_manual(values = col_vals, guide = "none") +
  scale_fill_manual(values  = col_vals, guide = "none") +
  scale_x_continuous(
    name   = "TP53 effect on LAG3 (beta, 95% CI)\n[LAG3 ~ tp53 * source + IDH]",
    breaks = seq(-0.1, 0.7, by = 0.1)
  ) +
  coord_cartesian(xlim = c(
    min(forest_data$ci_lower) - 0.15,
    max(forest_data$ci_upper) + 0.35
  )) +
  labs(
    y        = NULL,
    title    = "A  TP53 effect by source: interaction model (M_int)",
    subtitle = paste0(
      "GDC Grade4 (n=442). ",
      "All (adjusted): M_base estimate. ",
      "TCGA/CPTAC: source-specific effects from M_int."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 10, face = "bold"),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, color = "grey40"),
    plot.margin        = margin(10, 130, 10, 10)
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS24a_interaction_forest.pdf"), width = 11, height = 4.5)
print(fig_a)
dev.off()
agg_png(file.path(OUT_DIR, "figS24a_interaction_forest_450dpi.png"),
        width = 11, height = 4.5, units = "in", res = 450)
print(fig_a)
dev.off()
log_msg("保存: figS24a_interaction_forest.pdf/png")

# =============================================================================
# 7. Panel B: EMM（source x TP53の4点）
# =============================================================================

log_msg("--- Panel B: EMM 4点プロット 作成 ---")

# M_intから4点の予測値を手動計算
# 参照: LAG3 = intercept + beta_tp53*tp53 + beta_source*source + beta_idh*idh_mean
#              + beta_interaction*tp53*source
# idh_binの平均値で固定（marginal mean）
idh_mean <- mean(gdc_g4$idh_bin, na.rm = TRUE)
b_int    <- coef(fit_int)

pred_grid <- expand.grid(
  tp53_bin   = c(0, 1),
  source_bin = c(0, 1)
) %>%
  mutate(
    pred = b_int["(Intercept)"] +
      b_int["tp53_bin"]   * tp53_bin +
      b_int["source_bin"] * source_bin +
      b_int["idh_bin"]    * idh_mean +
      b_int["tp53_bin:source_bin"] * tp53_bin * source_bin,
    tp53_label   = factor(ifelse(tp53_bin == 1, "TP53 Mut", "TP53 WT"),
                          levels = c("TP53 WT", "TP53 Mut")),
    source_label = factor(ifelse(source_bin == 0, "TCGA", "CPTAC/HCMI"),
                          levels = c("TCGA", "CPTAC/HCMI"))
  )

# 実データのraw点も重ねる
raw_summary <- gdc_g4 %>%
  group_by(source_label, tp53_label) %>%
  summarise(
    median_val = median(LAG3_log2tpm, na.rm = TRUE),
    q1         = quantile(LAG3_log2tpm, 0.25, na.rm = TRUE),
    q3         = quantile(LAG3_log2tpm, 0.75, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  )

fig_b <- ggplot() +
  # 実データ（IQR棒）
  geom_linerange(data = raw_summary,
                 aes(x = tp53_label, ymin = q1, ymax = q3,
                     color = source_label),
                 linewidth = 4, alpha = 0.25,
                 position = position_dodge(width = 0.5)) +
  # 実データ中央値
  geom_point(data = raw_summary,
             aes(x = tp53_label, y = median_val,
                 color = source_label, shape = source_label),
             size = 4, stroke = 1.2,
             position = position_dodge(width = 0.5)) +
  # モデル予測値（EMM）
  geom_line(data = pred_grid,
            aes(x = tp53_label, y = pred,
                group = source_label, color = source_label,
                linetype = source_label),
            linewidth = 1.0,
            position = position_dodge(width = 0.5)) +
  geom_point(data = pred_grid,
             aes(x = tp53_label, y = pred,
                 color = source_label),
             shape = 18, size = 5,
             position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("TCGA" = COL_TCGA, "CPTAC/HCMI" = COL_CPTAC),
                     name = "Source") +
  scale_shape_manual(values = c("TCGA" = 16, "CPTAC/HCMI" = 17),
                     name = "Source") +
  scale_linetype_manual(values = c("TCGA" = "solid", "CPTAC/HCMI" = "dashed"),
                        name = "Source") +
  labs(
    x        = NULL,
    y        = "LAG3 expression [log2(TPM+1)]",
    title    = "B  Estimated marginal means: TP53 x source",
    subtitle = paste0(
      "Diamonds: model-predicted EMM (IDH at mean). ",
      "Circles/triangles: observed median. ",
      "Bars: IQR."
    )
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8.5, color = "grey40")
  )

# 単体出力
pdf(file.path(OUT_DIR, "figS24b_emm.pdf"), width = 7, height = 5)
print(fig_b)
dev.off()
agg_png(file.path(OUT_DIR, "figS24b_emm_450dpi.png"),
        width = 7, height = 5, units = "in", res = 450)
print(fig_b)
dev.off()
log_msg("保存: figS24b_emm.pdf/png")

# =============================================================================
# 8. Combined: A + B
# =============================================================================

log_msg("--- Combined 上下配置 作成 ---")

# 交互作用検定の結論文
int_conclusion <- ifelse(
  int_p < 0.05,
  sprintf("Interaction significant (p=%.3f): source-specific effect confirmed.", int_p),
  sprintf("Interaction non-significant (p=%.3f): source difference not statistically clear.", int_p)
)
log_msg(sprintf("  結論: %s", int_conclusion))

fig_combined <- wrap_plots(fig_a) / wrap_plots(fig_b) +
  plot_layout(heights = c(1.1, 1.2)) +
  plot_annotation(
    title   = "TP53 x source interaction analysis: LAG3 effect in GDC Grade4",
    caption = paste0(
      "GDC Grade4 (n=442). ",
      "M_int: LAG3 ~ tp53*source + IDH. ",
      "M_int_immune: + T-cell + APM scores. ",
      "EMM computed at mean IDH. ",
      int_conclusion
    ),
    theme = theme(
      plot.title   = element_text(size = 11, face = "bold"),
      plot.caption = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

pdf(file.path(OUT_DIR, "figS24_combined.pdf"), width = 11, height = 10)
print(fig_combined)
dev.off()
agg_png(file.path(OUT_DIR, "figS24_combined_450dpi.png"),
        width = 11, height = 10, units = "in", res = 450)
print(fig_combined)
dev.off()
log_msg("保存: figS24_combined.pdf/png")

# =============================================================================
# 9. 完了サマリー
# =============================================================================

log_msg("=== Step 24: 完了 ===")
log_msg(sprintf("交互作用検定結果: p=%.4f (%s)", int_p,
                ifelse(int_p < 0.05, "有意", "非有意")))
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 24 完了\n")
cat(sprintf("交互作用 p=%.4f (%s)\n", int_p,
            ifelse(int_p < 0.05, "有意", "非有意")))
cat("出力ファイル:\n")
cat(sprintf("  %s/figS24a_interaction_forest.pdf/png\n", OUT_DIR))
cat(sprintf("  %s/figS24b_emm.pdf/png\n",               OUT_DIR))
cat(sprintf("  %s/figS24_combined.pdf/png\n",           OUT_DIR))
cat(sprintf("  %s/step24_interaction_results.csv\n",    OUT_DIR))
cat(sprintf("  %s/step24_emm_results.csv\n",            OUT_DIR))
cat(sprintf("  %s/step24_log.txt\n",                    OUT_DIR))
cat("============================\n")
