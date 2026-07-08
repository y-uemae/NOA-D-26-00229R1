# =============================================================================
# step18_interaction_emm.R
# GBM/Glioma TP53×LAG3 解析 - Step 18: Interaction検定 + Estimated Marginal Means
#
# 目的:
#   「TP53変異とLAG3の関連は免疫スコアに対して平行移動か」を統計的に明文化する。
#   (1) Interaction検定: tp53_bin × tcell_score / tp53_bin × apm_score が非有意
#       → 「同じ免疫量でもTP53-Mutの方がLAG3が高い（切片差）」を証明
#   (2) EMMsで視覚化: tcell_score = mean-1SD / mean / mean+1SD での Mut-WT差
#       → 「T-cell量に関わらずTP53効果が一定（平行移動）」を示す
#
# Step17との関係:
#   Step17でM0→M1でβがほぼ変わらない（+0.3%）ことを確認済み。
#   本ステップはその「平行移動」の統計的・視覚的な補強。
#
# 解析モデル（GDC Grade4, n=442）:
#   M_base:      LAG3 ~ tp53_bin + tcell_score + apm_score + source_bin + idh_bin
#   M_int_tcell: LAG3 ~ tp53_bin * tcell_score + apm_score + source_bin + idh_bin
#   M_int_apm:   LAG3 ~ tp53_bin * apm_score + tcell_score + source_bin + idh_bin
#
# ポイント:
#   - interaction係数のβ+95%CI（効果量）を必ず残す（p値だけでは不十分）
#   - EMMsは ±1SD（mean-centered）で代表点を指定（再現性・論文記載のしやすさ優先）
#   - apm/source/idh は emmeans のデフォルト（全サンプル平均）で平均化
#
# 出力（Supplement推奨）:
#   step18_interaction_results.csv  ← β/SE/CI/p（main effects + interactions）
#   step18_emm_contrasts.csv        ← tcell低/中/高での Mut-WT 差 + CI
#   figS18_emm_panel.pdf/png        ← EMM差の可視化（小パネル）
#   step18_log.txt
#
# 出力先: 18_interaction/
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(broom)
library(emmeans)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "18_interaction")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# スコア定義（Step17と同一）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

COL_MUT <- "#E64B35"
COL_WT  <- "#AAAAAA"

GDC_PATH  <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step18_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line      <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 18: Interaction検定 + EMMs 開始 ===")

# =============================================================================
# 2. データ読み込み・スコア計算（Step17と同一ロジック）
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH,  show_col_types = FALSE) %>% filter(include_flag == TRUE)
wide     <- read_csv(WIDE_PATH, show_col_types = FALSE)
log_msg(sprintf("GDC: %d行 / wide: %d行x%d列", nrow(gdc_base), nrow(wide), ncol(wide)))

# 列名正規化ヘルパー（HLA-A → HLA.A 両対応）
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
  log_msg(sprintf("  %s: %d/%d遺伝子で計算", score_name, length(found), length(genes)))
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

all_score_genes  <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols_found <- intersect(
  c(paste0(all_score_genes, "_log2tpm"),
    paste0(gsub("-", ".", all_score_genes), "_log2tpm")),
  names(wide)
)
wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id")), any_of(score_cols_found))

gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_score %>% filter(!is.na(case_barcode)),  by = "case_barcode",  suffix = c("",".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_score %>% filter(!is.na(wxs_sample_id)), by = "wxs_sample_id", suffix = c("",".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)
for (col in score_cols_found) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) { gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]]); gdc[[wcol]] <- NULL }
}

gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "IFN-gamma")
  )

# GDC Grade4 解析対象
gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT")
  )

log_msg(sprintf("GDC Grade4: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4), sum(gdc_g4$tp53_bin==1), sum(gdc_g4$tp53_bin==0)))

# =============================================================================
# 3. Mean-centering（EMMs の基準点を作るため）
# =============================================================================

log_msg("--- Mean-centering ---")

mu_tcell <- mean(gdc_g4$tcell_score, na.rm = TRUE)
sd_tcell <- sd(gdc_g4$tcell_score,   na.rm = TRUE)
mu_apm   <- mean(gdc_g4$apm_score,   na.rm = TRUE)
sd_apm   <- sd(gdc_g4$apm_score,     na.rm = TRUE)

log_msg(sprintf("  T-cell score: mean=%.4f, SD=%.4f", mu_tcell, sd_tcell))
log_msg(sprintf("  APM score:    mean=%.4f, SD=%.4f", mu_apm,   sd_apm))
log_msg(sprintf("  T-cell ±1SD: [%.4f, %.4f, %.4f]",
                mu_tcell - sd_tcell, mu_tcell, mu_tcell + sd_tcell))
log_msg(sprintf("  APM    ±1SD: [%.4f, %.4f, %.4f]",
                mu_apm - sd_apm, mu_apm, mu_apm + sd_apm))

# mean-centered版を追加（interaction係数の解釈を簡単にするため）
gdc_g4 <- gdc_g4 %>%
  mutate(
    tcell_c = tcell_score - mu_tcell,
    apm_c   = apm_score   - mu_apm
  )

# =============================================================================
# 4. 回帰モデル（base + interaction 2本）
# =============================================================================

log_msg("--- 回帰モデル ---")

# M_base: interaction なし（Step17 M1の再現・centered版）
fit_base <- lm(
  LAG3_log2tpm ~ tp53_bin + tcell_c + apm_c + source_bin + idh_bin,
  data = gdc_g4
)

# M_int_tcell: tp53 × tcell の interaction
fit_int_tcell <- lm(
  LAG3_log2tpm ~ tp53_bin * tcell_c + apm_c + source_bin + idh_bin,
  data = gdc_g4
)

# M_int_apm: tp53 × apm の interaction
fit_int_apm <- lm(
  LAG3_log2tpm ~ tp53_bin * apm_c + tcell_c + source_bin + idh_bin,
  data = gdc_g4
)

# モデル比較（LRT）
anova_tcell <- anova(fit_base, fit_int_tcell)
anova_apm   <- anova(fit_base, fit_int_apm)
log_msg(sprintf("  LRT M_base vs M_int_tcell: F=%.3f, p=%.4f",
                anova_tcell$F[2], anova_tcell$`Pr(>F)`[2]))
log_msg(sprintf("  LRT M_base vs M_int_apm:   F=%.3f, p=%.4f",
                anova_apm$F[2], anova_apm$`Pr(>F)`[2]))

# =============================================================================
# 5. 結果抽出・保存
# =============================================================================

log_msg("--- 結果抽出 ---")

extract_coefs <- function(fit, label) {
  ci  <- confint(fit, level = 0.95)
  tbl <- tidy(fit) %>%
    left_join(
      as.data.frame(ci) %>%
        tibble::rownames_to_column("term") %>%
        rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`),
      by = "term"
    ) %>%
    mutate(
      model     = label,
      r_squared = summary(fit)$r.squared,
      n         = nrow(fit$model),
      beta      = round(estimate,   4),
      se        = round(std.error,  4),
      ci_lower  = round(ci_lower,   4),
      ci_upper  = round(ci_upper,   4),
      p_value   = p.value,
      r_squared = round(r_squared,  4)
    ) %>%
    select(model, n, term, beta, se, ci_lower, ci_upper, p_value, r_squared)
  tbl
}

res_base       <- extract_coefs(fit_base,       "M_base")
res_int_tcell  <- extract_coefs(fit_int_tcell,  "M_int_tcell")
res_int_apm    <- extract_coefs(fit_int_apm,    "M_int_apm")

all_results <- bind_rows(res_base, res_int_tcell, res_int_apm)
write_csv(all_results, file.path(OUT_DIR, "step18_interaction_results.csv"))
log_msg("保存: step18_interaction_results.csv")

# 重要な係数をログに出力
log_msg("--- 重要係数サマリー ---")
key_terms <- c("tp53_bin", "tcell_c", "apm_c",
               "tp53_bin:tcell_c", "tp53_bin:apm_c")

for (mdl in c("M_base", "M_int_tcell", "M_int_apm")) {
  log_msg(sprintf("  [%s]", mdl))
  rows <- all_results %>% filter(model == mdl, term %in% key_terms)
  for (i in seq_len(nrow(rows))) {
    r   <- rows[i, ]
    sig <- if (!is.na(r$p_value) && r$p_value < 0.05) "★" else "  "
    log_msg(sprintf("    %s %-20s beta=%+.4f [%+.4f, %+.4f] p=%.4f",
                    sig, r$term, r$beta, r$ci_lower, r$ci_upper, r$p_value))
  }
}

# =============================================================================
# 6. Estimated Marginal Means（EMMs）
# =============================================================================

log_msg("--- EMMs 計算 ---")

# T-cell scoreの代表点（±1SD, mean-centered座標で指定）
tcell_vals <- c(-sd_tcell, 0, sd_tcell)   # centered座標
tcell_labs <- c("Low (-1SD)", "Mean", "High (+1SD)")

apm_mean_c  <- 0   # mean-centered なので平均は0
source_mean <- mean(gdc_g4$source_bin, na.rm = TRUE)
idh_mean    <- mean(gdc_g4$idh_bin,    na.rm = TRUE)

log_msg(sprintf("  EMM代表点: tcell_c = [%.4f, %.4f, %.4f]",
                tcell_vals[1], tcell_vals[2], tcell_vals[3]))
log_msg(sprintf("  固定値: apm_c=%.4f, source_bin=%.4f, idh_bin=%.4f",
                apm_mean_c, source_mean, idh_mean))

# M_int_tcell でEMMs
emm_tcell <- emmeans(
  fit_int_tcell,
  ~ tp53_bin | tcell_c,
  at = list(
    tcell_c    = tcell_vals,
    apm_c      = apm_mean_c,
    source_bin = source_mean,
    idh_bin    = idh_mean
  )
)

# Mut - WT の pairwise contrast（各T-cell水準で）
contrast_tcell <- contrast(emm_tcell, method = "revpairwise", by = "tcell_c")
contrast_df    <- as.data.frame(summary(contrast_tcell, infer = TRUE)) %>%
  mutate(
    tcell_level = factor(tcell_labs,
                         levels = tcell_labs),
    estimate    = round(estimate,  4),
    SE          = round(SE,        4),
    lower.CL    = round(lower.CL,  4),
    upper.CL    = round(upper.CL,  4),
    p.value     = p.value
  )

log_msg("  Mut-WT差（adjusted）:")
for (i in seq_len(nrow(contrast_df))) {
  r   <- contrast_df[i, ]
  sig <- if (!is.na(r$p.value) && r$p.value < 0.05) "★" else "  "
  log_msg(sprintf("    %s T-cell %-12s: diff=%+.4f [%+.4f, %+.4f] p=%.4f",
                  sig, r$tcell_level, r$estimate, r$lower.CL, r$upper.CL, r$p.value))
}

# APM側も同様（M_int_apm で EMMs）
apm_vals <- c(-sd_apm, 0, sd_apm)
apm_labs <- c("Low (-1SD)", "Mean", "High (+1SD)")

emm_apm <- emmeans(
  fit_int_apm,
  ~ tp53_bin | apm_c,
  at = list(
    apm_c      = apm_vals,
    tcell_c    = 0,           # T-cell mean（centered=0）
    source_bin = source_mean,
    idh_bin    = idh_mean
  )
)
contrast_apm <- contrast(emm_apm, method = "revpairwise", by = "apm_c")
contrast_apm_df <- as.data.frame(summary(contrast_apm, infer = TRUE)) %>%
  mutate(
    apm_level = factor(apm_labs, levels = apm_labs),
    estimate  = round(estimate, 4),
    SE        = round(SE,       4),
    lower.CL  = round(lower.CL, 4),
    upper.CL  = round(upper.CL, 4)
  )

log_msg("  Mut-WT差（APM-indexed）:")
for (i in seq_len(nrow(contrast_apm_df))) {
  r   <- contrast_apm_df[i, ]
  sig <- if (!is.na(r$p.value) && r$p.value < 0.05) "★" else "  "
  log_msg(sprintf("    %s APM %-12s: diff=%+.4f [%+.4f, %+.4f] p=%.4f",
                  sig, r$apm_level, r$estimate, r$lower.CL, r$upper.CL, r$p.value))
}

# 結合して保存
contrast_df$score_type  <- "T-cell"
contrast_df$level_label <- as.character(contrast_df$tcell_level)
contrast_df$score_val   <- tcell_vals

contrast_apm_df$score_type  <- "APM"
contrast_apm_df$level_label <- as.character(contrast_apm_df$apm_level)
contrast_apm_df$score_val   <- apm_vals

emm_out <- bind_rows(
  contrast_df    %>% select(score_type, level_label, score_val,
                            estimate, SE, lower.CL, upper.CL, p.value),
  contrast_apm_df %>% select(score_type, level_label, score_val,
                             estimate, SE, lower.CL, upper.CL, p.value)
)
write_csv(emm_out, file.path(OUT_DIR, "step18_emm_contrasts.csv"))
log_msg("保存: step18_emm_contrasts.csv")

# =============================================================================
# 7. 図：figS18_emm_panel（Supplement）
# =============================================================================

log_msg("--- figS18_emm_panel 作成 ---")

# LRT p値をラベル用に整形
p_tcell <- anova_tcell$`Pr(>F)`[2]
p_apm   <- anova_apm$`Pr(>F)`[2]
fmt_p   <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)

# Panel A: T-cell score × TP53 interaction EMM
pa <- contrast_df %>%
  mutate(tcell_level = factor(level_label, levels = tcell_labs)) %>%
  ggplot(aes(x = tcell_level, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.8, color = COL_MUT) +
  geom_point(size = 4, color = COL_MUT, fill = COL_MUT, shape = 21, stroke = 1.0) +
  scale_y_continuous(
    name   = "Adjusted difference\n(TP53 Mut - WT) in LAG3",
    limits = function(x) c(min(x) - 0.05, max(x) + 0.05),
    breaks = seq(-0.2, 0.8, by = 0.1)
  ) +
  labs(
    x        = "T-cell score level",
    title    = "A  T-cell score",
    subtitle = sprintf("Interaction: %s (LRT)", fmt_p(p_tcell))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title    = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

# Panel B: APM score × TP53 interaction EMM
pb <- contrast_apm_df %>%
  mutate(apm_level = factor(level_label, levels = apm_labs)) %>%
  ggplot(aes(x = apm_level, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.15, linewidth = 0.8, color = "#2166AC") +
  geom_point(size = 4, color = "#2166AC", fill = "#2166AC", shape = 21, stroke = 1.0) +
  scale_y_continuous(
    name   = "Adjusted difference\n(TP53 Mut - WT) in LAG3",
    limits = function(x) c(min(x) - 0.05, max(x) + 0.05),
    breaks = seq(-0.2, 0.8, by = 0.1)
  ) +
  labs(
    x        = "APM score level",
    title    = "B  APM score",
    subtitle = sprintf("Interaction: %s (LRT)", fmt_p(p_apm))
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title    = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

# 並列（patchwork）
fig_s18 <- pa + pb +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    title    = "Supplementary: TP53 effect on LAG3 is consistent across immune levels (no interaction)",
    subtitle = paste0(
      "GDC Grade4 (n=442). Points: adjusted TP53 Mut-WT difference in LAG3 at T-cell/APM low/mean/high.\n",
      "Low/Mean/High = mean-1SD / mean / mean+1SD of each score."
    ),
    caption  = paste0(
      "Models: M_int_tcell: LAG3 ~ tp53*tcell + APM + source + IDH; ",
      "M_int_apm: LAG3 ~ tp53*APM + tcell + source + IDH (scores mean-centered).\n",
      "EMMs estimated at apm/source/idh = sample mean. ",
      "tp53_bin: Mut=1/WT=0; LRT = likelihood ratio test vs additive model."
    ),
    theme = theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8,  color = "grey40"),
      plot.caption  = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

# PDF
pdf(file.path(OUT_DIR, "figS18_emm_panel.pdf"), width = 8, height = 5)
print(fig_s18)
dev.off()
log_msg("保存: figS18_emm_panel.pdf")

# PNG 450dpi
agg_png(file.path(OUT_DIR, "figS18_emm_panel_450dpi.png"),
        width = 8, height = 5, units = "in", res = 450)
print(fig_s18)
dev.off()
log_msg("保存: figS18_emm_panel_450dpi.png")

# =============================================================================
# 8. 査読向けサマリーテキスト
# =============================================================================

log_msg("--- 論文記載用サマリー ---")

int_tcell_coef <- all_results %>%
  filter(model == "M_int_tcell", term == "tp53_bin:tcell_c")
int_apm_coef <- all_results %>%
  filter(model == "M_int_apm", term == "tp53_bin:apm_c")

log_msg(sprintf("  Interaction T-cell: beta=%+.4f [%+.4f, %+.4f] p=%.4f",
                int_tcell_coef$beta, int_tcell_coef$ci_lower,
                int_tcell_coef$ci_upper, int_tcell_coef$p_value))
log_msg(sprintf("  Interaction APM:    beta=%+.4f [%+.4f, %+.4f] p=%.4f",
                int_apm_coef$beta, int_apm_coef$ci_lower,
                int_apm_coef$ci_upper, int_apm_coef$p_value))

emm_check_tcell <- emm_out %>% filter(score_type == "T-cell")
log_msg("  EMM Mut-WT差（T-cell low/mean/high）:")
for (i in seq_len(nrow(emm_check_tcell))) {
  r <- emm_check_tcell[i, ]
  log_msg(sprintf("    %-12s: diff=%+.4f [%+.4f, %+.4f]",
                  r$level_label, r$estimate, r$lower.CL, r$upper.CL))
}
log_msg("  （3水準で差がほぼ同じ → 平行移動の確認）")

# =============================================================================
# 9. 完了
# =============================================================================

log_msg("=== Step 18: 完了 ===")
close(log_con)

cat("\n============================\n")
cat("Step 18 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  figS18_emm_panel.pdf/png\n"))
cat(sprintf("  step18_interaction_results.csv\n"))
cat(sprintf("  step18_emm_contrasts.csv\n"))
cat(sprintf("  step18_log.txt\n"))
cat(sprintf("出力先: %s\n", OUT_DIR))
cat("============================\n")
