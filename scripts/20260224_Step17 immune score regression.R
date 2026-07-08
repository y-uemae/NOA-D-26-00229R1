# =============================================================================
# step17_immune_score_regression.R
# GBM/Glioma TP53×LAG3 解析 - Step 17: 免疫スコア作成 + 調整回帰 + Fig3
#
# 目的:
#   T-cell / APM / IFNγ スコアを既存28遺伝子から作成し、
#   「免疫状態で調整してもTP53効果が残るか」を回帰で示す。
#   Fig3（GDC Grade4メイン） + Supplement（GLASS同方向確認）を出力。
#
# スコア定義（事前固定）:
#   T-cell score : CD3D, CD3E, CD3G, CD8A, CD8B, GZMA, GZMB, PRF1  (8遺伝子)
#   APM score    : B2M, TAP1, TAP2, TAPBP, HLA-A, HLA-B, HLA-C, NLRC5  (8遺伝子)
#   IFNγ score   : STAT1, IRF1, IRF9, CXCL9, CXCL10, CXCL11,
#                  GBP1, GBP2, GBP4, GBP5, IDO1  (11遺伝子)
#   算出: 各遺伝子の log2(TPM+1) の算術平均（scale前）
#
# 回帰モデル（GDC Grade4）:
#   M0: LAG3 ~ tp53_bin + source_bin + idh_bin               (基準・Step09b再現)
#   M1: LAG3 ~ tp53_bin + source_bin + idh_bin + tcell_score + apm_score
#   M2: LAG3 ~ tp53_bin + source_bin + idh_bin + tcell_score + apm_score + ifng_score
#
# Fig3 構成:
#   Panel A-left : LAG3 vs Tcell score（TP53 Mut/WT で色分け、GDC Grade4）
#   Panel A-right: LAG3 vs APM score（同上）
#   Panel B      : TP53係数（β）の調整前→調整後 変化（Forest風 dot+CI）
#
# Supplement:
#   figS17_glass_immune_score.pdf/png
#   GLASS WXS（n=79）での同方向確認（点を薄く）
#
# 出力先: 17_immune_score/
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(broom)
library(ragg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "17_immune_score")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# スコア定義（事前固定）
TCELL_GENES <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1")
APM_GENES   <- c("B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5")
IFNG_GENES  <- c("STAT1", "IRF1", "IRF9", "CXCL9", "CXCL10", "CXCL11",
                 "GBP1", "GBP2", "GBP4", "GBP5", "IDO1")

# 色設定（引継書10章に準拠）
COL_MUT  <- "#E64B35"   # TP53 Mut
COL_WT   <- "#AAAAAA"   # TP53 WT
COL_TCGA <- "#3C5488"
COL_CPTAC <- "#E07B54"

# 入力ファイル
GDC_PATH   <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")
WIDE_PATH  <- file.path(RESULT_DIR, "06_gene_expression/gene_expression_wide.csv")
GLASS_PATH <- file.path(RESULT_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step17_log.txt")
log_con  <- file(log_file, open = "wt")

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  writeLines(line, con = log_con)
  if (also_print) message(line)
}

log_msg("=== Step 17: 免疫スコア作成 + 調整回帰 + Fig3 開始 ===")
log_msg(sprintf("T-cell genes (%d): %s", length(TCELL_GENES), paste(TCELL_GENES, collapse=", ")))
log_msg(sprintf("APM genes    (%d): %s", length(APM_GENES),   paste(APM_GENES,   collapse=", ")))
log_msg(sprintf("IFNg genes   (%d): %s", length(IFNG_GENES),  paste(IFNG_GENES,  collapse=", ")))

# =============================================================================
# 2. データ読み込み・結合
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)
wide <- read_csv(WIDE_PATH, show_col_types = FALSE)

log_msg(sprintf("GDC: %d行 / wide: %d行×%d列", nrow(gdc_base), nrow(wide), ncol(wide)))

# スコア計算に必要な列名（log2tpm）
all_score_genes <- c(TCELL_GENES, APM_GENES, IFNG_GENES)
score_cols      <- paste0(gsub("-", ".", all_score_genes), "_log2tpm")
# HLA-A等のハイフンはRの列名でドットになっている場合があるため確認
actual_cols <- names(wide)
score_cols_found <- intersect(
  c(paste0(all_score_genes, "_log2tpm"),
    paste0(gsub("-", ".", all_score_genes), "_log2tpm")),
  actual_cols
)
log_msg(sprintf("スコア遺伝子列（wide内に存在）: %d / %d",
                length(score_cols_found), length(all_score_genes)))

# 不足列の確認
genes_from_cols <- gsub("_log2tpm$", "", score_cols_found)
genes_from_cols <- gsub("\\.", "-", genes_from_cols)  # ドット→ハイフンに戻す
missing_genes <- setdiff(all_score_genes, genes_from_cols)
if (length(missing_genes) > 0) {
  log_msg(sprintf("WARNING: 列が見つからない遺伝子: %s", paste(missing_genes, collapse=", ")))
}

# wideからスコア用列を抽出
wide_score <- wide %>%
  select(any_of(c("case_barcode", "wxs_sample_id", "rna_file_id")),
         any_of(score_cols_found))

# GDC結合（Step16aと同じロジック）
gdc_tcga  <- gdc_base %>% filter(source == "TCGA") %>%
  left_join(wide_score %>% filter(!is.na(case_barcode)),
            by = "case_barcode", suffix = c("", ".w"))
gdc_cptac <- gdc_base %>% filter(source == "CPTAC_HCMI") %>%
  left_join(wide_score %>% filter(!is.na(wxs_sample_id)),
            by = "wxs_sample_id", suffix = c("", ".w"))
gdc <- bind_rows(gdc_tcga, gdc_cptac)

# .w列のマージ
for (col in score_cols_found) {
  wcol <- paste0(col, ".w")
  if (wcol %in% names(gdc)) {
    gdc[[col]] <- coalesce(gdc[[wcol]], gdc[[col]])
    gdc[[wcol]] <- NULL
  }
}
log_msg(sprintf("GDC結合後: %d行", nrow(gdc)))

# =============================================================================
# 3. 免疫スコア計算（各遺伝子log2tpmの算術平均）
# =============================================================================

log_msg("--- 免疫スコア計算 ---")

# 列名正規化ヘルパー（HLA-A → HLA-A or HLA.A 両方対応）
get_log2_col <- function(gene, df) {
  c1 <- paste0(gene, "_log2tpm")
  c2 <- paste0(gsub("-", ".", gene), "_log2tpm")
  if (c1 %in% names(df)) return(c1)
  if (c2 %in% names(df)) return(c2)
  return(NA_character_)
}

calc_score <- function(df, genes, score_name) {
  cols <- sapply(genes, get_log2_col, df = df)
  found <- cols[!is.na(cols)]
  missing <- genes[is.na(cols)]
  if (length(missing) > 0)
    log_msg(sprintf("  %s: 不足遺伝子 %s", score_name, paste(missing, collapse=",")))
  log_msg(sprintf("  %s: %d/%d遺伝子で計算", score_name, length(found), length(genes)))
  
  mat <- df %>% select(all_of(found)) %>% mutate(across(everything(), as.numeric))
  rowMeans(mat, na.rm = TRUE)
}

# GDC スコア
gdc <- gdc %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "T-cell"),
    apm_score   = calc_score(., APM_GENES,   "APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "IFNg")
  )

# GLASS スコア
glass_base <- read_csv(GLASS_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE, tp53_status %in% c("Mut", "WT"))

glass_base <- glass_base %>%
  mutate(
    tcell_score = calc_score(., TCELL_GENES, "GLASS T-cell"),
    apm_score   = calc_score(., APM_GENES,   "GLASS APM"),
    ifng_score  = calc_score(., IFNG_GENES,  "GLASS IFNg")
  )

log_msg(sprintf("スコアNA率確認（GDC Grade4）:"))
gdc_g4 <- gdc %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant"),
    tp53_label = ifelse(tp53_status == "mutant", "TP53 Mut", "TP53 WT")
  )

for (sc in c("tcell_score", "apm_score", "ifng_score")) {
  log_msg(sprintf("  %s: NA=%d, median=%.3f",
                  sc, sum(is.na(gdc_g4[[sc]])),
                  median(gdc_g4[[sc]], na.rm = TRUE)))
}
log_msg(sprintf("GDC Grade4 解析対象: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$tp53_bin==1), sum(gdc_g4$tp53_bin==0)))

# =============================================================================
# 4. スコア・相関確認ログ
# =============================================================================

log_msg("--- スコア相関確認（GDC Grade4） ---")
cor_lag3_tcell <- cor(gdc_g4$LAG3_log2tpm, gdc_g4$tcell_score, use="complete.obs")
cor_lag3_apm   <- cor(gdc_g4$LAG3_log2tpm, gdc_g4$apm_score,   use="complete.obs")
cor_lag3_ifng  <- cor(gdc_g4$LAG3_log2tpm, gdc_g4$ifng_score,  use="complete.obs")
cor_tcell_apm  <- cor(gdc_g4$tcell_score,  gdc_g4$apm_score,   use="complete.obs")
log_msg(sprintf("  LAG3 vs Tcell: r=%.3f", cor_lag3_tcell))
log_msg(sprintf("  LAG3 vs APM:   r=%.3f", cor_lag3_apm))
log_msg(sprintf("  LAG3 vs IFNg:  r=%.3f", cor_lag3_ifng))
log_msg(sprintf("  Tcell vs APM:  r=%.3f (共線性確認)", cor_tcell_apm))

# =============================================================================
# 5. 回帰解析（M0/M1/M2）
# =============================================================================

log_msg("--- 回帰解析 ---")

run_model <- function(df, formula_str, label) {
  df_use <- df %>%
    select(all_of(all.vars(as.formula(formula_str)))) %>%
    na.omit()
  fit <- lm(as.formula(formula_str), data = df_use)
  ci  <- confint(fit, "tp53_bin", level = 0.95)
  tbl <- tidy(fit) %>% filter(term == "tp53_bin")
  tibble(
    model     = label,
    formula   = formula_str,
    n         = nrow(df_use),
    beta      = round(tbl$estimate,   4),
    ci_lower  = round(ci[1],          4),
    ci_upper  = round(ci[2],          4),
    se        = round(tbl$std.error,  4),
    p_value   = tbl$p.value,
    r_squared = round(summary(fit)$r.squared, 4)
  )
}

models <- list(
  M0 = "LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin",
  M1 = "LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin + tcell_score + apm_score",
  M2 = "LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin + tcell_score + apm_score + ifng_score"
)

reg_results <- bind_rows(mapply(run_model,
                                formula_str = models,
                                label       = names(models),
                                MoreArgs    = list(df = gdc_g4),
                                SIMPLIFY    = FALSE))

# β変化率の計算
beta_m0 <- reg_results$beta[reg_results$model == "M0"]
reg_results <- reg_results %>%
  mutate(
    beta_change_pct = round((beta - beta_m0) / abs(beta_m0) * 100, 1),
    model_label = case_when(
      model == "M0" ~ "M0: TP53 + source + IDH\n(base)",
      model == "M1" ~ "M1: + T-cell + APM",
      model == "M2" ~ "M2: + T-cell + APM + IFN\u03b3"
    )
  )

log_msg("回帰結果:")
for (i in seq_len(nrow(reg_results))) {
  r <- reg_results[i,]
  log_msg(sprintf("  %s: beta=%+.4f [%+.4f, %+.4f] p=%.2e R2=%.3f (beta_change=%+.1f%%)",
                  r$model, r$beta, r$ci_lower, r$ci_upper,
                  r$p_value, r$r_squared, r$beta_change_pct))
}

# 保存
write_csv(reg_results, file.path(OUT_DIR, "step17_regression_results.csv"))
log_msg("保存: step17_regression_results.csv")

# =============================================================================
# 6. Fig3 作成
# =============================================================================

log_msg("--- Fig3 作成 ---")

# --- Panel A-left: LAG3 vs T-cell score ---
pa_left <- ggplot(gdc_g4, aes(x = tcell_score, y = LAG3_log2tpm,
                              color = tp53_label, fill = tp53_label)) +
  geom_point(alpha = 0.45, size = 1.2, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
  scale_color_manual(values = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT),
                     name = NULL) +
  scale_fill_manual(values  = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT),
                    name = NULL) +
  labs(x = "T-cell score [mean log2(TPM+1)]",
       y = "LAG3 expression\n[log2(TPM+1)]",
       title = "A") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.2f", cor_lag3_tcell),
           hjust = -0.2, vjust = 1.4, size = 3, color = "grey30") +
  theme_bw(base_size = 10) +
  theme(legend.position   = "bottom",
        legend.text       = element_text(size = 8.5),
        panel.grid.minor  = element_blank(),
        plot.title        = element_text(face = "bold", size = 11))

# --- Panel A-right: LAG3 vs APM score ---
pa_right <- ggplot(gdc_g4, aes(x = apm_score, y = LAG3_log2tpm,
                               color = tp53_label, fill = tp53_label)) +
  geom_point(alpha = 0.45, size = 1.2, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
  scale_color_manual(values = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT),
                     name = NULL) +
  scale_fill_manual(values  = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT),
                    name = NULL) +
  labs(x = "APM score [mean log2(TPM+1)]",
       y = NULL,
       title = " ") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.2f", cor_lag3_apm),
           hjust = -0.2, vjust = 1.4, size = 3, color = "grey30") +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 11))

# --- Panel B: TP53係数の調整前→調整後 Forest風 ---
pb_data <- reg_results %>%
  mutate(
    model_label = factor(model_label,
                         levels = rev(c(
                           "M0: TP53 + source + IDH\n(base)",
                           "M1: + T-cell + APM",
                           paste0("M2: + T-cell + APM + IFN\u03b3")
                         ))),
    beta_pct_label = sprintf("%+.1f%%", beta_change_pct),
    is_base = (model == "M0")
  )

pb <- ggplot(pb_data, aes(x = beta, y = model_label,
                          color = is_base, fill = is_base)) +
  geom_vline(xintercept = 0,       linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_vline(xintercept = beta_m0, linetype = "dotted", color = COL_MUT,   linewidth = 0.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.2, linewidth = 0.7) +
  geom_point(shape = 21, size = 4, stroke = 1.0) +
  geom_text(aes(label = beta_pct_label, x = ci_upper),
            hjust = -0.15, size = 3, color = "grey30") +
  scale_color_manual(values = c("TRUE" = COL_MUT, "FALSE" = "#2166AC"), guide = "none") +
  scale_fill_manual(values  = c("TRUE" = COL_MUT, "FALSE" = "#2166AC"), guide = "none") +
  scale_x_continuous(
    name   = expression(paste("TP53 coefficient (", beta[TP53], ") with 95% CI")),
    limits = c(NA, max(pb_data$ci_upper, na.rm = TRUE) * 1.25),
    breaks = seq(0, 0.5, by = 0.1)
  ) +
  labs(y = NULL, title = "B") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8.5),
        plot.title  = element_text(face = "bold", size = 11))

# --- 結合 ---
fig3 <- (pa_left + pa_right) / pb +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title    = "TP53 mutation is independently associated with LAG3 expression after immune adjustment",
    subtitle = "GDC Grade4 (n=442). Scores: mean log2(TPM+1) of component genes.",
    caption  = paste0(
      "T-cell score: CD3D/E/G, CD8A/B, GZMA/B, PRF1 (n=8). ",
      "APM score: B2M, TAP1/2, TAPBP, HLA-A/B/C, NLRC5 (n=8). ",
      "IFN\u03b3 score: STAT1, IRF1/9, CXCL9/10/11, GBP1/2/4/5, IDO1 (n=11).\n",
      "M0: LAG3 ~ tp53 + source + IDH. ",
      "M1: + T-cell + APM. M2: + T-cell + APM + IFN\u03b3. ",
      "Red dotted line: M0 (base) \u03b2. % = change from M0."
    ),
    theme = theme(
      plot.title    = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "grey40"),
      plot.caption  = element_text(size = 7,  color = "grey50", hjust = 0)
    )
  )

# --- PDF出力 ---
pdf_path <- file.path(OUT_DIR, "fig3_immune_adjustment.pdf")
pdf(pdf_path, width = 10, height = 8)
print(fig3)
dev.off()
log_msg("保存: fig3_immune_adjustment.pdf")

# --- PNG出力 ---
png_path <- file.path(OUT_DIR, "fig3_immune_adjustment_450dpi.png")
agg_png(png_path, width = 10, height = 8, units = "in", res = 450)
print(fig3)
dev.off()
log_msg("保存: fig3_immune_adjustment_450dpi.png")

# =============================================================================
# 7. Supplement図（GLASS 同方向確認）
# =============================================================================

log_msg("--- Supplement: GLASS 同方向確認 ---")

glass_g4 <- glass_base %>%
  mutate(
    tp53_label = ifelse(tp53_status == "Mut", "TP53 Mut", "TP53 WT"),
    tp53_bin   = as.integer(tp53_status == "Mut"),
    idh_bin    = as.integer(idh_status  == "Mut")
  )

cor_glass_tcell <- cor(glass_g4$LAG3_log2tpm, glass_g4$tcell_score, use="complete.obs")
cor_glass_apm   <- cor(glass_g4$LAG3_log2tpm, glass_g4$apm_score,   use="complete.obs")
log_msg(sprintf("GLASS: LAG3 vs Tcell r=%.3f, LAG3 vs APM r=%.3f",
                cor_glass_tcell, cor_glass_apm))

# GLASS 回帰
glass_reg <- bind_rows(
  run_model(glass_g4, "LAG3_log2tpm ~ tp53_bin", "M0_GLASS"),
  run_model(glass_g4, "LAG3_log2tpm ~ tp53_bin + tcell_score + apm_score", "M1_GLASS"),
  run_model(glass_g4, "LAG3_log2tpm ~ tp53_bin + tcell_score + apm_score + ifng_score", "M2_GLASS")
) %>%
  mutate(
    beta_change_pct = round((beta - beta[model=="M0_GLASS"]) /
                              abs(beta[model=="M0_GLASS"]) * 100, 1)
  )

log_msg("GLASS 回帰結果:")
for (i in seq_len(nrow(glass_reg))) {
  r <- glass_reg[i,]
  log_msg(sprintf("  %s: beta=%+.4f [%+.4f, %+.4f] p=%.3f (change=%+.1f%%)",
                  r$model, r$beta, r$ci_lower, r$ci_upper,
                  r$p_value, r$beta_change_pct))
}

write_csv(glass_reg, file.path(OUT_DIR, "step17_glass_regression.csv"))
log_msg("保存: step17_glass_regression.csv")

# GLASS 散布図（薄め・参考用）
ps_left <- ggplot(glass_g4, aes(x = tcell_score, y = LAG3_log2tpm,
                                color = tp53_label)) +
  geom_point(alpha = 0.5, size = 1.5, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15,
              aes(fill = tp53_label)) +
  scale_color_manual(values = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT), name = NULL) +
  scale_fill_manual(values  = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT), name = NULL) +
  labs(x = "T-cell score", y = "LAG3 [log2(TPM+1)]",
       title = sprintf("GLASS WXS (n=%d)", nrow(glass_g4)),
       subtitle = sprintf("r = %.2f", cor_glass_tcell)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ps_right <- ggplot(glass_g4, aes(x = apm_score, y = LAG3_log2tpm,
                                 color = tp53_label)) +
  geom_point(alpha = 0.5, size = 1.5, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15,
              aes(fill = tp53_label)) +
  scale_color_manual(values = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT), name = NULL) +
  scale_fill_manual(values  = c("TP53 Mut" = COL_MUT, "TP53 WT" = COL_WT), name = NULL) +
  labs(x = "APM score", y = NULL,
       subtitle = sprintf("r = %.2f", cor_glass_apm)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# GLASS forest風
beta_glass_m0 <- glass_reg$beta[glass_reg$model == "M0_GLASS"]
ps_forest <- glass_reg %>%
  mutate(
    model_label = factor(
      c("M0: TP53 only", "M1: + T-cell + APM",
        paste0("M2: + T-cell + APM + IFN\u03b3")),
      levels = rev(c("M0: TP53 only", "M1: + T-cell + APM",
                     paste0("M2: + T-cell + APM + IFN\u03b3")))
    ),
    is_base = (model == "M0_GLASS")
  ) %>%
  ggplot(aes(x = beta, y = model_label, color = is_base, fill = is_base)) +
  geom_vline(xintercept = 0,             linetype = "dashed", color="#888888", linewidth=0.4) +
  geom_vline(xintercept = beta_glass_m0, linetype = "dotted", color=COL_MUT,   linewidth=0.5) +
  geom_errorbar(aes(xmin=ci_lower, xmax=ci_upper),
                orientation="y", width=0.2, linewidth=0.7) +
  geom_point(shape=21, size=4, stroke=1.0) +
  scale_color_manual(values=c("TRUE"=COL_MUT, "FALSE"="#2166AC"), guide="none") +
  scale_fill_manual(values =c("TRUE"=COL_MUT, "FALSE"="#2166AC"), guide="none") +
  labs(x = expression(paste(beta[TP53], " with 95% CI")), y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank())

fig_s17 <- (ps_left + ps_right) / ps_forest +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title   = "Supplementary: GLASS validation - immune score adjustment",
    caption = "GLASS WXS-only, TCGA-excluded (n=79). Same direction as GDC.",
    theme   = theme(plot.title   = element_text(size=10, face="bold"),
                    plot.caption = element_text(size=7, color="grey50", hjust=0))
  )

pdf(file.path(OUT_DIR, "figS17_glass_immune_score.pdf"), width=10, height=8)
print(fig_s17)
dev.off()
agg_png(file.path(OUT_DIR, "figS17_glass_immune_score_450dpi.png"),
        width=10, height=8, units="in", res=450)
print(fig_s17)
dev.off()
log_msg("保存: figS17_glass_immune_score.pdf/png")

# =============================================================================
# 8. スコア記述統計の保存
# =============================================================================

score_stats <- gdc_g4 %>%
  group_by(tp53_label) %>%
  summarise(
    n            = n(),
    tcell_median = round(median(tcell_score, na.rm=TRUE), 3),
    tcell_Q1     = round(quantile(tcell_score, .25, na.rm=TRUE), 3),
    tcell_Q3     = round(quantile(tcell_score, .75, na.rm=TRUE), 3),
    apm_median   = round(median(apm_score,   na.rm=TRUE), 3),
    apm_Q1       = round(quantile(apm_score,   .25, na.rm=TRUE), 3),
    apm_Q3       = round(quantile(apm_score,   .75, na.rm=TRUE), 3),
    ifng_median  = round(median(ifng_score,  na.rm=TRUE), 3),
    ifng_Q1      = round(quantile(ifng_score,  .25, na.rm=TRUE), 3),
    ifng_Q3      = round(quantile(ifng_score,  .75, na.rm=TRUE), 3),
    lag3_median  = round(median(LAG3_log2tpm, na.rm=TRUE), 3),
    .groups = "drop"
  )

write_csv(score_stats, file.path(OUT_DIR, "step17_score_descriptives.csv"))
log_msg("保存: step17_score_descriptives.csv")

log_msg("--- スコア記述統計（GDC Grade4） ---")
for (i in seq_len(nrow(score_stats))) {
  r <- score_stats[i,]
  log_msg(sprintf("  [%s] n=%d | Tcell=%.3f [%.3f-%.3f] | APM=%.3f [%.3f-%.3f] | IFNg=%.3f [%.3f-%.3f]",
                  r$tp53_label, r$n,
                  r$tcell_median, r$tcell_Q1, r$tcell_Q3,
                  r$apm_median,  r$apm_Q1,   r$apm_Q3,
                  r$ifng_median, r$ifng_Q1,  r$ifng_Q3))
}

# =============================================================================
# 9. 完了
# =============================================================================

log_msg("=== Step 17: 完了 ===")
log_msg(sprintf("出力: %s", OUT_DIR))
close(log_con)

cat("\n============================\n")
cat("Step 17 完了\n")
cat("出力ファイル:\n")
cat(sprintf("  %s/fig3_immune_adjustment.pdf\n",           OUT_DIR))
cat(sprintf("  %s/fig3_immune_adjustment_450dpi.png\n",    OUT_DIR))
cat(sprintf("  %s/figS17_glass_immune_score.pdf\n",        OUT_DIR))
cat(sprintf("  %s/figS17_glass_immune_score_450dpi.png\n", OUT_DIR))
cat(sprintf("  %s/step17_regression_results.csv\n",        OUT_DIR))
cat(sprintf("  %s/step17_glass_regression.csv\n",          OUT_DIR))
cat(sprintf("  %s/step17_score_descriptives.csv\n",        OUT_DIR))
cat(sprintf("  %s/step17_log.txt\n",                       OUT_DIR))
cat("============================\n")
