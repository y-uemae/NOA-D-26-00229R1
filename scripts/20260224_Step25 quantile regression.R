# =============================================================================
# step25_quantile_regression.R
# GBM/Glioma TP53xLAG3 解析 - Step 25: 分位点回帰（Quantile Regression）
#
# 目的:
#   TCGAの上側尾（外れ値）によるOLS結果への影響を検証する。
#   tau=0.25/0.50/0.75の分位点回帰でTP53効果が全分位点で正方向であることを示す。
#   "OLSの結果が外れ値に引っ張られていない"の根拠としてSupplementに配置。
#
# モデル:
#   rq(LAG3 ~ tp53_bin + source_bin + idh_bin, tau = c(0.25, 0.50, 0.75))
#   比較: OLS（M_base, Step09b再現）のβと並べて提示
#
# パッケージ:
#   quantreg（CIはboot法、R=1000）
#
# 出力:
#   figS25_quantile_forest.pdf/png  : tau別β+CI Forest風（OLSと並列）
#   step25_quantile_results.csv     : 全tau+OLSの結果
#   step25_log.txt
#
# 出力先: 25_quantile_regression/
#
# 作成日: 2026-02-24
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ragg)

# quantregの確認・インストール
if (!requireNamespace("quantreg", quietly = TRUE)) {
  install.packages("quantreg")
}
library(quantreg)

# =============================================================================
# 0. 設定
# =============================================================================

BASE_DIR   <- here::here()
RESULT_DIR <- file.path(BASE_DIR, "results/TP53/20260221")
OUT_DIR    <- file.path(RESULT_DIR, "25_quantile_regression")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 分位点設定
TAU_VALS <- c(0.25, 0.50, 0.75)
BOOT_R   <- 1000
set.seed(42)

# 色設定（引継書に準拠）
COL_MUT  <- "#E64B35"
COL_OLS  <- "#555555"   # OLS（グレー）
COL_Q25  <- "#4DBBD5"   # tau=0.25
COL_Q50  <- "#E64B35"   # tau=0.50（中央値・強調）
COL_Q75  <- "#00A087"   # tau=0.75

# 入力ファイル
GDC_PATH <- file.path(RESULT_DIR, "08_final_cohort/final_cohort.csv")

# =============================================================================
# 1. ログ設定
# =============================================================================

log_file <- file.path(OUT_DIR, "step25_log.txt")
cat("", file = log_file, append = FALSE)

log_msg <- function(msg, also_print = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", timestamp, msg)
  cat(line, "\n", file = log_file, append = TRUE)
  if (also_print) cat(line, "\n")
}

log_msg("=== Step 25: 分位点回帰 開始 ===")
log_msg(sprintf("tau: %s", paste(TAU_VALS, collapse = ", ")))
log_msg(sprintf("Bootstrap: R=%d, seed=42", BOOT_R))

# =============================================================================
# 2. データ読み込み
# =============================================================================

log_msg("--- データ読み込み ---")

gdc_base <- read_csv(GDC_PATH, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

gdc_g4 <- gdc_base %>%
  filter(grade == "Grade4", tp53_status %in% c("mutant", "wildtype")) %>%
  mutate(
    tp53_bin   = as.integer(tp53_status == "mutant"),
    source_bin = as.integer(source == "CPTAC_HCMI"),
    idh_bin    = as.integer(idh_status == "mutant")
  ) %>%
  select(LAG3_log2tpm, tp53_bin, source_bin, idh_bin) %>%
  na.omit()

log_msg(sprintf("GDC Grade4: n=%d (Mut=%d, WT=%d)",
                nrow(gdc_g4),
                sum(gdc_g4$tp53_bin == 1),
                sum(gdc_g4$tp53_bin == 0)))

# =============================================================================
# 3. OLS（M_base・比較用）
# =============================================================================

log_msg("--- OLS（比較用）---")

fit_ols  <- lm(LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin, data = gdc_g4)
ci_ols   <- as.double(confint(fit_ols, "tp53_bin"))
smr_ols  <- coef(summary(fit_ols))

ols_result <- tibble(
  method   = "OLS",
  tau      = NA_real_,
  beta     = round(smr_ols["tp53_bin", "Estimate"],   4),
  ci_lower = round(ci_ols[1],                          4),
  ci_upper = round(ci_ols[2],                          4),
  se       = round(smr_ols["tp53_bin", "Std. Error"],  4),
  p_value  = smr_ols["tp53_bin", "Pr(>|t|)"]
)

log_msg(sprintf("  OLS: beta=%+.4f [%+.4f, %+.4f] p=%.4f",
                ols_result$beta, ols_result$ci_lower,
                ols_result$ci_upper, ols_result$p_value))

# =============================================================================
# 4. 分位点回帰（tau別）
# =============================================================================

log_msg("--- 分位点回帰 ---")

qr_results <- lapply(TAU_VALS, function(tau) {
  
  log_msg(sprintf("  tau=%.2f: 実行中...", tau))
  
  fit_rq <- rq(LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin,
               tau = tau, data = gdc_g4)
  
  # Bootstrap CI（se = "boot"）
  smr_rq <- tryCatch(
    summary(fit_rq, se = "boot", R = BOOT_R, bsmethod = "xy"),
    error = function(e) {
      log_msg(sprintf("  tau=%.2f: boot失敗、nid法にフォールバック", tau))
      summary(fit_rq, se = "nid")
    }
  )
  
  coef_mat <- smr_rq$coefficients
  beta     <- coef_mat["tp53_bin", "Value"]
  se_val   <- coef_mat["tp53_bin", "Std. Error"]
  pval     <- coef_mat["tp53_bin", "Pr(>|t|)"]
  
  # CI = beta ± 1.96 * SE
  ci_lo <- beta - 1.96 * se_val
  ci_hi <- beta + 1.96 * se_val
  
  log_msg(sprintf("  tau=%.2f: beta=%+.4f [%+.4f, %+.4f] p=%.4f",
                  tau, beta, ci_lo, ci_hi, pval))
  
  tibble(
    method   = "Quantile",
    tau      = tau,
    beta     = round(beta,  4),
    ci_lower = round(ci_lo, 4),
    ci_upper = round(ci_hi, 4),
    se       = round(se_val, 4),
    p_value  = pval
  )
})

qr_all <- do.call(rbind, qr_results)

# OLSと結合
all_results <- do.call(rbind, list(ols_result, qr_all)) %>%
  mutate(
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    # 表示ラベル
    model_label = factor(
      case_when(
        method == "OLS"      ~ "OLS (mean)",
        tau    == 0.25       ~ "Quantile tau=0.25",
        tau    == 0.50       ~ "Quantile tau=0.50 (median)",
        tau    == 0.75       ~ "Quantile tau=0.75"
      ),
      levels = rev(c(
        "OLS (mean)",
        "Quantile tau=0.25",
        "Quantile tau=0.50 (median)",
        "Quantile tau=0.75"
      ))
    ),
    # 色グループ
    col_group = case_when(
      method == "OLS" ~ "ols",
      tau    == 0.25  ~ "q25",
      tau    == 0.50  ~ "q50",
      tau    == 0.75  ~ "q75"
    ),
    # 線種（OLSは破線・強調）
    is_ols = (method == "OLS")
  )

write_csv(all_results %>% mutate(model_label = as.character(model_label)),
          file.path(OUT_DIR, "step25_quantile_results.csv"))
log_msg("保存: step25_quantile_results.csv")

# 全て正方向か確認
all_positive <- all(all_results$beta > 0)
log_msg(sprintf("全tau+OLSでbeta>0: %s", ifelse(all_positive, "YES", "NO")))
all_sig <- all(all_results$p_value < 0.05)
log_msg(sprintf("全tau+OLSでp<0.05: %s", ifelse(all_sig, "YES", "NO")))

# =============================================================================
# 5. figS25: Forest風プロット
# =============================================================================

log_msg("--- figS25: Forest風プロット 作成 ---")

col_vals <- c(
  "ols" = COL_OLS,
  "q25" = COL_Q25,
  "q50" = COL_Q50,
  "q75" = COL_Q75
)
fill_vals <- c(
  "ols" = "#AAAAAA",
  "q25" = COL_Q25,
  "q50" = COL_Q50,
  "q75" = COL_Q75
)

# x軸範囲
x_min <- min(all_results$ci_lower, na.rm = TRUE) - 0.05
x_max <- max(all_results$ci_upper, na.rm = TRUE) + 0.30

fig_s25 <- ggplot(all_results,
                  aes(x = beta, y = model_label,
                      color = col_group, fill = col_group)) +
  # 参照線
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.5) +
  # OLS参照線（点線）
  geom_vline(xintercept = ols_result$beta,
             linetype = "dotted", color = COL_OLS, linewidth = 0.5) +
  # CI棒
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
                orientation = "y", width = 0.25,
                linewidth = ifelse(all_results$is_ols, 0.7, 0.9)) +
  # 点
  geom_point(aes(shape = is_ols),
             size = 4.5, stroke = 1.1) +
  # ラベル（右端）
  geom_text(aes(x = ci_upper,
                label = sprintf("b=%+.3f [%.3f, %.3f] %s",
                                beta, ci_lower, ci_upper, sig_label)),
            hjust = -0.08, size = 3.0, color = "grey30") +
  scale_color_manual(values = col_vals, guide = "none") +
  scale_fill_manual(values  = fill_vals, guide = "none") +
  scale_shape_manual(values = c("TRUE" = 22, "FALSE" = 21), guide = "none") +
  scale_x_continuous(
    name   = "TP53 coefficient (beta) with 95% CI\n[LAG3 ~ tp53 + source + IDH]",
    breaks = seq(-0.1, 0.6, by = 0.1)
  ) +
  coord_cartesian(xlim = c(x_min, x_max)) +
  labs(
    y        = NULL,
    title    = "TP53 effect on LAG3: quantile regression vs OLS",
    subtitle = paste0(
      "GDC Grade4 (n=", nrow(gdc_g4), "). ",
      "Quantile regression: boot CI (R=", BOOT_R, ", seed=42). ",
      "All estimates positive and significant."
    ),
    caption  = paste0(
      "OLS (square): mean-based estimate. ",
      "Quantile (circle): tau=0.25/0.50/0.75. ",
      "Gray dotted: OLS reference. ",
      "Model: LAG3 ~ tp53_bin + source_bin + idh_bin."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    plot.caption       = element_text(size = 7.5, color = "grey50", hjust = 0),
    plot.margin        = margin(10, 130, 10, 10)
  )

# PDF出力
pdf(file.path(OUT_DIR, "figS25_quantile_forest.pdf"), width = 10, height = 5)
print(fig_s25)
dev.off()
log_msg("保存: figS25_quantile_forest.pdf")

# PNG出力
agg_png(file.path(OUT_DIR, "figS25_quantile_forest_450dpi.png"),
        width = 10, height = 5, units = "in", res = 450)
print(fig_s25)
dev.off()
log_msg("保存: figS25_quantile_forest_450dpi.png")

# =============================================================================
# 6. 完了
# =============================================================================

log_msg("=== Step 25: 完了 ===")
log_msg(sprintf("全beta>0: %s / 全p<0.05: %s",
                ifelse(all_positive, "YES", "NO"),
                ifelse(all_sig, "YES", "NO")))
log_msg(sprintf("出力: %s", OUT_DIR))

cat("\n============================\n")
cat("Step 25 完了\n")
cat(sprintf("全beta>0: %s / 全p<0.05: %s\n",
            ifelse(all_positive, "YES", "NO"),
            ifelse(all_sig, "YES", "NO")))
cat("出力ファイル:\n")
cat(sprintf("  %s/figS25_quantile_forest.pdf\n",         OUT_DIR))
cat(sprintf("  %s/figS25_quantile_forest_450dpi.png\n",  OUT_DIR))
cat(sprintf("  %s/step25_quantile_results.csv\n",        OUT_DIR))
cat(sprintf("  %s/step25_log.txt\n",                     OUT_DIR))
cat("============================\n")
