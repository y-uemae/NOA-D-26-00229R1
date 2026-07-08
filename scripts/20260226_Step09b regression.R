# ============================================================
# step09b_regression.R  Grade4 回帰解析（共変量調整）
# 目的: TP53効果がsource/IDH調整後も残るか確認
# 入力: 08_final_cohort/final_cohort.csv
# 出力: 09_statistics/step09b_regression_results.csv      （TP53行のみ・既存）
#       09_statistics/step09b_regression_results_full.csv （全係数・新規追加）
# ============================================================

library(tidyverse)

BASE_DIR <- here::here("results", "TP53", "20260221")
IN_FILE  <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
OUT_DIR  <- file.path(BASE_DIR, "09_statistics")

cohort <- read_csv(IN_FILE, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

g4 <- cohort %>%
  filter(grade == "Grade4") %>%
  mutate(
    tp53_bin   = if_else(tp53_status == "mutant",    1L, 0L),
    idh_bin    = if_else(idh_status  == "mutant",    1L, 0L),
    source_bin = if_else(source      == "TCGA",      1L, 0L)  # TCGA=1, CPTAC_HCMI=0
  )

cat("Grade4 n=", nrow(g4), "\n")
cat("TP53 mutant:", sum(g4$tp53_bin), "/ wildtype:", sum(1 - g4$tp53_bin), "\n")
cat("IDH mutant:", sum(g4$idh_bin), "\n")
cat("Source TCGA:", sum(g4$source_bin), "/ CPTAC_HCMI:", sum(1 - g4$source_bin), "\n\n")

# ============================================================
# モデル1: 単変量（TP53のみ）
# モデル2: TP53 + source
# モデル3: TP53 + source + IDH（メインモデル M0）
# ============================================================

fit1      <- lm(LAG3_log2tpm ~ tp53_bin,                        data = g4)
fit2      <- lm(LAG3_log2tpm ~ tp53_bin + source_bin,           data = g4)
fit3      <- lm(LAG3_log2tpm ~ tp53_bin + source_bin + idh_bin, data = g4)

print_model <- function(fit, label) {
  s         <- summary(fit)
  coef_tp53 <- coef(s)["tp53_bin", ]
  cat(sprintf(
    "[%s] β_TP53=%.3f (95CI: %.3f–%.3f) | SE=%.3f | t=%.2f | p=%.4e | R²=%.3f\n",
    label,
    coef_tp53["Estimate"],
    confint(fit)["tp53_bin", 1],
    confint(fit)["tp53_bin", 2],
    coef_tp53["Std. Error"],
    coef_tp53["t value"],
    coef_tp53["Pr(>|t|)"],
    s$r.squared
  ))
}

cat("=== 回帰モデル比較（アウトカム: LAG3_log2tpm）===\n\n")
print_model(fit1, "Model1: TP53 only")
print_model(fit2, "Model2: TP53 + source")
print_model(fit3, "Model3: TP53 + source + IDH")

# Model3詳細
cat("\n=== Model3 全係数 ===\n")
print(summary(fit3)$coefficients)

cat("\n=== Model3 95%信頼区間 ===\n")
print(confint(fit3))

# ============================================================
# IDH WT限定でも確認（Model2相当: TP53 + source）
# ============================================================
cat("\n=== IDH wildtype限定（n=410）: TP53 + source ===\n")
g4_idhwt    <- g4 %>% filter(idh_status == "wildtype")
fit_idhwt   <- lm(LAG3_log2tpm ~ tp53_bin + source_bin, data = g4_idhwt)
print_model(fit_idhwt, "IDH_WT: TP53 + source")
cat("全係数:\n")
print(summary(fit_idhwt)$coefficients)
cat("95%CI:\n")
print(confint(fit_idhwt))

# ============================================================
# 残差正規性チェック（簡易）
# ============================================================
cat("\n=== Model3 残差確認 ===\n")
res3 <- residuals(fit3)
cat(sprintf(
  "残差: mean=%.4f, SD=%.3f, skewness=%.3f\n",
  mean(res3), sd(res3),
  mean((res3 - mean(res3))^3) / sd(res3)^3
))
sw <- shapiro.test(sample(res3, min(length(res3), 5000)))
cat(sprintf("Shapiro-Wilk（n=%d）: W=%.4f, p=%.4e\n", length(res3), sw$statistic, sw$p.value))
cat("※ n=442の線形回帰は中心極限定理により正規性への厳密な依存は低い\n")

# ============================================================
# 保存1: 既存形式（TP53行のみ）← 変更なし・論文照合済み数値を保護
# ============================================================
extract_model <- function(fit, model_label, subset_label = "Grade4_all") {
  s   <- summary(fit)
  ci  <- confint(fit)
  cf  <- coef(s)["tp53_bin", ]
  tibble(
    model         = model_label,
    subset        = subset_label,
    n             = nrow(fit$model),
    beta_tp53     = cf["Estimate"],
    ci_lower      = ci["tp53_bin", 1],
    ci_upper      = ci["tp53_bin", 2],
    se            = cf["Std. Error"],
    t_value       = cf["t value"],
    p_value       = cf["Pr(>|t|)"],
    r_squared     = s$r.squared,
    adj_r_squared = s$adj.r.squared
  )
}

results <- bind_rows(
  extract_model(fit1,     "TP53_only",       "Grade4_all"),
  extract_model(fit2,     "TP53+source",     "Grade4_all"),
  extract_model(fit3,     "TP53+source+IDH", "Grade4_all"),
  extract_model(fit_idhwt,"TP53+source",     "Grade4_IDH_wildtype")
)

write_csv(results, file.path(OUT_DIR, "step09b_regression_results.csv"))
cat("\n✅ 保存（既存形式）:", file.path(OUT_DIR, "step09b_regression_results.csv"), "\n")

# ============================================================
# 保存2: 全係数出力（新規）← Supplementary Table S1 用
# ============================================================
# term列の表示名を論文の変数定義に合わせてリネーム
TERM_LABELS <- c(
  "(Intercept)" = "Intercept",
  "tp53_bin"    = "TP53 (mutant vs WT)",
  "source_bin"  = "Source (TCGA vs CPTAC/HCMI)",
  "idh_bin"     = "IDH (mutant vs WT)"
)

extract_model_full <- function(fit, model_label, subset_label = "Grade4_all") {
  s  <- summary(fit)
  ci <- confint(fit)
  cf <- coef(s)
  
  tibble(
    model         = model_label,
    subset        = subset_label,
    n             = nrow(fit$model),
    term          = rownames(cf),
    term_label    = TERM_LABELS[rownames(cf)],
    beta          = cf[, "Estimate"],
    ci_lower      = ci[rownames(cf), 1],
    ci_upper      = ci[rownames(cf), 2],
    se            = cf[, "Std. Error"],
    t_value       = cf[, "t value"],
    p_value       = cf[, "Pr(>|t|)"],
    r_squared     = s$r.squared,
    adj_r_squared = s$adj.r.squared
  )
}

results_full <- bind_rows(
  extract_model_full(fit1,      "TP53_only",       "Grade4_all"),
  extract_model_full(fit2,      "TP53+source",     "Grade4_all"),
  extract_model_full(fit3,      "TP53+source+IDH", "Grade4_all"),
  extract_model_full(fit_idhwt, "TP53+source",     "Grade4_IDH_wildtype")
)

write_csv(results_full, file.path(OUT_DIR, "step09b_regression_results_full.csv"))
cat("✅ 保存（全係数）:", file.path(OUT_DIR, "step09b_regression_results_full.csv"), "\n")

# ============================================================
# コンソール確認: Model3（M0）全係数
# ============================================================
cat("\n=== Supplementary Table S1 用: M0全係数（TP53+source+IDH, Grade4_all）===\n")
results_full %>%
  filter(model == "TP53+source+IDH", subset == "Grade4_all") %>%
  select(term_label, beta, ci_lower, ci_upper, p_value, r_squared) %>%
  print()

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("【step09b 完了サマリ】\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\n■ TP53係数（既存ファイル確認用）\n")
print(results %>% select(model, subset, n, beta_tp53, ci_lower, ci_upper, p_value, r_squared))
cat("\n■ 全係数（新規ファイル確認用）\n")
print(results_full %>% select(model, subset, term_label, beta, ci_lower, ci_upper, p_value))
