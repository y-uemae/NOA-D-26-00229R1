# =============================================================================
# step30c_screening_stage1_stage2.R 【v1.0】
# 目的: 2段階スクリーニング
#   Stage1: LAG3 ~ score_k + source + IDH（相関探索）
#   Stage2: Δβ計算（TP53のβがスコア追加でどれだけ変わるか）
#     Base:  LAG3 ~ TP53 + source + IDH + Tcell + APM + IFNg + TumorPurity
#     Test:  LAG3 ~ TP53 + source + IDH + Tcell + APM + IFNg + TumorPurity + score_k
# 入力: 30_screening/analysis_dataset_grade4.rds
# 出力先: 30_screening/
#   stage1_results.rds / .csv
#   stage2_results.rds / .csv
# 作成日: 2026-02-27
# バージョン履歴:
#   v1.0 - 初版
# =============================================================================

library(tidyverse)

# =============================================================================
# 0. パス設定
# =============================================================================
base_dir <- here::here("results", "TP53", "20260221")
out_dir  <- file.path(base_dir, "30_screening")

dataset_path <- file.path(out_dir, "analysis_dataset_grade4.rds")

output_stage1_rds <- file.path(out_dir, "stage1_results.rds")
output_stage1_csv <- file.path(out_dir, "stage1_results.csv")
output_stage2_rds <- file.path(out_dir, "stage2_results.rds")
output_stage2_csv <- file.path(out_dir, "stage2_results.csv")

# =============================================================================
# 1. データ読み込み
# =============================================================================
cat("データ読み込み中...\n")
df <- readRDS(dataset_path)
cat("サンプル数:", nrow(df), "/ 列数:", ncol(df), "\n")

# スクリーニング対象スコア列を取得
score_cols <- names(df)[grepl("^HALLMARK_", names(df))]
score_cols <- c(score_cols, "TumorPurity")  # ESTIMATEも含める
cat("スクリーニング対象スコア数:", length(score_cols), "\n")

# 因子型への変換
df <- df %>%
  mutate(
    tp53_status = factor(tp53_status, levels = c("wildtype", "mutant")),
    source      = factor(source,      levels = c("TCGA", "CPTAC_HCMI")),
    idh_status  = factor(idh_status,  levels = c("wildtype", "mutant"))
  )

# =============================================================================
# 2. Stage 1: 相関スクリーニング
#    LAG3 ~ score_k + source + IDH
# =============================================================================
cat("\n========================================\n")
cat("Stage 1: 相関スクリーニング開始\n")
cat("========================================\n")

run_stage1 <- function(score_col, df) {
  fml <- as.formula(paste("LAG3 ~", score_col, "+ source + idh_status"))
  fit <- lm(fml, data = df)
  s   <- summary(fit)$coefficients
  # score_kの行を取得
  row_name <- score_col
  if (!row_name %in% rownames(s)) return(NULL)
  data.frame(
    score      = score_col,
    beta       = s[row_name, "Estimate"],
    se         = s[row_name, "Std. Error"],
    t_value    = s[row_name, "t value"],
    p_value    = s[row_name, "Pr(>|t|)"],
    r2_adj     = summary(fit)$adj.r.squared
  )
}

stage1_list <- map(score_cols, run_stage1, df = df)
stage1_df   <- bind_rows(stage1_list) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("Stage1完了。上位10件:\n")
print(stage1_df %>% select(score, beta, p_value, p_adj) %>% head(10),
      digits = 3)

saveRDS(stage1_df, output_stage1_rds)
write_csv(stage1_df, output_stage1_csv)
cat("Stage1保存完了:", output_stage1_csv, "\n")

# =============================================================================
# 3. Stage 2: Δβ計算
#    Base model: LAG3 ~ TP53 + source + IDH + Tcell + APM + IFNg + TumorPurity
#    Test model: Base + score_k（1つずつ追加）
# =============================================================================
cat("\n========================================\n")
cat("Stage 2: Δβ計算開始\n")
cat("========================================\n")

# --- Base model ---
base_fml <- LAG3 ~ tp53_status + source + idh_status + Tcell + APM + IFNg + TumorPurity
base_fit  <- lm(base_fml, data = df)
beta_base <- coef(summary(base_fit))["tp53_statusmutant", "Estimate"]
se_base   <- coef(summary(base_fit))["tp53_statusmutant", "Std. Error"]
p_base    <- coef(summary(base_fit))["tp53_statusmutant", "Pr(>|t|)"]
r2_base   <- summary(base_fit)$adj.r.squared

cat("Base model: β_TP53 =", round(beta_base, 4),
    "/ p =", formatC(p_base, format = "e", digits = 2),
    "/ adj.R² =", round(r2_base, 3), "\n")

# --- Test models（score_kを1つずつ追加）---
run_stage2 <- function(score_col, df, beta_base) {
  fml <- as.formula(paste(
    "LAG3 ~ tp53_status + source + idh_status + Tcell + APM + IFNg + TumorPurity +",
    score_col
  ))
  fit      <- lm(fml, data = df)
  s        <- coef(summary(fit))
  beta_new <- s["tp53_statusmutant", "Estimate"]
  se_new   <- s["tp53_statusmutant", "Std. Error"]
  p_new    <- s["tp53_statusmutant", "Pr(>|t|)"]
  r2_new   <- summary(fit)$adj.r.squared
  
  # score_k自体のβ・p
  row_score  <- score_col
  beta_score <- if (row_score %in% rownames(s)) s[row_score, "Estimate"] else NA
  p_score    <- if (row_score %in% rownames(s)) s[row_score, "Pr(>|t|)"] else NA
  
  delta_beta    <- beta_new - beta_base
  delta_beta_pct <- (delta_beta / abs(beta_base)) * 100
  
  data.frame(
    score          = score_col,
    beta_tp53_base = beta_base,
    beta_tp53_test = beta_new,
    delta_beta     = delta_beta,
    delta_beta_pct = delta_beta_pct,
    p_tp53_test    = p_new,
    beta_score     = beta_score,
    p_score        = p_score,
    r2_adj_test    = r2_new
  )
}

stage2_list <- map(score_cols, run_stage2, df = df, beta_base = beta_base)
stage2_df   <- bind_rows(stage2_list) %>%
  mutate(p_score_adj = p.adjust(p_score, method = "BH")) %>%
  arrange(desc(abs(delta_beta_pct)))

cat("Stage2完了。Δβ%上位10件:\n")
print(stage2_df %>%
        select(score, beta_tp53_base, beta_tp53_test, delta_beta_pct, p_score) %>%
        head(10),
      digits = 3)

cat("\n★ 最大Δβ%:", round(max(abs(stage2_df$delta_beta_pct)), 1), "%\n")
cat("★ 全スコアのΔβ%範囲:",
    round(min(stage2_df$delta_beta_pct), 1), "〜",
    round(max(stage2_df$delta_beta_pct), 1), "%\n")

saveRDS(stage2_df, output_stage2_rds)
write_csv(stage2_df, output_stage2_csv)
cat("Stage2保存完了:", output_stage2_csv, "\n")

# =============================================================================
# 4. 完了サマリー
# =============================================================================
cat("\n========================================\n")
cat("Step 30 完了\n")
cat("========================================\n")
cat("Stage1（相関）: ", nrow(stage1_df), "スコア →", output_stage1_csv, "\n")
cat("Stage2（Δβ）: ", nrow(stage2_df), "スコア →", output_stage2_csv, "\n")
cat("Base β_TP53:", round(beta_base, 4), "\n")
cat("最大Δβ%:", round(max(abs(stage2_df$delta_beta_pct)), 1), "%\n")
