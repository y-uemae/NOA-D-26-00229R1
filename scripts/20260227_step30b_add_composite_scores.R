# =============================================================================
# step30b_add_composite_scores.R 【v1.0】
# 目的: analysis_dataset_grade4にTcell/APM/IFNgの複合スコアを追加
# 入力: 30_screening/analysis_dataset_grade4.rds
# 出力: 30_screening/analysis_dataset_grade4.rds / .csv（上書き）
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

# =============================================================================
# 1. データ読み込み
# =============================================================================
cat("データ読み込み中...\n")
df <- readRDS(dataset_path)
cat("サンプル数:", nrow(df), "/ 列数:", ncol(df), "\n")

# =============================================================================
# 2. 複合スコアの遺伝子構成を定義
# Step09bで使用した構成遺伝子（引継書より）
# =============================================================================

# APM（Antigen Presentation Machinery）
apm_genes <- c("B2M_log2tpm", "TAP1_log2tpm", "TAP2_log2tpm",
               "TAPBP_log2tpm", "HLA-A_log2tpm", "HLA-B_log2tpm",
               "HLA-C_log2tpm", "NLRC5_log2tpm")

# IFNg シグネチャー
ifng_genes <- c("STAT1_log2tpm", "IRF1_log2tpm", "IRF9_log2tpm",
                "CXCL9_log2tpm", "CXCL10_log2tpm", "CXCL11_log2tpm",
                "GBP1_log2tpm", "GBP2_log2tpm", "GBP4_log2tpm",
                "GBP5_log2tpm", "IDO1_log2tpm")

# T細胞スコア
tcell_genes <- c("CD3D_log2tpm", "CD3E_log2tpm", "CD3G_log2tpm",
                 "CD8A_log2tpm", "CD8B_log2tpm",
                 "GZMA_log2tpm", "GZMB_log2tpm", "PRF1_log2tpm")

# =============================================================================
# 3. 遺伝子の存在確認
# =============================================================================
check_genes <- function(genes, score_name, df) {
  missing <- setdiff(genes, names(df))
  if (length(missing) > 0) {
    cat("⚠", score_name, "- 欠損列:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("✔", score_name, "- 全", length(genes), "遺伝子OK\n")
  }
  return(length(missing) == 0)
}

cat("\n遺伝子列の存在確認:\n")
ok_apm   <- check_genes(apm_genes,   "APM",   df)
ok_ifng  <- check_genes(ifng_genes,  "IFNg",  df)
ok_tcell <- check_genes(tcell_genes, "Tcell", df)

# =============================================================================
# 4. 複合スコアを行平均で計算
# =============================================================================
cat("\n複合スコア計算中...\n")

if (ok_apm) {
  df <- df %>%
    mutate(APM = rowMeans(select(., all_of(apm_genes)), na.rm = TRUE))
  cat("APM summary:\n")
  print(summary(df$APM))
}

if (ok_ifng) {
  df <- df %>%
    mutate(IFNg = rowMeans(select(., all_of(ifng_genes)), na.rm = TRUE))
  cat("IFNg summary:\n")
  print(summary(df$IFNg))
}

if (ok_tcell) {
  df <- df %>%
    mutate(Tcell = rowMeans(select(., all_of(tcell_genes)), na.rm = TRUE))
  cat("Tcell summary:\n")
  print(summary(df$Tcell))
}

# LAG3列名を統一
df <- df %>% rename(LAG3 = LAG3_log2tpm)
cat("LAG3 summary:\n")
print(summary(df$LAG3))

# =============================================================================
# 5. 必要列の最終確認
# =============================================================================
cat("\n最終確認:\n")
required_cols <- c("LAG3", "tp53_status", "source", "idh_status",
                   "Tcell", "APM", "IFNg", "TumorPurity")
missing_cols  <- setdiff(required_cols, names(df))

if (length(missing_cols) > 0) {
  cat("⚠ まだ欠損列があります:", paste(missing_cols, collapse = ", "), "\n")
} else {
  cat("✔ 必要列すべて揃いました\n")
  cat("\n回帰変数のsummary:\n")
  print(summary(df[, required_cols]))
}

cat("\nサンプル数:", nrow(df), "/ 列数:", ncol(df), "\n")

# =============================================================================
# 6. 上書き保存
# =============================================================================
saveRDS(df, dataset_path)
write_csv(df, sub("\\.rds$", ".csv", dataset_path))
cat("上書き保存完了:", dataset_path, "\n")
