# =============================================================================
# step29_ssgsea_hallmark.R 【v1.0】
# 目的: ssGSEA によるHallmark50遺伝子セットスコア計算（GDC + GLASS）
# 入力: expression_full_log2tpm_wide.csv（951×59428, pair_id=連番）
#       glass_expression_full_log2tpm_wide.csv（79×35434）
# 出力先: 29_ssgsea_hallmark/
#   ssgsea_hallmark_gdc.rds / .csv
#   ssgsea_hallmark_glass.rds / .csv
# 列構成: pair_id, HALLMARK_xxx（50列）
# 作成日: 2026-02-27
# バージョン履歴:
#   v1.0 - 初版
# =============================================================================

library(tidyverse)
library(GSVA)
library(msigdbr)

# =============================================================================
# 0. パス設定
# =============================================================================
base_dir <- here::here("results", "TP53", "20260221")
expr_dir <- file.path(base_dir, "27_expression_matrix")
out_dir  <- file.path(base_dir, "29_ssgsea_hallmark")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("フォルダ作成:", out_dir, "\n")
} else {
  cat("出力フォルダ確認済み:", out_dir, "\n")
}

input_gdc   <- file.path(expr_dir, "expression_full_log2tpm_wide.csv")
input_glass <- file.path(expr_dir, "glass_expression_full_log2tpm_wide.csv")

output_gdc_rds   <- file.path(out_dir, "ssgsea_hallmark_gdc.rds")
output_gdc_csv   <- file.path(out_dir, "ssgsea_hallmark_gdc.csv")
output_glass_rds <- file.path(out_dir, "ssgsea_hallmark_glass.rds")
output_glass_csv <- file.path(out_dir, "ssgsea_hallmark_glass.csv")

# =============================================================================
# 1. Hallmark遺伝子セット取得
# =============================================================================
cat("Hallmark遺伝子セット取得中...\n")
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat("遺伝子セット数:", length(hallmark_list), "\n")
cat("先頭3セットの遺伝子数:\n")
print(sapply(hallmark_list[1:3], length))

# =============================================================================
# 2. 関数定義: run_ssgsea()
# =============================================================================
run_ssgsea <- function(expr_wide_path, cohort_name, hallmark_list,
                       out_rds, out_csv) {
  cat("\n========================================\n")
  cat(cohort_name, "処理開始\n")
  cat("========================================\n")
  
  # --- 読み込み ---
  cat("発現行列読み込み中...\n")
  expr_wide <- read_csv(expr_wide_path, show_col_types = FALSE)
  cat("読み込み完了:", nrow(expr_wide), "サンプル×",
      ncol(expr_wide) - 1, "遺伝子\n")
  
  # --- 行列化: 行=遺伝子, 列=サンプル ---
  cat("行列変換中...\n")
  pair_ids <- as.character(expr_wide$pair_id)
  
  mat <- expr_wide %>%
    select(-pair_id) %>%
    as.matrix()
  rownames(mat) <- pair_ids
  mat <- t(mat)  # 行=遺伝子, 列=サンプル（pair_id）
  
  cat("行列サイズ:", nrow(mat), "遺伝子×", ncol(mat), "サンプル\n")
  
  # --- ssGSEA実行 ---
  # GSVA v2.x の新API: gsvaParam() でパラメータオブジェクトを作成
  cat("ssGSEA実行中（数分〜十数分かかります）...\n")
  param <- ssgseaParam(
    exprData  = mat,
    geneSets  = hallmark_list,
    normalize = TRUE   # スコアを[-1,1]に正規化
  )
  ssgsea_mat <- gsva(param, verbose = TRUE)
  
  cat("ssGSEA完了。出力サイズ:", nrow(ssgsea_mat), "セット×",
      ncol(ssgsea_mat), "サンプル\n")
  
  # --- 整形: wide形式（行=サンプル, 列=遺伝子セット）---
  result_df <- ssgsea_mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("pair_id")
  
  # GDCは pair_id を整数に変換
  if (all(grepl("^[0-9]+$", result_df$pair_id))) {
    result_df <- result_df %>% mutate(pair_id = as.integer(pair_id))
  }
  
  cat("サンプル数:", nrow(result_df), "\n")
  cat("列数（pair_id + 遺伝子セット）:", ncol(result_df), "\n")
  cat("スコアの範囲確認（HALLMARK_GLYCOLYSIS）:\n")
  if ("HALLMARK_GLYCOLYSIS" %in% names(result_df)) {
    print(summary(result_df$HALLMARK_GLYCOLYSIS))
  }
  
  # --- 上書き保存 ---
  saveRDS(result_df, out_rds)
  write_csv(result_df, out_csv)
  cat("保存完了（RDS）:", out_rds, "\n")
  cat("保存完了（CSV）:", out_csv, "\n")
  
  return(result_df)
}

# =============================================================================
# 3. GDC実行
# =============================================================================
ssgsea_gdc <- run_ssgsea(
  expr_wide_path = input_gdc,
  cohort_name    = "GDC",
  hallmark_list  = hallmark_list,
  out_rds        = output_gdc_rds,
  out_csv        = output_gdc_csv
)

# =============================================================================
# 4. GLASS実行
# =============================================================================
ssgsea_glass <- run_ssgsea(
  expr_wide_path = input_glass,
  cohort_name    = "GLASS",
  hallmark_list  = hallmark_list,
  out_rds        = output_glass_rds,
  out_csv        = output_glass_csv
)

# =============================================================================
# 5. 完了サマリー
# =============================================================================
cat("\n========================================\n")
cat("Step 29 完了\n")
cat("========================================\n")
cat("GDC:   ", nrow(ssgsea_gdc),   "サンプル×", ncol(ssgsea_gdc)-1,   "セット →", output_gdc_csv,   "\n")
cat("GLASS: ", nrow(ssgsea_glass), "サンプル×", ncol(ssgsea_glass)-1, "セット →", output_glass_csv, "\n")
cat("列: pair_id, HALLMARK_xxx（50列）\n")
