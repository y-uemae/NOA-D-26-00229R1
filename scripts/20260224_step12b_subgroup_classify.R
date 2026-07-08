# =============================================================================
# 20260224_step12b_subgroup_classify.R（修正版）
# Step 12b: TP53 4群分類確定
#
# 方針変更:
#   mutation_raw_variantsはTCGAのみ → CPTAC結合不可
#   → tp53_mutation_table.csvをベースに再分類（全source対応）
#   結合キー: cohort$wxs_sample_id ↔ tp53_table$wxs_sample_id
#
# 分類ルール（事前固定）:
#   優先1: WT          … tp53_status == "wildtype"
#   優先2: Hotspot     … Missense & tp53_HGVSpにR175/R248/R273/R282
#   優先3: Truncating  … Hotspot以外のNonsense/FrameShift/Splice/Start/Nonstop
#   優先4: Other_missense … 上記以外
#
# 入力:
#   08_final_cohort/final_cohort.csv
#   07_mutations/tp53_mutation_table.csv
#   07_mutations/mutation_raw_variants.csv  （TCGA Hotspot詳細確認用）
# 出力:
#   12_subgroup/step12b_subgroup_classified.csv  ★Step12c入力
#   12_subgroup/step12b_classification_log.csv
#   12_subgroup/step12b_crosstab.csv
# =============================================================================

library(tidyverse)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR  <- here::here("results", "TP53", "20260221")
COH_CSV   <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
TP53_CSV  <- file.path(BASE_DIR, "07_mutations/tp53_mutation_table.csv")
MUT_CSV   <- file.path(BASE_DIR, "07_mutations/mutation_raw_variants.csv")
OUT_DIR   <- file.path(BASE_DIR, "12_subgroup")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

HOTSPOT_PATTERN    <- "R175|R248|R273|R282"
TRUNCATING_CLASSES <- c(
  "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"
)

# ── 1. データ読み込み ─────────────────────────────────────────────────────────
cohort <- read_csv(COH_CSV, show_col_types = FALSE) %>%
  filter(include_flag == TRUE, grade == "Grade4")

tp53_table <- read_csv(TP53_CSV, show_col_types = FALSE)

mut_raw <- read_csv(MUT_CSV, show_col_types = FALSE) %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(case_barcode_12 = str_sub(Tumor_Sample_Barcode, 1, 12))

cat("=== 基本確認 ===\n")
cat("Grade4コホート n =", nrow(cohort), "\n")
cat("tp53_mutation_table 行数 =", nrow(tp53_table), "\n\n")

# ── 2. tp53_HGVSpの確認 ───────────────────────────────────────────────────────
cat("=== tp53_HGVSp サンプル（mutant上位10件）===\n")
tp53_table %>%
  filter(tp53_status == "mutant") %>%
  select(source, wxs_sample_id, tp53_subgroup, tp53_classifications, tp53_HGVSp) %>%
  head(10) %>% print()

# ── 3. 4群分類関数 ────────────────────────────────────────────────────────────
classify_tp53_class4 <- function(tp53_status, tp53_HGVSp,
                                 tp53_classifications, tp53_subgroup) {
  # 優先1: WT
  if (tp53_status == "wildtype") return(list(class = "WT", method = "tp53_status"))
  
  hgvsp <- tp53_HGVSp %||% ""
  cls   <- tp53_classifications %||% ""
  sub   <- tp53_subgroup %||% ""
  
  # 優先2: Hotspot（Missense & R175/R248/R273/R282）
  # tp53_HGVSp は "|" 区切りで複数あり得る
  hgvsp_parts <- str_split(hgvsp, "\\|")[[1]]
  is_hotspot <- any(
    str_detect(hgvsp_parts, HOTSPOT_PATTERN) &
      str_detect(cls, "Missense_Mutation"),
    na.rm = TRUE
  )
  # tp53_subgroupにR175/R248/R273/R282が直接入っているケースも拾う
  if (!is_hotspot && str_detect(sub, HOTSPOT_PATTERN) &&
      str_detect(cls, "Missense_Mutation")) {
    is_hotspot <- TRUE
  }
  if (is_hotspot) return(list(class = "Hotspot", method = "HGVSp"))
  
  # 優先3: Truncating（Variant_Classification文字列から判定）
  is_truncating <- any(str_detect(cls, paste(TRUNCATING_CLASSES, collapse = "|")))
  if (is_truncating) return(list(class = "Truncating", method = "tp53_classifications"))
  
  # 優先4: Other_missense
  return(list(class = "Other_missense", method = "tp53_classifications"))
}

# ── 4. Grade4コホートにtp53_tableを結合 ──────────────────────────────────────
# 結合キー: cohort$wxs_sample_id ↔ tp53_table$wxs_sample_id
cohort_with_tp53 <- cohort %>%
  select(wxs_sample_id, case_barcode, source, tp53_status, idh_status,
         LAG3_log2tpm) %>%
  left_join(
    tp53_table %>%
      select(wxs_sample_id, tp53_HGVSp, tp53_classifications, tp53_subgroup_orig = tp53_subgroup),
    by = "wxs_sample_id"
  )

# 結合確認
n_join_fail <- sum(is.na(cohort_with_tp53$tp53_HGVSp) &
                     cohort_with_tp53$tp53_status == "mutant")
cat(sprintf("\n=== 結合確認: mutantでtp53_HGVSp欠損 = %d件 ===\n", n_join_fail))

# ── 5. 分類実行 ───────────────────────────────────────────────────────────────
cohort_classified <- cohort_with_tp53 %>%
  mutate(
    result = pmap(
      list(tp53_status, tp53_HGVSp, tp53_classifications, tp53_subgroup_orig),
      classify_tp53_class4
    ),
    tp53_class4     = map_chr(result, "class"),
    classify_method = map_chr(result, "method"),
    tp53_class4     = factor(tp53_class4,
                             levels = c("WT", "Hotspot", "Truncating", "Other_missense"))
  ) %>%
  select(-result)

# ── 6. NAチェック ─────────────────────────────────────────────────────────────
n_na <- sum(is.na(cohort_classified$tp53_class4))
cat(sprintf("\n=== NAチェック: tp53_class4のNA = %d件 ===\n", n_na))

if (n_na > 0) {
  cat("⚠️  NA残存:\n")
  cohort_classified %>%
    filter(is.na(tp53_class4)) %>%
    select(wxs_sample_id, case_barcode, source, tp53_status) %>%
    print()
  stop("NAが残っています")
} else {
  cat("✅ NA ゼロ確認\n")
}

# ── 7. 分類結果確認 ───────────────────────────────────────────────────────────
cat("\n=== 新分類クロス表（source別）===\n")
crosstab <- cohort_classified %>%
  count(tp53_class4, source) %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))
print(crosstab)

cat("\n=== Mut合計の一致確認 ===\n")
tcga_mut  <- sum(cohort_classified$source == "TCGA" &
                   cohort_classified$tp53_status == "mutant")
cptac_mut <- sum(cohort_classified$source == "CPTAC_HCMI" &
                   cohort_classified$tp53_status == "mutant")
cat(sprintf("TCGA mutant:       %d（期待値79）%s\n",
            tcga_mut,  if (tcga_mut  == 79) "✅" else "⚠️"))
cat(sprintf("CPTAC_HCMI mutant: %d（期待値68）%s\n",
            cptac_mut, if (cptac_mut == 68) "✅" else "⚠️"))

cat(sprintf("\nCPTAC_HCMI Hotspot: %d（HGVSp情報による・0の場合は保守的推定）\n",
            sum(cohort_classified$source == "CPTAC_HCMI" &
                  cohort_classified$tp53_class4 == "Hotspot", na.rm = TRUE)))

cat("\n=== 分類方法の内訳 ===\n")
cohort_classified %>%
  filter(tp53_status == "mutant") %>%
  count(source, classify_method, tp53_class4) %>%
  print()

# ── 8. G245S の確認（Other_missenseに入っていること）────────────────────────
g245 <- cohort_classified %>%
  filter(str_detect(tp53_HGVSp %||% "", "G245"))
cat(sprintf("\n=== G245S確認: %d件 → tp53_class4 = %s（期待: Other_missense）===\n",
            nrow(g245),
            if (nrow(g245) > 0) as.character(g245$tp53_class4[1]) else "該当なし"))

# ── 9. CSV出力 ────────────────────────────────────────────────────────────────
write_csv(cohort_classified,
          file.path(OUT_DIR, "step12b_subgroup_classified.csv"))

write_csv(
  cohort_classified %>%
    select(wxs_sample_id, case_barcode, source, tp53_status,
           tp53_class4, classify_method,
           tp53_HGVSp, tp53_classifications, tp53_subgroup_orig),
  file.path(OUT_DIR, "step12b_classification_log.csv")
)

write_csv(crosstab, file.path(OUT_DIR, "step12b_crosstab.csv"))

cat("\n=== Step 12b 完了 ===\n")
cat("  step12b_subgroup_classified.csv  ★Step12c入力\n")
cat("  step12b_classification_log.csv\n")
cat("  step12b_crosstab.csv\n")
