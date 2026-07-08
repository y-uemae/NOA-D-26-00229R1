# ============================================================
# step05c_01_glass_mapping.R
# GLASS サンプルマッピング作成
# 出力先: results/TP53/20260221/05c_glass/
# ============================================================

library(tidyverse)
library(data.table)

# --- パス設定 ---
GLASS_DIR <- here::here("data", "raw", "external_validation", "difg_glass")
OUT_DIR   <- here::here("results", "TP53", "20260221", "05c_glass")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. ファイル読み込み
# ============================================================
cs <- read.table(file.path(GLASS_DIR, "data_clinical_sample.txt"),
                 sep = "\t", header = TRUE, comment.char = "#",
                 fill = TRUE, stringsAsFactors = FALSE, quote = "")

cp <- read.table(file.path(GLASS_DIR, "data_clinical_patient.txt"),
                 sep = "\t", header = TRUE, comment.char = "#",
                 fill = TRUE, stringsAsFactors = FALSE, quote = "")

# mutationはTumor_Sample_Barcode列のみ取得（152万行対策）
mut_ids_raw <- fread(file.path(GLASS_DIR, "data_mutations.txt"),
                     sep = "\t", select = "Tumor_Sample_Barcode",
                     skip = 0) %>%
  pull(Tumor_Sample_Barcode) %>%
  unique()

# RNAヘッダーのみ（列名取得）
rna_header <- read.table(file.path(GLASS_DIR, "data_mrna_seq_tpm.txt"),
                         sep = "\t", header = TRUE, nrows = 0,
                         comment.char = "#", stringsAsFactors = FALSE, quote = "")
rna_cols_dot  <- setdiff(colnames(rna_header), "Hugo_Symbol")
rna_cols_dash <- gsub("\\.", "-", rna_cols_dot)

cat("RNA列数:", length(rna_cols_dot), "\n")
cat("mutation ユニーク数:", length(mut_ids_raw), "\n")
cat("RNA→dash変換例:", paste(head(rna_cols_dash, 3), collapse = ", "), "\n")

n_rna_cs  <- sum(rna_cols_dash %in% cs$SAMPLE_ID)
n_rna_mut <- sum(rna_cols_dash %in% mut_ids_raw)
cat("RNA(dash) × SAMPLE_ID 一致:", n_rna_cs, "\n")
cat("RNA(dash) × mutation   一致:", n_rna_mut, "\n")

# ============================================================
# 2. RNA列名対応テーブル（dot ↔ dash）
# ============================================================
rna_id_map <- tibble(
  rna_col_dot  = rna_cols_dot,
  rna_col_dash = rna_cols_dash
)

# ============================================================
# 3. clinical_sample ベースでmapping構築
# ============================================================
mapping <- cs %>%
  transmute(
    patient_id   = PATIENT_ID,
    sample_id    = SAMPLE_ID,
    sample_type  = SAMPLE_TYPE,
    aliquot_type = ALIQUOT_ANALYSIS_TYPE,
    tumor_grade  = TUMOR_GRADE,
    idh_status   = IDH_STATUS,
    idh_codel    = IDH_CODEL_STATUS,
    dna_barcode  = DNA_ALIQUOT_BARCODE,
    rna_barcode  = RNA_ALIQUOT_BARCODE,
    has_rna_type = grepl("RNA", ALIQUOT_ANALYSIS_TYPE),
    has_rna_expr = SAMPLE_ID %in% rna_cols_dash,
    has_mutation = SAMPLE_ID %in% mut_ids_raw
  ) %>%
  mutate(has_rna = has_rna_expr)

cat("\n=== ALIQUOT_ANALYSIS_TYPE 集計 ===\n")
print(table(mapping$aliquot_type, useNA = "ifany"))

cat("\n=== has_rna_type vs has_rna_expr ===\n")
print(table(rna_type = mapping$has_rna_type, rna_expr = mapping$has_rna_expr, useNA = "ifany"))

cat("\n=== has_mutation 集計 ===\n")
print(table(mapping$has_mutation, useNA = "ifany"))

# ============================================================
# 4. Grade4判定
# ============================================================
cat("\nTUMOR_GRADE ユニーク値:", paste(unique(mapping$tumor_grade), collapse = ", "), "\n")

mapping <- mapping %>%
  mutate(
    is_grade4 = case_when(
      tumor_grade == "IV"                                        ~ TRUE,
      tumor_grade %in% c("II", "III")                           ~ FALSE,
      grepl("4|IV|GBM|Glioblastoma", tumor_grade,
            ignore.case = TRUE)                                  ~ TRUE,
      grepl("2|3|II|III|LGG", tumor_grade,
            ignore.case = TRUE)                                  ~ FALSE,
      TRUE                                                       ~ NA
    )
  )

cat("\nis_grade4 集計:\n")
print(table(mapping$is_grade4, useNA = "ifany"))

# ============================================================
# 5. SAMPLE_TYPE確認
# ============================================================
cat("\nSAMPLE_TYPE ユニーク値:", paste(unique(mapping$sample_type), collapse = ", "), "\n")
cat("\nSAMPLE_TYPE × is_grade4 クロス集計:\n")
print(table(sample_type = mapping$sample_type, grade4 = mapping$is_grade4, useNA = "ifany"))

# ============================================================
# 6. include_flag 付与
# ============================================================
mapping <- mapping %>%
  mutate(
    is_primary   = sample_type == "Tumor Primary",
    include_flag = has_rna & has_mutation & is_grade4 == TRUE & is_primary,
    exclude_reason = case_when(
      include_flag                                               ~ "included_ok",
      !is_primary & is_grade4 == TRUE & has_rna & has_mutation  ~ "not_primary",
      is.na(is_grade4) | is_grade4 == FALSE                     ~ "not_grade4",
      !has_rna                                                   ~ "no_rna",
      !has_mutation                                              ~ "no_mutation",
      TRUE                                                       ~ "other"
    )
  )

cat("\n=== exclude_reason 集計（is_grade4==TRUE のみ） ===\n")
print(table(mapping %>% filter(is_grade4 == TRUE) %>% pull(exclude_reason), useNA = "ifany"))

cat("\n=== 最終 include_flag=TRUE 件数:", sum(mapping$include_flag, na.rm = TRUE), "===\n")
cat("（引継書参照値: 95例）\n")

mapping <- mapping %>%
  mutate(
    seq_type = case_when(
      grepl("WXS", aliquot_type) & grepl("RNA", aliquot_type) ~ "WXS",
      grepl("WGS", aliquot_type) & grepl("RNA", aliquot_type) ~ "WGS",
      TRUE                                                     ~ NA_character_
    )
  )

cat("\n=== include_flag=TRUE の seq_type 内訳 ===\n")
print(table(mapping %>% filter(include_flag) %>% pull(seq_type), useNA = "ifany"))
cat("（引継書参照値: WXS=79, WGS=16）\n")

# ============================================================
# 7. RNA列名対応をmappingに追加
# ============================================================
mapping <- mapping %>%
  left_join(rna_id_map, by = c("sample_id" = "rna_col_dash"))

# ============================================================
# 8. 保存
# ============================================================
write.csv(mapping, file.path(OUT_DIR, "step05c_01_glass_mapping.csv"),
          row.names = FALSE)

included <- mapping %>% filter(include_flag)
write.csv(included, file.path(OUT_DIR, "step05c_01_glass_mapping_included.csv"),
          row.names = FALSE)

cat("\n✅ 保存完了\n")
cat("全mapping:", nrow(mapping), "行\n")
cat("included: ", nrow(included), "行\n")
