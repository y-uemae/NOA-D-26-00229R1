# ============================================================
# step05c_02_glass_final_cohort.R
# GLASS 28遺伝子抽出 + TP53/IDH変異抽出 + final_cohort構築
# mutation coverage必須版
# 出力先: results/TP53/20260221/05c_glass/
# ============================================================

library(tidyverse)
library(data.table)

# --- パス設定 ---
GLASS_DIR <- here::here("data", "raw", "external_validation", "difg_glass")
OUT_DIR   <- here::here("results", "TP53", "20260221", "05c_glass")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 0. 定数定義
# ============================================================
GENES_28 <- c(
  "B2M", "TAP1", "TAP2", "TAPBP", "HLA-A", "HLA-B", "HLA-C", "NLRC5",
  "STAT1", "IRF1", "IRF9",
  "CXCL9", "CXCL10", "CXCL11",
  "GBP1", "GBP2", "GBP4", "GBP5",
  "IDO1",
  "CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
  "GZMA", "GZMB", "PRF1",
  "LAG3"
)

TP53_NONSILENT <- c(
  "Missense_Mutation", "Nonsense_Mutation",
  "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
  "Translation_Start_Site", "Nonstop_Mutation"
)

# ============================================================
# 1. mapping読み込み（101例）
# ============================================================
mapping <- read.csv(file.path(OUT_DIR, "step05c_01_glass_mapping_included.csv"),
                    stringsAsFactors = FALSE)
stopifnot(nrow(mapping) == 101)
cat("✅ mapping: 101行\n")

# ============================================================
# 2. mutation全件読み込み → coverage確認
# ============================================================
cat("\nmutation全件読み込み中（coverage確認用）...\n")
mut_all <- fread(file.path(GLASS_DIR, "data_mutations.txt"), sep = "\t") %>%
  as_tibble()

# ssm2_pass_call == "t" のみ
mut_pass_all <- mut_all %>%
  filter(tolower(as.character(ssm2_pass_call)) == "t")
cat("pass行数（全遺伝子）:", nrow(mut_pass_all), "\n")

# 101サンプルのうちmutationファイルに存在するサンプル（全遺伝子）
mut_covered_samples <- mut_pass_all %>%
  filter(Tumor_Sample_Barcode %in% mapping$sample_id) %>%
  pull(Tumor_Sample_Barcode) %>%
  unique()
cat("101例中 mutation coverage あり:", length(mut_covered_samples), "例\n")
cat("mutation coverage なし:", 101 - length(mut_covered_samples), "例\n")

# coverage フラグをmappingに付与
mapping <- mapping %>%
  mutate(
    mutation_coverage = sample_id %in% mut_covered_samples
  )

cat("\nmutation_coverage 集計:\n")
print(table(mapping$mutation_coverage, useNA = "ifany"))

# ============================================================
# 3. TP53 / IDH 変異抽出（coverage ありサンプルのみ）
# ============================================================
mut_target <- mut_pass_all %>%
  filter(Hugo_Symbol %in% c("TP53", "IDH1", "IDH2"),
         Tumor_Sample_Barcode %in% mut_covered_samples)
cat("\nTP53/IDH1/IDH2（pass & covered）行数:", nrow(mut_target), "\n")

# --- TP53 ---
hotspot_tp53 <- c("p.R273C","p.R175H","p.R248Q","p.R273H",
                  "p.R282W","p.R248W","p.H179R","p.G245S")

tp53_mut <- mut_target %>%
  filter(Hugo_Symbol == "TP53",
         Variant_Classification %in% TP53_NONSILENT) %>%
  group_by(Tumor_Sample_Barcode) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    sample_id     = Tumor_Sample_Barcode,
    tp53_status   = "Mut",
    tp53_HGVSp    = HGVSp_Short,
    tp53_variant  = Variant_Classification,
    tp53_subgroup = case_when(
      HGVSp_Short %in% hotspot_tp53                         ~ HGVSp_Short,
      Variant_Classification == "Missense_Mutation"         ~ "Missense_other",
      Variant_Classification %in% c("Nonsense_Mutation",
                                    "Frame_Shift_Del","Frame_Shift_Ins",
                                    "Translation_Start_Site","Nonstop_Mutation")        ~ "Truncating",
      Variant_Classification == "Splice_Site"               ~ "Splice_Site",
      TRUE                                                  ~ "Other_nonsilent"
    )
  )
cat("TP53 Mut サンプル数:", nrow(tp53_mut), "\n")

# --- IDH1 / IDH2 ---
idh_mut <- mut_target %>%
  filter(Hugo_Symbol %in% c("IDH1", "IDH2"),
         Variant_Classification %in% TP53_NONSILENT) %>%
  group_by(Tumor_Sample_Barcode) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    sample_id    = Tumor_Sample_Barcode,
    idh_mut_gene = Hugo_Symbol,
    idh_HGVSp    = HGVSp_Short,
    idh_subgroup = case_when(
      HGVSp_Short == "p.R132H" & Hugo_Symbol == "IDH1" ~ "IDH1_R132H",
      Hugo_Symbol == "IDH1"                            ~ "IDH1_other",
      HGVSp_Short == "p.R172K" & Hugo_Symbol == "IDH2" ~ "IDH2_R172K",
      Hugo_Symbol == "IDH2"                            ~ "IDH2_other",
      TRUE                                             ~ NA_character_
    )
  )
cat("IDH Mut サンプル数:", nrow(idh_mut), "\n")

# ============================================================
# 4. RNA発現量抽出（28遺伝子 × 101サンプル）
# ============================================================
cat("\nRNAデータ読み込み中...\n")
rna_all <- read.table(file.path(GLASS_DIR, "data_mrna_seq_tpm.txt"),
                      sep = "\t", header = TRUE, comment.char = "#",
                      stringsAsFactors = FALSE, quote = "")

rna_28 <- rna_all %>%
  filter(Hugo_Symbol %in% GENES_28) %>%
  group_by(Hugo_Symbol) %>% slice(1) %>% ungroup()  # 重複は先頭採用

target_cols_dot <- mapping$rna_col_dot
stopifnot(length(target_cols_dot) == 101)

rna_long <- rna_28 %>%
  select(Hugo_Symbol, all_of(target_cols_dot)) %>%
  pivot_longer(-Hugo_Symbol, names_to = "rna_col_dot", values_to = "tpm") %>%
  mutate(log2tpm = log2(tpm + 1)) %>%
  left_join(mapping %>% select(sample_id, rna_col_dot), by = "rna_col_dot")

gene_expr <- rna_long %>%
  select(sample_id, Hugo_Symbol, tpm) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = tpm,
              names_glue = "{Hugo_Symbol}_tpm") %>%
  left_join(
    rna_long %>%
      select(sample_id, Hugo_Symbol, log2tpm) %>%
      pivot_wider(names_from = Hugo_Symbol, values_from = log2tpm,
                  names_glue = "{Hugo_Symbol}_log2tpm"),
    by = "sample_id"
  ) %>%
  select(all_of(intersect(
    c("sample_id", as.vector(rbind(paste0(GENES_28, "_tpm"),
                                   paste0(GENES_28, "_log2tpm")))),
    colnames(.)
  )))

# ============================================================
# 5. final_cohort構築（mutation coverage 必須）
# ============================================================
final_all <- mapping %>%
  select(patient_id, sample_id, sample_type, tumor_grade,
         seq_type, aliquot_type, mutation_coverage) %>%
  rename(
    case_barcode = patient_id,
    pair_id      = sample_id,
    grade        = tumor_grade
  ) %>%
  mutate(
    grade  = "Grade4",
    cohort = "GLASS",
    source = paste0("GLASS_", seq_type)
  ) %>%
  left_join(gene_expr,  by = c("pair_id" = "sample_id")) %>%
  left_join(tp53_mut %>% select(sample_id, tp53_status, tp53_subgroup,
                                tp53_HGVSp, tp53_variant),
            by = c("pair_id" = "sample_id")) %>%
  mutate(
    tp53_status   = case_when(
      tp53_status == "Mut"     ~ "Mut",
      mutation_coverage        ~ "WT",       # coverage あり & 変異なし → WT
      TRUE                     ~ "Unknown"   # coverage なし → Unknown
    ),
    tp53_subgroup = if_else(is.na(tp53_subgroup),
                            if_else(mutation_coverage, "WT", "Unknown"),
                            tp53_subgroup)
  ) %>%
  left_join(idh_mut, by = c("pair_id" = "sample_id")) %>%
  mutate(
    idh_status   = case_when(
      !is.na(idh_mut_gene) ~ "Mut",
      mutation_coverage    ~ "WT",
      TRUE                 ~ "Unknown"
    ),
    idh_subgroup = case_when(
      !is.na(idh_subgroup) ~ idh_subgroup,
      mutation_coverage    ~ "WT",
      TRUE                 ~ "Unknown"
    )
  ) %>%
  select(-idh_mut_gene) %>%
  mutate(
    include_flag   = mutation_coverage,
    exclude_reason = if_else(mutation_coverage, "included_ok",
                             "no_mutation_coverage"),
    lag3_status    = if_else(!is.na(LAG3_log2tpm), "ok", "missing"),
    maf_status     = if_else(mutation_coverage, "ok", "no_coverage")
  )

# ============================================================
# 6. サマリー
# ============================================================
cat("\n=== 全101例 サマリー ===\n")
cat("exclude_reason:\n")
print(table(final_all$exclude_reason, useNA = "ifany"))

final <- final_all %>% filter(include_flag)
cat("\n=== 主解析セット（mutation coverage あり）===\n")
cat("行数:", nrow(final), "（参照値: 95例）\n")

cat("\nTP53 WT/Mut:\n")
print(table(final$tp53_status, useNA = "ifany"))
cat("（参照値: WT=66, Mut=29）\n")

cat("\nIDH WT/Mut:\n")
print(table(final$idh_status, useNA = "ifany"))

cat("\nTP53 × IDH クロス集計:\n")
print(table(TP53 = final$tp53_status, IDH = final$idh_status, useNA = "ifany"))

cat("\nsource（WXS/WGS）:\n")
print(table(final$source, useNA = "ifany"))
cat("（参照値: WXS=79, WGS=16）\n")

cat("\nLAG3 log2(TPM+1):\n")
lag3 <- final$LAG3_log2tpm
cat("n=", sum(!is.na(lag3)),
    " median=", round(median(lag3, na.rm = TRUE), 3),
    " max=",    round(max(lag3, na.rm = TRUE), 3), "\n")

# ============================================================
# 7. 保存
# ============================================================
# 全101例（監査用）
write.csv(final_all, file.path(OUT_DIR, "glass_final_cohort_all101.csv"),
          row.names = FALSE)

# 主解析セット（coverage あり）
write.csv(final, file.path(OUT_DIR, "glass_final_cohort.csv"),
          row.names = FALSE)

cat("\n✅ 保存完了\n")
cat("監査用（全101例）:", file.path(OUT_DIR, "glass_final_cohort_all101.csv"), "\n")
cat("主解析セット:      ", file.path(OUT_DIR, "glass_final_cohort.csv"), "\n")
cat("主解析件数:", nrow(final), "\n")
