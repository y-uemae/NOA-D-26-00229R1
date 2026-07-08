# ============================================================
# step09_main.R  GDC 統計解析
# 修正点: include_flag==TRUE, tp53_status="mutant"/"wildtype"
#         idh_status="mutant"/"wildtype", LAG3_log2tpm列名確認済
# ============================================================

library(tidyverse)

BASE_DIR <- here::here("results", "TP53", "20260221")
IN_FILE  <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
OUT_DIR  <- file.path(BASE_DIR, "09_statistics")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c(
  "B2M","TAP1","TAP2","TAPBP","HLA-A","HLA-B","HLA-C","NLRC5",
  "STAT1","IRF1","IRF9",
  "CXCL9","CXCL10","CXCL11",
  "GBP1","GBP2","GBP4","GBP5",
  "IDO1",
  "CD3D","CD3E","CD3G","CD8A","CD8B",
  "GZMA","GZMB","PRF1",
  "LAG3"
)

log2tpm_col <- function(g) paste0(g, "_log2tpm")

# ---- 読込・フィルタ（実データ仕様に合わせて修正）----
cohort <- read_csv(IN_FILE, show_col_types = FALSE) %>%
  filter(include_flag == TRUE)

cat("include_flag==TRUE後:", nrow(cohort), "\n")
cat("Grade:\n");  print(table(cohort$grade))
cat("TP53:\n");   print(table(cohort$tp53_status))
cat("IDH:\n");    print(table(cohort$idh_status))
cat("Source:\n"); print(table(cohort$source))

# ============================================================
# 共通関数（ラベルを実データ仕様に修正）
# ============================================================
compare_two_groups <- function(df, expr_col,
                               group_col  = "tp53_status",
                               mut_label  = "mutant",    # ← 修正
                               wt_label   = "wildtype") { # ← 修正
  df2 <- df %>%
    filter(.data[[group_col]] %in% c(mut_label, wt_label)) %>%
    select(expr = all_of(expr_col), group = all_of(group_col)) %>%
    drop_na()
  
  mut_vals <- df2$expr[df2$group == mut_label]
  wt_vals  <- df2$expr[df2$group == wt_label]
  n_mut <- length(mut_vals)
  n_wt  <- length(wt_vals)
  
  if (n_mut < 3 || n_wt < 3) {
    return(tibble(
      n_mut=n_mut, n_wt=n_wt,
      median_mut=NA_real_, median_wt=NA_real_, median_diff=NA_real_,
      mean_mut=NA_real_, mean_wt=NA_real_,
      sd_mut=NA_real_, sd_wt=NA_real_,
      HL_estimate=NA_real_, p_wilcox=NA_real_,
      cliffs_delta=NA_real_, effect_size_label=NA_character_
    ))
  }
  
  wtest <- wilcox.test(mut_vals, wt_vals, exact = FALSE, conf.int = TRUE)
  
  dominance <- outer(mut_vals, wt_vals, FUN = function(a, b) sign(a - b))
  delta <- sum(dominance) / (length(mut_vals) * length(wt_vals))
  
  abs_d <- abs(delta)
  label <- case_when(
    abs_d < 0.147 ~ "negligible",
    abs_d < 0.330 ~ "small",
    abs_d < 0.474 ~ "medium",
    TRUE          ~ "large"
  )
  
  tibble(
    n_mut        = n_mut,
    n_wt         = n_wt,
    median_mut   = median(mut_vals),
    median_wt    = median(wt_vals),
    median_diff  = median(mut_vals) - median(wt_vals),
    mean_mut     = mean(mut_vals),
    mean_wt      = mean(wt_vals),
    sd_mut       = sd(mut_vals),
    sd_wt        = sd(wt_vals),
    HL_estimate  = as.numeric(wtest$estimate),
    p_wilcox     = wtest$p.value,
    cliffs_delta = delta,
    effect_size_label = label
  )
}

# ============================================================
# 1. Grade4 主解析
# ============================================================
g4 <- cohort %>% filter(grade == "Grade4")
cat("\n=== Grade4 主解析 n=", nrow(g4), "===\n")
cat("TP53:", table(g4$tp53_status), "\n")

res_g4_main <- compare_two_groups(g4, "LAG3_log2tpm") %>%
  mutate(analysis="Grade4_main", gene="LAG3", subset="All", p_BH=p_wilcox)
print(res_g4_main %>% select(n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label))

# ============================================================
# 2. Source別層別
# ============================================================
cat("\n=== Grade4 Source別 ===\n")
res_source <- map_dfr(c("TCGA","CPTAC_HCMI"), function(src) {
  df_src <- g4 %>% filter(source == src)
  cat(src, "n=", nrow(df_src), " TP53:", table(df_src$tp53_status), "\n")
  compare_two_groups(df_src, "LAG3_log2tpm") %>%
    mutate(analysis="Grade4_by_source", gene="LAG3", subset=src)
}) %>% mutate(p_BH = p.adjust(p_wilcox, method="BH"))
print(res_source %>% select(subset, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))

# ============================================================
# 3. 感度解析: IDH WT限定
# ============================================================
cat("\n=== Grade4 IDH WT限定 ===\n")
g4_idhwt <- g4 %>% filter(idh_status == "wildtype")  # ← 修正
cat("IDH wildtype n=", nrow(g4_idhwt), " TP53:", table(g4_idhwt$tp53_status), "\n")

res_idh_sens <- compare_two_groups(g4_idhwt, "LAG3_log2tpm") %>%
  mutate(analysis="Grade4_IDH_WT_sensitivity", gene="LAG3",
         subset="IDH_wildtype", p_BH=p_wilcox)
print(res_idh_sens %>% select(n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label))

# ============================================================
# 4. Grade2/3 サブ解析
# ============================================================
cat("\n=== Grade2/3 サブ解析 ===\n")
res_sub <- map_dfr(c("Grade2","Grade3"), function(gr) {
  df_gr <- cohort %>% filter(grade == gr)
  cat(gr, "n=", nrow(df_gr), " TP53:", table(df_gr$tp53_status), "\n")
  compare_two_groups(df_gr, "LAG3_log2tpm") %>%
    mutate(analysis="SubAnalysis", gene="LAG3", subset=gr)
}) %>% mutate(p_BH = p.adjust(p_wilcox, method="BH"))
print(res_sub %>% select(subset, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))

# ============================================================
# 5. 28遺伝子探索
# ============================================================
cat("\n=== 28遺伝子探索 ===\n")

run_28genes <- function(df, analysis_label, subset_label) {
  map_dfr(GENES, function(g) {
    col <- log2tpm_col(g)
    if (!col %in% names(df)) { warning("列なし: ", col); return(NULL) }
    compare_two_groups(df, col) %>%
      mutate(gene=g, analysis=analysis_label, subset=subset_label)
  })
}

res_28_g4 <- run_28genes(g4, "28gene_Grade4", "Grade4") %>%
  mutate(p_BH = p.adjust(p_wilcox, method="BH"))
res_28_g2 <- run_28genes(cohort %>% filter(grade=="Grade2"), "28gene_Grade2", "Grade2") %>%
  mutate(p_BH = p.adjust(p_wilcox, method="BH"))
res_28_g3 <- run_28genes(cohort %>% filter(grade=="Grade3"), "28gene_Grade3", "Grade3") %>%
  mutate(p_BH = p.adjust(p_wilcox, method="BH"))

cat("Grade4 BH<0.05:\n")
print(res_28_g4 %>% filter(p_BH<0.05) %>%
        select(gene, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))
cat("Grade2 BH<0.05:\n")
print(res_28_g2 %>% filter(p_BH<0.05) %>%
        select(gene, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))
cat("Grade3 BH<0.05:\n")
print(res_28_g3 %>% filter(p_BH<0.05) %>%
        select(gene, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))

# ============================================================
# 6. 保存
# ============================================================
res_lag3_summary <- bind_rows(res_g4_main, res_source, res_idh_sens, res_sub) %>%
  select(analysis, subset, gene, n_mut, n_wt,
         median_mut, median_wt, median_diff,
         mean_mut, mean_wt, sd_mut, sd_wt,
         HL_estimate, p_wilcox, p_BH, cliffs_delta, effect_size_label)

write_csv(res_lag3_summary, file.path(OUT_DIR, "step09_lag3_summary.csv"))

res_28_all <- bind_rows(res_28_g4, res_28_g2, res_28_g3) %>%
  select(analysis, subset, gene, n_mut, n_wt,
         median_diff, HL_estimate, p_wilcox, p_BH, cliffs_delta, effect_size_label)

write_csv(res_28_all, file.path(OUT_DIR, "step09_28genes_exploratory.csv"))

cat("\n✅ 保存完了\n")
cat(" ", file.path(OUT_DIR, "step09_lag3_summary.csv"), "\n")
cat(" ", file.path(OUT_DIR, "step09_28genes_exploratory.csv"), "\n")

# ============================================================
# 7. 最終サマリ
# ============================================================
cat("\n", paste(rep("=",60), collapse=""), "\n")
cat("【Step09 完了サマリ】\n")
cat(paste(rep("=",60), collapse=""), "\n\n")

cat("■ LAG3 Grade4 全体\n")
with(res_g4_main, cat(sprintf(
  "  mutant n=%d, wildtype n=%d | median diff=%.3f | p=%.4f | δ=%.3f (%s)\n",
  n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label)))

cat("\n■ Source別\n")
pwalk(res_source, function(subset, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label, ...) {
  cat(sprintf("  %s: mutant n=%d, wildtype n=%d | diff=%.3f | p=%.4f | p_BH=%.4f | δ=%.3f (%s)\n",
              subset, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label))
})

cat("\n■ IDH wildtype限定\n")
with(res_idh_sens, cat(sprintf(
  "  mutant n=%d, wildtype n=%d | diff=%.3f | p=%.4f | δ=%.3f (%s)\n",
  n_mut, n_wt, median_diff, p_wilcox, cliffs_delta, effect_size_label)))

cat("\n■ Grade2/3\n")
pwalk(res_sub, function(subset, n_mut, n_wt, median_diff, p_wilcox, p_BH, cliffs_delta, effect_size_label, ...) {
  cat(sprintf("  %s: mutant n=%d, wildtype n=%d | diff=%.3f | p_BH=%.4f | δ=%.3f (%s)\n",
              subset, n_mut, n_wt, median_diff, p_BH, cliffs_delta, effect_size_label))
})
