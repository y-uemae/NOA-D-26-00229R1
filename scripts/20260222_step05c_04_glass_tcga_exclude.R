# ============================================================
# step05c_04_glass_tcga_exclude.R
# GLASS内TCGAサンプル除外 + 最終セット確定
# 入力: glass_final_cohort_wxs.csv / glass_final_cohort_wxs_wgs.csv
# 出力: glass_final_cohort_wxs_notcga.csv（主解析）
#       glass_final_cohort_wxs_wgs_notcga.csv（感度解析）
# 出力先: results/TP53/20260221/05c_glass/
# ============================================================

library(tidyverse)

OUT_DIR <- here::here("results", "TP53", "20260221", "05c_glass")

# ============================================================
# 1. 読み込み
# ============================================================
wxs     <- read.csv(file.path(OUT_DIR, "glass_final_cohort_wxs.csv"),
                    stringsAsFactors = FALSE)
wxs_wgs <- read.csv(file.path(OUT_DIR, "glass_final_cohort_wxs_wgs.csv"),
                    stringsAsFactors = FALSE)

cat("WXS-only:  ", nrow(wxs), "行\n")
cat("WXS+WGS:   ", nrow(wxs_wgs), "行\n")

# ============================================================
# 2. TCGA混在確認
# ============================================================
# case_barcode（patient_id）の先頭でTCGA判定
# GLASSのIDは "GLSS-XX-XXXX" 形式
# TCGAは "TCGA-XX-XXXX" 形式
cat("\n=== case_barcode プレフィックス確認 ===\n")
wxs_wgs <- wxs_wgs %>%
  mutate(
    id_prefix = sub("^([A-Za-z]+)-.*", "\\1", case_barcode)
  )
cat("WXS+WGS id_prefix テーブル:\n")
print(table(wxs_wgs$id_prefix, useNA = "ifany"))

# TCGA該当サンプル
tcga_samples <- wxs_wgs %>% filter(id_prefix == "TCGA")
cat("\nTCGA該当サンプル数:", nrow(tcga_samples), "\n")
if (nrow(tcga_samples) > 0) {
  cat("TCGA case_barcode 例:", paste(head(tcga_samples$case_barcode, 5), collapse = ", "), "\n")
  cat("TCGA source:\n")
  print(table(tcga_samples$source))
}

# ============================================================
# 3. TCGA除外 → 最終セット
# ============================================================
wxs_notcga     <- wxs     %>%
  mutate(id_prefix = sub("^([A-Za-z]+)-.*", "\\1", case_barcode)) %>%
  filter(id_prefix != "TCGA") %>%
  select(-id_prefix)

wxs_wgs_notcga <- wxs_wgs %>%
  filter(id_prefix != "TCGA") %>%
  select(-id_prefix)

cat("\n=== TCGA除外後 サマリー ===\n")
for (label in c("WXS-only（主解析）", "WXS+WGS（感度解析）")) {
  df <- if (grepl("WXS-only", label)) wxs_notcga else wxs_wgs_notcga
  cat("\n---", label, "---\n")
  cat("n =", nrow(df), "\n")
  cat("TP53:\n"); print(table(df$tp53_status, useNA = "ifany"))
  cat("IDH:\n");  print(table(df$idh_status,  useNA = "ifany"))
  cat("source:\n"); print(table(df$source, useNA = "ifany"))
  lag3 <- df$LAG3_log2tpm
  cat("LAG3: median=", round(median(lag3, na.rm=TRUE), 3),
      " max=", round(max(lag3, na.rm=TRUE), 3), "\n")
}

cat("\n参照値: n=95, TP53_Mut=29, TP53_WT=66, WXS=79, WGS=16\n")

# TP53別LAG3中央値（主解析セット）
cat("\nTP53別 LAG3中央値（WXS-only, TCGA除外）:\n")
wxs_notcga %>%
  group_by(tp53_status) %>%
  summarise(n      = n(),
            median = round(median(LAG3_log2tpm, na.rm=TRUE), 3),
            .groups = "drop") %>%
  print()
cat("（引継書参照値: GLASS Grade4 median_diff=+0.418, p=0.0037, Cliff's δ=0.376）\n")

# ============================================================
# 4. 保存
# ============================================================
write.csv(wxs_notcga,
          file.path(OUT_DIR, "glass_final_cohort_wxs_notcga.csv"),
          row.names = FALSE)
write.csv(wxs_wgs_notcga,
          file.path(OUT_DIR, "glass_final_cohort_wxs_wgs_notcga.csv"),
          row.names = FALSE)

cat("\n✅ 保存完了\n")
cat("主解析（WXS, TCGA除外）:      ",
    file.path(OUT_DIR, "glass_final_cohort_wxs_notcga.csv"), "\n")
cat("感度解析（WXS+WGS, TCGA除外）:",
    file.path(OUT_DIR, "glass_final_cohort_wxs_wgs_notcga.csv"), "\n")
