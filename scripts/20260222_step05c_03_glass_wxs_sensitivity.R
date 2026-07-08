# ============================================================
# step05c_03_glass_wxs_sensitivity.R
# GLASS WXS-only主解析セット作成
# 入力: glass_final_cohort.csv（101例）
# 出力: glass_final_cohort_wxs.csv（WXSのみ）
#       glass_final_cohort_wxs_wgs.csv（WXS+WGS=101例・既存の別名保存）
# 出力先: results/TP53/20260221/05c_glass/
# ============================================================

library(tidyverse)

OUT_DIR <- here::here("results", "TP53", "20260221", "05c_glass")

# ============================================================
# 1. 読み込み
# ============================================================
cohort_all <- read.csv(file.path(OUT_DIR, "glass_final_cohort.csv"),
                       stringsAsFactors = FALSE)
cat("読み込み:", nrow(cohort_all), "行\n")
cat("source 内訳:\n")
print(table(cohort_all$source, useNA = "ifany"))

# ============================================================
# 2. WXS-only 抽出
# ============================================================
cohort_wxs <- cohort_all %>% filter(source == "GLASS_WXS")
cat("\nWXS-only 件数:", nrow(cohort_wxs), "（参照値: 79例）\n")

# TP53
cat("\nTP53 WT/Mut（WXS-only）:\n")
print(table(cohort_wxs$tp53_status, useNA = "ifany"))
cat("（参照値: WT=66, Mut=29）\n")

# IDH
cat("\nIDH WT/Mut（WXS-only）:\n")
print(table(cohort_wxs$idh_status, useNA = "ifany"))

# TP53 × IDH
cat("\nTP53 × IDH（WXS-only）:\n")
print(table(TP53 = cohort_wxs$tp53_status, IDH = cohort_wxs$idh_status, useNA = "ifany"))

# LAG3
lag3 <- cohort_wxs$LAG3_log2tpm
cat("\nLAG3 log2(TPM+1)（WXS-only）:\n")
cat("n=", sum(!is.na(lag3)),
    " median=", round(median(lag3, na.rm = TRUE), 3),
    " max=",    round(max(lag3, na.rm = TRUE), 3), "\n")

# TP53別LAG3中央値（参照値確認用）
cat("\nTP53別 LAG3中央値（WXS-only）:\n")
cohort_wxs %>%
  group_by(tp53_status) %>%
  summarise(n       = n(),
            median  = round(median(LAG3_log2tpm, na.rm = TRUE), 3),
            .groups = "drop") %>%
  print()

# ============================================================
# 3. WXS+WGS（101例）との比較サマリー
# ============================================================
cat("\n=== WXS-only vs WXS+WGS 比較 ===\n")
summary_tbl <- bind_rows(
  cohort_wxs %>%
    summarise(set       = "WXS-only",
              n         = n(),
              tp53_mut  = sum(tp53_status == "Mut"),
              tp53_wt   = sum(tp53_status == "WT"),
              idh_mut   = sum(idh_status  == "Mut"),
              lag3_med  = round(median(LAG3_log2tpm, na.rm = TRUE), 3)),
  cohort_all %>%
    summarise(set       = "WXS+WGS",
              n         = n(),
              tp53_mut  = sum(tp53_status == "Mut"),
              tp53_wt   = sum(tp53_status == "WT"),
              idh_mut   = sum(idh_status  == "Mut"),
              lag3_med  = round(median(LAG3_log2tpm, na.rm = TRUE), 3))
)
print(summary_tbl)
cat("\n参照値: n=95, TP53_Mut=29, TP53_WT=66, IDH_Mut=不明\n")

# ============================================================
# 4. 保存
# ============================================================
write.csv(cohort_wxs,
          file.path(OUT_DIR, "glass_final_cohort_wxs.csv"),
          row.names = FALSE)
write.csv(cohort_all,
          file.path(OUT_DIR, "glass_final_cohort_wxs_wgs.csv"),
          row.names = FALSE)

cat("\n✅ 保存完了\n")
cat("WXS-only:  ", file.path(OUT_DIR, "glass_final_cohort_wxs.csv"), "\n")
cat("WXS+WGS:   ", file.path(OUT_DIR, "glass_final_cohort_wxs_wgs.csv"), "\n")
