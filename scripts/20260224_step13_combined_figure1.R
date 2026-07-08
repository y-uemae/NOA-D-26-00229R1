# =============================================================================
# 20260224_step13_combined_figure1.R
# Main Figure 1: 4パネル結合
#
# Layout:
#   [ A: GDC WT/Mut (2facet) ] [ B: GLASS WT/Mut ] [ C: Forest ]
#   [          D: TP53 4群 boxplot (full width)                 ]
#
# 入力:
#   08_final_cohort/final_cohort.csv
#   05c_glass/glass_final_cohort_wxs_notcga.csv
#   12_subgroup/step12b_subgroup_classified.csv
#   11_meta_analysis/step11_meta_forest_summary.csv
#   11_meta_analysis/step11_meta_results.csv
# 出力:
#   13_visualization/fig1_combined.pdf
#   13_visualization/fig1_combined_450dpi.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ragg)

# ── 0. 設定 ──────────────────────────────────────────────────────────────────
BASE_DIR   <- here::here("results", "TP53", "20260221")
GDC_CSV    <- file.path(BASE_DIR, "08_final_cohort/final_cohort.csv")
GLASS_CSV  <- file.path(BASE_DIR, "05c_glass/glass_final_cohort_wxs_notcga.csv")
SUB_CSV    <- file.path(BASE_DIR, "12_subgroup/step12b_subgroup_classified.csv")
FOREST_CSV <- file.path(BASE_DIR, "11_meta_analysis/step11_meta_forest_summary.csv")
RESULT_CSV <- file.path(BASE_DIR, "11_meta_analysis/step11_meta_results.csv")
OUT_DIR    <- file.path(BASE_DIR, "13_visualization")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

DPI       <- 450
WIDTH_IN  <- 14.0
HEIGHT_IN <- 10.0

# ── 共通色（引継書固定）──────────────────────────────────────────────────────
COL_WT    <- "#AAAAAA"
COL_MUT   <- "#E64B35"
COL_4 <- c(
  "WT"             = "#AAAAAA",
  "Hotspot"        = "#E64B35",
  "Truncating"     = "#4DBBD5",
  "Other_missense" = "#00A087"
)
COL_GDC_TCGA  <- "#3C5488"
COL_GDC_CPTAC <- "#E07B54"
COL_GLASS     <- "#8491B4"
COL_RE        <- "#E64B35"
COL_FE        <- "#777777"

# ── 共通テーマ ────────────────────────────────────────────────────────────────
theme_paper <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, face = "bold"),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    axis.text.x      = element_text(size = 9),
    legend.position  = "none",
    plot.title       = element_text(size = 11, face = "bold"),
    panel.spacing    = unit(3, "mm"),
    plot.margin      = margin(4, 4, 4, 4, unit = "mm")
  )

# =============================================================================
# Panel A: GDC Grade4 WT/Mut（source別facet）
# =============================================================================
gdc_raw <- read_csv(GDC_CSV, show_col_types = FALSE) %>%
  filter(include_flag == TRUE, grade == "Grade4") %>%
  mutate(
    tp53_label   = factor(if_else(tp53_status == "mutant", "Mut", "WT"),
                          levels = c("WT", "Mut")),
    source_label = factor(if_else(source == "CPTAC_HCMI", "CPTAC/HCMI", source),
                          levels = c("TCGA", "CPTAC/HCMI"))
  )

gdc_stats <- gdc_raw %>%
  group_by(source_label) %>%
  summarise(
    n_wt     = sum(tp53_label == "WT"),
    n_mut    = sum(tp53_label == "Mut"),
    med_diff = median(LAG3_log2tpm[tp53_label == "Mut"], na.rm = TRUE) -
      median(LAG3_log2tpm[tp53_label == "WT"],  na.rm = TRUE),
    p_val    = wilcox.test(LAG3_log2tpm[tp53_label == "Mut"],
                           LAG3_log2tpm[tp53_label == "WT"],
                           exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = if_else(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val)),
    n_label = sprintf("WT n=%d\nMut n=%d", n_wt, n_mut)
  )

y_max_ab <- max(gdc_raw$LAG3_log2tpm,
                read_csv(GLASS_CSV, show_col_types = FALSE)$LAG3_log2tpm,
                na.rm = TRUE)
y_top_ab <- ceiling(y_max_ab * 10) / 10 + 0.2
y_lim_ab <- c(-0.15, y_top_ab + 0.75)
sig_y_ab <- y_top_ab + 0.10
dlt_y_ab <- y_top_ab + 0.42

panelA <- ggplot(gdc_raw,
                 aes(x = tp53_label, y = LAG3_log2tpm,
                     color = tp53_label)) +
  geom_jitter(width = 0.17, alpha = 0.28, size = 0.7, shape = 16) +
  geom_boxplot(aes(fill = tp53_label), width = 0.42, outlier.shape = NA,
               alpha = 0.15, color = "black", linewidth = 0.45) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.38, linewidth = 0.7, color = "black", fatten = 1) +
  # ブラケット
  geom_segment(data = gdc_stats,
               aes(x=1.07, xend=1.93, y=sig_y_ab+0.08, yend=sig_y_ab+0.08),
               inherit.aes=FALSE, color="black", linewidth=0.35) +
  geom_segment(data = gdc_stats,
               aes(x=1.07, xend=1.07, y=sig_y_ab+0.03, yend=sig_y_ab+0.08),
               inherit.aes=FALSE, color="black", linewidth=0.35) +
  geom_segment(data = gdc_stats,
               aes(x=1.93, xend=1.93, y=sig_y_ab+0.03, yend=sig_y_ab+0.08),
               inherit.aes=FALSE, color="black", linewidth=0.35) +
  geom_text(data = gdc_stats,
            aes(x=1.5, y=sig_y_ab, label=p_label),
            inherit.aes=FALSE, size=2.8, fontface="italic", color="black") +
  geom_text(data = gdc_stats,
            aes(x=1.5, y=dlt_y_ab,
                label=sprintf("\u0394med=+%.3f", med_diff)),
            inherit.aes=FALSE, size=2.5, color="gray35") +
  geom_text(data = gdc_stats %>%
              pivot_longer(c(n_wt,n_mut), names_to="grp", values_to="n") %>%
              mutate(x=if_else(grp=="n_wt",1,2),
                     lbl=sprintf("n=%d",n)),
            aes(x=x, y=-0.08, label=lbl),
            inherit.aes=FALSE, size=2.4, color="gray45") +
  facet_wrap(~source_label, nrow=1) +
  scale_color_manual(values=c("WT"=COL_WT,"Mut"=COL_MUT)) +
  scale_fill_manual( values=c("WT"=COL_WT,"Mut"=COL_MUT)) +
  scale_y_continuous(limits=y_lim_ab, expand=c(0,0),
                     breaks=seq(0,floor(y_top_ab)+1,by=1)) +
  labs(title="A   GDC Grade 4", x=NULL,
       y="LAG3 expression\n[log2(TPM+1)]") +
  theme_paper

# =============================================================================
# Panel B: GLASS WT/Mut
# =============================================================================
glass_raw <- read_csv(GLASS_CSV, show_col_types = FALSE) %>%
  mutate(tp53_label = factor(if_else(tp53_status=="Mut","Mut","WT"),
                             levels=c("WT","Mut")))

glass_stats <- glass_raw %>%
  summarise(
    n_wt=sum(tp53_label=="WT"), n_mut=sum(tp53_label=="Mut"),
    med_diff=median(LAG3_log2tpm[tp53_label=="Mut"],na.rm=TRUE)-
      median(LAG3_log2tpm[tp53_label=="WT"],na.rm=TRUE),
    p_val=wilcox.test(LAG3_log2tpm[tp53_label=="Mut"],
                      LAG3_log2tpm[tp53_label=="WT"],exact=FALSE)$p.value
  ) %>%
  mutate(p_label=if_else(p_val<0.001,"p < 0.001",sprintf("p = %.3f",p_val)))

panelB <- ggplot(glass_raw,
                 aes(x=tp53_label, y=LAG3_log2tpm, color=tp53_label)) +
  geom_jitter(width=0.17, alpha=0.38, size=0.9, shape=16) +
  geom_boxplot(aes(fill=tp53_label), width=0.42, outlier.shape=NA,
               alpha=0.15, color="black", linewidth=0.45) +
  stat_summary(fun=median, geom="crossbar",
               width=0.38, linewidth=0.7, color="black", fatten=1) +
  annotate("segment", x=1.07,xend=1.93, y=sig_y_ab+0.08,yend=sig_y_ab+0.08,
           color="black",linewidth=0.35) +
  annotate("segment", x=1.07,xend=1.07, y=sig_y_ab+0.03,yend=sig_y_ab+0.08,
           color="black",linewidth=0.35) +
  annotate("segment", x=1.93,xend=1.93, y=sig_y_ab+0.03,yend=sig_y_ab+0.08,
           color="black",linewidth=0.35) +
  annotate("text", x=1.5, y=sig_y_ab, label=glass_stats$p_label,
           size=2.8, fontface="italic", color="black") +
  annotate("text", x=1.5, y=dlt_y_ab,
           label=sprintf("\u0394med=+%.3f",glass_stats$med_diff),
           size=2.5, color="gray35") +
  annotate("text", x=c(1,2), y=-0.08,
           label=c(sprintf("n=%d",glass_stats$n_wt),
                   sprintf("n=%d",glass_stats$n_mut)),
           size=2.4, color="gray45") +
  scale_color_manual(values=c("WT"=COL_WT,"Mut"=COL_MUT)) +
  scale_fill_manual( values=c("WT"=COL_WT,"Mut"=COL_MUT)) +
  scale_y_continuous(limits=y_lim_ab, expand=c(0,0),
                     breaks=seq(0,floor(y_top_ab)+1,by=1)) +
  labs(title="B   GLASS (WXS, non-TCGA)", x=NULL, y=NULL) +
  theme_paper

# =============================================================================
# Panel C: Forest plot（theme_void版・軽量）
# =============================================================================
forest_raw <- read_csv(FOREST_CSV, show_col_types = FALSE)
result_raw <- read_csv(RESULT_CSV, show_col_types = FALSE)
re_row <- filter(result_raw, model=="Random Effect")
fe_row <- filter(result_raw, model=="Fixed Effect")

# ── パネルC座標定義（forest図1.5倍拡大・右上固定版）──────────────────────────
# 文字位置は変更なし、forest描画ゾーンのみ拡大
# 右上固定: XCF_MAX / YC_MAX を基準に左下方向へ拡張

XC_MIN=-2.60; XC_MAX=2.20
XCF_MIN=-0.70   # 左に拡張（従来-0.35 → -0.70）
XCF_MAX= 0.95   # 右端固定
XC_LBL=-2.58; XC_BETA=1.00; XC_P=1.82
YC_HDR=7.3; YC_MIN=0.2; YC_MAX=7.8; YC_SEP=3.0

# forest描画ゾーンのy範囲も拡張（上端固定・下方向へ）
# 行間隔を広げることで縦方向1.5倍
YC_ROWS <- c(         # 各studyのy座標を再定義
  GDC_TCGA       = 7.0,
  GDC_CPTAC_HCMI = 5.8,
  GLASS_WXS      = 4.6,
  Pooled_RE      = 2.8,
  Pooled_FE      = 1.6
)

pC_df <- tibble(
  y   = as.numeric(YC_ROWS),
  beta=c(0.375,0.135,0.497,re_row$beta_pooled,fe_row$beta_pooled),
  ci_lo=c(0.188,-0.009,0.192,re_row$ci_lo,fe_row$ci_lo),
  ci_hi=c(0.563,0.279,0.802,re_row$ci_hi,fe_row$ci_hi),
  col=c(COL_GDC_TCGA,COL_GDC_CPTAC,COL_GLASS,COL_RE,COL_FE),
  rt=c("ind","ind","ind","re","fe"),
  lbl=c("GDC – TCGA","GDC – CPTAC/HCMI","GLASS – WXS","Pooled (RE)","Pooled (FE)"),
  n_lbl=c("Mut=79/WT=166","Mut=68/WT=129","Mut=24/WT=55","",""),
  bv=c("0.375 (0.188, 0.563)","0.135 (−0.009, 0.279)","0.497 (0.192, 0.802)",
       sprintf("%.3f (%.3f, %.3f)",re_row$beta_pooled,re_row$ci_lo,re_row$ci_hi),
       sprintf("%.3f (%.3f, %.3f)",fe_row$beta_pooled,fe_row$ci_lo,fe_row$ci_hi)),
  pv=c("< 0.001","0.067","0.002","0.004","< 0.001"),
  fc=c("plain","plain","plain","bold","plain")
)

# x軸・区切り線のy座標も更新
YC_AXIS_Y  <- 1.0    # x軸ラインのy位置
YC_SEP     <- 3.7    # 区切り破線
YC_HDR     <- 8.2    # ヘッダー行

SZ <- 3.2  # forest内文字サイズ

panelC <- ggplot() +
  annotate("segment",x=0,xend=0,y=YC_AXIS_Y+0.05,yend=YC_HDR-0.2,
           color="black",linewidth=0.4) +
  annotate("segment",x=XCF_MIN,xend=XCF_MAX,
           y=YC_AXIS_Y,yend=YC_AXIS_Y,
           color="black",linewidth=0.5) +
  { tks=seq(-0.4,0.8,0.2)   # 拡張されたxゾーン用に目盛り追加
  list(annotate("segment",x=tks,xend=tks,
                y=YC_AXIS_Y,yend=YC_AXIS_Y-0.18,
                color="black",linewidth=0.4),
       annotate("text",x=tks,y=YC_AXIS_Y-0.42,
                label=sprintf("%.1f",tks),size=2.8,hjust=0.5)) } +
  annotate("text",x=(XCF_MIN+XCF_MAX)/2,y=YC_AXIS_Y-0.88,
           label=expression(beta~"coeff. (TP53 Mut vs WT)"),
           size=3.0,hjust=0.5) +
  annotate("segment",x=XC_MIN,xend=XC_MAX,
           y=YC_SEP,yend=YC_SEP,
           color="gray70",linewidth=0.35,linetype="dashed") +
  geom_segment(data=filter(pC_df,rt=="ind"),
               aes(x=ci_lo,xend=ci_hi,y=y,yend=y,color=col),linewidth=1.3) +
  geom_segment(data=filter(pC_df,rt=="re"),
               aes(x=ci_lo,xend=ci_hi,y=y,yend=y),
               color=COL_RE,linewidth=2.6) +
  geom_segment(data=filter(pC_df,rt=="fe"),
               aes(x=ci_lo,xend=ci_hi,y=y,yend=y),
               color=COL_FE,linewidth=1.2,linetype="dashed") +
  geom_point(data=filter(pC_df,rt=="ind"),
             aes(x=beta,y=y,color=col),shape=15,size=3.8) +
  geom_point(data=filter(pC_df,rt=="re"),aes(x=beta,y=y),
             shape=23,size=6.0,fill=COL_RE,color=COL_RE) +
  geom_point(data=filter(pC_df,rt=="fe"),aes(x=beta,y=y),
             shape=23,size=4.2,fill=COL_FE,color=COL_FE) +
  geom_text(data=pC_df,aes(x=XC_LBL,y=y,label=lbl,
                           color=col,fontface=fc),hjust=0,size=SZ) +
  geom_text(data=filter(pC_df,n_lbl!=""),
            aes(x=XC_LBL,y=y-0.55,label=n_lbl),
            inherit.aes=FALSE,hjust=0,size=SZ-0.5,color="gray55") +
  annotate("text",x=XC_LBL,y=YC_HDR,label="Study",
           hjust=0,size=SZ+0.3,fontface="bold",color="black") +
  annotate("text",x=XC_BETA,y=YC_HDR,label="\u03b2 (95% CI)",
           hjust=0,size=SZ+0.3,fontface="bold",color="black") +
  annotate("text",x=XC_P,y=YC_HDR,label="p",
           hjust=0,size=SZ+0.3,fontface="bold",color="black") +
  geom_text(data=pC_df,aes(x=XC_BETA,y=y,label=bv,
                           color=col,fontface=fc),hjust=0,size=SZ-0.2) +
  geom_text(data=pC_df,aes(x=XC_P,y=y,label=pv,
                           color=col,fontface=fc),hjust=0,size=SZ-0.2) +
  annotate("text",x=XC_MIN,y=YC_AXIS_Y-1.5,
           label=sprintf("I\u00b2=%.1f%%  \u03c4\u00b2=%.4f  Qp=%.3f",
                         re_row$I2,re_row$tau2,re_row$Q_p),
           hjust=0,size=2.6,color="gray40") +
  scale_color_identity() +
  scale_x_continuous(limits=c(XC_MIN,XC_MAX),expand=c(0,0)) +
  scale_y_continuous(limits=c(YC_AXIS_Y-1.8, YC_HDR+0.3),expand=c(0,0)) +
  labs(title="C   Meta-analysis (GDC + GLASS)") +
  theme_void() +
  theme(plot.title=element_text(size=11,face="bold"),
        plot.margin=margin(4,4,4,4,unit="mm"))

# =============================================================================
# Panel D: TP53 4群 boxplot
# =============================================================================
sub_df <- read_csv(SUB_CSV, show_col_types = FALSE) %>%
  mutate(tp53_class4=factor(tp53_class4,
                            levels=c("WT","Hotspot","Truncating","Other_missense")))

n_grp <- sub_df %>% count(tp53_class4)
x_lbl <- setNames(
  sprintf("%s\n(n=%d)", n_grp$tp53_class4, n_grp$n),
  n_grp$tp53_class4
)

wt_v   <- sub_df$LAG3_log2tpm[sub_df$tp53_class4=="WT"]
p_bh   <- p.adjust(c(
  wilcox.test(sub_df$LAG3_log2tpm[sub_df$tp53_class4=="Hotspot"],   wt_v,exact=FALSE)$p.value,
  wilcox.test(sub_df$LAG3_log2tpm[sub_df$tp53_class4=="Truncating"],wt_v,exact=FALSE)$p.value,
  wilcox.test(sub_df$LAG3_log2tpm[sub_df$tp53_class4=="Other_missense"],wt_v,exact=FALSE)$p.value
), method="BH")

y_top_d <- ceiling(max(sub_df$LAG3_log2tpm,na.rm=TRUE)*10)/10+0.2
y_lim_d <- c(-0.15, y_top_d+0.95)
b1y <- y_top_d+0.15; b2y <- y_top_d+0.52
fmt_p <- function(p) if(p<0.001) "p[BH] < 0.001" else sprintf("p[BH] = %.3f",p)

brk <- tibble(
  x1=c(1,1), x2=c(2,3),
  by=c(b1y,b2y), py=c(b1y+0.12,b2y+0.12),
  lab=c(fmt_p(p_bh[1]),fmt_p(p_bh[2]))
)

panelD <- ggplot(sub_df,aes(x=tp53_class4,y=LAG3_log2tpm,
                            color=tp53_class4,fill=tp53_class4)) +
  geom_jitter(width=0.17,alpha=0.28,size=0.7,shape=16) +
  geom_boxplot(width=0.42,outlier.shape=NA,
               alpha=0.15,color="black",linewidth=0.45) +
  stat_summary(fun=median,geom="crossbar",
               width=0.38,linewidth=0.7,color="black",fatten=1) +
  geom_segment(data=brk,aes(x=x1+0.08,xend=x2-0.08,y=by,yend=by),
               inherit.aes=FALSE,color="black",linewidth=0.4) +
  geom_segment(data=brk,aes(x=x1+0.08,xend=x1+0.08,y=by-0.07,yend=by),
               inherit.aes=FALSE,color="black",linewidth=0.4) +
  geom_segment(data=brk,aes(x=x2-0.08,xend=x2-0.08,y=by-0.07,yend=by),
               inherit.aes=FALSE,color="black",linewidth=0.4) +
  geom_text(data=brk,aes(x=(x1+x2)/2,y=py,label=lab),
            inherit.aes=FALSE,size=3.0,color="black",fontface="italic",
            parse=TRUE) +
  scale_color_manual(values=COL_4) +
  scale_fill_manual( values=COL_4) +
  scale_x_discrete(labels=x_lbl) +
  scale_y_continuous(limits=y_lim_d,expand=c(0,0),
                     breaks=seq(0,floor(y_top_d)+1,by=1)) +
  labs(title="D   TP53 subgroup analysis  –  GDC Grade 4 (pre-specified)",
       x=NULL, y="LAG3 expression\n[log2(TPM+1)]") +
  theme_paper

# =============================================================================
# 結合・出力
# =============================================================================
# 上段: A(2facet) + B(1)  → widths=c(2,1)
# 下段: C(forest) + D(4群) → widths=c(1.4,1)

top_row <- panelA + panelB +
  plot_layout(widths = c(2, 1))

bot_row <- panelC + panelD +
  plot_layout(widths = c(1.4, 1))

fig1 <- top_row / bot_row +
  plot_layout(heights = c(1, 1))

pdf_path <- file.path(OUT_DIR, "fig1_combined.pdf")
pdf(pdf_path, width = WIDTH_IN, height = HEIGHT_IN)
print(fig1)
dev.off()
cat("PDF:", pdf_path, "\n")

png_path <- file.path(OUT_DIR, sprintf("fig1_combined_%ddpi.png", DPI))
agg_png(png_path, width = WIDTH_IN, height = HEIGHT_IN,
        units = "in", res = DPI, scaling = 1.0)
print(fig1)
dev.off()
cat("PNG:", png_path, "\n")

cat("\n=== Step 13 combined Figure 1 完了 ===\n")
cat("レイアウト: 上段(A+B) / 下段(C+D)\n")
cat("サイズ:", WIDTH_IN, "×", HEIGHT_IN, "inch /", DPI, "dpi\n")
