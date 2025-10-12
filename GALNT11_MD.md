# Galnt11 数据整理
## Fig1A
```R

cd ~/home/syc/Data_Prepare

setwd("~/home/syc/Data_Prepare/Figure1/")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
 library(org.Hs.eg.db)


gene_txdb <- as.data.frame(GenomicFeatures::genes(txdb))
gene_txdb$SYMBOL <- mapIds(org.Hs.eg.db, 
  keys = gene_txdb$gene_id, 
  column = "SYMBOL", 
  keytype = "ENTREZID", 
  multiVals = "first")

gene_7q_all = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 148200001 & gene_txdb$end <= 159345973),]
#write.csv(gene_7q_all,"gene_7q36.1_36.3_all.csv")
gene_7q_all = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 148200001 & gene_txdb$end <= 159345973),"SYMBOL"] #96个基因

gene_list <- gene_txdb[gene_txdb$SYMBOL%in%c(na.omit(gene_7q_all)),] #167
table(gene_list$seqnames)
# 1. 确保 CNV 数据列名正确
cnv_data <- read.table("~/home/syc/Public_Data/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg", stringsAsFactors = FALSE,header=TRUE)

#提取LAML
library(stringr)
library(data.table)
Alter_tcga <- read.table("~/home/syc/Galnt11_DXT/alterations_across_samples.tsv",sep="\t",header=T)

cnv_data <- cnv_data %>%
  mutate(Sample_short = str_sub(Sample, 1, 15))
clinic_data <- fread("~/home/syc/Public_Data/laml_tcga_pan_can_atlas_2018_clinical_data.tsv")
cnv_data_laml <- cnv_data[cnv_data$Sample_short%in%c(clinic_data$'Sample ID'),]

library(AnnotationHub)
library(biovizBase)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)

chr7_cnv <- cnv_data_laml %>%
  filter(Chromosome == 7)

chr7_cnv <- chr7_cnv %>%
  mutate(CNV_Type = case_when(
    Segment_Mean < -0.1 ~ "Loss",
    Segment_Mean > 0.1  ~ "Gain",
    TRUE ~ "Neutral"
  ))
 
p1 <- ggplot(chr7_cnv %>% filter(CNV_Type != "Neutral"),  # 只显示 Loss 或 Gain
       aes(x = Start, xend = End, y = Sample_short, yend = Sample_short, color = CNV_Type)) +
  geom_segment(size = 1) +
  scale_color_manual(values = c("Loss" = "blue", "Gain" = "red")) +
  labs(title = "CNV Distribution on Chromosome 7 in LAML Samples",
       x = "Chromosome 7 Position (bp)",
       y = "Sample",
       color = "CNV Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))


sample <- c("TCGA-AB-2952-03","TCGA-AB-2944-03","TCGA-AB-2941-03","TCGA-AB-2928-03","TCGA-AB-2929-03","TCGA-AB-2917-03","TCGA-AB-2920","TCGA-AB-2894-03",
  "TCGA-AB-2885-03","TCGA-AB-2882-03","TCGA-AB-2883-03","TCGA-AB-2874-03","TCGA-AB-2857-03","TCGA-AB-2838-03","TCGA-AB-2820-03","TCGA-AB-2817-03",
  "TCGA-AB-2887-03","TCGA-AB-2935-03","TCGA-AB-2915-03","TCGA-AB-2878-03",
  "TCGA-AB-3007-03","TCGA-AB-2923-03",
  "TCGA-AB-2805-03","TCGA-AB-2939-03","TCGA-AB-2949-03","TCGA-AB-2938-03","TCGA-AB-2904-03","TCGA-AB-2806-03")

library(ggplot2)
library(dplyr)

chr7_cnv2 <- chr7_cnv %>%
  filter(Sample_short %in% sample)

chr7_cnv2$Sample_short <- factor(chr7_cnv2$Sample_short, levels = sample)
p1 <- ggplot(chr7_cnv2 %>% filter(CNV_Type != "Neutral"),
             aes(x = Start, xend = End, y = Sample_short, yend = Sample_short, color = CNV_Type)) +
  geom_segment(size = 1) +
  scale_color_manual(values = c("Loss" = "blue", "Gain" = "red")) +
  labs(title = "CNV Distribution on Chromosome 7 in LAML Samples",
       x = "Chromosome 7 Position (bp)",
       y = "Sample",
       color = "CNV Type") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6))

print(p1)

#查看区域
chr7_cnv[chr7_cnv$Sample_short%in%c("TCGA-AB-2938-03"),]
chr7_cnv[chr7_cnv$Sample_short%in%c("TCGA-AB-2939-03"),]

library(ggplot2)
library(dplyr)

# 假设 chr7_cnv2 已经准备好，Sample_short 是因子顺序
# gene_7q_all = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 127110765 & gene_txdb$end <= 129810467),]
# gene_7q_all = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 144105941 & gene_txdb$end <= 158385118),]


library(ggplot2)
library(dplyr)


highlight_regions <- data.frame(
  xmin = c(127110765, 144105941, 140212419),
  xmax = c(129810467, 158385118, 144105940)
)

genes_to_mark <- data.frame(
  gene = c("EZH2", "KMT2C", "GALNT11"),
  pos  = c(148807257, 152134922, 152025674)
)

library(dplyr)
library(tidyr)
library(dplyr)

#取 Loss 片段，并计算中点与长度
loss_df <- chr7_cnv2 %>%
  dplyr::filter(CNV_Type == "Loss") %>%
  dplyr::mutate(
    mid   = (Start + End) / 2,
    width = pmax(End - Start, 1)  # 作为权重，避免0
  ) %>%
  dplyr::select(Sample_short, Start, End, mid, width)


p <- ggplot() +
  geom_segment(data = chr7_cnv2 %>% filter(CNV_Type != "Neutral"),
               aes(x = Start, xend = End, y = Sample_short, yend = Sample_short, color = CNV_Type),
               size = 1) +
  scale_color_manual(values = c("Loss" = "blue", "Gain" = "red")) +
  geom_density(data = loss_df, 
              aes(x = (Start + End)/2, y = ..scaled.. * length(unique(chr7_cnv2$Sample_short)) + max(as.numeric(chr7_cnv2$Sample_short)) * 0.05),
               color = "black", fill = "lightblue", alpha = 0.4) +
  geom_rect(data = highlight_regions,
            inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = NA, color = "red", size = 1.2) +
  geom_vline(data = genes_to_mark,
             aes(xintercept = pos),
             color = "red", linetype = "dashed", size = 0.8) +
  geom_text(data = genes_to_mark,
            aes(x = pos, y = max(as.numeric(chr7_cnv2$Sample_short)) + 1, label = gene),
            color = "red", angle = 90, vjust = -0.5, size = 3) +
  labs(title = "Chromosome 7 CNV Heatmap with Loss Density",
       x = "Chromosome 7 Position (bp)",
       y = "Sample",
       color = "CNV Type") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6))

print(p)

# 5. 保存图
ggsave(p, file="Fig1A-TCGA-LAML_CNV_subset.png", height=4, width=6)
```
<img src="/Github_Fig/Figure1/Fig1A-TCGA-LAML_CNV_subset.png">



## Supplementary Fig1A
```R
/usr/local/R4.2/bin/R

#TCGA-LAML
degs <- read.csv("~/home/syc/Galnt11_DXT/TCGA/TCGA-Bulk_7Q_DEG.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
res$symbol <- mapIds(x = anno_data,keys = rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res$entrez <- mapIds(x = anno_data,keys = rownames(res),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res$symbol
AA <- as.character(AA)
res$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res$GENENAME <- as.character(res$GENENAME)

#GSE6891
DEG_limma_GSE6891 <- read.csv("~/home/syc/Public_Data/GSE6891_DEGs_7q_vs_NN.csv")

#GSE10358
DEG_limma_GSE10358 <- read.csv("~/home/syc/Public_Data/GSE10358_DEGs_7q_vs_NN.csv")

DEG_limma_GSE6891[DEG_limma_GSE6891$X=="GALNT11",]



#读取CDRs上的编码基因
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
 library(org.Hs.eg.db)


gene_txdb <- as.data.frame(GenomicFeatures::genes(txdb))
gene_txdb$SYMBOL <- mapIds(org.Hs.eg.db, 
  keys = gene_txdb$gene_id, 
  column = "SYMBOL", 
  keytype = "ENTREZID", 
  multiVals = "first")

gene_7q_CDR1 = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 127110765 & gene_txdb$end <= 129810467),"SYMBOL"] #50个基因
gene_7q_CDR2 = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 140212419 & gene_txdb$end <= 144105940),"SYMBOL"] #42个基因
gene_7q_CDR3 = gene_txdb[which(gene_txdb$seqnames == "chr7" & gene_txdb$start >= 144105941 & gene_txdb$end <= 158385118),"SYMBOL"] #170个基因

gene_7q_CDRs <- c(gene_7q_CDR1,gene_7q_CDR2,gene_7q_CDR3)
library(tidyverse)
library(GEOquery)
library(tinyarray)
wide_means_TCGA  <- res[res$symbol%in%c(na.omit(unique(gene_7q_CDRs))),] #180
wide_means_TCGA <- wide_means_TCGA[wide_means_TCGA$log2FoldChange<0,]
combined <- wide_means_TCGA
combined <- combined[is.finite(combined$log2FoldChange), ]
# 1. 按 log2FoldChange 排序并创建 rank（纵轴）
combined <- combined %>%
  arrange(log2FoldChange) %>%
  mutate(rank_ordered = row_number())

top5_genes <- combined %>%
  arrange(log2FoldChange) %>%   # log2FoldChange 从小到大
  slice(1:10) %>%        # 取最小的5个
  pull(symbol)
 [1] "LEP"      "TMEM176B" "OR9A4"    "BLACE"    "CLEC5A"   "TMEM176A"
 [7] "TAS2R39"  "FLNC"     "TPK1"     "RNY1"

# 3. 创建标记列
combined <- combined %>%
  mutate(Highlight = case_when(
    symbol == "KMT2C" ~ "KMT2C",
    symbol == "GALNT11" ~ "GALNT11",
    TRUE ~ "Other"

  ))

p_TCGA <- ggplot(combined, aes(x = rank_ordered, y = log2FoldChange)) +
  geom_point(aes(color = Highlight), size = 2) +
  geom_text(data = subset(combined, Highlight %in% c("KMT2C","GALNT11")),
            aes(label = symbol),
            vjust = -1, size = 4) +
  scale_color_manual(values = c("KMT2C" = "red", "GALNT11" = "red", "Other" = "blue")) +
  labs(title = "Scatter plot of log2FC vs rank",
       x = "Rank (sorted by log2FoldChange)",
       y = "log2FoldChange") +
  theme_bw()


wide_means_GSE6891  <- DEG_limma_GSE6891[DEG_limma_GSE6891$X%in%c(na.omit(unique(gene_7q_CDRs))),] #71
wide_means_GSE6891 <- wide_means_GSE6891[wide_means_GSE6891$logFC<0,]
combined <- wide_means_GSE6891
combined <- combined[is.finite(combined$logFC), ]
combined <- combined %>%
  arrange(logFC) %>%
  mutate(rank_ordered = row_number())

top5_genes <- combined %>%
  arrange(logFC) %>%   # logFC 从小到大
  slice(1:10) %>%        # 取最小的5个
  pull(X)
[1] "GALNT11" "GSTK1"   "TPK1"    "MRPS33"  "FASTK"   "SLC4A2"  "MKRN1"
 [8] "CDK5"    "CHPF2"   "CALU"

combined <- combined %>%
  mutate(Highlight = case_when(
    X == "KMT2C" ~ "KMT2C",
    X == "GALNT11" ~ "GALNT11",
    TRUE ~ "Other"
  ))

p_GSE1 <- ggplot(combined, aes(x = rank_ordered, y = logFC)) +
  geom_point(aes(color = Highlight), size = 2) +
  geom_text(data = subset(combined, Highlight %in% c("KMT2C","GALNT11")),
            aes(label = X),
            vjust = -1, size = 4) +
  scale_color_manual(values = c("KMT2C" = "red", "GALNT11" = "red", "Other" = "blue")) +
  labs(title = "Scatter plot of logFC vs rank",
       x = "Rank (sorted by logFC)",
       y = "logFC") +
  theme_bw()



wide_means_GSE10358  <- DEG_limma_GSE10358[DEG_limma_GSE10358$X%in%c(na.omit(unique(gene_7q_CDRs))),] #71
wide_means_GSE10358 <- wide_means_GSE10358[wide_means_GSE10358$logFC<0,]
combined <- wide_means_GSE10358
combined <- combined[is.finite(combined$logFC), ]
combined <- combined %>%
  arrange(logFC) %>%
  mutate(rank_ordered = row_number())

top5_genes <- combined %>%
  arrange(logFC) %>%   # logFC 从小到大
  slice(1:10) %>%        # 取最小的5个
  pull(X)
 [1] "C7orf33"      "NDUFB2-AS1"   "CNPY1"        "SMARCD3"      "LINC00996"
 [6] "MKRN1"        "STRIP2"       "ATP6V0E2-AS1" "TMEM176A"     "ZNF746"

combined <- combined %>%
  mutate(Highlight = case_when(
    X == "KMT2C" ~ "KMT2C",
    X == "GALNT11" ~ "GALNT11",
    TRUE ~ "Other"
  ))

p_GSE2 <- ggplot(combined, aes(x = rank_ordered, y = logFC)) +
  geom_point(aes(color = Highlight), size = 2) +
  geom_text(data = subset(combined, Highlight %in% c("KMT2C","GALNT11")),
            aes(label = X),
            vjust = -1, size = 4) +
  scale_color_manual(values = c("KMT2C" = "red", "GALNT11" = "red", "Other" = "blue")) +
  labs(title = "Scatter plot of logFC vs rank",
       x = "Rank (sorted by logFC)",
       y = "logFC") +
  theme_bw()

library(cowplot)
p <- plot_grid(p_TCGA,p_GSE1,p_GSE2,ncol=3)

ggsave(p,file="Supp_Figure1A_scatter_plot_3-Cohort_CDRs_proteincodingGene.png",height=3,width=13)

```
<img src="/Github_Fig/Figure1/Supp_Figure1A_scatter_plot_3-Cohort_CDRs_proteincodingGene.png">

## Fig1B-C
```R


library(ggvenn)
# 样例数据
a <- list(`TCGA` = wide_means_TCGA$symbol,
          `GSE68` = wide_means_GSE6891$X,
          `GSE10` = wide_means_GSE10358$X)
#可视化绘制
library(ggvenn)

# 假设 a 是你的集合列表
# a <- list(Cohort1 = ..., Cohort2 = ..., Cohort3 = ...)

# 设置字体
opar <- par(family = "Roboto Condensed")  

# 绘图
p1 <- ggvenn(
  a,
  fill_color = c("#1b9e77", "#d95f02", "#7570b3"),   # 自定义填充颜色
  fill_alpha = 0.6,                                  # 透明度
  stroke_color = "black",                             # 边框颜色
  stroke_linetype = "longdash",                       # 边框线型
  stroke_size = 1,                                    # 边框宽度
  set_name_size = 10,                                 # 集合名称字体大小
  text_size = 6                                       # 重叠数字字体大小
)

ggsave(p1,file="Figure1B-Venn_3Cohort_CDRs_80_Gene.png",height=6,width=6)

gene_use <- intersect(wide_means_TCGA$symbol,wide_means_GSE6891$X)
gene_use <- intersect(gene_use,wide_means_GSE10358$X)
gene_use


wide_means_TCGA1 <-  wide_means_TCGA[wide_means_TCGA$symbol%in%c(gene_use),]
wide_means_GSE11 <-  wide_means_GSE6891[wide_means_GSE6891$X%in%c(gene_use),]
wide_means_GSE21 <-  wide_means_GSE10358[wide_means_GSE10358$X%in%c(gene_use),]

wide_means_TCGA1$log2FoldChange <- abs(wide_means_TCGA1$log2FoldChange)
wide_means_GSE11$logFC <- abs(wide_means_GSE11$logFC)
wide_means_GSE21$logFC <- abs(wide_means_GSE21$logFC)

wide_means_TCGA1$z11 <- (wide_means_TCGA1$log2FoldChange - min(wide_means_TCGA1$log2FoldChange))/(max(wide_means_TCGA1$log2FoldChange)-min(wide_means_TCGA1$log2FoldChange))
wide_means_GSE11$z21 <- (wide_means_GSE11$logFC - min(wide_means_GSE11$logFC))/(max(wide_means_GSE11$logFC)-min(wide_means_GSE11$logFC))
wide_means_GSE21$z31 <- (wide_means_GSE21$logFC - min(wide_means_GSE21$logFC))/(max(wide_means_GSE21$logFC)-min(wide_means_GSE21$logFC))


library(dplyr)
library(tidyr) 
colnames(wide_means_TCGA1)[8] <- "gene"
colnames(wide_means_GSE11)[1] <- "gene"
colnames(wide_means_GSE21)[1] <- "gene"
combined <- full_join(wide_means_TCGA1, wide_means_GSE11, by="gene") %>%
            full_join(wide_means_GSE21, by="gene") %>%
            mutate(z11=replace_na(z11,0), z21=replace_na(z21,0), z31=replace_na(z31,0),
                   combined_score = rowMeans(cbind(z11,z21,z31)))

combined[order(combined$combined_score),]



combined$combined_score <- combined$combined_score*(-1)
# 1. 按 log2FC 排序并创建 rank（纵轴）
combined <- combined %>%
  arrange(combined_score) %>%
  mutate(rank_ordered = row_number())

top5_genes <- combined %>%
  top_n(5, wt = abs(combined_score)) %>%  # 可以根据绝对 log2FC 排前5
  pull(gene)

# 3. 创建标记列
combined <- combined %>%
  mutate(Highlight = case_when(
    gene == "GALNT11" ~ "GALNT11",
    gene == "KMT2C" ~ "KMT2C",
    gene == "EZH2" ~ "EZH2",
    TRUE ~ "Other"
  ))

p <- ggplot(combined, aes(x = rank_ordered, y = combined_score)) +
  geom_point(aes(color = Highlight), size = 2) +
  geom_text(data = subset(combined, Highlight %in% c("KMT2C","GALNT11","EZH2")),
            aes(label = gene),
            vjust = -1, size = 4) +
  scale_color_manual(values = c("KMT2C" = "red", "GALNT11" = "red", "EZH2" = "red", "Other" = "blue")) +
  labs(title = "",
       x = "Rank (sorted by combined_score)",
       y = "combined_score") +
  theme_bw()

# 5. 显示图
ggsave(p,file="Figure1C_scatter_plot_3Cohort_CDRs.png",height=2.8,width=4)
```
<img src="/Github_Fig/Figure1/Figure1B-Venn_3Cohort_CDRs_80_Gene.png">
<img src="/Github_Fig/Figure1/Figure1C_scatter_plot_3Cohort_CDRs.png">


## Fig1D and Supplemental Fig1B
```R
setwd("~/home/syc/Data_Prepare/Figure1/")
Count <- read.csv("~/DATAbase/ALL_TCGA_DATA/RNA/Counts/TCGA/TCGA-LAML_RNA_expr_Counts.csv",check.names=F)
Count_tumor <- colnames(Count)[as.integer(substr(colnames(Count),14,15)) < 10]

Count_tmp1 <- Count[,colnames(Count)%in%c(Count_tumor[2:152])]
extracted_ids <- gsub("([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[0-9]+)[A-Z]+.*", "\\1", colnames(Count_tmp1))

colnames(Count_tmp1) <- extracted_ids

rownames(Count_tmp1) <- Count[,1]

library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
Count_tmp1_tmp <- Count_tmp1
Count_tmp1$symbol <- mapIds(x = anno_data,keys = rownames(Count_tmp1),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
Count_tmp2 <- Count_tmp1[Count_tmp1$symbol%in%c("GALNT11"),]
rownames(Count_tmp2) <- Count_tmp2$symbol
Count_tmp2 <- subset(Count_tmp2,select=-c(symbol))


Count_tmp3 <- data.frame(t(Count_tmp2))

df_trns <- Count_tmp3
mean(df_trns$GALNT11)
[1] 1308.722
median(df_trns$GALNT11)
[1] 1301

library("survival")
library("survminer")

df_trns <- data.frame(df_trns)
df_trns$Group_mean <- "GALNT11_low"
df_trns[which(df_trns$GALNT11 > 1308.722),'Group_mean'] <- 'GALNT11_high'
df_trns$Group_median <- "GALNT11_low"
df_trns[which(df_trns$GALNT11 > 1301),'Group_median'] <- 'GALNT11_high'
df_trns$Sample <- rownames(df_trns)
cilinc <- read.csv("~/DATAbase/ALL_TCGA_DATA/clinical_info/ALL_info_includ_RNA_DNA/TCGA/TCGA-LAML_clinical.csv")
cilinc[which(is.na(cilinc$days_to_death)),"days_to_death"] <- cilinc[which(is.na(cilinc$days_to_death)),"days_to_last_follow_up"]
cilinic <- subset(cilinc,select=c(submitter_id,vital_status,days_to_death))
cilinic$submitter_id <- paste(cilinic$submitter_id,"-03",sep="")
colnames(cilinic) <- c("Sample","Death","OS")
dat_new <- merge(df_trns,cilinic,by="Sample",all.x=TRUE)

#根据基因
dat_new$Death <- as.character(dat_new$Death)
dat_new$status <- ifelse(dat_new$Death=="Alive",0,1)
dat_new$OS <- as.numeric(dat_new$OS)
coxph_result <- coxph(formula = Surv(OS, status) ~ GALNT11, data = dat_new)
summary(coxph_result,data=dat_new)
ggforest(coxph_result, data =dat_new, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)

dat_new.cut <- surv_cutpoint(
   dat_new,
   time = "OS",
   event = "status", 
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat_new.cut)
p1 <- plot(dat_new.cut, "GALNT11",palette = "npg")

dat_new.cut.cat <- surv_categorize(dat_new.cut) 
fit <- survfit(Surv(OS, status) ~ GALNT11, data = dat_new.cut.cat)
p2 <- ggsurvplot(fit, data = dat_new.cut.cat,
surv.median.line = "hv",
color = NULL, 
#palette = "lancet", 
pval = TRUE,
risk.table=TRUE,
xlim=c(0,1500),
ggtheme = theme_pubr(),
palette = c("#DE5F4C","#74B5DE"),
break.x.by=500)
png(width=600,height=600,"GLNT11_SURVIVE_gene_TCGA.png")
p2
dev.off()

Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
studyID = "GSE37642"
library(GEOquery)
eSet <- getGEO(studyID, destdir = "./", getGPL = F)
pdata = pData(eSet[[1]])#表型信息，里面包含了需要的分组信息
pdata <- data.frame(GEO_ID=pdata$geo_accession,Sample_ID=pdata$title)

library(tidyverse)
library(GEOquery)
library(tinyarray)

Rawdata1=read.table('GSE37642-GPL570_series_matrix.txt.gz',
                   sep = '\t',quote ="",fill = T,
                   comment.char = "!",header = T)
Rawdata2=read.table('GSE37642-GPL96_series_matrix.txt.gz',
                   sep = '\t',quote ="",fill = T,
                   comment.char = "!",header = T)
Rawdata3=read.table('GSE37642-GPL97_series_matrix.txt.gz',
                   sep = '\t',quote ="",fill = T,
                   comment.char = "!",header = T)

clinic=read.table('GSE37642_Survival_data.txt',header=TRUE)

colnames(Rawdata1) <- gsub("X.","",colnames(Rawdata1))
colnames(Rawdata1) <- gsub("\\.","",colnames(Rawdata1))
colnames(Rawdata2) <- gsub("X.","",colnames(Rawdata2))
colnames(Rawdata2) <- gsub("\\.","",colnames(Rawdata2))
colnames(Rawdata3) <- gsub("X.","",colnames(Rawdata3))
colnames(Rawdata3) <- gsub("\\.","",colnames(Rawdata3))

library(data.table)
ids_GPL570 <- fread("~/home/syc/Regression/GPL570.annot")
head(ids_GPL570)

ids_GPL96 <- fread("~/home/syc/Regression/GPL96.annot")
head(ids_GPL96)

ids_GPL97 <- fread("~/home/syc/Regression/GPL97.annot")
head(ids_GPL97)

ids_GPL570[ids_GPL570$'Gene symbol'=="GALNT11",]
ids_GPL96[ids_GPL96$'Gene symbol'=="GALNT11",]
ids_GPL97[ids_GPL97$'Gene symbol'=="GALNT11",]

Rawdata1$ID_REF <- gsub("\"","",Rawdata1$ID_REF)
Rawdata2$ID_REF <- gsub("\"","",Rawdata2$ID_REF)
dat1 <- Rawdata1[Rawdata1$ID_REF=="219013_at",]
dat2 <- Rawdata2[Rawdata2$ID_REF=="219013_at",]
dat2 <- subset(dat2,select=-c(ID_REF))
dat_use <- cbind(dat1,dat2)
dat_use <- subset(dat_use,select=-c(ID_REF))
rownames(dat_use) <- "GALNT11"
dat_use <- data.frame(t(dat_use))
dat_use$SAMPLE <- rownames(dat_use)

data_final <- merge(dat_use,clinic,by="SAMPLE")
dat <- na.omit(data_final)

library(survival)
library(survminer)
library(dplyr)

dat <- dat %>%
  mutate(
    OS = as.numeric(OS),
    event = case_when(
      Death %in% c("dead", "Dead", "DEAD", 1, "1") ~ 1,
      Death %in% c("alive", "Alive", "ALIVE", 0, "0") ~ 0,
      TRUE ~ NA_real_
    ),
    GALNT11 = as.numeric(GALNT11)
  ) %>%
  filter(!is.na(OS), !is.na(event), !is.na(GALNT11), OS > 0)


library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS, event) ~ GALNT11, data = dat)
summary(coxph_result,data=dat)
ggforest(coxph_result, data =dat, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat.cut <- surv_cutpoint(
   dat,
   time = "OS",
   event = "event",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)

summary(dat.cut)
plot(dat.cut, "GALNT11")
dat.cut.cat <- surv_categorize(dat.cut) 
library(survival)
fit <- survfit(Surv(OS, event) ~ GALNT11, data = dat.cut.cat)
p1 <- ggsurvplot(fit, data = dat.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500, 
palette = c("#A52A2A","#4682B4"))

png(width=600,height=600,"GLNT11_SURVIVE_GSE37642.png")
p1
dev.off()

```
<img src="/Github_Fig/Figure1/GLNT11_SURVIVE_gene_TCGA.png">
<img src="/Github_Fig/Figure1/GLNT11_SURVIVE_GSE37642.png">



## Fig3A-B
```R

**********************7q LOSS TMB****************************


setwd("~/home/syc/Data_Prepare/Figure3")
/usr/local/R4.2/bin/R
library(ggplot2)
library(ggpubr)
library(trqwe)
library(maftools)
PAAD_maf <- mcreadRDS("~/home/pym/human_data/gdc_download_20230213_095040.200153/all_maf/TCGA_LAML_MAF.rds",mc.cores=20)
PAAD_maf <- read.maf(PAAD_maf)

#TMB计算
aml_tmb <-  tmb(maf = PAAD_maf) #131
aml_tmb <- as.data.frame(aml_tmb)
aml_tmb$Sample <- gsub("-",".",aml_tmb$Tumor_Sample_Barcode)
aml_tmb$Sample <- substr(aml_tmb$Sample,1,15)
aml_tmb <- aml_tmb[!duplicated(aml_tmb$Sample),]  #124

Alter_tcga <- read.table("~/home/syc/Galnt11_DXT/alterations_across_samples.tsv",sep="\t",header=T)
#定义7q
Status_7q <- subset(Alter_tcga,select=c(Altered,Sample.ID,GALNT11..HETLOSS.HOMDEL.MUT.AMP.GAIN.FUSION,ABCB8..HETLOSS.HOMDEL.MUT.AMP.GAIN.FUSION))
Status_7q$Group <- "Diploid"
Status_7q[which(Status_7q$Altered==1),"Group"] <- "Loss"
Status_7q[which(Status_7q$Altered==0),"Group"] <- "Diploid"

rownames(Status_7q) <- Status_7q$Sample.ID

table(Status_7q$Group)

Diploid    Loss
    165      25
aml_tmb$Sample <- gsub("\\.","-",aml_tmb$Sample)
dat <- merge(Status_7q,aml_tmb,by.x="Sample.ID",by.y="Sample",all=TRUE)
dat_new <- subset(dat,select=c(Group,total_perMB))

dat_new <- na.omit(dat_new)

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(Seurat)
dat_new$Group <- factor(dat_new$Group,levels=c("Diploid","Loss"))
dat_new1 <- dat_new[which(dat_new$total_perMB<1),]

p  <-  ggplot(dat_new1, aes(x=Group, y=total_perMB,fill = Group)) +
  geom_boxplot()+RotatedAxis() + theme_classic()
p


dat_new1$Group <- as.factor(dat_new1$Group)
bino_model <- glm(dat_new1$Group ~ dat_new1$total_perMB, family="binomial")

anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 0.01294755

**********************7q LOSS ANEUPLOIDY_SCORE****************************


tcga_clinical_sample1 <- read.table("~/home/pym/human_data/laml_tcga_pan_can_atlas_2018/data_clinical_sample.txt",sep="\t",header=T)
dat <- merge(Status_7q,tcga_clinical_sample1,by.x="Sample.ID",by.y="SAMPLE_ID",all=TRUE)
dat_new <- subset(dat,select=c(Group,ANEUPLOIDY_SCORE))

dat_new <- na.omit(dat_new)

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(Seurat)
dat_new$Group <- factor(dat_new$Group,levels=c("Diploid","Loss"))
dat_new1 <- dat_new[which(dat_new$ANEUPLOIDY_SCORE<10),]

p1  <-  ggplot(dat_new1, aes(x=Group, y=ANEUPLOIDY_SCORE,fill = Group)) +
  geom_boxplot()+RotatedAxis() + theme_classic()
p1

g <- p+p1
ggsave(g,file="Figure3A-B-total_perMB_in7q_ANEUPLOID_SCORE.png",height=3,width=5)


dat_new1$Group <- as.factor(dat_new1$Group)
bino_model <- glm(dat_new1$Group ~ dat_new1$ANEUPLOIDY_SCORE, family="binomial")
 anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 3.843946e-05
```
<img src="/Github_Fig/Figure3/Figure3A-B-total_perMB_in7q_ANEUPLOID_SCORE.png">

## Fig3C
```R
degs <- read.csv("~/home/syc/Galnt11_DXT/TCGA/TCGA-Bulk_7Q_DEG.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
res$symbol <- mapIds(x = anno_data,keys = rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res$entrez <- mapIds(x = anno_data,keys = rownames(res),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res$symbol
AA <- as.character(AA)
res$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res$GENENAME <- as.character(res$GENENAME)


###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- res
P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

geneset <- c("REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="Figure3C-GSEA_TCGA-AML_7q_vs_Diploid2_DNA_DAMAGEV1007_2.png",height=4,width=4.5)
```
<img src="/Github_Fig/Figure3/Figure3C-GSEA_TCGA-AML_7q_vs_Diploid2_DNA_DAMAGEV1007_2.png">
<img src="/Github_Fig/Figure3/Figure3C-GSEA_TCGA-AML_7q_vs_Diploid2_DNA_DAMAGEV1007_2.png">


## Fig3D-E
```R
**********************GALNT11 Expr TMB****************************

$$$$$$$$ GALNT11 Exp vs TMB $$$$$$$$ 分成high low
$$$$$$$$ GALNT11 Exp vs TMB $$$$$$$$ 分成high low

library(ggplot2)
library(ggpubr)
library(tidyr)
library(trqwe)
library(maftools)
Count <- read.csv("~/DATAbase/ALL_TCGA_DATA/RNA/Counts/TCGA/TCGA-LAML_RNA_expr_Counts.csv",check.names=F)
Count_tumor <- colnames(Count)[as.integer(substr(colnames(Count),14,15)) < 10]

Count_tmp1 <- Count[,colnames(Count)%in%c(Count_tumor[2:152])]
extracted_ids <- gsub("([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[0-9]+)[A-Z]+.*", "\\1", colnames(Count_tmp1))

colnames(Count_tmp1) <- extracted_ids

rownames(Count_tmp1) <- Count[,1]

library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
Count_tmp1_tmp <- Count_tmp1
Count_tmp1$symbol <- mapIds(x = anno_data,keys = rownames(Count_tmp1),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
Count_tmp2 <- Count_tmp1[Count_tmp1$symbol%in%c("GALNT11"),]
Count_tmp2 <- subset(Count_tmp2,select=-c(symbol))
Count_tmp3 <- data.frame(t(Count_tmp2))
colnames(Count_tmp3) <- "GALNT11"
Count_tmp3$sample <- rownames(Count_tmp3)

library(ggplot2)
library(ggpubr)
library(tidyr)
library(trqwe)
library(maftools)
PAAD_maf <- mcreadRDS("~/home/pym/human_data/gdc_download_20230213_095040.200153/all_maf/TCGA_LAML_MAF.rds",mc.cores=20)
PAAD_maf <- read.maf(PAAD_maf)

#TMB计算
aml_tmb <-  tmb(maf = PAAD_maf) #131
aml_tmb <- as.data.frame(aml_tmb)
aml_tmb$sample <- substr(aml_tmb$Tumor_Sample_Barcode,1,15)
aml_tmb <- aml_tmb[!duplicated(aml_tmb$sample),]  #124

intersect_gene <- intersect(Count_tmp3$sample,aml_tmb$sample)


Count_tmp4 <- Count_tmp3[Count_tmp3$sample%in%c(intersect_gene),]
aml_tmb <- aml_tmb[aml_tmb$sample%in%c(intersect_gene),]

Count_tmp4 <- Count_tmp4[order(Count_tmp4$GALNT11),]
Count_tmp4_lo <- head(Count_tmp4,10)#0.1
Count_tmp4_high <- tail(Count_tmp4,10)#0.1
Count_tmp4_lo$Group <- "low"
Count_tmp4_high$Group <- "high"

galnt11_tmb_exp1 <- rbind(Count_tmp4_lo,Count_tmp4_high)

galnt11_tmb_exp2 <- merge(aml_tmb,galnt11_tmb_exp1,by='sample')

library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(Seurat)
galnt11_tmb_exp2$Group <- factor(galnt11_tmb_exp2$Group,levels=c("high","low"))
dat1 <- galnt11_tmb_exp2[galnt11_tmb_exp2$total_perMB<1,]
p1  <-  ggplot(dat1, aes(x=Group, y=total_perMB,fill = Group)) +
  geom_boxplot()+RotatedAxis() + theme_classic()

dat1$Group <- as.factor(dat1$Group)
bino_model <- glm(dat1$Group ~ dat1$total_perMB, family="binomial")
anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 0.05420457


**********************GALNT11 Expr ANEUPLOIDY_SCORE****************************

tcga_clinical_sample1 <- read.table("~/home/pym/human_data/laml_tcga_pan_can_atlas_2018/data_clinical_sample.txt",sep="\t",header=T)
dat11 <- merge(dat1,tcga_clinical_sample1,by.x="sample",by.y="SAMPLE_ID")

p2  <-  ggplot(dat11, aes(x=Group, y=ANEUPLOIDY_SCORE,fill = Group)) +
  geom_boxplot()+ RotatedAxis() + theme_classic()
p2

g <- p1+p2
ggsave(g,file="Figure3E-F-TCGA-GALNT11_ANEUPLOID_SCORE.png",height=3,width=5)


dat11$Group <- as.factor(dat11$Group)
bino_model <- glm(dat11$Group ~ dat11$ANEUPLOIDY_SCORE, family="binomial")

anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 0.0284016

```
<img src="/Github_Fig/Figure3/Figure3E-F-TCGA-GALNT11_ANEUPLOID_SCORE.png">


# Figure3F
```R

########分成high和low进行分析
Count <- read.csv("~/DATAbase/ALL_TCGA_DATA/RNA/Counts/TCGA/TCGA-LAML_RNA_expr_Counts.csv",check.names=F)
Count_tumor <- colnames(Count)[as.integer(substr(colnames(Count),14,15)) < 10]

Count_tmp1 <- Count[,colnames(Count)%in%c(Count_tumor[2:152])]
extracted_ids <- gsub("([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[0-9]+)[A-Z]+.*", "\\1", colnames(Count_tmp1))

colnames(Count_tmp1) <- extracted_ids

rownames(Count_tmp1) <- Count[,1]

library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
Count_tmp1_tmp <- Count_tmp1
all_expr <- Count_tmp1_tmp
all_expr <- all_expr[,colnames(all_expr)%in%c(galnt11_tmb_exp1$sample)]

dat_new1 <- all_expr[,match(galnt11_tmb_exp1$sample,colnames(all_expr))]

idx <- match(galnt11_tmb_exp1$sample, colnames(all_expr))  # 用样本表去匹配列名
ok <- !is.na(idx)
dat_new1 <- all_expr[, idx[ok], drop = FALSE]              # 重排/筛选列
meta1 <- galnt11_tmb_exp1[ok, , drop = FALSE]           # 同步重排临床表

library(DESeq2)
colData <- data.frame(row.names = colnames(dat_new1),group_list= galnt11_tmb_exp1$Group)
colData


dds <- DESeqDataSetFromMatrix(countData = dat_new1,colData = colData,design = ~ group_list)
vsd <- vst(dds, blind = FALSE)

table(galnt11_tmb_exp1$Group)

high  low
  10   10


count <- assay(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
count_and_nor <- cbind(count,normalized_counts)
write.csv(count_and_nor,file="7Q_count_and_normalized_data_GALNT11_high_1007.csv")



dds <- DESeq(dds)
#进行差异基因分析
resultsNames(dds)
res <-  results(dds, contrast=c("group_list","low","high"))
head(res)
write.csv(res, file="Bulk_7Q_DEG_Galnt11_1007.csv")

degs <- read.csv("~/home/syc/Galnt11_DXT/GALNT11_expr/Bulk_7Q_DEG_Galnt11_1007.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
res$symbol <- mapIds(x = anno_data,keys = rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res$entrez <- mapIds(x = anno_data,keys = rownames(res),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res$symbol
AA <- as.character(AA)
res$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res$GENENAME <- as.character(res$GENENAME)


###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- res
P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-TCGA-AML7GALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

geneset <- c("HALLMARK_DNA_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="Figure3G-GSEA_TCGA-galnt11_low_vs_high_DNA_DAMAGEV1007.png",height=4,width=4.5)

# geneset <- c("ATM_TARGET_GENES")
# p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
# p
# p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
# p
# gsea_df[gsea_df$ID%in%c("ATM_TARGET_GENES"),]
# ggsave(p,file="Figure4D-GSEA_TCGA-galnt11_low_vs_high_ATM.png",height=4,width=4.5)

```
<img src="/Github_Fig/Figure3/Figure3G-GSEA_TCGA-galnt11_low_vs_high_DNA_DAMAGEV1007.png">



# Figure3G and Supplemental Figure4A
```R
library(tidyverse)
library(GEOquery)
library(tinyarray)
library(limma)
library(edgeR)

beat.AML.RNA.raw.counts <- read.csv(row.names=1,'~/home/pym/human_data/aml_ohsu_2022/BeatAML/beataml_waves1to4_counts_dbgap.AML.samples.csv')
exprSet2_GALNT11 <- data.frame(t(beat.AML.RNA.raw.counts[rownames(beat.AML.RNA.raw.counts)%in%c("GALNT11"),]))
exprSet2_GALNT11$sample <- rownames(exprSet2_GALNT11)
exprSet2_GALNT11 <- exprSet2_GALNT11[order(exprSet2_GALNT11$GALNT11,decreasing=FALSE),]
nrow(exprSet2_GALNT11) #671
exprSet2_GALNT11_lo <- head(exprSet2_GALNT11,67) #0.1
exprSet2_GALNT11_hi <- tail(exprSet2_GALNT11,67) #0.1

exprSet2_GALNT11_lo$Group <- "low"
exprSet2_GALNT11_hi$Group <- "high"

group_list <- rbind(exprSet2_GALNT11_lo,exprSet2_GALNT11_hi)
group_list <- group_list[order(group_list$Group),]
exprSet1_use <- beat.AML.RNA.raw.counts[,colnames(beat.AML.RNA.raw.counts)%in%c(group_list$sample)]
exprSet1_use_log <- log2(exprSet1_use+1)
expr_matrix <- as.matrix(exprSet1_use_log[, group_list$sample])

meta_info <- data.frame(sample=c(colnames(expr_matrix)),
  group=c(group_list$Group)
  )
design <- model.matrix(~0 + factor(meta_info$group))
colnames(design) <- c("high", "low")

# limma 分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(SevenQ_vs_NN = low - high, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 输出差异表达基因
DEG_limma <- topTable(fit2, number = Inf, adjust.method = "BH")
head(DEG_limma)

write.csv(DEG_limma,file="beatAML-DEGs_GALNT11_low_vs_high.csv")
###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-beatAML-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

geneset <- c("REACTOME_DNA_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="beatAML-galnt11_low_vs_high_DNA_DAMAGEV1007.png",height=4,width=4.5)

geneset <- c("ATM_DN.V1_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="beatAML-galnt11_low_vs_high_ATM.png",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]

```
<img src="/Github_Fig/Figure3/beatAML-galnt11_low_vs_high_DNA_DAMAGEV1007.png">
<img src="/Github_Fig/Figure4/beatAML-galnt11_low_vs_high_ATM.png">


# Figure3H and Figure4B
```R

GSE10358_exprSet1 <- read.csv("~/home/syc/Public_Data/GSE10358_expr_tumor_anno.csv")
library(edgeR)
rownames(GSE10358_exprSet1) <- GSE10358_exprSet1$X
exprSet2 <- subset(GSE10358_exprSet1,select=-c(X))
exprSet2_GALNT11 <- data.frame(t(exprSet2[rownames(exprSet2)%in%c("GALNT11"),]))
exprSet2_GALNT11$sample <- rownames(exprSet2_GALNT11)
exprSet2_GALNT11 <- exprSet2_GALNT11[order(exprSet2_GALNT11$GALNT11,decreasing=FALSE),]
nrow(exprSet2_GALNT11) #98
exprSet2_GALNT11_lo <- head(exprSet2_GALNT11,10) #0.1
exprSet2_GALNT11_hi <- tail(exprSet2_GALNT11,10) #0.1

exprSet2_GALNT11_lo$Group <- "low"
exprSet2_GALNT11_hi$Group <- "high"

group_list <- rbind(exprSet2_GALNT11_lo,exprSet2_GALNT11_hi)
group_list <- group_list[order(group_list$Group),]

exprSet1_use <- exprSet2[,colnames(exprSet2)%in%c(group_list$sample)]
exprSet1_use_log <- log2(exprSet1_use+1)
expr_matrix <- as.matrix(exprSet1_use_log[, group_list$sample])

meta_info <- data.frame(sample=c(colnames(expr_matrix)),
  group=c(group_list$Group)
  )
design <- model.matrix(~0 + factor(meta_info$group))
colnames(design) <- c("high", "low")

# limma 分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(SevenQ_vs_NN = low - high, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 输出差异表达基因
DEG_limma <- topTable(fit2, number = Inf, adjust.method = "BH")
head(DEG_limma)
write.csv(DEG_limma,file="GSE10358_DEGs_GALNT11_low_vs_high.csv")

###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-GSE10358-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
geneset <- c("REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="GSE10358-galnt11_low_vs_high_DNA_DAMAGEV2.png",height=4,width=4.5)

geneset <- c("PID_ATM_PATHWAY")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
gsea_df[gsea_df$ID%in%c("PID_ATM_PATHWAY"),]
ggsave(p,file="GSE10358-galnt11_low_vs_high_ATM.png",height=4,width=4.5)

```
<img src="/Github_Fig/Figure3/GSE10358-galnt11_low_vs_high_DNA_DAMAGEV2.png">
<img src="/Github_Fig/Figure4/GSE10358-galnt11_low_vs_high_ATM.png">


# Figure3I
```R
library(tidyverse)
library(GEOquery)
library(tinyarray)
DEG_limma <- read.csv("GSE37642-_DEGs_GALNT11_low_vs_high.csv")
###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)
custom_GSEA_GMT <- read.gmt("~/2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH')  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-GSE37642-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
geneset <- c("HALLMARK_DNA_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="GSE37642-galnt11_low_vs_high_DNA_DAMAGE.png",height=4,width=4.5)

# geneset <- c("PID_ATM_PATHWAY")
# p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
# p
# p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
# p
# ggsave(p,file="GSE37642-galnt11_low_vs_high_ATM.png",height=4,width=4.5)

```
<img src="/Github_Fig/Figure3/GSE37642-galnt11_low_vs_high_DNA_DAMAGE.png">


# Figure4A
```R
set.seed(100)
DEG_limma <- read.csv("~/home/syc/Public_Data/GSE10358_DEGs_7q_vs_NN.csv")
DEG_limma$group <- ifelse(abs(DEG_limma$logFC) > 0 ,ifelse(DEG_limma$logFC > 0,"up","down"),"ns")


#GSEA富
library(clusterProfiler)
library(trqwe)
library(enrichplot)

options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
DEG_limma$symbol <- DEG_limma$X
P_subcluster <- DEG_limma
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
geneset <- c("ATM_DN.V1_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="FIG4E-GSEA_GSE10358_7q_vs_Diploid2_ATM_0916.png",height=4,width=4.5)

```
<img src="/Github_Fig/Figure4/GSEA_GSE10358_7q_vs_Diploid2_ATM_0916">


# Figure4C
```R
#


# 需要的包
library(data.table)
library(stringr)

files <- Sys.glob("GSM49234*_Rawcounts.B16_*.txt.gz")
stopifnot(length(files) == 6)

read_one_counts <- function(f){
  dt <- fread(f)  # 自动解压 .gz
  # 选基因列（优先找字符列），选计数列（优先找数值列）
  char_cols <- which(vapply(dt, is.character, logical(1)))
  num_cols  <- which(vapply(dt, is.numeric,   logical(1)))
  if (length(char_cols) == 0) char_cols <- 1L
  if (length(num_cols)  == 0) num_cols  <- setdiff(seq_along(dt), char_cols)[1]
  gene_col  <- char_cols[1]
  count_col <- num_cols[1]
  dt <- dt[, .(gene = as.character(.SD[[1]]),
               count = as.numeric(.SD[[2]])),
           .SDcols = c(gene_col, count_col)]
  dt <- dt[!is.na(gene) & gene != "" & !is.na(count)]
  dt <- dt[, .(count = sum(count, na.rm = TRUE)), by = gene]
  sample <- sub(".*Rawcounts\\.(.*)\\.txt.*$", "\\1", basename(f))
  setnames(dt, "count", sample)
  dt[]
}

# 3) 读取并合并到一个矩阵
lst <- lapply(files, read_one_counts)
counts_dt <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), lst)
counts_dt[is.na(counts_dt)] <- 0

# 转成基因×样本矩阵
counts_mat <- as.data.frame(counts_dt)
rownames(counts_mat) <- counts_mat$gene
counts_mat$gene <- NULL
counts_mat <- as.matrix(counts_mat)
storage.mode(counts_mat) <- "integer"

# 4) 构建分组信息（VC vs ATMKO）
samples <- colnames(counts_mat)
group   <- sub("^B16_([^_]+)_\\d+$", "\\1", samples)   # 取 VC / ATMKO
coldata <- data.frame(sample = samples,
                      group  = factor(group, levels = c("VC","ATMKO")),
                      row.names = samples,
                      check.names = FALSE)

# 5) 快速检查与导出
cat("Matrix dimension: ", paste(dim(counts_mat), collapse = " x "), "\n")
print(head(rownames(counts_mat), 5))
print(colSums(counts_mat))

write.table(counts_mat, "counts_matrix.tsv", sep = "\t",
            quote = FALSE, col.names = NA)
write.csv(coldata, "sample_metadata.csv", row.names = TRUE)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ] #低表达基因过滤
dds <- DESeq(dds)
res <- results(dds, contrast = c("group","ATMKO","VC"))

#基因注释
library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_id_version") %>%
  mutate(ENSEMBL = sub("\\..*$", "", ensembl_id_version))  # 去掉 .版本号

map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = res_df$ENSEMBL,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME", "ENTREZID")
)

map_dedup <- map %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL   = dplyr::first(na.omit(SYMBOL)),
    GENENAME = dplyr::first(na.omit(GENENAME)),
    ENTREZID = dplyr::first(na.omit(ENTREZID)),
    .groups = "drop"
  )

res_anno <- res_df %>%
  left_join(map_dedup, by = "ENSEMBL") %>%
  relocate(ENSEMBL, SYMBOL, GENENAME, ENTREZID)

write.csv(res_anno, "DESeq2_ATMKO_vs_VC_annotated.csv", row.names = FALSE)
res_anno2 <- res_anno
res_anno2 <- res_anno2[,c("log2FoldChange","SYMBOL")]
res_anno2 <- na.omit(res_anno2)
custom_GSEA_GMT <- data.frame(gene=head(res_anno2[order(res_anno2$log2FoldChange),"SYMBOL"],100))
custom_GSEA_GMT$term <-"ATMKO_DN"
setwd("~/home/syc/Data_Prepare/Figure4")
library(clusterProfiler)
library(trqwe)
library(enrichplot)
res_1 <- read.csv("~/home/syc/Galnt11_DXT/Glant11_bulk_ckit/PXY_3_shG494_VS_shRen_result.csv")
colnames(res_1) <- c("symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")

options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- res_1
P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- custom_GSEA_GMT[,c("term","gene")]
write.csv(custom_GSEA_GMT,"ATM_KO_DN_GSE161922.csv")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_AM
 $ enrichmentScore: num -0.467
 $ NES            : num -1.4
 $ pvalue         : num 0.0281
 $ p.adjust       : num 0.0281


library(ggpubr)
library(ggplot2)

geneset <- c("ATMKO_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3','#459943'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="Figure4D-GSEA_shGalnt11_vs_shRen_ATM.png",height=4,width=4.5)
```
<img src="/Github_Fig/Figure4/Figure4D-GSEA_shGalnt11_vs_shRen_ATM.png">


# Supplemental Fig4B
```R

library(tidyverse)
library(GEOquery)
library(tinyarray)
library(limma)
library(edgeR)

exprSet1 <- read.csv("~/home/syc/Public_Data/GSE6891expr_tumor_anno.csv")
rownames(exprSet1) <- exprSet1$X
exprSet2 <- subset(exprSet1,select=-c(X))
exprSet2_GALNT11 <- data.frame(t(exprSet2[rownames(exprSet2)%in%c("GALNT11"),]))
exprSet2_GALNT11$sample <- rownames(exprSet2_GALNT11)
exprSet2_GALNT11 <- exprSet2_GALNT11[order(exprSet2_GALNT11$GALNT11,decreasing=FALSE),]
nrow(exprSet2_GALNT11) #217
exprSet2_GALNT11_lo <- head(exprSet2_GALNT11,12) #0.1
exprSet2_GALNT11_hi <- tail(exprSet2_GALNT11,12) #0.1

exprSet2_GALNT11_lo$Group <- "low"
exprSet2_GALNT11_hi$Group <- "high"

group_list <- rbind(exprSet2_GALNT11_lo,exprSet2_GALNT11_hi)
group_list <- group_list[order(group_list$Group),]

exprSet1_use <- exprSet2[,colnames(exprSet2)%in%c(group_list$sample)]
expr_matrix <- as.matrix(exprSet1_use[, group_list$sample])

meta_info <- data.frame(sample=c(colnames(expr_matrix)),
  group=c(group_list$Group)
  )
design <- model.matrix(~0 + factor(meta_info$group))
colnames(design) <- c("high", "low")

# limma 分析
fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(SevenQ_vs_NN = low - high, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 输出差异表达基因
DEG_limma <- topTable(fit2, number = Inf, adjust.method = "BH")
head(DEG_limma)
write.csv(DEG_limma,file="GSE6891_DEGs_GALNT11_low_vs_high.csv")


###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("~/2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-GSE6891-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

# geneset <- c("GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE")
# p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
# p
# p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
# p
# ggsave(p,file="GSE6891-galnt11_low_vs_high_DNA_DAMAGEV1007.png",height=4,width=4.5)

geneset <- c("ATM_DN.V1_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
gsea_df[gsea_df$ID%in%c("ATM_DN.V1_DN"),]
ggsave(p,file="GSE6891-galnt11_low_vs_high_ATM.png",height=4,width=4.5)
```

<img src="/Github_Fig/Figure4/GSE6891-galnt11_low_vs_high_ATM.png">



# Figure5

```R
DEG_limma <- read.csv("~/home/syc/Public_Data/GSE10358_DEGs_7q_vs_NN.csv")
DEG_limma$group <- ifelse(abs(DEG_limma$logFC) > 0 ,ifelse(DEG_limma$logFC > 0,"up","down"),"ns")

#GSEA富
library(clusterProfiler)
library(trqwe)
library(enrichplot)

options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
DEG_limma$symbol <- DEG_limma$X
P_subcluster <- DEG_limma
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)
res_sig <- read.csv("~/home/syc/Regression/GSE103424_DEGs_CR_vs_NonCR.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=FALSE),"X"],200))
custom_GSEA_GMT$term <- "resistance_signature"
custom_GSEA_GMT <- custom_GSEA_GMT[,c("term","gene")]
#write.csv(custom_GSEA_GMT,file="custom_GSEA_GMT_GSE103424_Resistance_sig.csv")

gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM
geneset <- c("resistance_signature")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="SuppleFigure5A-GSEA_GSE10358_7q_vs_Diploid2_Resistance_gse103424_0916.png",height=4,width=4.5)


DEG_limma <- read.csv("~/home/syc/Public_Data/GSE6891_DEGs_7q_vs_NN.csv")
DEG_limma$group <- ifelse(abs(DEG_limma$logFC) > 0 ,ifelse(DEG_limma$logFC > 0,"up","down"),"ns")

#GSEA富
library(clusterProfiler)
library(trqwe)
library(enrichplot)

options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
DEG_limma$symbol <- DEG_limma$X
P_subcluster <- DEG_limma
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)
res_sig <- read.csv("~/home/syc//Regression/GSE103424_DEGs_CR_vs_NonCR.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=FALSE),"X"],200))
custom_GSEA_GMT$term <- "resistance_signature"
library(homologene)
custom_GSEA_GMT <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM
geneset <- c("resistance_signature")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="SuppleFig5B-GSEA_GSE6891_7q_vs_Diploid2_Resistance_GSE103424_0916.png",height=4,width=4.5)


degs <- read.csv("~/home/syc/Galnt11_DXT/TCGA/Bulk_7Q_DEG.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
res$symbol <- mapIds(x = anno_data,keys = rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res$entrez <- mapIds(x = anno_data,keys = rownames(res),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res$symbol
AA <- as.character(AA)
res$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res$GENENAME <- as.character(res$GENENAME)
###GSEA
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- res
P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

res_sig <- read.csv("~/home/syc//Regression/GSE83533_DEGs_relapse_vs_diagnosis.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Relapse_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM

geneset <- c("Relapse_signature")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3','#459943'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="SuppleFig5C-GSEA_TCGA-AML_7Q_VS_DIPLOID_Resistance_signature0916_GSE83533.png",height=4,width=4.5)


degs <- read.csv("~/home/syc/Galnt11_DXT/GALNT11_expr/Bulk_7Q_DEG_Galnt11_1007.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
res$symbol <- mapIds(x = anno_data,keys = rownames(res),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res$entrez <- mapIds(x = anno_data,keys = rownames(res),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res$symbol
AA <- as.character(AA)
res$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res$GENENAME <- as.character(res$GENENAME)


###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- res
P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)



@@@@@@@@@@@@@@@@@@@@Resistance

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
DEG_limma <- read.csv("~/home/syc/GALNT11_expr_Public_Data/GSE10358_DEGs_GALNT11_low_vs_high.csv")
colnames(DEG_limma)[1] <- "symbol"
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

res_sig <- read.csv("~/home/syc/Regression/GSE103424_DEGs_CR_vs_NonCR.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=FALSE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM
gsea_df <- gsea_AM@result
geneset <- c("Resistance_signature")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
gsea_df[gsea_df$ID%in%c("Resistance_signature"),]
ggsave(p,file="TCGA-AML-galnt11_low_vs_high_Resistance.png",height=4,width=4.5)



setwd("~/home/syc/Data_Prepare/Figure5")
DEG_limma <- read.csv("~/home/syc/Public_Data/GSE106291_DEGs_res_vs_ses.csv")

#GSEA富
library(clusterProfiler)
library(trqwe)
library(enrichplot)

options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
P_subcluster <- P_subcluster[,c("logFC","X")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$X
geneList = sort(aa,decreasing = T)

res_sig <- read.csv("~/home/syc/Galnt11_DXT/Glant11_bulk_ckit/PXY_3_shG494_VS_shRen_result.csv")
colnames(res_sig) <- c("symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
res_sig$ccscore <- res_sig$baseMean*res_sig$log2FoldChange^3
res_sig1 <- na.omit(res_sig[res_sig$pvalue<0.05,])
res_sig1 <- res_sig1[res_sig1$log2FoldChange>0,]
#custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$pvalue),"symbol"],100))
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$ccscore,decreasing=TRUE),"symbol"],100))

custom_GSEA_GMT$term <- "shGalnt11_vs_shRen_Up"

library(homologene)
gene1 <- mouse2human(custom_GSEA_GMT$gene)$humanGene
custom_GSEA_GMT1 <- data.frame(gene=gene1)
custom_GSEA_GMT1$term <- "shGalnt11_vs_shRen_Up"
custom_GSEA_GMT1 <- custom_GSEA_GMT1[,c("term","gene")]
#write.csv(custom_GSEA_GMT1,file="custom_GSEA_GMT1_shGALANT11_UP_new.csv")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM

geneset <- c("shGalnt11_vs_shRen_Up")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3','#459943'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="Figure5B-GSEA_GSEA_GSE106291.png",height=4,width=4.5)

```
<img src="/Github_Fig/Figure5/SuppleFigure5A-GSEA_GSE10358_7q_vs_Diploid2_Resistance_gse103424_0916.png">
<img src="/Github_Fig/Figure5/SuppleFig5B-GSEA_GSE6891_7q_vs_Diploid2_Resistance_GSE103424_0916.png">
<img src="/Github_Fig/Figure5/SuppleFig5C-GSEA_TCGA-AML_7Q_VS_DIPLOID_Resistance_signature0916_GSE103424.png">
<img src="/Github_Fig/Figure5/TCGA-AML-galnt11_low_vs_high_Resistance.png">
<img src="/Github_Fig/Figure5/Figure5B-GSEA_GSEA_GSE106291.png">

