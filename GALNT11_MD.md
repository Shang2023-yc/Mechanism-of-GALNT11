# Galnt11 数据整理,更新20260306
## Fig1A
```R

cd ~/home/syc/Data_Prepare

setwd("./Data_Prepare_update/Figure1/")
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
<img src="/Data_Prepare_update/Figure1/Fig1A-TCGA-LAML_CNV_subset.png">



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
<img src="/Data_Prepare_update/Figure1/Supp_Figure1A_scatter_plot_3-Cohort_CDRs_proteincodingGene.png">

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
<img src="/Data_Prepare_update/Figure1/Figure1B-Venn_3Cohort_CDRs_80_Gene.png">
<img src="/Data_Prepare_update/Figure1/Figure1C_scatter_plot_3Cohort_CDRs.png">


## Fig1D-G
```R
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


###TCGA-AML
setwd("./Data_Prepare_update/Figure1")

library(clusterProfiler)
library(org.Hs.eg.db)
dlbcl_rna <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/RNA/FPKM/TCGA/TCGA-LAML_RNA_expr.csv")
rownames(dlbcl_rna) <- dlbcl_rna$X
dlbcl_rna$symbol <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
dlbcl_rna$entrez <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- dlbcl_rna$symbol
AA <- as.character(AA)
dlbcl_rna$GENENAME <- mapIds(x = org.Hs.eg.db,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
dlbcl_rna$GENENAME <- as.character(dlbcl_rna$GENENAME)
dlbcl_rna <- dlbcl_rna[,2:(ncol(dlbcl_rna)-2)]
colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)] <- substr(colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)],1,16)
gene_list <- c("GALNT11")
length(unique(gene_list))

dlbcl_rna_chr17 <- dlbcl_rna[dlbcl_rna$symbol %in% gene_list,]
rownames(dlbcl_rna_chr17)<-dlbcl_rna_chr17$symbol
dlbcl_rna_chr17 <- dlbcl_rna_chr17[,1:(ncol(dlbcl_rna)-1)]
dat_new <- data.frame(t(dlbcl_rna_chr17))
dat_new$sample <- substr(rownames(dat_new), 1, 15)


cilinc <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/clinical_info/ALL_info_includ_RNA_DNA/TCGA/TCGA-LAML_clinical.csv")
cilinc[which(is.na(cilinc$days_to_death)),"days_to_death"] <- cilinc[which(is.na(cilinc$days_to_death)),"days_to_last_follow_up"]

cilinic <- subset(cilinc,select=c(submitter_id,vital_status,days_to_death))
cilinic$submitter_id <- paste(cilinic$submitter_id,"-03",sep="")
colnames(cilinic) <- c("Sample","Death","OS")
cilinic$sample <- gsub("-","\\.",cilinic$Sample)

df_trns <- merge(dat_new,cilinic,by="sample",all.x=TRUE)

library("survival")
library("survminer")
df_trns$Death <- as.character(df_trns$Death)
df_trns$status <- ifelse(df_trns$Death=="Alive",0,1)
df_trns$OS <- as.numeric(df_trns$OS)

coxph_result <- coxph(formula = Surv(OS, status) ~ GALNT11, data = df_trns)
summary(coxph_result,data=df_trns)
ggforest(coxph_result, data =df_trns, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)

df_trns.cut <- surv_cutpoint(
   df_trns,
   time = "OS",
   event = "status", 
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(df_trns.cut)
p1 <- plot(df_trns.cut, "GALNT11",palette = "npg")

df_trns.cut.cat <- surv_categorize(df_trns.cut) 
fit <- survfit(Surv(OS, status) ~ GALNT11, data = df_trns.cut.cat)
p2 <- ggsurvplot(fit, data = df_trns.cut.cat,
surv.median.line = "hv",
color = NULL, 
palette = "lancet", 
pval = TRUE,
risk.table=TRUE,
xlim=c(0,1500),
break.x.by=500)

png("Figure1E_TCGA_GLNT11_SURVIVE_gene_0227.png")
p2
dev.off()


<<<<<<<<<<<<<<<<<<<<提取非7q样本

alter_tcga <- read.table("./human_data/AML_samples_TCGA.tsv",sep="\t",header=T)
ABCB8_deletion_sample <- alter_tcga[which(alter_tcga$ABCB8..HETLOSS.HOMDEL.AMP.GAIN %in% c("HETLOSS","HOMDEL")),]
ABCB8_intact_sample <- alter_tcga[which(alter_tcga$ABCB8..HETLOSS.HOMDEL.AMP.GAIN == "no alteration"),]
ABCB8_deletion_sample <- ABCB8_deletion_sample$Sample.ID  #25
ABCB8_intact_sample <- ABCB8_intact_sample$Sample.ID  #166
GALNT11_deletion_sample <- alter_tcga[which(alter_tcga$GALNT11..HETLOSS.HOMDEL.AMP.GAIN %in% c("HETLOSS","HOMDEL")),]
GALNT11_intact_sample <- alter_tcga[which(alter_tcga$GALNT11..HETLOSS.HOMDEL.AMP.GAIN == "no alteration"),]
GALNT11_deletion_sample <- GALNT11_deletion_sample$Sample.ID  #25
GALNT11_intact_sample <- GALNT11_intact_sample$Sample.ID  #166
length(intersect(ABCB8_deletion_sample,GALNT11_deletion_sample))
ABCB8_deletion_sample <- as.character(ABCB8_deletion_sample)
ABCB8_intact_sample <- as.character(ABCB8_intact_sample)
ABCB8_deletion_patinet <- substr(ABCB8_deletion_sample,1,12)
ABCB8_intact_patient <- substr(ABCB8_intact_sample,1,12)
ABCB8_deletion_patinet <- gsub("-",".",ABCB8_deletion_patinet)
ABCB8_intact_patient <- gsub("-",".",ABCB8_intact_patient)
ABCB8_deletion_patinet <- unique(ABCB8_deletion_patinet)
ABCB8_intact_patient <- unique(ABCB8_intact_patient)
length(ABCB8_intact_patient)
length(ABCB8_deletion_patinet)
sampleTable <- data.frame(sample=c(ABCB8_intact_patient,ABCB8_deletion_patinet),
  group=c(rep("intact_7q",166),rep("del_7q",25)))

sampleTable_intact <- sampleTable[sampleTable$group%in%c("intact_7q"),]

df_trns$sample <- substr(df_trns$sample, 1, 12)


df_trns_no7q <- df_trns[df_trns$sample%in%c(sampleTable_intact$sample),]
coxph_result <- coxph(formula = Surv(OS, status) ~ GALNT11, data = df_trns_no7q)
summary(coxph_result,data=df_trns_no7q)
ggforest(coxph_result, data =df_trns_no7q, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)

df_trns_no7q.cut <- surv_cutpoint(
   df_trns_no7q,
   time = "OS",
   event = "status", 
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(df_trns_no7q.cut)
p1 <- plot(df_trns_no7q.cut, "GALNT11",palette = "npg")

df_trns_no7q.cut.cat <- surv_categorize(df_trns_no7q.cut) 
fit <- survfit(Surv(OS, status) ~ GALNT11, data = df_trns_no7q.cut.cat)
p2 <- ggsurvplot(fit, data = df_trns_no7q.cut.cat,
surv.median.line = "hv",
color = NULL, 
palette = "lancet", 
pval = TRUE,
risk.table=TRUE,
xlim=c(0,1500),
break.x.by=500)

png("Figure1E_TCGA-GLNT11_SURVIVE_gene_0227_no7q.png")
p2
dev.off()


<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cancer discovery Cohort <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< cancer discovery Cohort <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

setwd("./Data_Prepare_update/Figure1")
library(readxl)
data11 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_1.1_Clinical_summary_186_Patients.xlsx")
data12 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_1.2_Clinical_summary_186_Patients_Description.xlsx")
data13 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_3_Drug_response_sDSS_164S_17Healthy.xlsx")
data14 <- read.csv("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/RNA_seq_AML_norm.counts.csv")


data11 <- data.frame(data11)
data14_GALNT11 <- data14[data14$X%in%c("GALNT11"),]
rownames(data14_GALNT11) <- data14_GALNT11$X
data14_GALNT11 <- subset(data14_GALNT11,select=-c(X))

data14_GALNT11 <- data.frame(t(data14_GALNT11))
data14_GALNT11$Patient_ID <- substr(rownames(data14_GALNT11), 1, 7)


dat_use <- merge(data14_GALNT11,data11,by="Patient_ID")
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(os_time, os_event) ~ GALNT11, data = dat_use)
summary(coxph_result,data=dat_use)
ggforest(coxph_result, data =dat_use, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat.cut <- surv_cutpoint(
   dat_use,
   time = "os_time",
   event = "os_event",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat.cut)
plot(dat.cut, "GALNT11")
dat.cut.cat <- surv_categorize(dat.cut) 
library(survival)
fit <- survfit(Surv(os_time, os_event) ~ GALNT11, data = dat.cut.cat)
p1 <- ggsurvplot(fit, data = dat.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500, 
palette = c("#A52A2A","#4682B4"))
png(file="Figure1F_GALNT11_SURVIVE_cancerdis.png")
p1
dev.off()



chr7q_sample <- c("AML_008","AML_020","AML_035","AML_041","AML_068","AML_076","AML_100","AML_104","AML_110","AML_138","AML_139","AML_157","AML_092")
dat_use <- merge(data14_GALNT11,data11,by="Patient_ID")
dat_use_no7q <- dat_use[!(dat_use$Patient_ID %in% chr7q_sample), ]
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(os_time, os_event) ~ GALNT11, data = dat_use_no7q)
summary(coxph_result,data=dat_use_no7q)
ggforest(coxph_result, data =dat_use_no7q, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat.cut <- surv_cutpoint(
   dat_use_no7q,
   time = "os_time",
   event = "os_event",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat.cut)
plot(dat.cut, "GALNT11")
dat.cut.cat <- surv_categorize(dat.cut) 
library(survival)
fit <- survfit(Surv(os_time, os_event) ~ GALNT11, data = dat.cut.cat)
p1 <- ggsurvplot(fit, data = dat.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500, 
palette = c("#A52A2A","#4682B4"))
png(file="Figure1F_GALNT11_SURVIVE_cancerdis_no7q.png")
p1
dev.off()


<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 陈赛娟 Cohort <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 剔除PML::RARA 患者
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 陈赛娟 Cohort <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 剔除PML::RARA 患者

setwd("./Galnt11_Pancancer/blood_csj")
library(readxl)
blood_path <- "./Galnt11_Pancancer/blood_csj/blood_bld-2024-027692-mmc2.xlsx"
excel_sheets(blood_path)
blood_df <- read_excel(blood_path, sheet = 1)

dim(blood_df)
head(blood_df)
library(openxlsx)
expr_path <- "Gene expression matrix of 361 AML patients with RNA-seq data in transcript per million (TPM) format.xlsx"
getSheetNames(expr_path)
expr_df <- read.xlsx(expr_path, sheet = 1)

dim(expr_df)
expr_df[1:5, 1:5]

library(data.table)
prot_path <- "Protein abundance of 374 AML patients after adjusting batch effects.csv"
protein_df <- fread(prot_path, data.table = FALSE)
dim(protein_df)
protein_df[1:5, 1:5]

setwd("./Data_Prepare_update/Figure1")

expr_df[1:5, 1:5]
expr_df_GALNT11 <- expr_df[expr_df$symbol%in%c("GALNT11","ATM"),]
rownames(expr_df_GALNT11) <- expr_df_GALNT11$symbol

expr_df_GALNT11 <- subset(expr_df_GALNT11,select=-c(geneid,symbol))
expr_df_GALNT11 <- data.frame(t(expr_df_GALNT11))
expr_df_GALNT11$Sample_ID <- rownames(expr_df_GALNT11)


blood_df_no_RARA <- blood_df[!(blood_df$WHO_classification_2022%in%c("PML::RARA")),]

dat_use <- merge(expr_df_GALNT11,blood_df_no_RARA,by="Sample_ID")
dat_use_no_7q <- dat_use[dat_use$Minus7_7q==0,]
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS, OS_status) ~ GALNT11, data = dat_use_no_7q)
summary(coxph_result,data=dat_use_no_7q)
ggforest(coxph_result, data =dat_use_no_7q, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat.cut <- surv_cutpoint(
   dat_use_no_7q,
   time = "OS",
   event = "OS_status",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat.cut)
plot(dat.cut, "GALNT11")
dat.cut.cat <- surv_categorize(dat.cut) 
library(survival)
fit <- survfit(Surv(OS, OS_status) ~ GALNT11, data = dat.cut.cat)
p1 <- ggsurvplot(fit, data = dat.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500, 
palette = c("#A52A2A","#4682B4"))


png(file="Figure1G_GALNT11_SURVIVE_chen_no_7q_new_v0303.png")
p1
dev.off()


dat_use_7q <- dat_use[dat_use$Minus7_7q==0|dat_use$Minus7_7q==1,]
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS, OS_status) ~ GALNT11, data = dat_use_7q)
summary(coxph_result,data=dat_use_7q)
ggforest(coxph_result, data =dat_use_7q, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat.cut <- surv_cutpoint(
   dat_use_7q,
   time = "OS",
   event = "OS_status",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat.cut)
plot(dat.cut, "GALNT11")
dat.cut.cat <- surv_categorize(dat.cut) 
library(survival)
fit <- survfit(Surv(OS, OS_status) ~ GALNT11, data = dat.cut.cat)
p1 <- ggsurvplot(fit, data = dat.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500, 
palette = c("#A52A2A","#4682B4"))
png(file="Figure1G_GALNT11_SURVIVE_chen_7q_new_v0303.png")
p1
dev.off()

```
| <img src="/Data_Prepare_update/Figure1/GLNT11_SURVIVE_GSE37642.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1E_TCGA_GLNT11_SURVIVE_gene_0227.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1E_TCGA-GLNT11_SURVIVE_gene_0227_no7q.png" width="200">
| <img src="/Data_Prepare_update/Figure1/Figure1F_GALNT11_SURVIVE_cancerdis.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1F_GALNT11_SURVIVE_cancerdis_no7q.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1G_GALNT11_SURVIVE_chen_no_7q_new_v0303.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1G_GALNT11_SURVIVE_chen_7q_new_v0303.png" width="200"> |  |


## Fig1H-I
```R

/usr/local/R4.2/bin/R

setwd("./Data_Prepare_update/Figure1")
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(topGO)
library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(topGO)
library(ggplot2)
gene_gmt_list <- c("GALNT11")
tcga_surv <- readRDS("./workshop/DATABASE/ALL_TCGA_DATA/xena/All_TCGA_survival.rds")
type_code <- substr(sapply(strsplit(tcga_surv$sample, "-"), `[`, 4), 1, 2)
tcga_surv_tumor <- tcga_surv[type_code %in% sprintf("%02d", 1:9), ]
table(substr(sapply(strsplit(tcga_surv_tumor$sample, "-"), `[`, 4), 1, 2))
head(tcga_surv_tumor)
tcga_surv_tumor <- tcga_surv_tumor %>%
  dplyr::select(sample, tumor_type, OS.time, OS)

tcga_surv_tumor$sample <- substr(tcga_surv_tumor$sample, 1, 16)

TCGA_gene_expr <- readRDS("./workshop/DATABASE/ALL_TCGA_DATA/RNA/FPKM/TCGA/All_TCGA_FPKM.rds")
TCGA_gene_expr$symbol <- mapIds(x = org.Hs.eg.db,keys = rownames(TCGA_gene_expr),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
TCGA_gene_expr_sub <- TCGA_gene_expr[TCGA_gene_expr$symbol%in%c(gene_gmt_list),]
rownames(TCGA_gene_expr_sub) <- TCGA_gene_expr_sub$symbol
TCGA_gene_expr_sub <- subset(TCGA_gene_expr_sub,select=-c(symbol))
TCGA_gene_expr_sub_t <- data.frame(t(TCGA_gene_expr_sub))
# TCGA_gene_expr_sub_t <- TCGA_gene_expr_sub_t[which(rowSums(TCGA_gene_expr_sub_t) > 0),]
# TCGA_gene_expr_sub_t[,ncol(TCGA_gene_expr_sub_t)+1] <- apply(TCGA_gene_expr_sub_t,1,mean)
colnames(TCGA_gene_expr_sub_t)[ncol(TCGA_gene_expr_sub_t)] <- "GALNT11"
# TCGA_gene_expr_sub_t <- subset(TCGA_gene_expr_sub_t,select=c("GALNT11"))
library(dplyr)
library(survival)

expr_df <- as.data.frame(TCGA_gene_expr_sub_t)
expr_df$sample <- substr(rownames(expr_df), 1, 16)  # 前 16 位：TCGA-XX-XXXX-01A
head(expr_df)


# TCGA样本类型编码（第14-15位）：
# 01 = Primary Solid Tumor
# 02 = Recurrent Solid Tumor
# 03 = Primary Blood Derived Cancer - Peripheral Blood
# 05 = Additional - New Primary
# 06 = Metastatic
# 07 = Additional Metastatic
# 11 = Solid Tissue Normal
# 12 = Blood Derived Normal
# 13 = Bone Marrow Normal
# tcga_surv_tumor <- tcga_surv_tumor[substr(tcga_surv_tumor$sample, 14, 15) == "01", ]
tumor_codes <- c("01", "02", "03", "05", "06", "07")
tcga_surv_tumor <- tcga_surv_tumor[
  substr(tcga_surv_tumor$sample, 14, 15) %in% tumor_codes, 
]

dat_all <- inner_join(
  tcga_surv_tumor,
  expr_df,
  by = "sample"
)

head(dat_all)

dat_all1 <- dat_all
# 2.1 按 tumor_type 拆成 list
cancer_list <- split(dat_all1, dat_all1$tumor_type)
# 2.2 预先建一个空结果表
cox_results <- data.frame(
  tumor_type = character(),
  n         = integer(),
  HR        = numeric(),
  lower95   = numeric(),
  upper95   = numeric(),
  p         = numeric(),
  stringsAsFactors = FALSE
)

for (cancer in names(cancer_list)) {
  df <- cancer_list[[cancer]]
  df <- df[!is.na(df$OS.time) & !is.na(df$OS) & !is.na(df$GALNT11), ]
  if (nrow(df) < 10 || sd(df$GALNT11) == 0) {
    cox_results <- rbind(
      cox_results,
      data.frame(
        tumor_type = cancer,
        n         = nrow(df),
        HR        = NA,
        lower95   = NA,
        upper95   = NA,
        p         = NA,
        stringsAsFactors = FALSE
      )
    )
    next
  }
  # 核心：Cox 回归
  fit <- coxph(Surv(OS.time, OS) ~ GALNT11, data = df)
  s   <- summary(fit)
  
  hr      <- as.numeric(s$coef[, "exp(coef)"])
  lower   <- as.numeric(s$conf.int[, "lower .95"])
  upper   <- as.numeric(s$conf.int[, "upper .95"])
  p.value <- as.numeric(s$coef[, "Pr(>|z|)"])
  
  cox_results <- rbind(
    cox_results,
    data.frame(
      tumor_type = cancer,
      n         = nrow(df),
      HR        = hr,
      lower95   = lower,
      upper95   = upper,
      p         = p.value,
      stringsAsFactors = FALSE
    )
  )
}

head(cox_results)

cox_results <- cox_results %>%
  mutate(
    FDR    = p.adjust(p, method = "BH"),
    log2HR = log2(HR),
    direction = case_when(
      !is.na(p) & p < 0.05 & HR > 1 ~ "Risk (HR>1)",
      !is.na(p) & p < 0.05 & HR < 1 ~ "Protective (HR<1)",
      TRUE                              ~ "NS"
    )
  ) %>%
  arrange(p)
head(cox_results, 10)
write.csv(cox_results,"1_cox_results_GANLT11_gene_new_use.csv")



cox_results <- read.csv("./Galnt11_Pancancer/1_cox_results_GANLT11_gene_new_use.csv")
library(ggplot2)
library(RColorBrewer)
# 先确保 cancer_group 是 factor
cox_results$tumor_type <- factor(cox_results$tumor_type)

# 为每个 tumor_type 生成一个唯一颜色
n_grp   <- nlevels(cox_results$tumor_type)
my_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n_grp)  # 自动插值生成 n_grp 个颜色

library(dplyr)
library(ggplot2)

cox_plot <- cox_results %>%
  mutate(
    tumor_type = reorder(tumor_type, -log2HR),
    neglog10p  = -log10(p)
  )


cox_plot$tumor_type <- factor(
  cox_plot$tumor_type,
  levels = cox_plot$tumor_type[order(cox_plot$log2HR)]
)

cox_plot$color_value <- cox_plot$neglog10p
# cox_plot$color_value[cox_plot$p > 0.05] <- 0.3

p2 <- ggplot(cox_plot, aes(x = tumor_type, y = log2HR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    aes(size = n, color = color_value),
    alpha = 0.95
  ) +
  scale_size_continuous(
    name  = "Cohort size (n)",
    range = c(3, 12)
  ) +scale_color_gradientn(
  name = expression(-log[10](p)),
  colours = c(
    "#ffffff",  # 白
    "#fff0f3",  # 极浅粉
    "#fdd0d8",  # 浅粉
    "#fca5b5",  # 粉红
    "#fb6a4a",  # 红（更显著才开始明显红）
    "#cb181d"   # 深红（最显著）
  ),
  values = scales::rescale(c(
    0,                 # p=1
    0.6,               # 很不显著
    1.0,               # 不显著
    1.3,               # p≈0.05（仍然只是粉红，不会突然变红）
    2.0,               # 更显著开始红
    max(cox_plot$neglog10p)
  )),
  limits = c(0, max(cox_plot$neglog10p))
)+
  theme_bw() +
  theme(
    axis.text.x = element_text(
      size = 8, color = "black",
      angle = 90, hjust = 1, vjust = 0.5
    ),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  labs(
    y     = "log2(HR) of GALNT11 expression",
    title = "Pan-cancer GALNT11"
  )

p2
ggsave(p2,file="Figure1H-HR_of_Pancaner_dotplot_v3_neworder.png",height=2.7,width=5.5)



tcga_survival <- readRDS("./workshop/DATABASE/ALL_TCGA_DATA/xena/All_TCGA_survival.rds")
dlbcl_survival <- tcga_survival[tcga_survival$tumor_type=="TCGA-KIRP",]
dlbcl_survival$sample <- gsub("-", ".", dlbcl_survival$sample)
dlbcl_rna <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/RNA/FPKM/TCGA/TCGA-KIRP_RNA_expr.csv")
rownames(dlbcl_rna) <- dlbcl_rna$X
dlbcl_rna$symbol <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
dlbcl_rna$entrez <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- dlbcl_rna$symbol
AA <- as.character(AA)
dlbcl_rna$GENENAME <- mapIds(x = org.Hs.eg.db,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
dlbcl_rna$GENENAME <- as.character(dlbcl_rna$GENENAME)
dlbcl_rna <- dlbcl_rna[,2:(ncol(dlbcl_rna)-2)]
colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)] <- substr(colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)],1,16)
gene_list <- c("GALNT11")
length(unique(gene_list))
dlbcl_rna_chr17 <- dlbcl_rna[dlbcl_rna$symbol %in% gene_list,]
rownames(dlbcl_rna_chr17)<-dlbcl_rna_chr17$symbol
dlbcl_rna_chr17 <- dlbcl_rna_chr17[,1:(ncol(dlbcl_rna)-1)]
dat_new <- data.frame(t(dlbcl_rna_chr17))
dat_new$sample <- rownames(dat_new)

dat_new <- merge(dlbcl_survival,dat_new,by="sample",all.x=TRUE)
dat_new <- na.omit(dat_new)
df_trns <- dat_new
library(survival)
library(survminer)
library(dplyr)
dat_new <- df_trns
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS.time, OS) ~ GALNT11, data = dat_new)
summary(coxph_result,data=dat_new)
ggforest(coxph_result, data =dat_new, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat_new.cut <- surv_cutpoint(
   dat_new,
   time = "OS.time",
   event = "OS",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat_new.cut)
plot(dat_new.cut, "GALNT11")
dat_new.cut.cat <- surv_categorize(dat_new.cut) 
library(survival)
fit <- survfit(Surv(OS.time, OS) ~ GALNT11, data = dat_new.cut.cat)
p1 <- ggsurvplot(fit, data = dat_new.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500,   
palette = c("#A52A2A","#4682B4"))
png(file="Figure1I-V_Survival_plot_TCGA-KIRP.png")
p1
dev.off()

tcga_survival <- readRDS("./workshop/DATABASE/ALL_TCGA_DATA/xena/All_TCGA_survival.rds")
dlbcl_survival <- tcga_survival[tcga_survival$tumor_type=="TCGA-LUAD",]
dlbcl_survival$sample <- gsub("-", ".", dlbcl_survival$sample)
dlbcl_rna <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/RNA/FPKM/TCGA/TCGA-LUAD_RNA_expr.csv")
rownames(dlbcl_rna) <- dlbcl_rna$X
dlbcl_rna$symbol <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
dlbcl_rna$entrez <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- dlbcl_rna$symbol
AA <- as.character(AA)
dlbcl_rna$GENENAME <- mapIds(x = org.Hs.eg.db,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
dlbcl_rna$GENENAME <- as.character(dlbcl_rna$GENENAME)
dlbcl_rna <- dlbcl_rna[,2:(ncol(dlbcl_rna)-2)]
colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)] <- substr(colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)],1,16)
gene_list <- c("GALNT11")
length(unique(gene_list))
dlbcl_rna_chr17 <- dlbcl_rna[dlbcl_rna$symbol %in% gene_list,]
rownames(dlbcl_rna_chr17)<-dlbcl_rna_chr17$symbol
dlbcl_rna_chr17 <- dlbcl_rna_chr17[,1:(ncol(dlbcl_rna)-1)]
dat_new <- data.frame(t(dlbcl_rna_chr17))
dat_new$sample <- rownames(dat_new)

dat_new <- merge(dlbcl_survival,dat_new,by="sample",all.x=TRUE)
dat_new <- na.omit(dat_new)
df_trns <- dat_new
library(survival)
library(survminer)
library(dplyr)
dat_new <- df_trns
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS.time, OS) ~ GALNT11, data = dat_new)
summary(coxph_result,data=dat_new)
ggforest(coxph_result, data =dat_new, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat_new.cut <- surv_cutpoint(
   dat_new,
   time = "OS.time",
   event = "OS",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat_new.cut)
plot(dat_new.cut, "GALNT11")
dat_new.cut.cat <- surv_categorize(dat_new.cut) 
library(survival)
fit <- survfit(Surv(OS.time, OS) ~ GALNT11, data = dat_new.cut.cat)
p1 <- ggsurvplot(fit, data = dat_new.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500,   
palette = c("#A52A2A","#4682B4"))
png(file="Figure1I-V_Survival_plot_TCGA-LUAD.png")
p1
dev.off()



tcga_survival <- readRDS("./workshop/DATABASE/ALL_TCGA_DATA/xena/All_TCGA_survival.rds")
dlbcl_survival <- tcga_survival[tcga_survival$tumor_type=="TCGA-PAAD",]
dlbcl_survival$sample <- gsub("-", ".", dlbcl_survival$sample)
dlbcl_rna <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/RNA/FPKM/TCGA/TCGA-PAAD_RNA_expr.csv")
rownames(dlbcl_rna) <- dlbcl_rna$X
dlbcl_rna$symbol <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
dlbcl_rna$entrez <- mapIds(x = org.Hs.eg.db,keys = rownames(dlbcl_rna),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- dlbcl_rna$symbol
AA <- as.character(AA)
dlbcl_rna$GENENAME <- mapIds(x = org.Hs.eg.db,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
dlbcl_rna$GENENAME <- as.character(dlbcl_rna$GENENAME)
dlbcl_rna <- dlbcl_rna[,2:(ncol(dlbcl_rna)-2)]
colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)] <- substr(colnames(dlbcl_rna)[1:(ncol(dlbcl_rna)-1)],1,16)
gene_list <- c("GALNT11")
length(unique(gene_list))
dlbcl_rna_chr17 <- dlbcl_rna[dlbcl_rna$symbol %in% gene_list,]
rownames(dlbcl_rna_chr17)<-dlbcl_rna_chr17$symbol
dlbcl_rna_chr17 <- dlbcl_rna_chr17[,1:(ncol(dlbcl_rna)-1)]
dat_new <- data.frame(t(dlbcl_rna_chr17))
dat_new$sample <- rownames(dat_new)

dat_new <- merge(dlbcl_survival,dat_new,by="sample",all.x=TRUE)
dat_new <- na.omit(dat_new)
df_trns <- dat_new
library(survival)
library(survminer)
library(dplyr)
dat_new <- df_trns
library("survival")
library("survminer")
coxph_result <- coxph(formula = Surv(OS.time, OS) ~ GALNT11, data = dat_new)
summary(coxph_result,data=dat_new)
ggforest(coxph_result, data =dat_new, 
         main = "Hazard ratio", 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 3)
dat_new.cut <- surv_cutpoint(
   dat_new,
   time = "OS.time",
   event = "OS",
   variables = c("GALNT11"),
   progressbar=TRUE,
   minprop=0.1
)
summary(dat_new.cut)
plot(dat_new.cut, "GALNT11")
dat_new.cut.cat <- surv_categorize(dat_new.cut) 
library(survival)
fit <- survfit(Surv(OS.time, OS) ~ GALNT11, data = dat_new.cut.cat)
p1 <- ggsurvplot(fit, data = dat_new.cut.cat,
surv.median.line = "hv",
pval = TRUE,
pval.method = TRUE,ggtheme = theme_pubr(),
risk.table=TRUE,xlim=c(0,1500),break.time.by = 500,   
palette = c("#A52A2A","#4682B4"))
png(file="Figure1I-V_Survival_plot_TCGA-PAAD.png")
p1
dev.off()

```
| <img src="/Data_Prepare_update/Figure1/Figure1H-HR_of_Pancaner_dotplot_v3_neworder.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1I-V_Survival_plot_TCGA-LUAD.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1I-V_Survival_plot_TCGA-PAAD.png" width="200"> | <img src="/Data_Prepare_update/Figure1/Figure1I-V_Survival_plot_TCGA-KIRP.png" width="200">

## Supplementary1B
```R
setwd("./Galnt11_Pancancer/blood_csj")
library(readxl)
blood_path <- "./Galnt11_Pancancer/blood_csj/blood_bld-2024-027692-mmc2.xlsx"
excel_sheets(blood_path)
blood_df <- read_excel(blood_path, sheet = 1)
dim(blood_df)
head(blood_df)
library(openxlsx)
expr_path <- "Gene expression matrix of 361 AML patients with RNA-seq data in transcript per million (TPM) format.xlsx"
getSheetNames(expr_path)
expr_df <- read.xlsx(expr_path, sheet = 1)
dim(expr_df)
expr_df[1:5, 1:5]
library(data.table)
prot_path <- "Protein abundance of 374 AML patients after adjusting batch effects.csv"
protein_df <- fread(prot_path, data.table = FALSE)
dim(protein_df)
protein_df[1:5, 1:5]
setwd("./Data_Prepare_update/Figure1")


expr_df[1:5, 1:5]
expr_df_GALNT11 <- expr_df[expr_df$symbol%in%c("GALNT11","ATM"),]
rownames(expr_df_GALNT11) <- expr_df_GALNT11$symbol
expr_df_GALNT11 <- subset(expr_df_GALNT11,select=-c(geneid,symbol))
expr_df_GALNT11 <- data.frame(t(expr_df_GALNT11))
expr_df_GALNT11$Sample_ID <- rownames(expr_df_GALNT11)
dat_use <- merge(expr_df_GALNT11,blood_df,by="Sample_ID")
library(dplyr)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

df <- dat_use %>%
  mutate(
    across(everything(), ~ na_if(as.character(.), ".")),
    Age       = parse_number(Age),
    BM_blasts = parse_number(BM_blasts),
    WBC       = parse_number(WBC),
    HGB       = parse_number(HGB),
    PLT       = parse_number(PLT),
    GALNT11   = as.numeric(GALNT11),
    OS_status  = as.numeric(OS_status),
    EFS_status = as.numeric(EFS_status),
    HSCT       = as.numeric(HSCT)
  ) %>%
  filter(!is.na(GALNT11))

q <- quantile(df$GALNT11, probs = c(0.1, 0.9), na.rm = TRUE)

df <- df %>%
  mutate(
    GALNT11_HL = case_when(
      GALNT11 <= q[1] ~ "Low",
      GALNT11 >= q[2] ~ "High",
      TRUE ~ NA_character_
    ),
    GALNT11_HL = factor(GALNT11_HL, levels = c("High","Low"))
  ) %>%
  filter(!is.na(GALNT11_HL)) %>%      # 只保留 High/Low 两端
  arrange(desc(GALNT11))              # 样本按 GALNT11 高→低

samps <- df$Sample_ID

col_age  <- colorRamp2(range(df$Age, na.rm=TRUE), c("white","red"))
col_hgb  <- colorRamp2(range(df$HGB, na.rm=TRUE), c("white","orange"))
col_plt  <- colorRamp2(range(df$PLT, na.rm=TRUE), c("white","green"))
col_gal  <- colorRamp2(range(df$GALNT11, na.rm=TRUE), c("white","black"))

ha_top <- HeatmapAnnotation(
  Group = df$GALNT11_HL,
  GALNT11 = df$GALNT11,
  Sex = df$Sex,
  Age = df$Age,
  FAB = df$FAB,
  ELN = df$ELN_risk_2017,
  HGB = df$HGB,
  PLT = df$PLT,
  HSCT = df$HSCT,
  OS_status = df$OS_status,
  EFS_status = df$EFS_status,
  col = list(
    GALNT11 = col_gal,
    Age = col_age,
    HGB = col_hgb,
    PLT = col_plt
  ),
  annotation_name_side = "left"
)

karyo_cols <- c(
  "Normal_karyotype","Complex_karyotype","Monokaryotype",
  "Trisomy8","Minus5_5q","Minus7_7q","Minus17_abn17p"
)

karyo_mat <- df %>%
  select(all_of(karyo_cols)) %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%
  as.matrix() %>% t()
colnames(karyo_mat) <- samps
rownames(karyo_mat) <- karyo_cols

who <- df$WHO_classification_2022
who_levels <- sort(unique(who[!is.na(who)]))
who_mat <- sapply(who_levels, function(lv) as.numeric(who == lv)) %>% t()
colnames(who_mat) <- samps
rownames(who_mat) <- paste0("WHO: ", who_levels)

mat <- rbind(who_mat, karyo_mat)
col_bin <- c("0"="white", "1"="#2c7fb8")
ht <- Heatmap(
  mat,
  name = "Status",
  col = col_bin,
  na_col = "grey90",
  top_annotation = ha_top,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  column_split = df$GALNT11_HL,
  column_gap = unit(2, "mm")
)

png("Supplementary1B-heatmap_AML_clinical_01.png")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

```
<img src="/Data_Prepare_update/Figure1/Supplementary1B-heatmap_AML_clinical_01.png" width="300">

## Fig4D
```R
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

tab <- read.delim("./workshop/RNAseq/RNAseq_83_ZJN_20210629_12samples/sortedByCoord.out.bam/shG494_vs_shRen_3_vs_6_work_file/GSEA/shG494_versus_shRen_h/my_analysis.Gsea.1625373361910/ranked_gene_list_shG494_versus_shRen_1625373361910.xls",
                  check.names = FALSE, stringsAsFactors = FALSE)
P_subcluster <- tab
P_subcluster <- P_subcluster[,c("SCORE","NAME")]
P_subcluster <- P_subcluster[order(P_subcluster$SCORE,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$SCORE
names(aa) <- P_subcluster$NAME
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("./msigdb.v2022.1.Mm.symbols.gmt")
custom_GSEA_GMT <- read.gmt("./programme/gsea/msigdb_v7.1/msigdb_v7.1_GMTs/indivi_gmt/h.all.v7.1.symbols.gmt")
custom_GSEA_GMT <- custom_GSEA_GMT[custom_GSEA_GMT$term=="HALLMARK_DNA_REPAIR",]

gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',minGSSize = 15, maxGSSize = 500,
    nPerm = 1000)  #指定 p 值校正方法

geneset <- c("HALLMARK_DNA_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3','#459943'),base_size=12,subplots = 1:2)
p
pdf("GSEA_shGALNT11_DNA_DAMAGE.pdf",height=4,width=4.5)
p
dev.off()

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3','#459943'),base_size=12,subplots = 1:2)
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="GSEA_TCGA-GSEA_shGALNT11_DNA_DAMAGE.png",height=4,width=4.5)
 $ enrichmentScore: num -0.367
 $ NES            : num -1.46
 $ pvalue         : num 0.00926
 $ p.adjust       : num 0.00926
```

## Fig4E-H
```R
<<<<<<<<<<<<<<<<<<GSE37642
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- read.csv("5_GSE37642-_DEGs_GALNT11_low_vs_high.csv")
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)
custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH')  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"5-GSE37642-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

gsea_df[grep("REPAIR",gsea_df$ID),]
gsea_df[grep("DAMAGE",gsea_df$ID),]
gsea_df[grep("RECOM",gsea_df$ID),]

library(ggplot2)
library(ggpubr)
geneset <- c("HALLMARK_DNA_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p 
ggsave(p,file="5_GSE37642-galnt11_low_vs_high_DNA_DAMAGE.png",height=4,width=4.5)
ggsave(p,file="5_GSE37642-galnt11_low_vs_high_DNA_DAMAGE.pdf",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]
                    enrichmentScore       NES       pvalue     p.adjust
HALLMARK_DNA_REPAIR      -0.5173262 -1.711447 6.825839e-06 0.0001246911

<<<<<<<<<<<<<<<<<<TCGA
degs <- read.csv("./New_analysis/1_TCGALAML_DEGs_data_GALNT11_expr.csv")
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

custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-TCGALAML_GALNT11-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

gsea_df[grep("REPAIR",gsea_df$ID),]
geneset <- c("KEGG_NUCLEOTIDE_EXCISION_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="1_TCGALAML-GSEA_DNA_pathway.png",height=4,width=4.5)
ggsave(p,file="1_TCGALAML-GSEA_DNA_pathway.pdf",height=4,width=4.5)

gsea_df[gsea_df$ID%in%c("KEGG_NUCLEOTIDE_EXCISION_REPAIR"),]
                               enrichmentScore       NES      pvalue
KEGG_NUCLEOTIDE_EXCISION_REPAIR      -0.4826737 -1.787912 0.002442002

<<<<<<<<<<<<<<<beatAML
degs <- read.csv("./New_analysis/2_beatAML_DEGs_data_GALNT11_expr.csv")
rownames(degs) <- degs$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res<- degs
colnames(res)[1] <- "symbol"
res[res$symbol%in%c("GALNT11"),]

      symbol baseMean log2FoldChange     lfcSE      stat       pvalue
GALNT11 GALNT11 1206.079       -1.63727 0.1169566 -13.99896 1.581756e-44
                padj
GALNT11 5.411504e-40
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
custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"2_beatAML_GALNT11-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

geneset <- c("KEGG_MISMATCH_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="2_beatAML-GSEA_DNA_pathway.png",height=4,width=4.5)
ggsave(p,file="2_beatAML-GSEA_DNA_pathway.pdf",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]
                enrichmentScore       NES     pvalue   p.adjust     qvalue
KEGG_MISMATCH_REPAIR      -0.6845746 -2.840519 0.00101833 0.01387207 0.01049712


<<<<<<<<<<<<<<<<<<<< GSE10358

###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- read.csv("4—GSE10358_DEGs_GALNT11_low_vs_high.csv")
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"4-GSE10358-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

gsea_df[grep("DNA",gsea_df$ID),]
gsea_df[grep("DAMAGE",gsea_df$ID),]
gsea_df[grep("RECOM",gsea_df$ID),]

library(ggplot2)
library(ggpubr)
geneset <- c("REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV2.png",height=4,width=4.5)
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV2.pdf",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]

                                                        NES      pvalue
REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR -1.801405 0.000400641
                                                    p.adjust     qvalue rank
REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR 0.01284095 0.01053379 5420

```

| <img src="/Data_Prepare_update/Figure4/GSEA_TCGA-GSEA_shGALNT11_DNA_DAMAGE.png" width="200"> | <img src="/Data_Prepare_update/Figure4/5_GSE37642-galnt11_low_vs_high_DNA_DAMAGE.png" width="200"> | <img src="/Data_Prepare_update/Figure4/1_TCGALAML-GSEA_DNA_pathway.png" width="200"> | <img src="/Data_Prepare_update/Figure4/2_beatAML-GSEA_DNA_pathway.png" width="200"> | <img src="/Data_Prepare_update/Figure4/4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV2.png" width="200">

## Fig4I
```R
/usr/local/R4.2/bin/R
setwd("./Data_Prepare_update/Figure4/")
library(ggplot2)
library(ggpubr)
library(tidyr)
library(trqwe)
library(maftools)
Count <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/RNA/Counts/TCGA/TCGA-LAML_RNA_expr_Counts.csv",check.names=F)
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

Count_tmp4 <- Count_tmp3[order(Count_tmp3$GALNT11),]
Count_tmp4_lo <- head(Count_tmp4,16)#0.1
Count_tmp4_high <- tail(Count_tmp4,16)#0.1
Count_tmp4_lo$Group <- "low"
Count_tmp4_high$Group <- "high"
galnt11_tmb_exp1 <- rbind(Count_tmp4_lo,Count_tmp4_high)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(trqwe)
library(maftools)
PAAD_maf <- read.csv("./workshop/DATABASE/ALL_TCGA_DATA/Mut_Maf/maf_files/TCGA_LAML_Mutation.csv")
PAAD_maf <- read.maf(PAAD_maf)

#TMB计算
aml_tmb <-  tmb(maf = PAAD_maf) #131
aml_tmb <- as.data.frame(aml_tmb)
aml_tmb$sample <- substr(aml_tmb$Tumor_Sample_Barcode,1,15)
# aml_tmb <- aml_tmb[!duplicated(aml_tmb$sample),]  
galnt11_tmb_exp1$sample <- substr(galnt11_tmb_exp1$sample,1,12)
galnt11_tmb_exp2 <- merge(aml_tmb,galnt11_tmb_exp1,by='sample')
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(Seurat)
galnt11_tmb_exp2$Group <- factor(galnt11_tmb_exp2$Group,levels=c("low","high"))
dat1 <- galnt11_tmb_exp2[galnt11_tmb_exp2$total_perMB<1,]
p1  <-  ggplot(dat1, aes(x=Group, y=total_perMB,fill = Group)) +
  geom_boxplot()+RotatedAxis() + theme_classic()
dat1$Group <- as.factor(dat1$Group)
bino_model <- glm(dat1$Group ~ dat1$total_perMB, family="binomial")
anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 0.04376428


tcga_clinical_sample1 <- read.table("./human_data/laml_tcga_pan_can_atlas_2018/data_clinical_sample.txt",sep="\t",header=T)
dat11 <- merge(dat1,tcga_clinical_sample1,by.x="sample",by.y="PATIENT_ID")
p2  <-  ggplot(dat11, aes(x=Group, y=ANEUPLOIDY_SCORE,fill = Group)) +
  geom_boxplot()+ RotatedAxis() + theme_classic()
p2
g <- p1+p2
ggsave(g,file="Figure4I_1_TCGALAML_boxplot_TMB_ANEUPOLID.png")
ggsave(g,file="Figure4I_1_TCGALAML_boxplot_TMB_ANEUPOLID.pdf")

dat11$Group <- as.factor(dat11$Group)
bino_model <- glm(dat11$Group ~ dat11$ANEUPLOIDY_SCORE, family="binomial")
anova(bino_model, test="F")$"Pr(>F)"[2]
[1] 0.0106542

```
<img src="/Data_Prepare_update/Figure4/Figure4I_1_TCGALAML_boxplot_TMB_ANEUPOLID.png" width="300">

## Fig5B-C
```R

<<<<<<<<<<<<<<<GSE10358
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- read.csv("4—GSE10358_DEGs_GALNT11_low_vs_high.csv")
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"4-GSE10358-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
geneset <- c("PID_ATM_PATHWAY")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
gsea_df[gsea_df$ID%in%c("PID_ATM_PATHWAY"),]
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_ATM.png",height=4,width=4.5)
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_ATM.png",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]

               NES    pvalue  p.adjust   qvalue rank
PID_ATM_PATHWAY -1.482635 0.0384127 0.1913271 0.156951 6813


<<<<<<<<<<<<<<<GSE6891
###GSEA
library(clusterProfiler)
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- reda.csv("3_GSE6891_DEGs_GALNT11_low_vs_high.csv")
P_subcluster$symbol <- rownames(P_subcluster)
P_subcluster <- P_subcluster[,c("logFC","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$logFC,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$logFC
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"3-GSE6891AML-GALNT11-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
geneset <- c("ATM_DN.V1_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
gsea_df[gsea_df$ID%in%c("ATM_DN.V1_DN"),]
ggsave(p,file="3_GSE6891-galnt11_low_vs_high_ATM.png",height=4,width=4.5)
ggsave(p,file="GSE6891-galnt11_low_vs_high_ATM.pdf",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]

                     ID  Description setSize enrichmentScore       NES
ATM_DN.V1_DN ATM_DN.V1_DN ATM_DN.V1_DN     142      -0.4770668 -1.681594
                   pvalue    p.adjust      qvalue rank
ATM_DN.V1_DN 0.0006793478 0.008214701 0.006172544 2186
```

| <img src="/Data_Prepare_update/Figure5/4-GSE10358-galnt11_low_vs_high_ATM.png" width="300"> | <img src="/Data_Prepare_update/Figure5/GSE6891-galnt11_low_vs_high_ATM.png" width="300">


## Fig5J-M
```R
<<<<<<<<<<<<<<<<<<<<CCLE
setwd("./Data_Prepare_update/Figure5/")

library(trqwe)
sample_info <- read.csv("./workshop/DATABASE/CCLE/sample_info.csv",row.names=1) #1775
sample_info$cellline_id <- rownames(sample_info)
RPPA <- read.csv("./workshop/DATABASE/CCLE/CCLE_RPPA_20181003.csv",row.names=1)
RPPA$cell_line_name <- rownames(RPPA)
RPPA <- merge(RPPA,sample_info[,c(2,ncol(sample_info))],by.x="cell_line_name",by.y = "CCLE_Name")
rownames(RPPA) <- RPPA$cellline_id
RPPA_info <- read.csv("./workshop/DATABASE/CCLE/CCLE_RPPA_Ab_info_20181226.csv",row.names=1)
expression <- read.csv("./workshop/DATABASE/CCLE/CCLE_expression.csv",row.names=1)
mutation <- read.csv("./workshop/DATABASE/CCLE/CCLE_mutations.csv",row.names=1)
CCLE_gene_cn <- read.csv("./workshop/DATABASE/CCLE/CCLE_gene_cn.csv",row.names=1)
CCLE_RNA <- mcreadRDS("./workshop/DATABASE/CCLE/CCLE_RNA_normalised.rds",mc.cores=20)
all_RPPA_GA <- RPPA[,c(grep("ATM",colnames(RPPA),value=TRUE),grep("Galnt",colnames(RPPA),value=TRUE),"cell_line_name")]
aml_sample <- sample_info[which(sample_info$lineage_subtype=="AML"),]
aml_mutation <- mutation[which(mutation$DepMap_ID %in% aml_sample$cellline_id),]
chr7_del_aml_mutation <- aml_mutation[which(aml_mutation$Variant_Type=="DEL" & aml_mutation$Chromosome==7),]
chr7q_del_aml_mutation <- chr7_del_aml_mutation[which(chr7_del_aml_mutation$Start_position > 60000000),]
RPPA_aml <- na.omit(RPPA_aml) #32

tmp1 <- mutation[which(mutation$Hugo_Symbol=="GALNT11" & mutation$Variant_Type=="DEL"),]
tmp2 <- mutation[which(mutation$Hugo_Symbol=="ABCB8" & mutation$Variant_Type=="DEL"),]
cellline <- c(as.character(tmp1$DepMap_ID),as.character(tmp2$DepMap_ID))
sample_info[cellline,]
tmp <- mutation[which(mutation$Hugo_Symbol=="GALNT11"),]
********** protein *********
********** protein *********
#GALNT11表达量分组
tmp <- CCLE_RNA[c("GALNT11"),]
tmp <- t(tmp)
tmp <- as.data.frame(tmp)
rownames(tmp) <- gsub("\\.","-",rownames(tmp))
rna_galnt11_protein_atm <- merge(tmp,RPPA[,c(grep("ATM",colnames(RPPA),value=TRUE),"cellline_id")],
  by.x="row.names",by.y="cellline_id")
rownames(rna_galnt11_protein_atm) <- rna_galnt11_protein_atm[,1]
rna_galnt11_protein_atm <- rna_galnt11_protein_atm[,-1]

library("survival")
library("survminer")
library(maxstat)
# mod <- maxstat.test(ATM ~ GALNT11, data=rna_galnt11_protein_atm, smethod="Wilcoxon", pmethod="HL")
# mod
# plot(mod)

summary(rna_galnt11_protein_atm$GALNT11)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.09761 3.77426 4.46336 4.15636 5.00518 7.30889
#Wilcoxon模型 cutpoint：2.286881
rna_galnt11_protein_atm$group <- ifelse(rna_galnt11_protein_atm$GALNT11 > 4.46336,"GALNT11_high","GALNT11_low")
rna_galnt11_protein_atm$group <- factor(rna_galnt11_protein_atm$group,levels=c("GALNT11_high","GALNT11_low"))
ggboxplot(rna_galnt11_protein_atm,x="group",y="ATM",color="group",add="jitter",title="GALNT11 RNA-ATM protein in all cellline",palette="npg")+
  stat_compare_means(label = "p.format")

ggscatter(rna_galnt11_protein_atm, x = "GALNT11", y = "ATM", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"), 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "GALNT11_norm_RNA", ylab = "ATM_protein",title="GALNT11-ATM in all cellline")

aml_sample <- sample_info[which(sample_info$lineage_subtype=="AML"),]
aml_sample <- rownames(aml_sample)
aml_rna_galnt11_protein_atm <- rna_galnt11_protein_atm[aml_sample,]
aml_rna_galnt11_protein_atm <- na.omit(aml_rna_galnt11_protein_atm)
summary(aml_rna_galnt11_protein_atm$GALNT11)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3219  3.1745  3.8689  3.5647  4.0811  4.5416
# mod <- maxstat.test(ATM ~ GALNT11, data=aml_rna_galnt11_protein_atm, smethod="Wilcoxon", pmethod="HL")
#2.920293

aml_rna_galnt11_protein_atm$group <- ifelse(aml_rna_galnt11_protein_atm$GALNT11 > 3.5647,"GALNT11_high","GALNT11_low")
aml_rna_galnt11_protein_atm$group <- factor(aml_rna_galnt11_protein_atm$group,levels=c("GALNT11_high","GALNT11_low"))
p1 <- ggboxplot(aml_rna_galnt11_protein_atm,x="group",y="ATM",color="group",add="none",title="GALNT11 RNA-ATM protein in AML",palette="npg")+
  stat_compare_means(label = "p.format",method="t.test")
ggsave(p1,file="Fig5K-ATM_protein_GALNT11.png")

**********RNA**********
**********RNA**********

tmp <- CCLE_RNA[c("ATM","GALNT11"),]
tmp <- t(tmp)
tmp <- as.data.frame(tmp)
rownames(tmp) <- gsub("\\.","-",rownames(tmp))
aml_rna <- tmp[which(rownames(tmp) %in% aml_sample),]
chr7q_del_aml_cellline <- unique(as.character(chr7q_del_aml_mutation$DepMap_ID))
chr7q_del_aml_cellline <- sample_info[chr7q_del_aml_cellline,]
aml_rna$group <- ifelse(rownames(aml_rna) %in% rownames(chr7q_del_aml_cellline),"del_7q","7q_intact")
aml_rna$group <- ifelse(aml_rna$GALNT11 > 3.5647,"GALNT11_high","GALNT11_low")
aml_rna$group <- factor(aml_rna$group,levels=c("GALNT11_high","GALNT11_low"))
p1 <- ggboxplot(aml_rna,x="group",y="ATM",color="group",add="none",title="GALNT11 RNA-ATM protein in AML",palette="npg")+
  stat_compare_means(label = "p.format",method="t.test")
ggsave(p1,file="Fig5J-ATM_RNA_GALNT11.png")

<<<<<<<<<<<<<<<<<<<<CCLE
setwd("./Data_Prepare_update/Figure5/")
library(trqwe)
sample_info <- read.csv("./workshop/DATABASE/CCLE/sample_info.csv",row.names=1) #1775
sample_info$cellline_id <- rownames(sample_info)
RPPA <- read.csv("./workshop/DATABASE/CCLE/CCLE_RPPA_20181003.csv",row.names=1)
RPPA$cell_line_name <- rownames(RPPA)
RPPA <- merge(RPPA,sample_info[,c(2,ncol(sample_info))],by.x="cell_line_name",by.y = "CCLE_Name")
rownames(RPPA) <- RPPA$cellline_id
RPPA_info <- read.csv("./workshop/DATABASE/CCLE/CCLE_RPPA_Ab_info_20181226.csv",row.names=1)
expression <- read.csv("./workshop/DATABASE/CCLE/CCLE_expression.csv",row.names=1)
mutation <- read.csv("./workshop/DATABASE/CCLE/CCLE_mutations.csv",row.names=1)
CCLE_gene_cn <- read.csv("./workshop/DATABASE/CCLE/CCLE_gene_cn.csv",row.names=1)
CCLE_RNA <- mcreadRDS("./workshop/DATABASE/CCLE/CCLE_RNA_normalised.rds",mc.cores=20)
all_RPPA_GA <- RPPA[,c(grep("ATM",colnames(RPPA),value=TRUE),grep("GALNT",colnames(RPPA),value=TRUE),"cell_line_name")]
tmp1 <- mutation[which(mutation$Hugo_Symbol=="GALNT11" & mutation$Variant_Type=="DEL"),]
tmp2 <- mutation[which(mutation$Hugo_Symbol=="ABCB8" & mutation$Variant_Type=="DEL"),]
#没有GALNT11 和 ABCB8 都del的细胞系，且都不是AML的
cellline <- c(as.character(tmp1$DepMap_ID),as.character(tmp2$DepMap_ID))
sample_info[cellline,]
tmp <- mutation[which(mutation$Hugo_Symbol=="GALNT11"),]
aml_sample <- sample_info[which(sample_info$lineage_subtype=="AML"),]
aml_mutation <- mutation[which(mutation$DepMap_ID %in% aml_sample$cellline_id),]
chr7_del_aml_mutation <- aml_mutation[which(aml_mutation$Variant_Type=="DEL" & aml_mutation$Chromosome==7),]
chr7q_del_aml_mutation <- chr7_del_aml_mutation[which(chr7_del_aml_mutation$Start_position > 60000000),]

chr7q_del_aml_cellline <- unique(as.character(chr7q_del_aml_mutation$DepMap_ID))
chr7q_del_aml_cellline <- sample_info[chr7q_del_aml_cellline,]
chr7q_del_aml_cellline <- as.character(chr7q_del_aml_cellline$CCLE_Name) #7个
RPPA_aml <- RPPA[which(RPPA$cell_line_name %in% aml_sample$CCLE_Name),]
RPPA_aml <- na.omit(RPPA_aml) #32
RPPA_aml$group <- ifelse(RPPA_aml$cell_line_name %in% chr7q_del_aml_cellline,"del_7q","7q_intact")
RPPA_aml$group <- factor(RPPA_aml$group,levels=c("7q_intact","del_7q"))
RPPA_aml$group <- as.factor(RPPA_aml$group)
bino_model <- glm(RPPA_aml$group ~ RPPA_aml$ATM, family="binomial")
anova(bino_model, test="F")$"Pr(>F)"[2]
anova(bino_model, test="Rao")$"Pr(>Chi)"[2]
anova(bino_model, test="LRT")$"Pr(>Chi)"[2]
anova(bino_model, test="Chisq")$"Pr(>Chi)"[2]
library(ggpubr)
ggboxplot(RPPA_aml,x="group",y="ATM",fill="group",add="jitter",title="RPPA_ATM in AML",palette="npg")+
  stat_compare_means(label = "p.format")

#所有7q-细胞系
chr7_del_mutation <- mutation[which(mutation$Variant_Type=="DEL" & mutation$Chromosome==7),]
chr7q_del_mutation <- chr7_del_mutation[which(chr7_del_mutation$Start_position > 60000000),]
chr7q_del_cellline <- unique(as.character(chr7q_del_mutation$DepMap_ID))
RPPA_tmp <- RPPA
RPPA_tmp <- RPPA_tmp[,-1] #897
RPPA_tmp <- RPPA_tmp[which(rownames(RPPA_tmp) %in% unique(mutation$DepMap_ID)),] #895
RPPA_tmp$group <- ifelse(rownames(RPPA_tmp) %in% chr7q_del_cellline,"del_7q","7q_intact")
RPPA_tmp$group <- factor(RPPA_tmp$group,levels=c("7q_intact","del_7q"))

ggboxplot(RPPA_tmp,x="group",y="ATM",color="group",add="jitter",title="RPPA_ATM in all cellline",palette="npg")+
  stat_compare_means(label = "p.format")
library(Seurat)
p1 <- ggviolin(RPPA_tmp, "group", "ATM", fill = "group",
   palette = "npg",add = c("boxplot"), add.params = list(fill = "white",alpha=0.8)) +
   stat_compare_means(comparisons=list(c("7q_intact","del_7q")),method="wilcox") + RotatedAxis()
ggsave(p1,file="Fig5L-ATM_RPPA_Expr_in7q.png")

<<<<<<<<<<<<<<<<<<<<GSE10358

DEG_limma <- read.csv("./Public_Data/GSE10358_DEGs_7q_vs_NN.csv")
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
custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH')  #指定 p 值校正方法
gsea_df <- gsea_AM@result
write.csv(gsea_df,"03_GSE10358-7q-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
geneset <- c("ATM_DN.V1_DN")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="03_GSE10358-7q_ATM.png",height=4,width=4.5)
ggsave(p,file="03_GSE10358-7q_ATM.pdf",height=4,width=4.5)
gsea_df[gsea_df$ID%in%c(geneset),]

                       ID  Description setSize enrichmentScore       NES
ATM_DN.V1_DN ATM_DN.V1_DN ATM_DN.V1_DN     142      -0.4441916 -1.780684
                   pvalue     p.adjust       qvalue rank
ATM_DN.V1_DN 2.649351e-05 0.0004539717 0.0003613243 5162
```

| <img src="/Data_Prepare_update/Figure5/Fig5K-ATM_protein_GALNT11.png" width="300"> | <img src="/Data_Prepare_update/Figure5/Fig5J-ATM_RNA_GALNT11.png" width="300"> | <img src="/Data_Prepare_update/Figure5/Fig5L-ATM_RPPA_Expr_in7q.png" width="300"> | <img src="/Data_Prepare_update/Figure5/03_GSE10358-7q_ATM.png" width="300">

## Fig6C-D
```R
setwd("./Data_Prepare_update/Figure6/")
#TCGA
degs <- read.csv("./TCGA/TCGA-Bulk_7Q_DEG.csv")
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

res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]

write.csv(custom_GSEA_GMT,file="custom_GSEA_GMT.csv")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
 $ NES            : num 2.95
 $ pvalue         : num 0.000121
 $ p.adjust       : num 0.000121
 $ qvalue         : logi NA

gsea_df <- gsea_AM@result
write.csv(gsea_df,"00_TCGA-LAML-7q_GSEA_resist_nonCR.csv")
geneset <- c("Resistance_signature")

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="00_TCGA-LAML-7q-GSEA_resist_pathway.png",height=4,width=4.5)
ggsave(p,file="00_TCGA-LAML-7q-GSEA_resist_pathway.pdf",height=4,width=4.5)

#GSE10358
DEG_limma <- read.csv("./Public_Data/GSE10358_DEGs_7q_vs_NN.csv")
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
res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

$ NES            : num 1.52
 $ pvalue         : num 0.00174
 $ p.adjust       : num 0.00174

gsea_df <- gsea_AM@result
write.csv(gsea_df,"00_GSE10358-7q_GSEA_resist_nonCR.csv")
geneset <- c("Resistance_signature")

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="00_GSE10358-7q-GSEA_resist_pathway.png",height=4,width=4.5)
ggsave(p,file="00_GSE10358-7q-GSEA_resist_pathway.pdf",height=4,width=4.5)

#GSE6891
DEG_limma <- read.csv("./Public_Data/GSE6891_DEGs_7q_vs_NN.csv")
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
res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

 $ NES            : num 2.38
 $ pvalue         : num 0.000167
 $ p.adjust       : num 0.000167

gsea_df <- gsea_AM@result
write.csv(gsea_df,"00_GSE6891-7q_GSEA_resist_nonCR.csv")
geneset <- c("Resistance_signature")

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="00_GSE6891-7q-GSEA_resist_pathway.png",height=4,width=4.5)
ggsave(p,file="00_GSE6891-7q-GSEA_resist_pathway.pdf",height=4,width=4.5)

#TCGA-AML
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
DEG_limma <- read.csv("./New_analysis/1_TCGALAML_DEGs_data_GALNT11_expr.csv")
options(future.globals.maxSize= 300*1024^3)
future::plan("multicore", workers = 30)
P_subcluster <- DEG_limma
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
P_subcluster$symbol <- mapIds(x = anno_data,keys = P_subcluster$ensg,keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")

P_subcluster <- P_subcluster[,c("log2FoldChange","symbol")]
P_subcluster <- P_subcluster[order(P_subcluster$log2FoldChange,decreasing=TRUE),]
P_subcluster <- na.omit(P_subcluster)
aa <- P_subcluster$log2FoldChange
names(aa) <- P_subcluster$symbol
geneList = sort(aa,decreasing = T)

res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法
gsea_AM
$ NES            : num 2.52
 $ pvalue         : num 0.000103
 $ p.adjust       : num 0.000103


gsea_df <- gsea_AM@result
write.csv(gsea_df,"1-TCGALAML_GALNT11-GSEA_resist_nonCR.csv")
geneset <- c("Resistance_signature")

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="1_TCGALAML-GSEA_resist_pathway.png",height=4,width=4.5)
ggsave(p,file="1_TCGALAML-GSEA_resist_pathway.pdf",height=4,width=4.5)

#GSER10358

DEG_limma <- read.csv("./New_analysis/4—GSE10358_DEGs_GALNT11_low_vs_high.csv")
colnames(DEG_limma)[1] <- "symbol"
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

custom_GSEA_GMT <- read.gmt("./2_reference/GSVA_7.1/msigdb.v7.1.symbols.gmt")
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

gsea_df <- gsea_AM@result
write.csv(gsea_df,"4-GSE10358-AMLGALNT11_LOW_HIGH-GSEA.csv")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

gsea_df[grep("DNA",gsea_df$ID),]
gsea_df[grep("DAMAGE",gsea_df$ID),]
gsea_df[grep("RECOM",gsea_df$ID),]

library(ggplot2)
library(ggpubr)
geneset <- c("GO_RECOMBINATIONAL_REPAIR")
p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c("#f47720", "#0074b3", "#2a9d8f"),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV1007.png",height=4,width=4.5)
ggsave(p,file="4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV1007.pdf",height=4,width=4.5)

#GSER6891
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
DEG_limma <- read.csv("./New_analysis/3_GSE6891_DEGs_GALNT11_low_vs_high.csv")

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

res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
gsea_AM <- GSEA(geneList,  #待富集的基因列表
    TERM2GENE = custom_GSEA_GMT1,  #基因集
    pvalueCutoff = 1,  #指定 p.adjust 值阈值（可指定 1 以输出全部）
    pAdjustMethod = 'BH',
    nPerm = 10000)  #指定 p 值校正方法

 $ NES            : num 1.95
 $ pvalue         : num 0.00018
 $ p.adjust       : num 0.00018

gsea_df <- gsea_AM@result
write.csv(gsea_df,"3_GSE6891_GALNT11-GSEA_resist_nonCR.csv")
geneset <- c("Resistance_signature")

p <- gseaplot2(gsea_AM,geneSetID=geneset,color = c('#f47720','#0074b3'),base_size=12,subplots = 1:2)
p
p[[1]] <- p[[1]]+theme(legend.position = "top",legend.direction = "vertical")
p
ggsave(p,file="3_GSE6891-GSEA_resist_pathway.png",height=4,width=4.5)
ggsave(p,file="3_GSE6891-GSEA_resist_pathway.pdf",height=4,width=4.5)

```
| <img src="/Data_Prepare_update/Figure6/00_TCGA-LAML-7q-GSEA_resist_pathway.png" width="300"> | <img src="/Data_Prepare_update/Figure6/00_GSE10358-7q-GSEA_resist_pathway.png" width="300"> | <img src="/Data_Prepare_update/Figure6/00_GSE6891-7q-GSEA_resist_pathway.png" width="300">
| <img src="/Data_Prepare_update/Figure6/1_TCGALAML-GSEA_resist_pathway.png" width="300"> | <img src="/Data_Prepare_update/Figure6/4-GSE10358-galnt11_low_vs_high_DNA_DAMAGEV1007.png" width="300"> | <img src="/Data_Prepare_update/Figure6/3_GSE6891-GSEA_resist_pathway.png" width="300">


## Fig6E and Supplementary Fig5A-B
```R
/usr/local/R4.2/bin/R
setwd("./Data_Prepare_update/Figure6/")
degs_7q <- read.csv("./TCGA/TCGA-Bulk_7Q_DEG.csv")
rownames(degs_7q) <- degs_7q$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res_7q<- degs_7q
res_7q$symbol <- mapIds(x = anno_data,keys = rownames(res_7q),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res_7q$entrez <- mapIds(x = anno_data,keys = rownames(res_7q),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res_7q$symbol
AA <- as.character(AA)
res_7q$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res_7q$GENENAME <- as.character(res_7q$GENENAME)
degs_GALNT11 <- read.csv("./New_analysis/1_TCGALAML_DEGs_data_GALNT11_expr.csv")
rownames(degs_GALNT11) <- degs_GALNT11$X
library(org.Hs.eg.db)
organism <-"hsa"
anno_data=org.Hs.eg.db
res_GALNT11<- degs_GALNT11
res_GALNT11$symbol <- mapIds(x = anno_data,keys = rownames(res_GALNT11),keytype = "ENSEMBL",column = "SYMBOL",multiVals = "first")
res_GALNT11$entrez <- mapIds(x = anno_data,keys = rownames(res_GALNT11),keytype ="ENSEMBL",column ="ENTREZID",multiVals="first")
AA <- res_GALNT11$symbol
AA <- as.character(AA)
res_GALNT11$GENENAME <- mapIds(x = anno_data,keys = AA,keytype ="SYMBOL",column ="GENENAME",multiVals="first")
res_GALNT11$GENENAME <- as.character(res_GALNT11$GENENAME)
res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]


library(dplyr)
library(ggplot2)
library(ggrepel)

highlight_genes1 <- c("VIP","CYYR1","PRG1LP1","FGD5","MEIS1","SALL4","HMGA2","CYYR1-AS1","NKAIN2")
highlight_genes <- custom_GSEA_GMT1$gene
df_q <- res_7q %>%
  as.data.frame() %>%
  select(symbol, lfc_7q = log2FoldChange, padj_7q = pvalue) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = toupper(symbol)) %>%
  inner_join(
    res_GALNT11 %>%
      as.data.frame() %>%
      select(symbol, lfc_GALNT11 = log2FoldChange, padj_GALNT11 = pvalue) %>%
      filter(!is.na(symbol), symbol != "") %>%
      mutate(symbol = toupper(symbol)),
    by = "symbol"
  )

# 4) 标记四象限 & 是否高亮
df_q <- df_q %>%
  mutate(
    is_highlight = symbol %in% highlight_genes,
    is_lable = symbol %in% highlight_genes1
  )

p <- ggplot(df_q, aes(x = lfc_7q, y = lfc_GALNT11)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(
    data = df_q %>% filter(!is_highlight),
    color = "grey75", alpha = 0.6, size = 1
  ) +
  geom_point(
    data = df_q %>% filter(is_highlight),
    color = "#d73027", alpha = 0.95, size = 1.8
  ) +
  ggrepel::geom_text_repel(
    data = df_q %>% filter(is_lable),
    aes(label = symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  labs(
    x = "log2FC (7q group)",
    y = "log2FC (GALNT11 Low vs High)",
    title = "Quadrant plot: 7q vs GALNT11 differential expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+xlim(-7,7)+ylim(-7,7)

p
ggsave(p,file="Figure6E_TCAG_AML_vocan.png",height=4,width=4)
**********************************************************
************************GSE10358**************************
**********************************************************

setwd("./Data_Prepare_update/Figure5/")
degs_7q <- read.csv("./Public_Data/GSE10358_DEGs_7q_vs_NN.csv")
res_7q<- degs_7q
colnames(res_7q)[1] <- "symbol"
degs_GALNT11 <- read.csv("./New_analysis/4—GSE10358_DEGs_GALNT11_low_vs_high.csv")
colnames(degs_GALNT11)[1] <- "symbol"
res_GALNT11 <- degs_GALNT11
res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]
library(dplyr)
library(ggplot2)
library(ggrepel)
highlight_genes1 <- c("VIP","CYYR1","PRG1LP1","FGD5","MEIS1","SALL4","HMGA2","CYYR1-AS1","NKAIN2","TANC1","EFNA1","CNN3","PLSCR4")
highlight_genes <- custom_GSEA_GMT1$gene
df_q <- res_7q %>%
  as.data.frame() %>%
  select(symbol, lfc_7q = logFC, padj_7q = P.Value) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = toupper(symbol)) %>%
  inner_join(
    res_GALNT11 %>%
      as.data.frame() %>%
      select(symbol, lfc_GALNT11 = logFC, padj_GALNT11 = P.Value) %>%
      filter(!is.na(symbol), symbol != "") %>%
      mutate(symbol = toupper(symbol)),
    by = "symbol"
  )
df_q <- df_q %>%
  mutate(
    is_highlight = symbol %in% highlight_genes,
    is_lable = symbol %in% highlight_genes1
  )
p <- ggplot(df_q, aes(x = lfc_7q, y = lfc_GALNT11)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(
    data = df_q %>% filter(!is_highlight),
    color = "grey75", alpha = 0.6, size = 1
  ) +
  geom_point(
    data = df_q %>% filter(is_highlight),
    color = "#d73027", alpha = 0.95, size = 1.8
  ) +
  ggrepel::geom_text_repel(
    data = df_q %>% filter(is_lable),
    aes(label = symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  labs(
    x = "log2FC (7q group)",
    y = "log2FC (GALNT11 Low vs High)",
    title = "Quadrant plot: 7q vs GALNT11 differential expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )+xlim(-4,4)+ylim(-4,4)
p
ggsave(p,file="Supplementary_Figure5A_GSE10358_vocan.png",height=4,width=4)

**********************************************************
************************GSE6891**************************
**********************************************************

degs_7q <- read.csv("./Public_Data/GSE6891_DEGs_7q_vs_NN.csv")
res_7q<- degs_7q
colnames(res_7q)[1] <- "symbol"
degs_GALNT11 <- read.csv("./New_analysis/3_GSE6891_DEGs_GALNT11_low_vs_high.csv")
colnames(degs_GALNT11)[1] <- "symbol"
res_GALNT11 <- degs_GALNT11
res_sig <- read.csv("./New_analysis/7_chemo_nonCR_vs_CR_GSE103424_DEGs.csv")
res_sig1 <- res_sig[res_sig$P.Value<0.05,]
custom_GSEA_GMT <- data.frame(gene=head(res_sig1[order(res_sig1$logFC,decreasing=TRUE),"X"],200))
custom_GSEA_GMT$term <- "Resistance_signature"
custom_GSEA_GMT1 <- custom_GSEA_GMT[,c("term","gene")]

library(dplyr)
library(ggplot2)
library(ggrepel)

highlight_genes1 <- c("VIP","CYYR1","PRG1LP1","FGD5","MEIS1","SALL4","HMGA2","CYYR1-AS1","NKAIN2","TANC1","EFNA1","CNN3","PLSCR4","HEMGN","FAM30A","NPDC1","MMRN1","WASIR2")
highlight_genes <- custom_GSEA_GMT1$gene
df_q <- res_7q %>%
  as.data.frame() %>%
  select(symbol, lfc_7q = logFC, padj_7q = P.Value) %>%
  filter(!is.na(symbol), symbol != "") %>%
  mutate(symbol = toupper(symbol)) %>%
  inner_join(
    res_GALNT11 %>%
      as.data.frame() %>%
      select(symbol, lfc_GALNT11 = logFC, padj_GALNT11 = P.Value) %>%
      filter(!is.na(symbol), symbol != "") %>%
      mutate(symbol = toupper(symbol)),
    by = "symbol"
  )
df_q <- df_q %>%
  mutate(
    is_highlight = symbol %in% highlight_genes,
    is_lable = symbol %in% highlight_genes1
  )
p <- ggplot(df_q, aes(x = lfc_7q, y = lfc_GALNT11)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
  geom_point(
    data = df_q %>% filter(!is_highlight),
    color = "grey75", alpha = 0.6, size = 1
  ) +
  geom_point(
    data = df_q %>% filter(is_highlight),
    color = "#d73027", alpha = 0.95, size = 1.8
  ) +
  ggrepel::geom_text_repel(
    data = df_q %>% filter(is_lable),
    aes(label = symbol),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey60"
  ) +
  labs(
    x = "log2FC (7q group)",
    y = "log2FC (GALNT11 Low vs High)",
    title = "Quadrant plot: 7q vs GALNT11 differential expression"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
ggsave(p,file="Supplementary_Figure5BGSE6891_vocan.png",height=4,width=4)

```
| <img src="/Data_Prepare_update/Figure6/Figure6E_TCAG_AML_vocan.png" width="300"> | <img src="/Data_Prepare_update/Figure6/Supplementary_Figure5A_GSE10358_vocan.png" width="300"> | <img src="/Data_Prepare_update/Figure6/Supplementary_Figure5BGSE6891_vocan.png" width="300">


## Supplementary Fig5C-D
```R

/usr/local/R4.2/bin/R
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
setwd("./Data_Prepare_update/Figure5")
data11 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_1.1_Clinical_summary_186_Patients.xlsx")
data12 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_1.2_Clinical_summary_186_Patients_Description.xlsx")
data13_1 <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_3_Drug_response_sDSS_164S_17Healthy.xlsx")
# data13 <- read_excel("./AML_Drug_response/File_3.2_Drug_response_DSS_sDSS_164S_17Healthy.xlsx")
data13_ID <- read_excel("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/File_2_Drug_library_515D.xlsx")
data14 <- read.csv("./human_data/GSE12417/Functional_Precision_Medicine_Tumor_Board_AML/RNA_seq_AML_norm.counts.csv")
gal <- data14 %>%
  as.data.frame() %>%
  column_to_rownames("X") %>%
  { .["GALNT11", , drop = TRUE] } %>%
  as.numeric()

names(gal) <- colnames(data14)[colnames(data14) != "X"]

# 分位数阈值（你要 0.25 和 0.75）
q25 <- quantile(gal, 0.25, na.rm = TRUE)
q75 <- quantile(gal, 0.75, na.rm = TRUE)

gal_group <- tibble(
  Sample = names(gal),
  GALNT11 = gal,
  Group = case_when(
    GALNT11 <= q25 ~ "Low",
    GALNT11 >= q75 ~ "High",
    TRUE ~ NA_character_
  )
) %>% filter(!is.na(Group))

table(gal_group$Group)

library(tidyr)
library(dplyr)
drug_long <- data13_1 %>%
  pivot_longer(
    cols = -c(Drug_ID, Drug_name),
    names_to = "Sample",
    values_to = "Response"
  ) %>%
  left_join(gal_group, by = "Sample") %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("Low", "High")))

drug_long <- data.frame(drug_long)
drug_long <- drug_long %>%
  mutate(Response = ifelse(is.na(Response), 0, Response))
drug_stats <- drug_long %>%
  group_by(Drug_ID, Drug_name) %>%
  summarise(
    n_low  = sum(Group == "Low"),
    n_high = sum(Group == "High"),
    p = suppressWarnings(wilcox.test(Response ~ Group)$p.value),
    diff_mean = mean(Response[Group=="Low"], na.rm=TRUE) - mean(Response[Group=="High"], na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  arrange(p)

data <- drug_stats[drug_stats$p<0.05,]
data <- data %>%
  mutate(log_p = -log10(p))
data <- data %>%
  mutate(Drug_rank = rank(p))

library(dplyr)
library(ggplot2)

top_drugs <- drug_stats[drug_stats$p<0.05&drug_stats$diff_mean<0,]
data13_ID <- data.frame(data13_ID)
data13_ID[data13_ID$Preferred.name%in%c(top_drugs$Drug_name),]

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

mat <- drug_long %>%
  semi_join(top_drugs, by = c("Drug_ID","Drug_name")) %>%
  select(Drug_name, Sample, Response) %>%
  pivot_wider(names_from = Sample, values_from = Response) %>%
  as.data.frame()

rownames(mat) <- mat$Drug_name
mat$Drug_name <- NULL
mat <- as.matrix(mat)

sample_order <- gal_group %>%
  arrange(factor(Group, levels = c("High","Low")), GALNT11) %>%
  pull(Sample)
sample_order <- intersect(sample_order, colnames(mat))
mat <- mat[, sample_order, drop = FALSE]

mat0 <- mat
mat0[is.na(mat0)] <- 0
# 2) 行 z-score（每个药物内部标准化）
mat_z <- t(scale(t(mat0)))
# mat_z <- mat0
mat_z[is.na(mat_z) | is.nan(mat_z)] <- 0

cap <- 2
mat_z[mat_z >  cap] <-  cap
mat_z[mat_z < -cap] <- -cap

grp_vec <- gal_group$Group[match(colnames(mat_z), gal_group$Sample)]
gal_vec <- gal_group$GALNT11[match(colnames(mat_z), gal_group$Sample)]

ha <- HeatmapAnnotation(
  Group = grp_vec,
  GALNT11 = gal_vec,
  col = list(
    Group = c(Low = "#2b8cbe", High = "#de2d26"),
    GALNT11 = colorRamp2(
      quantile(gal_vec, c(0.05, 0.5, 0.95), na.rm = TRUE),
      c("#f7fbff", "#6baed6", "#08306b")
    )
  ),
  annotation_name_side = "left"
)

col_fun <- colorRamp2(c(-cap, 0, cap), c("#2166ac", "white", "#b2182b"))

# 让 Low / High 视觉上更分明：按组 split，并在组间留空
col_split <- factor(grp_vec, levels = c("High","Low"))

ht <- Heatmap(
  mat_z,
  name = paste0("Row z-score\n(capped ±", cap, ")"),
  top_annotation = ha,
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = col_split,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_title = "AML samples (GALNT11 Low vs High)",
  row_title = paste0("Top ", nrow(mat_z), " drugs")
)

png("Supplementary_5C_Heatmap_0227_GALNT11_response.png")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

library(ggplot2)
library(dplyr)
library(ggpubr)
sel <- c("Pixantrone", "Docetaxel","ABT-751")  # 可以修改为你感兴趣的药名
drug_long$Group <- factor(drug_long$Group,levels=c("High","Low"))
p_violin <- drug_long %>%
  filter(Drug_name %in% sel) %>%
  ggplot(aes(x = Group, y = Response, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  # 小提琴图，trim=FALSE保留完整形状
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.7) +  # 添加散点展示原始数据
  facet_wrap(~ Drug_name, scales = "free_y") +  # 按药物分面，y轴自由缩放
  stat_compare_means(method = "wilcox.test", label = "p.format") +  # 添加Wilcoxon检验p值
  scale_fill_manual(values = c(Low = "#2b8cbe", High = "#de2d26")) +  # 自定义颜色
  theme_bw() +  # 黑白主题
  labs(
    x = NULL, 
    y = "Drug response", 
    title = "Selected drugs: GALNT11 High vs Low"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),  # 标题居中
    legend.position = "bottom"  # 图例放在底部
  )
print(p_violin)
library(ggpubr)
sel <- c("Pixantrone", "Docetaxel","ABT-751")
df_plot <- drug_long %>%
  filter(Drug_name %in% sel) %>%
  filter(!is.na(Response), !is.na(Group)) %>%
  mutate(
    Group = factor(Group, levels = c("High","Low"))  # Low 左，High 右（推荐）
  )
p_violin_pub <- ggplot(df_plot, aes(x = Group, y = Response, fill = Group)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.55,
    linewidth = 0.25
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.9,
    linewidth = 0.3
  ) +
  geom_point(
    aes(color = Group),
    position = position_jitter(width = 0.12, height = 0),
    size = 0.75,
    alpha = 0.25
  ) +
  facet_wrap(~ Drug_name, scales = "free_y", nrow = 1) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    size = 4
  ) +
  scale_fill_manual(values = c(Low = "#4575b4", High = "#d73027")) +
  scale_color_manual(values = c(Low = "#4575b4", High = "#d73027"), guide = "none") +
  labs(
    x = NULL,
    y = "Drug response (DSS / sDSS)",
    title = "Selected drugs: GALNT11 Low vs High"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "bold")
  )
ggsave("Supplementary5D_violin_boxplot_GALNT11_drugs.png", p_violin_pub, width = 6, height = 3)

```

| <img src="/Data_Prepare_update/Figure5/Supplementary_5C_Heatmap_0227_GALNT11_response.png" width="300"> | <img src="/Data_Prepare_update/Figure5/Supplementary5D_violin_boxplot_GALNT11_drugs.png" width="300">
