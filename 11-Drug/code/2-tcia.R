rm(list = ls())
load("./data/rs_risk.RData")
load("./data/TCGA_Data.RData")
group_tcga <- rs_risk[[1]] %>% 
  dplyr::mutate(id = rownames(.)) %>% 
  dplyr::select(id, Risk) %>% 
  `names<-`(c("id", "group"))     
TCIA<-read.delim("data/TCIA-ClinicalData.tsv")
rownames(TCIA)<-TCIA$barcode
TCIA <- TCIA[stringr::str_sub(group_tcga$id,1,12),]
TCIA$level <- group_tcga$group#具体分组信息

library(ggplot2)
library(ggpubr)
library(patchwork)
#绘图
plot1<-function(x){
  ggplot(TCIA, aes(x=level, y=get(x),fill=level)) +
    geom_violin()+
    geom_boxplot(width=0.2,outlier.size=0)+
    theme_classic(base_size = 15)+
    geom_signif(comparisons = list(c("HighRisk","LowRisk")),step_increase = 0,
                map_signif_level =T,test = wilcox.test,textsize = 6)+
    theme(legend.position = "none",
          legend.title=element_blank())+
    scale_fill_manual(values=c("#A6B5C8", "#E4A6BD"))+
    scale_color_manual(values=c("#A6B5C8", "#E4A6BD"))+
    labs(x="",y=x)
}
#拼图
p1<-plot1("ips_ctla4_neg_pd1_neg")
p2<-plot1("ips_ctla4_neg_pd1_pos")
p3<-plot1("ips_ctla4_pos_pd1_neg")
p4<-plot1("ips_ctla4_pos_pd1_pos")
p1+p2+p3+p4
ggsave("output/TCIA.pdf",height = 10,width = 10)

TCIA1 <- TCIA[,c('level','ips_ctla4_neg_pd1_neg')]
TCIA2 <- TCIA[,c('level','ips_ctla4_neg_pd1_pos')]
TCIA3 <- TCIA[,c('level','ips_ctla4_pos_pd1_neg')]
TCIA4 <- TCIA[,c('level','ips_ctla4_pos_pd1_pos')]

write.csv(TCIA1,'TCIA1.csv',row.names = F)
write.csv(TCIA2,'TCIA2.csv',row.names = F)
write.csv(TCIA3,'TCIA3.csv',row.names = F)
write.csv(TCIA4,'TCIA4.csv',row.names = F)

TCIA_all <- TCIA[,c('level','ips_ctla4_neg_pd1_neg','ips_ctla4_neg_pd1_pos',
                 'ips_ctla4_pos_pd1_neg','ips_ctla4_pos_pd1_pos')]
write.csv(TCIA_all,'TCIA_all.csv',row.names = F)



#####全代码自动实现
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("IOBR", quietly = TRUE)){
  devtools::install_github("IOBR/IOBR")
}
  

library(IOBR)
library(tidyverse)
library(ggpubr)

ips<-deconvo_tme(eset = dat_fpkm, method = "ips", plot= FALSE)
head(ips)
ips <- ips[ips$ID %in% group_tcga$id,]
group_tcga <- group_tcga[ips$ID,]
ips$Group <- group_tcga$group

write.csv(ips,'output/IPS.csv',row.names = F)
