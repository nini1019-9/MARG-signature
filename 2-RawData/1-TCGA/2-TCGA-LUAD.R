rm(list = ls())
library(dplyr)

###READ PD INFORMATION
PD <- data.table::fread("PD.csv",data.table = F)
colnames(PD)[1] <- "sample_id"

###FPKM DATA######
dat_fpkm <- data.table::fread("TCGA-LUAD_FPKM.txt",check.names = F,data.table = F)
dat_fpkm <- aggregate(dat_fpkm, by = dat_fpkm$gene_id %>% list(),FUN = mean)
rownames(dat_fpkm) <- dat_fpkm$Group.1
dat_fpkm <- dat_fpkm[,-(1:2)]

range(dat_fpkm)
dat_fpkm <- log2(dat_fpkm+1)
# dat_fpkm_back <- dat_fpkm
# colnames(dat_fpkm) <-  stringr::str_sub(colnames(dat_fpkm),1,12)

dat_fpkm <- dat_fpkm[,PD$sample_id]
write.csv(dat_fpkm[,PD$Status == "LUAD"],"Disease_Matrix_FPKM.csv")
write.csv(dat_fpkm,"Matrix_FPKM.csv")

###COUNTS DATA######
dat_counts <- data.table::fread("TCGA-LUAD_Count.txt",check.names = F,data.table = F)
dat_counts <- aggregate(dat_counts, by = dat_counts$gene_id %>% list(),FUN = mean)
rownames(dat_counts) <- dat_counts$Group.1
dat_counts <- dat_counts[,-(1:2)]
# colnames(dat_counts) <-  stringr::str_sub(colnames(dat_counts),1,12)

range(dat_counts)
# dat_counts <- log2(dat_counts+1)
dat_counts <- dat_counts[,PD$sample_id]
# dat_counts <- dat_counts[,group_filter$status == "Tumor"]
write.csv(dat_counts[,PD$Status == "LUAD"],"Disease_Matrix_Counts.csv")
write.csv(dat_counts,"Matrix_Counts.csv")

save(PD,dat_counts,dat_fpkm,file = 'TCGA_Data.RData')

