rm(list = ls())
library(dplyr)

#确认样本数量 alive = 0, dead = 1
id <- readxl::read_xlsx("TCGA-LUAD_sampleid.xlsx")
clinical_old <- readxl::read_xlsx("TCGA-LUAD_clinicalold.xlsx")
#去掉没有分组信息肿瘤和正常分组信息的样本
clinical_old <- clinical_old[!is.na(clinical_old$status),]
table(clinical_old$status)
LUAD <- clinical_old[clinical_old$status == "Tumor",]
LUAD <- LUAD[,c(1,2,6,183,184,180,182,160:162)]
LUAD <- LUAD[LUAD$RNAseq样本编号 %in% id$sample_id,]

setdiff(id$sample_id,LUAD$RNAseq样本编号)
setdiff(LUAD$RNAseq样本编号,id$sample_id)
# #[1] "TCGA-A8-A09C-01A-11R-A00Z-07" "TCGA-BH-A0BS-01A-11R-A12P-07"
# #[3] "TCGA-E2-A14S-01A-11R-A12D-07" "TCGA-LD-A7W5-01A-22R-A352-07"
# 
# grep("TCGA-A8-A09C-01A-11R-A00Z-07",id$sample_id)
# # [1] 233
# grep("TCGA-BH-A0BS-01A-11R-A12P-07",id$sample_id)
# # [1] 542
# grep("TCGA-E2-A14S-01A-11R-A12D-07",id$sample_id)
# # [1] 806
# grep("TCGA-LD-A7W5-01A-22R-A352-07",id$sample_id)
# # [1] 1014

# id <- id[c(233,542,806,1014),]
# id$sample <- stringr::str_sub(id$sample_id,1,12)
# 
# clinical_new <- readxl::read_xlsx("TCGA-GBM_rnaseq_clinical_new.xlsx")
# clinical_new <- clinical_new[clinical_new$sample_id %in% id$sample,]
# clinical_new_stage <- clinical_new[,c(1,)]
# write.csv(clinical_new,"new sample.csv")

Control <- clinical_old[clinical_old$status == "Normal",]
Control <- Control[,c(1,2,6,183,184,180,182,160:162)]

LUAD$status <- gsub("Tumor","LUAD",LUAD$status)

PD <- rbind(Control,LUAD)
PD <- na.omit(PD)
write.csv(PD,"PD.csv",row.names = F)
