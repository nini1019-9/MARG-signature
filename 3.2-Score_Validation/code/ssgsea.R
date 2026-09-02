rm(list=ls())
options(stringsAsFactors = F)

library(GSVA)
library(GSEABase)
library(tidyverse)

## 数据导入
load('data/GEOdata.RData')
phen_gene <- data.table::fread("data/MARGs.csv", data.table = F)
phen_name <- "MARGs"

## 基因集构建
geneSet <- phen_gene[,1] %>% 
  data.frame(check.names = F) %>% 
  `names<-`(phen_name)

## ssgsea----
rs1 <- gsva(expr = as.matrix(GSE13213_Matrix),
           gset.idx.list = geneSet,
           method = "ssgsea", 
           kcdf = 'Gaussian', 
           abs.ranking = TRUE)
score_rs1 <- rs1 %>% t() %>% data.frame(check.names = F)

## km----
#!!!确保OS.time是以日为单位
surv_score1 <- score_rs1 %>% 
  dplyr::mutate(Sample = rownames(score_rs1)) %>%
  dplyr::inner_join(GSE13213_Group) %>% 
  dplyr::mutate(OS.time = OS.time/365)


library(survminer)
res_cut <- surv_cutpoint(surv_score1, time = "OS.time", event = "OS", variables =phen_name)
res_cut[["cutpoint"]]
# res_cat <- surv_categorize(res_cut)
surv_score1 <- surv_score1 %>% 
  dplyr::mutate(group = dplyr::if_else(get(phen_name) > res_cut[["cutpoint"]]$cutpoint, "HighScore", "LowScore") %>% 
                  factor(levels = c("LowScore", "HighScore"))) %>% 
  dplyr::arrange(MARGs)
colnames(surv_score1)[1] <- 'MA.Score'
write.csv(surv_score1, "output/0-ssgsea_GSE13213.csv", row.names = F)
write.csv(surv_score1[,c(4,5,1)], "output/1-KM_GSE13213.csv", row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_score1)

p <- ggsurvplot(fit,
                data = surv_score1, 
                conf.int = FALSE,
                pval.size = 5,
                risk.table=F,
                legend.title = "Group",
                legend.labs = c("LowScore", "HighScore"), # 指定图例分组标签
                xlab = "Follow up time(y)", # 指定x轴标签
                palette = c("#D3E2B7","#eca680"),
                ggtheme = theme_classic(base_size = 17),
                pval = T)
p
ggsave("output/1-km_ssgsea_GSE13213.pdf", width = 5.3, height = 5)



## ssgsea----
rs2 <- gsva(expr = as.matrix(GSE30219_Matrix),
            gset.idx.list = geneSet,
            method = "ssgsea", 
            kcdf = 'Gaussian', 
            abs.ranking = TRUE)
score_rs2 <- rs2 %>% t() %>% data.frame(check.names = F)

## km----

surv_score2 <- score_rs2 %>% 
  dplyr::mutate(Sample = rownames(score_rs2)) %>%
  dplyr::inner_join(GSE30219_Group) %>% 
  dplyr::mutate(OS.time = OS.time/365)


library(survminer)
res_cut <- surv_cutpoint(surv_score2, time = "OS.time", event = "OS", variables =phen_name)
res_cut[["cutpoint"]]
# res_cat <- surv_categorize(res_cut)
surv_score2 <- surv_score2 %>% 
  dplyr::mutate(group = dplyr::if_else(get(phen_name) > res_cut[["cutpoint"]]$cutpoint, "HighScore", "LowScore") %>% 
                  factor(levels = c("LowScore", "HighScore"))) %>% 
  dplyr::arrange(MARGs)
colnames(surv_score2)[1] <- 'MA.Score'
write.csv(surv_score2, "output/0-ssgsea_GSE30219.csv", row.names = F)
write.csv(surv_score2[,c(4,5,1)], "output/2-KM_GSE30219.csv", row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_score2)

p <- ggsurvplot(fit,
                data = surv_score2, 
                conf.int = FALSE,
                pval.size = 5,
                risk.table=F,
                legend.title = "Group",
                legend.labs = c("LowScore", "HighScore"), # 指定图例分组标签
                xlab = "Follow up time(y)", # 指定x轴标签
                palette = c("#D3E2B7","#eca680"),
                ggtheme = theme_classic(base_size = 17),
                pval = T)
p
ggsave("output/2-km_ssgsea_GSE30219.pdf", width = 5.3, height = 5)


######################
rs3 <- gsva(expr = as.matrix(GSE11969_Matrix),
            gset.idx.list = geneSet,
            method = "ssgsea", 
            kcdf = 'Gaussian', 
            abs.ranking = TRUE)
score_rs3 <- rs3 %>% t() %>% data.frame(check.names = F)

## km----

surv_score3 <- score_rs3 %>% 
  dplyr::mutate(Sample = rownames(score_rs3)) %>%
  dplyr::inner_join(GSE11969_Group) %>% 
  dplyr::mutate(OS.time = OS.time/365)


library(survminer)
res_cut <- surv_cutpoint(surv_score3, time = "OS.time", event = "OS", variables =phen_name)
res_cut[["cutpoint"]]
# res_cat <- surv_categorize(res_cut)
surv_score3 <- surv_score3 %>% 
  dplyr::mutate(group = dplyr::if_else(get(phen_name) > res_cut[["cutpoint"]]$cutpoint, "HighScore", "LowScore") %>% 
                  factor(levels = c("LowScore", "HighScore"))) %>% 
  dplyr::arrange(MARGs)
colnames(surv_score3)[1] <- 'MA.Score'
write.csv(surv_score3, "output/0-ssgsea_GSE11969.csv", row.names = F)
write.csv(surv_score3[,c(4,5,1)], "output/2-KM_GSE11969.csv", row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_score3)

p <- ggsurvplot(fit,
                data = surv_score3, 
                conf.int = FALSE,
                pval.size = 5,
                risk.table=F,
                legend.title = "Group",
                legend.labs = c("LowScore", "HighScore"), # 指定图例分组标签
                xlab = "Follow up time(y)", # 指定x轴标签
                palette = c("#D3E2B7","#eca680"),
                ggtheme = theme_classic(base_size = 17),
                pval = T)
p
ggsave("output/3-km_ssgsea_GSE11969.pdf", width = 5.3, height = 5)
