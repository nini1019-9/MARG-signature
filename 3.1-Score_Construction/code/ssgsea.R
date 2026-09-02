rm(list=ls())
options(stringsAsFactors = F)

library(GSVA)
library(GSEABase)
library(tidyverse)

## 数据导入
load('data/TCGA_Data.RData')
load('data/GEOdata.RData')
phen_gene <- read.csv('data/MARGs.csv')
phen_name <- "MARGS"


## 基因集构建
geneSet <- phen_gene[,1] %>% 
  data.frame(check.names = F) %>% 
  `names<-`(phen_name)

## ssgsea----
rs <- gsva(expr = as.matrix(dat_fpkm),
           gset.idx.list = geneSet,
           method = "ssgsea", 
           kcdf = 'Gaussian', 
           abs.ranking = TRUE)
score_rs <- rs %>% t() %>% data.frame(check.names = F)

## km----
#！！！！！注意TCGA和GEO清洗的OS.time是不是单位为天
surv_score <- score_rs %>% 
  dplyr::mutate(sample_id = rownames(score_rs)) %>%
  dplyr::inner_join(PD) %>% 
  dplyr::mutate(OS.time = OS.time/365)


library(survminer)
res_cut <- surv_cutpoint(surv_score, time = "OS.time", event = "OS", variables =phen_name)
res_cut[["cutpoint"]]
# res_cat <- surv_categorize(res_cut)
surv_score <- surv_score %>% 
  dplyr::mutate(group = dplyr::if_else(get(phen_name) > res_cut[["cutpoint"]]$cutpoint, "HighScore", "LowScore") %>% 
                  factor(levels = c("LowScore", "HighScore"))) %>% 
  dplyr::arrange(MARGS)
colnames(surv_score)[1] <- "MA.Score"
write.csv(surv_score, "output/0-ssgsea_tcga.csv", row.names = F)
write.csv(surv_score[,c(5,6,1)], "output/1-KM.csv", row.names = F)

library(survival)
fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_score)

p <- ggsurvplot(fit,
                data = surv_score, 
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
ggsave("output/1-km_ssgsea_tcga.pdf", width = 5.3, height = 5)

## time_roc----

library(timeROC)
time_roc_res <- timeROC(
  T = surv_score$OS.time*365,
  delta = surv_score$OS,
  marker = surv_score$MA.Score,
  cause = 1,
  times = c(1 * 365, 2 * 365, 3 * 365),
  ROC = TRUE,
  iid = TRUE
)
time_ROC_df <- data.frame(
  TP_1year = time_roc_res$TP[, 1],
  FP_1year = time_roc_res$FP[, 1],
  TP_2year = time_roc_res$TP[, 2],
  FP_2year = time_roc_res$FP[, 2],
  TP_3year = time_roc_res$TP[, 3],
  FP_3year = time_roc_res$FP[, 3]
)
k1<-ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#67ADB7") +
  geom_line(aes(x = FP_2year, y = TP_2year), size = 1, color = "#F5E1D8") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#E6E0B0") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw(base_size = 17) +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#67ADB7"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 2 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#F5E1D8"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E6E0B0"
  ) +
  labs(x = "1-specificity", y = "Sensitivity") +
  theme(
    panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black")
    
  )
k1
ggsave("output/timeroc_ssgsea_tcga.pdf", width = 5.3, height = 5)

## genome----
library(cBioPortalData)
cbio <- cBioPortal()
studies = getStudies(cbio)
head(studies$studyId)
table(studies$cancerTypeId)
id = c("luad_tcga_pan_can_atlas_2018") # 改为自己的
clinical = clinicalData(cbio, id)
cli_df <- clinical[,c("patientId","sampleId","MSI_SENSOR_SCORE","TMB_NONSYNONYMOUS", "MUTATION_COUNT", "FRACTION_GENOME_ALTERED")]
colnames(cli_df)[3] = "MSI_score" 
cli_df[,c(3:6)] <- cli_df[,c(3:6)] %>% 
  apply(2, function(x){as.numeric(x)})
cli_df <- cli_df %>% 
  dplyr::mutate(TMB_NONSYNONYMOUS = log10(TMB_NONSYNONYMOUS))
write.csv(cli_df, "output/cli_df.csv", row.names = F)

library(ggplot2)
library(ggpubr)
surv_score <- data.table::fread("output/0-ssgsea_tcga.csv", data.table = F)
dat_genome <- surv_score[,c(2, 1, 12)] %>% 
  dplyr::mutate(sample = stringr::str_sub(sample_id, 1, 15)) %>% 
  dplyr::inner_join(cli_df[,c(2:6)], by = c("sample" = "sampleId")) %>% 
  dplyr::select(-2)
write.csv(dat_genome, "output/dat_genome.csv", row.names = F)

id <- c('MSI','TMB','Mutation','FGA')
dat_genome$group <- factor(dat_genome$group, levels = as.vector(unique(dat_genome$group)))
dat_genome$MUTATION_COUNT <- log2(dat_genome$MUTATION_COUNT + 1)
colnames(dat_genome)[4:7] <- c('MSI','TMB','Mutation','FGA')
dat_genome$MSI <- log2(dat_genome$MSI + 1)



for (i in id) {
  # i = "MSI"
  p <- ggplot(dat_genome[,c("group",i)] %>% na.omit(), aes(x = group,y = get(i), fill = group)) +
    geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) +
    scale_fill_manual(values = c("#D3E2B7","#eca680"))  + 
    # geom_jitter(size = 0.5) +
    stat_compare_means(method = "wilcox.test", size = 6 * 0.35, label.x.npc = 0.5,
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                      symbols = c("***", "**", "*", "ns")), label = "p.signif") +
    labs( y = paste0(i,' Score')) +
    # geom_hline(yintercept = 0.4, lty = 4) +
    theme_bw() + 
    theme(panel.grid=element_blank(),legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.text = element_text(size = 10, color = "black"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
        # axis.title.y = element_blank(), # 隐藏y轴标签
        # axis.text.x = element_blank(), # 隐藏x轴内容
        # axis.ticks.x = element_blank())
    )
 # 隐藏x轴刻度线
  filename <- paste0("output/genome_", i, ".pdf")
  ggsave(filename, p, width = 5.3, height = 5)
}

write.csv(dat_genome[,c(2,4:7)],'2-TMB-MSI.csv')
write.csv(dat_genome,'dat_genome.csv')

