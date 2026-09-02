rm(list = ls())
library(tidyverse)
library(survival) 
library(data.table)
library(glmnet)
exp=fread("data/IMvogpr210_tpm.csv")
exp=as.data.frame(exp)
rownames(exp)=exp$V1
exp=exp[,-1]
cli=fread("data/Imvigor210_pdata.csv")
cli=as.data.frame(cli)
rownames(cli)=cli$V1
cli=cli[,c(3,22,23)]
colnames(cli)=c("Response","OS.time","OS")
cli=na.omit(cli)
samesample=intersect(colnames(exp),rownames(cli))
exp=t(exp)
exp=exp[samesample,]
exp=as.data.frame(exp)
cli=cli[samesample,]
table(cli$Response)
load("./data/RSF + plsRcox.RData")
hub_genes <- read.csv('data/hub_genes.csv') #因为基因没有变少，所以直接用unicox筛选出的基因，一般要用ML筛选出的hub基因
list_sur_gene[[1]] <- gsub('-','_',list_sur_gene[[1]])
colnames(exp) <- gsub('-','_',colnames(exp))
exp=cbind(cli[,2:3],exp[,hub_genes$x])


riskScore = as.numeric(predict(fit, type = "lp", newdata = exp[, -c(1, 2)]))
Risk=as.vector(ifelse(riskScore>median(riskScore),"HighRisk","LowRisk"))
outTab=cbind(exp[,c(1,2)], riskScore=as.vector(riskScore),Risk)
write.table(cbind(id=rownames(outTab),outTab),file="output/risk.IMv.txt",sep="\t",quote=F,row.names=F)
# library(ggstatsplot)
library(ggplot2)
library(ggpubr)
imminput=cbind(outTab,cli[,1])
write.csv(imminput[,c(5,3)],'IMV_plot.csv',row.names = F)
colnames(imminput)[ncol(imminput)]="Response"
pdf("output/immunotherapy.pdf",height=6,width=6)
# ggplot2分组比较图比较Response组间的riskScore,添加颜色
ggplot(imminput, aes(x = Response, y = riskScore, fill = Response)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.tes", paired = FALSE) +
  stat_compare_means(label.y = 0.5) +
  # scale_fill_manual(c("blue", "red")) +
  theme_bw()
dev.off()
