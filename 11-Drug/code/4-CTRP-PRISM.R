rm(list = ls())
library(tidyverse) # 用于读取MAF文件
library(ISOpureR) # 用于纯化表达???
library(impute) # 用于KNN填补药敏数据
# install.packages("data/pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
library(pRRophetic) # 用于药敏预测
trace(calcPhenotype, edit = T) #在calcPhenotype函数的第6，8，10，12和66行后面加个[1]
trace(summarizeGenesByMean, edit = T)# summarizeGenesByMean函数的第19行后面加个[1]
library(SimDesign) # 用于禁止药敏预测过程输出的信???
library(ggplot2) # 绘图
library(cowplot) # 合并图像
library(data.table)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}
load('data/TCGA_Data.RData')
load('data/rs_risk.RData')
expr=dat_fpkm
class(expr)
runpure <- F # 如果想运行就把这个改为T
if(runpure) {
  set.seed(123)
  # Run ISOpureR Step 1 - Cancer Profile Estimation
  ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
  # For reproducible results, set the random seed
  set.seed(456);
  # Run ISOpureR Step 2 - Patient Profile Estimation
  ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
  pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

Sinfo <- rs_risk[[1]] %>% 
  dplyr::mutate(id = rownames(.)) %>% 
  dplyr::select(id, RS, Risk) 
Sinfo=as.data.frame(Sinfo)
rownames(Sinfo)=Sinfo$id
if(!runpure) {
  pure.tumoexpr <- expr
}

# 
load("data/train_ctrp.RData")
# 测试
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5 # 纯化的测试集取表达稳定的基因
# testExpr <- log2(pure.tumoexpr[keepgene,] + 1) # 表达谱对数化
testExpr <- pure.tumoexpr[keepgene,] 
# 取训练集和测试集共有的基???
trainExpr <- trainExpr_ctrp
trainPtype <- trainPtype_ctrp
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
testExpr=as.matrix(testExpr)
outTab <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能???0值，因此加一个较小的数值防止报???
  
  # 岭回归预测药物敏感???
  # ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
  #                                 trainingPtype = tmp,
  #                                 testExprData = testExpr,
  #                                 powerTransformPhenotype = F,
  #                                 selection = 1))
  ptypeOut <- calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1)
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对???
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab

## prism---
load("data/train_prism.RData")
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
# testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
testExpr <- pure.tumoexpr[keepgene,]
testExpr=as.matrix(testExpr)
trainExpr <- trainExpr_prism
trainPtype <- trainPtype_prism
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- testExpr[comgene,]
outTab <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainPtype)) { 
  display.progress(index = i,totalN = ncol(trainPtype))
  d <- colnames(trainPtype)[i]
  tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于PRISM的AUC可能???0值，因此加一个较小的数值防止报???
  # ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
  #                                 trainingPtype = tmp,
  #                                 testExprData = testExpr,
  #                                 powerTransformPhenotype = F,
  #                                 selection = 1))
  ptypeOut <- calcPhenotype(trainingExprData = trainExpr,
                                  trainingPtype = tmp,
                                  testExprData = testExpr,
                                  powerTransformPhenotype = F,
                                  selection = 1)
  ptypeOut <- 2^ptypeOut - 0.00001 # 反对???
  outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
prism.pred.auc <- outTab

##
Sinfo <- Sinfo[,c(1,2,3)]
top.pps <- Sinfo[Sinfo$Risk=="HighRisk",] # 定义上十分位的样???
bot.pps <- Sinfo[Sinfo$Risk=="LowRisk",] # 定义下十分位的样???

ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均???
  b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均???
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
range(ctrp.log2fc)
candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.1] # 这里我调整了阈值，控制结果数目

prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均???
  b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均???
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
range(prism.log2fc)
candidate.prism <- prism.log2fc[prism.log2fc > 0.1] #和前面最好一致

ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$RS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
range(ctrp.cor)
candidate.ctrp2.sig <- ctrp.cor.p[ctrp.cor.p < 0.05]
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.1]  # 这里我调整了阈值，控制结果数目
candidate.ctrp2 <- candidate.ctrp2.sig[names(candidate.ctrp2.sig) %in% names(candidate.ctrp2)]


prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$RS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
range(prism.cor)
candidate.prism2.sig <- prism.cor.p[prism.cor.p < 0.05]
candidate.prism2 <- prism.cor[prism.cor < -0.1]  #和前面最好一致
candidate.prism2 <- candidate.prism2.sig[names(candidate.prism2.sig) %in% names(candidate.prism2)]

prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
prism.candidate_cor <- prism.cor[prism.candidate]
prism.candidate_cor <- prism.candidate_cor[order(prism.candidate_cor)]
prism.candidate <- names(prism.candidate_cor)

ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))
ctrp.candidate_cor <- ctrp.cor[ctrp.candidate]
ctrp.candidate_cor <- ctrp.candidate_cor[order(ctrp.candidate_cor)]
ctrp.candidate <- names(ctrp.candidate_cor)

## plot----
darkblue <- "#eca680"
lightblue <- "#D3E2B7"
corcol <- "#F5E1D8"

ctrp.cor.p[ctrp.cor.p == 0] = min(ctrp.cor.p[ctrp.cor.p != 0])/2
cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
# rownames(cor.data) <- ctrp.candidate2
# cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = corcol) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(
                  expand = expansion(mult = c(0.01,0.01))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~" value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank(),
        text = element_text(size = 20))

# prism.candidate_top10 <- prism.candidate[1:10]
prism.cor.p[prism.candidate][prism.cor.p[prism.candidate] == 0] <- min(prism.cor.p[prism.candidate][prism.cor.p[prism.candidate] != 0])/2
cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = corcol) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(expand = expansion(mult = c(0.01,0.01))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~" value)")) + 
  theme(legend.position = "bottom",
        axis.line.y = element_blank(),
        text = element_text(size = 20)
        )

ctrp.pred.auc1 <- ctrp.pred.auc[ctrp.candidate,Sinfo$id]
ctrp.pred.auc1 <- t(ctrp.pred.auc1)
ctrp.pred.auc1 <- cbind(Sinfo$Risk,ctrp.pred.auc1)
write.csv(ctrp.pred.auc1,'ctrp.pred.auc.csv',row.names = F)

ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","ns")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = factor(rep(c("LowRisk","HighRisk"),c(nrow(bot.pps),nrow(top.pps))),levels = c("LowRisk","HighRisk")),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
# write.csv(ctrp.pred.auc,'ctrp_auc_boxplot.csv') #拼图太慢仙桃画图
# ctrp.boxdata$drug <- sapply(strsplit(ctrp.boxdata$drug," (",fixed = T), "[",1)

p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(lightblue,darkblue)) + 
  scale_color_manual(values = c(lightblue,darkblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 16),
        legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 20)) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.pred.auc1 <- prism.pred.auc[prism.candidate,Sinfo$id]
prism.pred.auc1 <- t(prism.pred.auc1)
prism.pred.auc1 <- cbind(Sinfo$Risk,prism.pred.auc1)
write.csv(prism.pred.auc1,'prism.pred.auc.csv',row.names = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(bot.pps)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(top.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","ns")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                    data.frame(drug = d,
                                               auc = c(a,b),
                                               p = p,
                                               s = s,
                                               group = factor(rep(c("LowRisk","HighRisk"),c(nrow(bot.pps),nrow(top.pps))),levels = c("LowRisk","HighRisk")),
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(lightblue,darkblue)) + 
  scale_color_manual(values = c(lightblue,darkblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 16),
        legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 20))
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
require(cowplot)
plot_grid(p1, p3, p2, p4, labels=c("A", "B", "C", "D"), 
          ncol=2, 
          rel_widths = c(2, 2)) #左右两列的宽度比???
ggsave(filename = "output/drug target2.pdf",width = 16.2,height = 10)

## cmap----
Cmap=fread("data/export.txt")
prism.candidate <- sapply(strsplit(prism.candidate," (",fixed = T), "[",1)
ctrp.candidate <- sapply(strsplit(ctrp.candidate," (",fixed = T), "[",1)

intersect(Cmap$Name, prism.candidate)
intersect(Cmap$Name,ctrp.candidate)
intersect(ctrp.candidate,prism.candidate)
target1 <- Cmap[Cmap$Name %in% intersect(Cmap$Name,ctrp.candidate),]
target2 <- Cmap[Cmap$Name %in% intersect(Cmap$Name,prism.candidate),]
target3 <- Cmap[Cmap$Name %in% intersect(ctrp.candidate,prism.candidate),]
target <- rbind(target1,target2)
write.csv(target, "output/Cmap_filter.csv")


write.csv(prism.candidate, "output/prism_candidate.csv",row.names = F)
write.csv(ctrp.candidate, "output/ctrp_candidate.csv",row.names = F)

cmap_gene <- read.csv('output/Cmap_gene.csv')
hub_gene <- read.csv('data/hub_genes.csv')
cmap_gene2 <- limma::alias2Symbol(cmap_gene$Gene,expand.symbols = T)
hub_gene2 <- limma::alias2Symbol(hub_gene$x,expand.symbols = T)
intersect(cmap_gene2,hub_gene2)
