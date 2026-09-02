

#进行小鼠和人的基因转换
rm(list = ls())
#安装homologene这个R包
#install.packages('homologene')
#加载homologene这个R包
library(homologene)


#读取数据集矩阵，矩阵里面的数据不重要，主要是小鼠ID列
#数据集矩阵是已经完成探针替换的表达矩阵
data <- data.table::fread("GEO-adjusted-mouse.csv")
#更多基因方法是一样的,先制作出一个genelist
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]
genelist<-rownames(data)
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
genelistid <- homologene(genelist, inTax = 10090, outTax = 9606)
genelistid2 <- genelistid[,c(1,2)]
colnames(genelistid2) <- c("ID","GENE")
genelistid2 <- as.data.frame(genelistid2)
data$ID <- row.names(data)

data <- merge(genelistid2,data, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
data <- na.omit(data)
data2 <- aggregate(x=data,by = data$GENE %>% list(),FUN = mean)#多个探针对应一个基因取平均
row.names(data2) <- data2$Group.1
data2 <- data2[,-c(1,2,3)]

write.csv(data2,"1-GEO-adjusted-human-id.csv")




















