rm(list=ls())
library(GEOquery)

###GSE13213######
##Download Datasets
aset <- getGEO("GSE13213",destdir = ".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)  ## 平台文件
# exprs()提取matrix；pData提取group info
GSE13213 <- exprs(aset$GSE13213_series_matrix.txt.gz)
GSE13213 <- as.data.frame(GSE13213)
PD13213_test <- pData(aset$GSE13213_series_matrix.txt.gz)

##GSE Platform
library(dplyr)
aset$GSE13213_series_matrix.txt.gz$platform_id #GPL570
GPL6480 <- getGEO("GPL6480")
anno = GPL6480@dataTable@table

##Group Information
#把数据集GSE13213矩阵和分组情况保存下来
write.csv(PD13213_test,'PD13213_test.csv',row.names = F)
PD13213 <- read.csv('PD13213_test2.csv')

#把pd里和dat里取交集 然后糅合到一起 删除dat多余的没pd那些样本的数据（因为pd只有一些数据，dat有很多。）
GSE13213 <- GSE13213[,match(PD13213$geo_accession,colnames(GSE13213))]
# GSE13213 <- log2(GSE13213 + 1)
colnames(anno)
range(GSE13213)

library(magrittr)
anno$`Gene Symbol` <- strsplit(anno$gene_assignment, split = " /// ") %>% lapply(function(x){
  # x = strsplit(dat$gene_assignment[70750:70753], split = " /// ")[[1]]
  x1 <- x
  if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
  if(length(x1) == 0) x1 <- x[1]
  vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
    unlist() %>% unique()
  vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
  vec_gene <- paste0(vec_gene, collapse = " /// ")
  return(vec_gene)
}) %>% unlist()
anno$GeneSymbol %>% table() %>% sort(decreasing = T) %>% .[1:10]

GPLanno <- anno[,c("ID", "Gene Symbol")]
GSE13213$ID = rownames(GSE13213) #新建了一个ID用于当做新的一个列 

GSE13213 <- merge(GPLanno,GSE13213, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
GSE13213 <- na.omit(GSE13213)
GSE13213 <- aggregate(x=GSE13213,by = GSE13213$`Gene Symbol` %>% list(),FUN = mean)#多个探针对应一个基因取平均
rownames(GSE13213) = GSE13213$Group.1 #把Symbol变成行名
GSE13213 <- GSE13213[-1,-(1:3)]
GSE13213 <- GSE13213[!grepl("///", row.names(GSE13213)),]
range(GSE13213)
# GSE13213 <- log2(GSE13213+1)

write.csv(GSE13213,"GSE13213_Matrix.csv",row.names = T)
write.csv(PD13213,"GSE13213_Group.csv",row.names = F)

# write.csv(row.names(GSE13213),'GSE13213_Gene.csv',row.names = F)


##Download Datasets
aset <- getGEO("GSE42568",destdir = ".",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)  ## 平台文件 
# exprs()提取matrix；pData提取group info
GSE42568 <- exprs(aset$GSE42568_series_matrix.txt.gz)
GSE42568 <- as.data.frame(GSE42568)
PD42568_test <- pData(aset$GSE42568_series_matrix.txt.gz)

##GSE Platform
library(dplyr)
aset$GSE42568_series_matrix.txt.gz$platform_id #GPL570
GPL570 <- getGEO("GPL570")
anno = GPL570@dataTable@table

##Group Information
#把数据集GSE13213矩阵和分组情况保存下来
write.csv(PD42568_test,'PD42568_test.csv',row.names = F)
PD42568 <- read.csv('PD42568_test2.csv')

#把pd里和dat里取交集 然后糅合到一起 删除dat多余的没pd那些样本的数据（因为pd只有一些数据，dat有很多。）
GSE42568 <- GSE42568[,match(PD42568$geo_accession,colnames(GSE42568))]
# GSE13213 <- log2(GSE13213 + 1)
colnames(anno)
range(GSE42568)

library(magrittr)
anno$`Gene Symbol` <- strsplit(anno$gene_assignment, split = " /// ") %>% lapply(function(x){
  # x = strsplit(dat$gene_assignment[70750:70753], split = " /// ")[[1]]
  x1 <- x
  if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
  if(length(x1) == 0) x1 <- x[1]
  vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
    unlist() %>% unique()
  vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
  vec_gene <- paste0(vec_gene, collapse = " /// ")
  return(vec_gene)
}) %>% unlist()
anno$GeneSymbol %>% table() %>% sort(decreasing = T) %>% .[1:10]

GPLanno <- anno[,c("ID", "Gene Symbol")]
GSE42568$ID = rownames(GSE42568) #新建了一个ID用于当做新的一个列 

GSE42568 <- merge(GPLanno,GSE42568, by='ID',all.x=T,all.y=T) #合并两个表 取有id的交集
GSE42568 <- na.omit(GSE42568)
GSE42568 <- aggregate(x=GSE42568,by = GSE42568$`Gene Symbol` %>% list(),FUN = mean)#多个探针对应一个基因取平均
rownames(GSE42568) = GSE42568$Group.1 #把Symbol变成行名
GSE42568 <- GSE42568[-1,-(1:3)]
GSE42568 <- GSE42568[!grepl("///", row.names(GSE42568)),]
range(GSE42568)
# GSE13213 <- log2(GSE13213+1)

write.csv(GSE42568,"GSE42568_Matrix.csv",row.names = T)
write.csv(PD42568,"GSE42568_Group.csv",row.names = F)

write.csv(row.names(GSE42568),'GSE42568_Gene.csv',row.names = F)
write.csv(row.names(GSE13213),'GSE13213_Gene.csv',row.names = F)

save(GSE42568,PD42568,GSE13213,PD13213,file = 'Valid_set.RData')















