rm(list = ls())
## 1.下载原始数据----
library(GEOquery)
library(stringr)
### 如果出现网络问题导致R中下载报错，可以网页上直接下载，存在相同的文件夹
getGEOSuppFiles("GSE11969", 
                baseDir ="data/raw", # 文件存放位置
                fetch_files = T,
                makeDirectory = T) # makeDirectory = T新建一个子文件夹
untar(tarfile = "data/raw/GSE11969/GSE11969_RAW.tar", exdir = "data/raw/GSE11969") # 解压
file.remove("data/raw/GSE11969/GSE11969_RAW.tar") # 删除压缩包

## 2.样本文件准备----
group_geo <- data.table::fread("data/raw/GSE11969_Group.csv", data.table = F) # 自己外部处理获得
filename <- list.files(path = "data/raw/GSE11969/") # 用于后续读入原始样本文件
filename
filename_filter <- filename[str_sub(filename, 1, 9) %in% group_geo$Sample] # 挑选需要的样本
file.remove(paste0("data/raw/GSE11969/", filename[!(str_sub(filename, 1, 9) %in% group_geo$Sample)])) # 删除不需要的样本
group_geo <- group_geo[match(str_sub(filename, 1, 9), group_geo$Sample),] 
Targets <- data.frame(FileName = filename_filter, # 此文件可用于后续质量分析
                      sample_id = group_geo$Sample,
                      group = group_geo$Group,
                      row.names = group_geo$Sample, 
                      stringsAsFactors = F)
write.table(x = Targets,
            file = "data/raw/Targets.txt",
            quote = F,
            sep = "\t",
            row.names = F)

## 3.读入原始文件----
library(limma)
Targets <- limma::readTargets(path = "data/raw/", row.names = "sample_id") # 读入Targets文件
library(affy)
GSE11969_raw <- read.maimages(files = Targets$FileName,
                              source = "agilent", # source代表是经过哪种程序得到的，有些Agilent芯片是通过genepix处理的
                              path = "data/raw/GSE11969/",
                              names = Targets$sample_id,
                              other.columns = "gIsWellAboveBG", # 读取进去是为了判断是否高于背景值
                              green.only = T) # 代表是个单色芯片（一般不用改），默认是false
### GSE11969_raw的E代表expression,EB代表background
boxplot(log2(GSE11969_raw$E + 1))
GSE11969_raw$targets <- Targets

## 4.标准化----
### 背景矫正
GSE11969_bgc <- backgroundCorrect(RG = GSE11969_raw, 
                                  method = "normexp", # 推荐normexp进行背景矫正
                                  offset = 50, # 补偿值50
                                  normexp.method = "mle") # limma包推荐的mle
GSE11969_norm <- normalizeBetweenArrays(GSE11969_bgc, # normalizeBetweenArrays适用于单色芯片
                                        method = "quantile") # normalizeinArrays适用于双色芯片
GSE11969_expr_norm <- GSE11969_norm$E
sum(is.na(GSE11969_expr_norm))
boxplot(GSE11969_expr_norm)
GSE11969_exp <- cbind(ProbeName = GSE11969_norm[["genes"]][["ProbeName"]], GSE11969_expr_norm)
data.table::fwrite(GSE11969_exp, "output/GSE11969_norm.csv", row.names = F) # 写出，还需探针注释

### 5.质检----
library(arrayQualityMetrics)
GSE11969_raw$targets <- Targets # 加入target信息
GSE11969_targets <- GSE11969_raw$targets
# arrayQualityMetrics不接受ElistRaw，所以要构建ExpressionSet
all(rownames(GSE11969_targets) == colnames(GSE11969_raw$E)) # 确认targets样本顺序与矩阵样本顺序一致
rownames(GSE11969_targets) <- colnames(GSE11969_raw$E)
eset <- ExpressionSet(assayData = GSE11969_raw$E, 
                      phenoData = AnnotatedDataFrame(data = GSE11969_targets))
### 质量检测（原始）
dir.create("data/raw/GSE11969_QC_raw")
arrayQualityMetrics(eset, 
                    outdir = "data/raw/GSE11969_QC_raw", 
                    do.logtransform = T, # 对intensity进行对数转换,NAs are not allowed in subscripted assignments这个报错我改成了F
                    force = TRUE, # 强制覆盖重名文件
                    intgroup = "group") # 需要用到分组
### 看index文件中的Array metadata and outlier detection overview，"x"多的样本，在样本量大的情况下可以考虑剔除

### 质量检测（处理后）
eset2 <- ExpressionSet(assayData = GSE11969_norm$E, 
                       phenoData = AnnotatedDataFrame(data = GSE11969_targets))
dir.create("data/raw/GSE11969_QC_norm")
arrayQualityMetrics(eset2, 
                    outdir = "data/raw/GSE11969_QC_norm", 
                    do.logtransform = T,
                    force = TRUE,
                    intgroup = "group")
