rm(list = ls())
#拟时序分析
library(monocle)
library(Seurat)
# 构建CDS数据
## counts
sc_filter <- readRDS('data/sc_integer.rds')
hubgenes <- read.csv('data/hub_genes.csv')

# DefaultAssay(sc_filter) <- "SCT"
data <- as(as.matrix(sc_filter@assays$RNA@data), "sparseMatrix")
pd <- new('AnnotatedDataFrame', data = sc_filter@meta.data) # 样本信息
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) # 基因名称
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 1) 
head(fData(monocle_cds))

#02.选择基因用于轨迹分析
disp_table = dispersionTable(monocle_cds)
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.1)
monocle_cds = setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
diff_test_res = differentialGeneTest(monocle_cds,fullModelFormulaStr = "~cell_type")
ordering_genes = row.names (subset(diff_test_res, qval < 0.01))

monocle_cds = setOrderingFilter(monocle_cds, ordering_genes)
monocle_cds = reduceDimension(monocle_cds, 
                              max_components = 2,
                              reduction_method = 'DDRTree')
rm(list = lsf.str(envir = .GlobalEnv), envir = .GlobalEnv)
source("code/orderCells.R")
library(igraph)
monocle_cds = orderCells(monocle_cds)

pdf('output/01.pseudotime.pdf',width = 8,height = 6)
plot_cell_trajectory(monocle_cds,color_by = "cell_type")
dev.off()
pdf('output/02.pseudotime.pdf',width = 8,height = 6)
plot_cell_trajectory(monocle_cds,color_by = "Pseudotime")
dev.off()
pdf('output/03.pseudotime.pdf',width = 8,height = 6)
plot_cell_trajectory(monocle_cds,color_by = "State")
dev.off()
pdf('output/04.pseudotime.pdf',width = 8,height = 6)
plot_cell_trajectory(monocle_cds,color_by = "cell_type")+facet_wrap(~cell_type,nrow=2)
dev.off()

inter_gene <- intersect(rownames(monocle_cds),hub_genes)
pdf('output/pseudotimeheatmap.pdf',width = 10,height = 7.5)
plot_pseudotime_heatmap(monocle_cds[inter_gene,],
                        # num_clusters = 7,
                        cores = 1,
                        show_rownames = T)
dev.off()

cds_subset <- monocle_cds[inter_gene, ]

pdf("output/pseudotime_genes_State.pdf", width = 7, height = 14)
plot_genes_in_pseudotime(cds_subset, color_by = "State")
dev.off()

