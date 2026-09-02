rm(list = ls())
gc()
library(Seurat)
library(magrittr)
library(ggplot2)
addr_files <- dir("./scRNA/", full.names = T)
list_sc <- data.table::fread('scRNA/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz')

annotation <- read.table('GSE131907_Lung_Cancer_cell_annotation.txt',sep = '\t',header = T)
annotation <- annotation[annotation$Sample %in% c("LUNG_N01","LUNG_N06","LUNG_N08","LUNG_N09","LUNG_N18","LUNG_N19","LUNG_N20",
                                                  "LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T34"),]

sce <- SingleCellExperiment(assays = list(counts = list_sc))

sc_umi = sce[,annotation$Index]
sc_umi = assay(sc_umi)
rownames(sc_umi) <- list_sc$Index
sc_merge <- CreateSeuratObject(count = sc_umi, min.cells = 3, min.feature = 200)

# # 添加细胞注释信息
sc_merge[["sample"]] <- substr(colnames(sc_merge), start=18, stop=27)
sc_merge[["tissue"]] <- ifelse(sc_merge$sample %in% c("LUNG_N01","LUNG_N06","LUNG_N08","LUNG_N09","LUNG_N18","LUNG_N19","LUNG_N20"),"Normal","LUAD")


# 计算sc_merge的线粒体基因表达并进行筛选
meta.data <- sc_merge@meta.data
# cellid <- stringr::str_extract(rownames(meta.data), ".*(?=-)")
cellid <- rownames(meta.data)
sc_merge$cellid <- !duplicated(rownames(meta.data))
sc_merge <- subset(sc_merge,subset = cellid)
sc_merge$log10GenesPerUMI <- log10(sc_merge$nFeature_RNA) / log10(sc_merge$nCount_RNA)

#人用MT，鼠是Mt
sc_merge[["percent.mt"]] = PercentageFeatureSet(sc_merge, pattern = "^MT-" )/100 

sc_merge <- subset(x = sc_merge, subset= (nCount_RNA >= 500) & 
                     (nFeature_RNA >= 250) & 
                     (log10GenesPerUMI > 0.80) & 
                     (percent.mt < 0.20))
counts <- GetAssayData(object = sc_merge, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
sc_filter <- CreateSeuratObject(filtered_counts, meta.data = sc_merge@meta.data)

vlnplot <- VlnPlot(sc_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   pt.size = 0, cols = '#eca680',ncol = 3)
vlnplot
ggsave(filename = "1-ViolinPlot.pdf", plot = vlnplot, width = 17, height = 12, units = "cm")

# PCA和细胞聚类分析
# Normalization and down to 2 dimenson
sc_filter <- NormalizeData(sc_filter) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
# Removal of batch effects
sc_integer <- sc_filter %>% 
  harmony::RunHarmony(group.by.vars = "sample")

normalization <- DimPlot(sc_filter, reduction = "pca", pt.size = 1, group.by = "sample")
normalization
ggsave(filename = "2-Normalization.pdf", plot = normalization, width = 16.2, height = 10.3, units = "cm")
DimPlot(sc_integer, reduction = "harmony", pt.size = 1, group.by = "sample")

# 细胞聚类分析
# Find clusters: default resolution = 0.8 --> resolution参数越大，得到的cluster数量越多
sc_integer <- sc_integer %>% 
  FindNeighbors(dims = 1:30, verbose = T, reduction = "harmony") %>%
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = 1:30, reduction = "harmony")

umap1 <- DimPlot(sc_integer, reduction = "umap", pt.size = 1, label=T) 
umap1
ggsave(filename = "3-UMAP.pdf",plot = umap1,width = 16.2, height =9.4, units = "cm")
# Assign group information
# sc_integer@meta.data$group <- ifelse(sc_integer@meta.data$sample %in% c("GSM4569780","GSM4569781","GSM4569782"),"ARDS","Control")
# Redraw the umap to make sure the batch effects is removed
# batcheffect <- DimPlot(sc_integer, reduction = "umap", label=T, group.by = "tissue") 
# batcheffect
# ggsave(filename = "4-BatchEffect.pdf", plot = batcheffect, width = 16.6, height = 12, units = "cm")

# SingleR: automatically cell type annotation
library(SingleR)
# BiocManager::install("celldex")
# BiocManager::install("SingleR")
hpca <- celldex::HumanPrimaryCellAtlasData()
# saveRDS(hpca,file = "hpca.rds")
# mrsd <- celldex::MouseRNAseqData()
igd <- celldex::ImmGenData()

testdata <- GetAssayData(sc_integer, slot="data")
# library(homologene)
# genelist <- homologene(testdata@Dimnames[[1]],inTax = 10116,outTax = 10090)
# genelist <- genelist[!duplicated(genelist$`10090`),]
# testdata <- testdata[genelist$`10116`,]
# dim(testdata)
# rownames(testdata) <- genelist$`10090`

clusters <- sc_integer@meta.data$seurat_clusters

cellpred <- SingleR::SingleR(test = testdata,
                             ref = hpca,
                             labels = hpca$label.main,
                             de.method = "classic",
                             clusters = clusters)
celltype = data.frame(ClusterID=rownames(cellpred),
                      cell=cellpred$labels,
                      stringsAsFactors = F)
# 人工注释：陈威代码
celltype2 <- celltype
celltype2$cell[which(celltype2$ClusterID == 9)] <- "cluster9"
celltype2$cell[which(celltype2$ClusterID == 13)] <- "cluster13"
celltype2$cell[which(celltype2$ClusterID == 12)] <- "cluster12"
sc_integer$Manual_anno <- factor(sc_integer$seurat_clusters, labels = as.character(celltype2$cell))

Idents(sc_integer) <- "Manual_anno"

markers_cluster9 <- FindMarkers(sc_integer,ident.1 = "cluster9")
markers_cluster9_1 <- markers_cluster9 %>%
  dplyr::filter(p_val_adj<0.05) %>%
  mutate(pct=pct.1-pct.2) %>%
  dplyr::arrange(desc(avg_log2FC),desc(pct))
head(markers_cluster9_1,5)
DotPlot(sc_integer, assay = "RNA",features = rownames(markers_cluster9_1)[1:10],group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
celltype2$cell[which(celltype2$ClusterID == 9)] <- "Mast_cell"

markers_cluster13 <- FindMarkers(sc_integer,ident.1 = "cluster13")
markers_cluster13_1 <- markers_cluster13 %>%
  dplyr::filter(p_val_adj<0.05) %>%
  mutate(pct=pct.1-pct.2) %>%
  dplyr::arrange(desc(avg_log2FC),desc(pct))
head(markers_cluster13_1,10)
DotPlot(sc_integer, assay = "RNA",features = rownames(markers_cluster13_1)[1:20],group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
celltype2$cell[which(celltype2$ClusterID == 13)] <- "Endothelial_cells"

markers_cluster12 <- FindMarkers(sc_integer,ident.1 = "cluster12")
markers_cluster12_1 <- markers_cluster12 %>%
  dplyr::filter(p_val_adj<0.05) %>%
  mutate(pct=pct.1-pct.2) %>%
  dplyr::arrange(desc(avg_log2FC),desc(pct))
head(markers_cluster12_1,10)
DotPlot(sc_integer, assay = "RNA",features = rownames(markers_cluster12_1)[1:20],group.by = "seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
celltype2$cell[which(celltype2$ClusterID == 12)] <- "Endothelial_cells"
# # 人工注释：mouse_marker.csv
# Significant bio-markers which can be used for cluster annotation and gene expression analysis
#"IGKC"-"B-cell"
#"CD3D"-"T-cell"
#"PLAUR"-"Monocyte"
#"GZMB"-"NK_cell"
#"MARCO"-"Macrophage"
#"KRT19"-"Epithelial_cells"
#"'PECAM1'-"Endothelial_cells"

# marker = read.csv('mouse_marker.csv',check.names = F)
# celltype = c('Mast','Macrophages','Epithelial','Fibroblasts','Endothelial','B',
#              'Natural killer T','Smooth muscle','T regulatory','Dendritic','T','Neutrophils')
# for (i in celltype){
#   markerlist =  paste0(toupper(substr(marker[,i], 1, 1)), tolower(substr(marker[,i], 2, nchar(marker[,i]))))
#   markerlist = markerlist[markerlist!='']
#   markerlist = unique(markerlist)
#   var = paste0(i,'marker')
#   assign(var,markerlist)
# }

# #注释T细胞
# DotPlot(sc_integer, features = Tmarker )#3 Ptprcap Cd3g Lat Cd3d
# #注释B细胞
# DotPlot(sc_integer, features = Bmarker )#8 Ms4a1 CD74
# #注释Mast细胞
# DotPlot(sc_integer, features = Mastmarker )#17 cpa3 cma1 Cyp11a1 Hdc
# #注释Neu细胞
# DotPlot(sc_integer, features = Neutrophilsmarker )#4,11 entpd1 sell
# #注释平滑肌细胞
# DotPlot(sc_integer, features = `Smooth musclemarker` )#13 ACTA2 Myl9 Myh11 hspb6
# #注释Macro细胞
# DotPlot(sc_integer, features = Macrophagesmarker)#0 cd68 cd40 mrc1 
# #注释成纤维细胞
# DotPlot(sc_integer, features = Fibroblastsmarker)#2,6, 12 ,14 col3a1 col5a2 fn1 gsn lrp1 
# #注释表皮细胞
# DotPlot(sc_integer, features = Epithelialmarker)#4 5 7 9 10 11  Epcam Cdh1 krt5 krt15
# #注释内皮细胞
# DotPlot(sc_integer, features = Endothelialmarker)#16 pecam1 egfl7 flt1 emcn
# 
# clustercelltype = c(
#   "0" = "Macrophages", 
#   "1" = "Neutrophils", 
#   "2" = "Fibroblasts",
#   "3" = "T cell", 
#   "4" = "Epithelial cell", 
#   "5" = "Epithelial cell", 
#   "6" = "Fibroblasts",
#   "7" = "Epithelial cell", 
#   "8" = "B cell", 
#   "9" = "Epithelial cell",
#   "10" = "Epithelial cell", 
#   "11" = "Epithelial cell", 
#   "12" = "Fibroblasts", 
#   "13" = "Smooth muscle cell",
#   "14" = "Fibroblasts",
#   "15" = "Smooth muscle cell",
#   "16" = "Endothelial cell",
#   "17" = "Mast cell")
# sc_integer[['cell_type']] = unname(clustercelltype[sc_integer@meta.data$seurat_clusters])
# 
# sc_integer = RenameIdents(sc_integer,
#                             `0` = "Macrophages", 
#                             `1` = "Neutrophils", 
#                             `2` = "Fibroblasts",
#                             `3` = "T cell", 
#                             `4` = "Epithelial cell", 
#                             `5` = "Epithelial cell", 
#                             `6` = "Fibroblasts",
#                             `7` = "Epithelial cell", 
#                             `8` = "B cell", 
#                             `9` = "Epithelial cell",
#                             `10` = "Epithelial cell", 
#                             `11` = "Epithelial cell", 
#                             `12` = "Fibroblasts", 
#                             `13` = "Smooth muscle cell",
#                             `14` = "Fibroblasts",
#                             `15` = "Smooth muscle cell",
#                             `16` = "Endothelial cell",
#                             `17` = "Mast cell")

# saveRDS(sc_integer, file = "sc_integer.rds")
Idents(sc_integer) <- "seurat_clusters"
sc_integer <- subset(sc_integer, idents =  celltype2$ClusterID)
sc_integer$cell_anno <- factor(sc_integer$seurat_clusters, labels = celltype2$cell)

# Annotated UMAP
Idents(sc_integer) <- "cell_anno"
umap2 <- DimPlot(sc_integer, group.by="cell_anno", label=T,pt.size = 1, reduction='umap',cols = c('#67ADB7','#F49869','#F5E1D8','#E4A6BD','#F3D8E1','#AFACB7','#A6B5C8','#C17688',"#D3E2B7"))+ labs(title = NULL)
umap2
ggsave(filename = "5-UMAP_Anno.pdf", plot = umap2, width = 16.2, height =9.4, units = "cm")

# Calculate cell ratio by sample info
Cellratio_sample <- prop.table(table(Idents(sc_integer), sc_integer$sample), margin = 2)#计算各组样本不同细胞群比例: 1 = row, 2 = column, default is NULL
Cellratio_sample

Cellratio_sample <- as.data.frame(Cellratio_sample)
colnames(Cellratio_sample)[1] <- "Cell"

library("ggsci")
ratioplot_sample <- ggplot(Cellratio_sample) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cell),
           stat = "identity", width = 0.7, size = 0.4)+ 
  theme_classic() + scale_fill_manual(values = c('#67ADB7','#F49869','#F5E1D8','#E4A6BD','#F3D8E1','#AFACB7','#A6B5C8','#C17688',"#D3E2B7")) +
  labs(x='Sample', y = 'Ratio')+
  # scale_fill_futurama()+
  coord_flip()
#theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
#legend.position = "bottom")
ratioplot_sample
ggsave("6-Cellratio_Sample.pdf", plot = ratioplot_sample, scale = 1, width = 16.2, height = 9, units =c("cm"))

# Calculate cell ratio by group info
Cellratio_group <- prop.table(table(Idents(sc_integer), sc_integer$tissue), margin = 2)#计算各组样本不同细胞群比例: 1 = row, 2 = column, default is NULL
Cellratio_group

Cellratio_group <- as.data.frame(Cellratio_group)
colnames(Cellratio_group)[1] <- "Cell"

ratioplot_group <- ggplot(Cellratio_group) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Cell),
           stat = "identity", width = 0.7, size = 0.4, colour = NA)+ 
  theme_classic() + scale_fill_manual(values = c('#67ADB7','#F49869','#F5E1D8','#E4A6BD','#F3D8E1','#AFACB7','#A6B5C8','#C17688',"#D3E2B7")) +
  labs(x='Group', y = 'Ratio')+
  # scale_fill_futurama()+
  coord_flip()
#theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"),
#legend.position = "bottom")
ratioplot_group
ggsave("7-Cellratio_Group.pdf", plot = ratioplot_group, scale = 1, width = 16.6, height = 10, units =c("cm"))

# Dotplot and Featureplot with interested biomarker
marker <-  data.table::fread("hub_genes.csv",data.table = F)
marker <- marker$x
marker <- intersect(marker,rownames(sc_integer))

dotplot1 <- DotPlot(sc_integer,features = marker,cols = c("#0072B5", "#BC3C29"),dot.scale = 10)+
  ggpubr::rotate_x_text(90)
dotplot1
ggsave("8-Dotplot.pdf", plot = dotplot1, width = 8.1*3, height = 5*3, units = "cm")

featureplot <- FeaturePlot(sc_integer,features = marker,cols = c("#0072B5", "#BC3C29"),pt.size = 1,label.size = 3)+
  ggpubr::rotate_x_text(90)
featureplot
ggsave("10-Featureplot1.pdf", plot = featureplot, width = 30, height = 30, units = "cm")

# express_data <- sc_integer@assays[["RNA"]]@data
# sample <- sc_integer@meta.data[["sample"]]
# 
# express_data <- as.data.frame(express_data)
# sample <- as.data.frame(sample)
# GSE_mat <- as.data.frame(t(express_data[marker,]))
# gene<- rep(colnames(GSE_mat),each=nrow(GSE_mat))
# gene <- factor(gene,levels = colnames(GSE_mat))
# a<- sample$sample
# group <- rep(a,ncol(GSE_mat))


# sc_integer_oec <- subset(sc_integer, tissue == "OEC")
saveRDS(sc_integer,file = "sc_integer.rds")
