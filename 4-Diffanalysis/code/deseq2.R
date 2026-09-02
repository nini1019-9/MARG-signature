rm(list = ls())

## 1. 矩阵+分组 ----
library(magrittr)
ssgsea_score <- data.table::fread('data/0-ssgsea_tcga.csv', data.table = F)
group <- ssgsea_score %>% 
  dplyr::select(2, 12) %>% 
  dplyr::arrange(group) ## 排序样本，tumor在前
head(group)
tail(group)
unique(group$group)
group_list <- factor(group$group, levels = c("HighScore", "LowScore"))

load('data/TCGA_Data.RData')

counts <- dat_counts %>% 
  dplyr::select(group$sample_id)
all(group$sample_id == colnames(counts)) ## TURE则排序一致

## 差异分析----
library(DESeq2)
colData <- data.frame(row.names = colnames(counts), group_list = group_list)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)
res <-  results(dds2, contrast = c("group_list", "HighScore", "LowScore"))
res1 <- res %>% 
  data.frame() %>% 
  dplyr::arrange(padj) ## 以padj排序

length(which((abs(res1$log2FoldChange) > 2) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 1) & (res1$padj < 0.05)))
length(which((abs(res1$log2FoldChange) > 0) & (res1$padj < 0.05)))
# phen_gene <- data.table::fread("data/NMRGs.csv", data.table = F)
res2 <- res1 %>% 
  dplyr::filter(abs(log2FoldChange) > 2 & padj < 0.05)
# vn_genes <- phen_gene[,1] %>% intersect(res2 %>% rownames())
write.csv(res1, "output/degs_tcga.csv", row.names = T)
# write.csv(data.frame(genes = vn_genes), "output/vn_genes.csv", row.names = F)

## vn----
# library(ggvenn)
# input_vn <- list(DEGs = rownames(res2), NMRGs = phengenes$NMRGs)
# input_vn <- lapply(input_vn, function(x) x[!is.na(x)])
# mycol <- c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280")
# p <- ggvenn(input_vn, 
#             show_percentage = F,
#             digits = 1,
#             fill_color = mycol,
#             fill_alpha = 0.5,
#             stroke_color = "white",
#             stroke_size = 0.75,
#             set_name_color = "black",
#             set_name_size = 6 * 0.35,
#             text_color = "black",
#             text_size = 6 * 0.35) +
#   theme(plot.margin = margin(-0.5, -0.5, -0.5, -0.5, "cm"))
# p
# ggsave("output/vn.pdf", p, width = 10, height = 10, units = "cm")

## volcano----
library(ggplot2)
mytheme <-
  theme(
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # 图像
    plot.title = element_text(size = 7,hjust = 0.5, vjust = 0.5)
  ) +
  theme(
    panel.background = element_blank(), # 面板
    panel.grid = element_blank(),
    panel.borde = element_rect(fill = NA, linewidth = 0.75 * 0.47)
  ) + # 添加外框
  theme(
    axis.line = element_line(size = 0.75 * 0.47), # 坐标轴
    axis.text = element_text(size = 6, color = "black"), # 坐标文字为黑色
    axis.title = element_text(size = 6),
    axis.ticks = element_line(size = 0.75 * 0.47)
  ) +
  theme(
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(c(0.3, 0.3), "cm"), # 图注
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.margin = margin(),
    legend.box.margin = margin(),
    legend.box.spacing = unit(0, "cm"),
    legend.background = element_blank(), legend.spacing = unit(0, "cm"),
    legend.box.background = element_blank()
  ) +
  theme(
    panel.grid = element_line(
      colour = "grey90",
      size = 0.75 * 0.47, linetype = 1
    ),
    plot.title = element_text(
      size = 7, hjust = 0, vjust = 0.5,
      margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    ),
    legend.position = "top", # 图注位置（上）
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    legend.title = element_blank(),
    legend.key.size = unit(c(0.15, 0.15), "cm")
  )

logfc <- 2
colnames(res1)
input_volcano <- res1 %>% 
  dplyr::rename(logFC = log2FoldChange) %>% 
  dplyr::mutate(threshold = ifelse(logFC >= logfc & padj < 0.05,"Up", 
                                   ifelse(logFC <= -logfc & padj < 0.05, "Down", "Not sig")))
write.csv(input_volcano[,c(2,6)],'output/1-input_volcano.csv')
library('ggthemes')

p <- ggplot(input_volcano, aes(x = logFC, y = -log10(padj))) +
  geom_point(shape = 21, aes(colour = factor(threshold, levels = c("Up", "Not sig", "Down")), 
                             fill = after_scale(alpha(colour, 0.5))), size = 1) +
  scale_colour_manual(values = c("Up"= "#A85774", "Down"="#CBD8EC", "Not sig"= "#7f7f7f")) +
  geom_hline(yintercept = -log10(0.05), size = 0.75 * 0.47, colour = "#999999", linetype = "longdash") +
  geom_vline(xintercept = c(logfc, -logfc), size = 0.75 * 0.47, colour = "#999999", linetype = "longdash") +
  xlab("Log2(Fold Change)") + ylab("-log10(P.adj)") +
  mytheme +
  theme(legend.title = element_blank(),
        legend.position = "top")
p
ggsave("output/volcano_ssgsea_tcga.pdf", p, width = 5, height = 5, units = "cm")
ggsave("output/volcano_sgsea_tcga.tiff", p, width = 5, height = 5, units = "cm")

## heatmap---
library(pheatmap)

res3 <- res1[order(res1$log2FoldChange),]
res3 <- na.omit(res3)
up_gene <- res3 %>% 
  dplyr::filter(padj < 0.05,log2FoldChange > 1) %>% 
  dplyr::top_n(-20, padj)
down_gene <- res3 %>% 
  dplyr::filter(padj < 0.05,log2FoldChange < -1) %>% 
  dplyr::top_n(-20, padj)
group2 <- group %>% 
  dplyr::arrange(desc(group))
input_heatmap <- dat_fpkm[c(rownames(up_gene), rownames(down_gene)),] %>% 
  dplyr::select(group2$sample_id)
# dplyr::filter(rownames(.) %in% vn_genes)
input_heatmap2 <- rbind(group2$group,colnames(input_heatmap),input_heatmap)
rownames(input_heatmap2)[1] <- '#Group'
rownames(input_heatmap2)[2] <- 'Gene'
write.csv(input_heatmap2,'output/2-input_heatmap.csv')

annotation_col <- group2 %>% 
  tibble::column_to_rownames("RNAseq")
ann_colors <- list(group = c(high= "#EE0000", low= "#3B4992"))
pdf("output/heatmap.pdf", width = 5, height = 5)
pheatmap(input_heatmap, scale = "row",
         cluster_cols= F,
         show_colnames = F,
         show_rownames = T,
         border_color = NA,
         fontsize = 6,
         annotation_col = annotation_col,
         color=colorRampPalette(c("#3B4992", "#ffffff","#EE0000"))(100),
         clustering_distance_cols = "manhattan",
         annotation_names_row = F,
         annotation_colors = ann_colors)
dev.off()

## pca----
library(scatterplot3d)
input_pca <- dat_fpkm[,group2$sample_id]
pca <- princomp(input_pca)
color <- factor(group2$group,labels = c("#A6B5C8","#E4A6BD"), levels = c("LowScore","HighScore"))
pdf('output/pca.pdf',onefile=TRUE,width=5.3,height=5)
diffangle <- function(ang){
  scatterplot3d(pca$loadings[,1:3],main='PCA',color=color,type='p',
                highlight.3d=F,angle=ang,grid=T,box=T,scale.y=1,
                cex.symbols=1.2,pch=16,col.grid='lightblue')
  legend("topright",c("LowScore", "HighScore"), fill=c("#A6B5C8","#E4A6BD"),box.col="grey")
}
sapply(seq(-360,360,5),diffangle)
dev.off()


table(group2$group)
library(Rtsne)
cla <- c(rep('LowScore',432),rep('HighScore',82))
tSNE_res <- Rtsne(t(input_pca),
                  dims=3,
                  perplexity=5,
                  verbose=F,
                  max_iter=500,
                  check_duplicates=F)
tsne <- data.frame(tSNE1 = tSNE_res[["Y"]][,1],
                   tSNE2 = tSNE_res[["Y"]][,2],
                   tSNE3 = tSNE_res[["Y"]][,3],
                   cluster = cla)

tsne$colors <- ifelse(tsne$cluster %in% "LowScore","#A6B5C8","#E4A6BD")
colnames(tsne) <- c("x","y","z","cluster","colors")
tsne$cluster <- as.factor(tsne$cluster)
pdf("output/3-Cluster_PCA.pdf",width = 8.1,height = 6)
library(scatterplot3d)
p3 <- scatterplot3d(tsne[,1:3],color = tsne$colors,main="Consensus Cluster PCA",pch = 21,bg = tsne$colors)
legend("bottom",col = "black", legend = levels(tsne$cluster),pt.bg = c("#A6B5C8","#E4A6BD"), pch = 21,
       inset = -0.2, xpd = TRUE, horiz = TRUE)
dev.off()

# legend("topright",paste("G",1:7,sep=""),fill=c(lblue,purple,orange,yellow,red2,green,"black"),box.col="grey")
