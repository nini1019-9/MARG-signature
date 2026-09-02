rm(list = ls())
library(magrittr)
library(ggpubr)
load('data/TCGA_Data.RData')

dat_tide_norm <- dat_fpkm %>% 
  apply(1, function(x)(x - mean(x))) %>% 
  t() %>% data.frame(check.names = F)
write.table(dat_tide_norm, file="output/TIDE_input.txt", sep = "\t", quote = F, row.names = T)

load("./data/rs_risk.RData")
group_tcga <- rs_risk[[1]] %>% 
  dplyr::mutate(id = rownames(.)) %>% 
  dplyr::select(id, Risk) %>% 
  `names<-`(c("id", "group"))
rs_tide <- data.table::fread("output/tide_output.csv", data.table = F)
input_tide <- rs_tide %>% 
  dplyr::select("Patient", "TIDE") %>% 
  dplyr::inner_join(group_tcga, by = c("Patient" = "id"))

## violin----
violin=ggviolin(input_tide, x="group", y="TIDE", fill = "group",
                xlab="Group", ylab="TIDE",
                legend.title="Risk Group",
                palette=c("#D3E2B7","#F18A69"),
                add="boxplot", add.params = list(fill="white"))+ 
  geom_signif(comparisons = list(c("high","low")),step_increase = 0,
              map_signif_level =T,test = wilcox.test,textsize = 6)
pdf(file=paste0("output/tide_vioplot.pdf"), width=8.1, height=5)
print(violin)
dev.off()
write.csv(input_tide[,c(3,2)],'tide_plot.csv',row.names = F)
## stacked_column----
library(plyr)
input_tide2 <- rs_tide %>% 
  dplyr::select("Patient", "Responder") %>% 
  dplyr::inner_join(group_tcga, by = c("Patient" = "id")) %>% 
  dplyr::select(2,3) %>% 
  table() %>% data.frame() %>% 
  ddply(.(group), transform, percent = Freq/sum(Freq) * 100) %>% 
  dplyr::mutate(label = paste0(sprintf("%.0f", percent), "%"))
input_tide2 <- input_tide2[order(input_tide2$group,decreasing = T),]
input_tide2$group <- factor(input_tide2$group,levels = c('LowRisk','HighRisk'))
bioCol=c("#D3E2B7","#eca680")
p=ggplot(input_tide2, aes(x = factor(group), y = percent, fill = Responder,stratum = Responder, alluvium = Responder)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  #scale_fill_manual(values=bioCol)+
  xlab("Risk Group")+ ylab("Percent Weight")+  guides(fill=guide_legend(title="Responder"))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  #coord_flip()+
  scale_fill_manual(values=bioCol) +
  theme_bw() + theme(panel.grid=element_blank(),text = element_text(size = 23))
pdf(file=paste0("output/1-tide_stacked_column.pdf"), width=5.3*2, height=5*2)
p
dev.off()
