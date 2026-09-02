require(data.table)

ADRGs <- read.csv('MARGs_raw.csv')

tcga <- read.csv('Matrix_FPKM.csv')

GSE13213 <- read.csv('GSE13213_Gene.csv')
GSE30219 <- read.csv('GSE30219_Gene.csv')
GSE11969 <- read.csv('GSE11969_Gene.csv')

ADRGs_filter <- intersect(ADRGs$MARGs,tcga$X)
ADRGs_filter <- intersect(ADRGs_filter,GSE13213$x)
ADRGs_filter <- intersect(ADRGs_filter,GSE30219$x)
ADRGs_filter <- intersect(ADRGs_filter,GSE11969$x)

write.csv(ADRGs_filter,'MARGS.csv',row.names = F)

