rm(list = ls())
library(tidyverse)
library(dplyr)


GSE13213_Group <- read.csv('GSE13213_Group.csv')
GSE13213_Matrix <- read.csv('GSE13213_Matrix.csv')
GSE13213_Matrix <- GSE13213_Matrix %>% column_to_rownames('X')


GSE30219_Group <- read.csv('GSE30219_Group.csv')
GSE30219_Matrix <- read.csv('GSE30219_Matrix.csv')
GSE30219_Matrix <- GSE30219_Matrix %>% column_to_rownames('X')

GSE11969_Group <- read.csv('GSE11969_Group.csv')
GSE11969_Matrix <- read.csv('GSE11969_Matrix.csv')
GSE11969_Matrix <- GSE11969_Matrix %>% column_to_rownames('X')
write.csv(rownames(GSE11969_Matrix),"GSE11969_Gene.csv")

save(GSE13213_Group,GSE13213_Matrix,GSE30219_Group,
        GSE30219_Matrix,GSE11969_Group,GSE11969_Matrix,file = "GEOdata.RData")
