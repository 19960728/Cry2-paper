library(limma)
library(dplyr)
library(magrittr)
library(tibble)
library(reshape2)
library(ggplot2)
library(GEOquery)
library(ggrepel)
library(stringr)
setwd("G:\\R\\GSE138260/data")
options(stringsAsFactors = F)
options('download.file.method.GEOquery'='libcurl')
rt <- getGEO('GSE138260',getGPL = F,destdir = '.',AnnotGPL = F,GSEMatrix = T)
exprSet <- exprs(rt[[1]])
pdata <- pData(rt[[1]])
exprSet <- as.data.frame(exprSet)

probe2symbol <- data.table::fread('GSE138260_family.soft.gz',skip = "ID")
probe2symbol2 <- probe2symbol[,c(1,7)]

probe2symbol3 <- probe2symbol2[probe2symbol2$GENE_SYMBOL!="",]

exprSet$ID <- data.frame(rownames(exprSet))

names(exprSet[,37]) <- c("ID")

library(dplyr)
library(tibble)

exprSet <- dplyr::inner_join(exprSet,probe2symbol3,by='ID')

exprSet <- dplyr::select(exprSet,GENE_SYMBOL,everything())

colnames(exprSet)

exprSet <- exprSet[,!names(exprSet) %in% c('ID')]

exprSet <- exprSet %>% mutate(rowMeans=rowMeans(exprSet[,-1])) %>% 
  arrange(desc(rowMeans)) %>% 
  distinct(GENE_SYMBOL,.keep_all = T) %>% 
  select(-rowMeans) %>% 
  column_to_rownames("GENE_SYMBOL")

group_list <- factor(c(rep("ND",19),rep("AD",17)))
group_list=factor(group_list,c("AD","ND"))

boxplot(exprSet,outline=FALSE,notch=T,col=group_list,las=2)
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE,notch=T,col=group_list,las=2)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC){ ex[which(ex<=0)] <- NaN
exprSet <- log2(ex)
print('log2 transform finished')}else{print('log2 transform not needed')}

save(exprSet,pdata,file = 'GSE138260_clean_data.Rdata')


