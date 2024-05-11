library(limma)
library(dplyr)
library(magrittr)
library(tibble)
library(reshape2)
library(ggplot2)
library(GEOquery)
library(ggrepel)
library(stringr)
options(stringsAsFactors = F)
options('download.file.method.GEOquery'='libcurl')
rt <- getGEO('GSE48350',getGPL = F,destdir = '.',AnnotGPL = F,GSEMatrix = T)
exprSet <- exprs(rt[[1]])
pdata <- pData(rt[[1]])
exprSet <- as.data.frame(exprSet)

library(hgu133plus2.db)
probe2symbol <- toTable(hgu133plus2SYMBOL)
new<- as.data.frame(rownames(exprSet))
colnames(new) <- c('probe_id')
exprSet <- data.frame(new,exprSet)

###
exprSet <- dplyr::inner_join(exprSet,probe2symbol,by='probe_id')
exprSet <- exprSet[,!names(exprSet) %in% c('probe_id')]
exprSet <- dplyr::select(exprSet,symbol,everything())
###

index <- order(rowMeans(exprSet[,-1]),decreasing = T)
exprSet_ordered=exprSet[index,]
keep=!duplicated(exprSet_ordered$symbol)
exprSet=exprSet_ordered[keep,]
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]


group <- sapply(strsplit(pdata$title,'_'),'[',4)

group_list <- factor(c(rep("ND",173),rep("AD",80)))
group_list=factor(group_list,c("ND","AD"))

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

group_list <- factor(c(rep("ND",173),rep("AD",80)))
group_list=factor(group_list,c("ND","AD"))
group_GSE48350 <- group
pdata_GSE48350 <- pdata
save(group_GSE48350,pdata_GSE48350,exprSet,file = 'GSE48350_Clean_data.Rdata')