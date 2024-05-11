
load("data/clean_data.Rdata")
GSE138260 <- as.data.frame(GSE138260)
GSE29378 <- as.data.frame(GSE29378)
GSE48350 <- as.data.frame(GSE48350)
GSE5281 <- as.data.frame(GSE5281)
#
GSE138260_clock <- GSE138260[c("CLOCK","PER1","PER2","PER3","CRY1","CRY2"),]
GSE138260_clock <- as.data.frame(t(GSE138260_clock))
GSE138260_clock$group <- c(rep("Ctrl",19),rep("AD",17))

#
GSE29378_clock <- GSE29378[c("CLOCK","PER1","PER2","PER3","CRY1","CRY2"),]
GSE29378_clock <- as.data.frame(t(GSE29378_clock))
GSE29378_clock$group <- c(rep("Ctrl",32),rep("AD",31))

#
GSE48350_clock <- GSE48350[c("CLOCK","PER1","PER2","PER3","CRY1","CRY2"),]
GSE48350_clock <- as.data.frame(t(GSE48350_clock))
GSE48350_clock$group <- GSE48350_group

##第三个是GSE5281
GSE5281_clock <- GSE5281[c("CLOCK","PER1","PER2","PER3","CRY1","CRY2"),]
GSE5281_clock <- as.data.frame(t(GSE5281_clock))
GSE5281_clock$group <- GSE5281_group

clock01 <- as.data.frame(t(GSE138260["CLOCK",]))
clock01$group <- GSE138260_group
#
clock02 <- as.data.frame(t(GSE29378["CLOCK",]))
clock02$GSM <- rownames(clock02)
clock02 <- dplyr::inner_join(clock02,GSE29378_group,"GSM")
rownames(clock02) <- clock02$GSM
clock02 <- clock02[,-c(2)]
clock02$group <- c(rep("Ctrl",32),rep("AD",31))
#
clock03 <- as.data.frame(t(GSE48350["CLOCK",]))
clock03$group <- GSE48350_group
#
clock04 <- as.data.frame(t(GSE5281["CLOCK",]))
clock04$group <- GSE5281_group

clock01$supp <- c(rep("1",36))
clock02$supp <- c(rep("2",63))
clock03$supp <- c(rep("3",253))
clock04$supp <- c(rep("4",161))
clock <- rbind(clock01,clock02,clock03,clock04)
clock$group <- factor(clock$group,levels=c("Ctrl","AD"))
library(ggforce)
library(ggpubr)
clock_map <- ggplot(clock,aes(x=supp,y=CLOCK))+geom_boxplot(
  aes(fill = group),
  position = position_dodge(0.9))+ 
  scale_fill_manual(values = c("#63B8FF", "#FA8072"))+
  theme_classic()
library(export)
graph2eps(clock_map,'clock expression',font='TT Times New Roman',aspectr=2,height=3,width=4)  

