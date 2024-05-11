load("data/GSE138260_clean_data.Rdata")

GSE138260 <- exprSet
GSE138260_group <- c(rep("Ctrl",19),rep("AD",17))

load("data/GSE5281_clean_data.Rdata")
GSE5281 <- exprSet
GSE5281_group <- c(rep("Ctrl",74),rep("AD",87))

load("data/GSE29378_Clean_data.Rdata")
GSE29378 <- exprSet
GSE29378_group <- pdata_match_GSE29378

load("data/GSE48350_Clean_data.Rdata")
GSE48350 <-exprSet
GSE48350_group <- c(rep("Ctrl",173),rep("AD",80))

save(GSE28146,GSE28146_group,GSE5281,GSE5281_group,
     GSE29378,GSE29378_group,GSE48350,GSE48350_group,file ="data/clean_data.Rdata")
