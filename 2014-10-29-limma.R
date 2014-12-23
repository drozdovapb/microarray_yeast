# for linux 
setwd() #put your path here

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

library(limma)

alldata <-read.table("./normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus") #don't forget to change the path if u need to
alldata.matrix <- as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, del-minus, plus-del, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

#everything above taken from limma user's guide

#extracting variables and writing tables

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="p")
del_minus<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="p")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="p")

#writing to files
#write.table(plus_minus, "../processed_tables_R/plus-minus.csv", sep="\t", dec=".", col.names=NA)
#write.table(del_minus, "../processed_tables_R/del-minus.csv", sep="\t", dec=".", col.names=NA)
#write.table(plus_del, "../processed_tables_R/minus-del.csv", sep="\t", dec=".", col.names=NA)

logdata<-cbind(plus_minus$logFC, del_minus$logFC, plus_del$logFC)
colnames(logdata)<-c("plus_min", "del_min", "plus_del")
row.names(logdata)<-row.names(plus_minus)
#write.table(logdata, "../processed_tables_R/logdata.csv", sep="\t", dec=".", col.names=NA, quote=FALSE)