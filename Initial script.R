# for linux 
setwd("/media/drozdovapb//441DB6652D7ED741/Work/Current work/Projects/bioinfo/raw_data_our/")

#for windows
#setwd("D:/Polina/Work/Ingenuity")

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

library(limma)

alldata <-read.table("./Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)


#everything above taken from limma user's guide



#extracting variables and writing tables

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="p")
minus_del<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="p")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="p")

#write.table(plus_minus, "./data/plus-minus.csv", sep="\t", dec=".", col.names=NA)
#write.table(minus_del, "./data/minus-del.csv", sep="\t", dec=".", col.names=NA)
#write.table(plus_del, "./data/minus-del.csv", sep="\t", dec=".", col.names=NA)
            

#subsetting

#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}

#applying function!!!



#creating a matrix
creatematrix <-function(comparison) {
  abslogFCs <- c(0, 0.5, 1, 2)
  p.values <-c(0.05, 0.01, 0.001)
mymatrix <-matrix(nro=length(abslogFCs), nco=length(p.values))
dimnames(mymatrix) <- list(abslogFCs=abslogFCs, p.values=p.values)
for (i in 1:length(abslogFCs)){
  for (j in 1:length(p.values)){
  mymatrix[i,j]<-nrow(cutoff(comparison, abslogFCs[i], p.values[j]))
  }
  }
print(mymatrix)

filename=paste("./data/",deparse(substitute(comparison)),"_","number of genes",".csv", sep="")
#mydata.frame<-as.data.frame(mymatrix)
write.table(mymatrix, file=filename, col.names=NA)
}

creatematrix(plus_minus)
creatematrix(minus_del)
creatematrix(plus_del)

printcutoffs <-function(table) {
  
  abslogFCs <- c(0,0.5, 1, 2)
  p.values <-c(0.05, 0.01, 0.001)
  
  for (l in 1:length(abslogFCs)){
    for (m in 1:length(p.values)){
      p.value.name <- deparse(p.values[m])
      abslogFC.name <-deparse(abslogFCs[l])
      write.table(cutoff(table, abslogFCs[l], p.values[m]), 
        file=paste("./data/",deparse(substitute(table)), 
        "abslogFC",abslogFC.name,"adj.p.val", p.value.name,".csv", sep="_"), col.names=NA)
    }
  }
}

printcutoffs(plus_minus)
printcutoffs(plus_del)
printcutoffs(minus_del)

plus_minus_cut <- cutoff(plus_minus, 0.5, 0.001)
plus_del_cut <- cutoff(plus_del, 0.5, 0.001)
minus_del_cut <- cutoff(minus_del, 0.5, 0.001)

write.table(plus_minus_cut, "plus_minus_cut.csv", col.names=NA)
write.table(plus_del_cut, "plus_del_cut.csv", col.names=NA)
write.table(minus_del_cut, "minus_del_cut.csv", col.names=NA)

# splitting into upregulated and downregulated
plus_minus_cut_up <- subset(plus_minus_cut, plus_minus_cut$logFC > 0)
write.table(plus_minus_cut_up, "plus_minus_cut_up.csv", col.names=NA, sep="\t")
plus_minus_cut_down <- subset(plus_minus_cut, plus_minus_cut$logFC < 0)
write.table(plus_minus_cut_down, "plus_minus_cut_down.csv", col.names=NA, sep="\t")
plus_del_cut_up <-subset(plus_del_cut, logFC > 0) 
write.table(plus_del_cut_up, "plus_del_cut_up.csv", col.names=NA, sep="\t")
plus_del_cut_down <- subset(plus_del_cut, logFC < 0)
write.table(plus_del_cut_down, "plus_del_cut_down.csv", col.names=NA, sep="\t")
minus_del_cut_up <- subset(minus_del_cut, logFC > 0)
write.table(minus_del_cut_up, "minus_del_cut_up.csv", col.names=NA, sep="\t")
minus_del_cut_down <- subset(minus_del_cut, logFC <0)
write.table(minus_del_cut_down, "minus_del_cut_down.csv", col.names=NA, sep="\t")
                         

#sorting by FC, by descending abs(!) logFC
plus_minus_cut_up <- plus_minus_cut_up[order(-plus_minus_cut_up$logFC),]
plus_minus_cut_down <- plus_minus_cut_down[order(plus_minus_cut_down$logFC),]
plus_del_cut_up <-plus_del_cut_up[order(-plus_del_cut_up$logFC),] 
plus_del_cut_down <- plus_del_cut_down[order (plus_del_cut_down$logFC),]
minus_del_cut_up <- minus_del_cut_up[order(-minus_del_cut_up$logFC),]
minus_del_cut_down <- minus_del_cut_down [order(minus_del_cut_down$logFC),]

#writing tables

#install.packages("WriteXLS")
#library(WriteXLS)
#I should finish this later

#venn diagrams
#install.packages("qpcR")
library(qpcR)
columnnames <- c("plus_minus", "plus_del", "minus_del")
plus_minus_cut_NA <- c(plus_minus_cut$ID, rep.int(NA,
        times=max(nrow(plus_minus_cut),nrow(plus_del_cut),nrow(minus_del_cut))
                                                  -nrow(plus_minus_cut)))
plus_del_cut_NA <- c(plus_del_cut$ID, rep.int(NA,
        times=max(nrow(plus_minus_cut),nrow(plus_del_cut),nrow(minus_del_cut))
                                                  -nrow(plus_del_cut)))
minus_del_cut_NA <- c(minus_del_cut$ID, rep.int(NA,
        times=max(nrow(plus_minus_cut),nrow(plus_del_cut),nrow(minus_del_cut))
                                                  -nrow(minus_del_cut)))

