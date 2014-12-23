#functions

setwd("/media/57C506FA65AE1944/Work/Projects/bioinfo/proto_scripts/")

#packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("VennDiagram")
#biocLite("seqinr")
#biocLite("limma")
library(VennDiagram)
library (limma)
library(seqinr)

#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}

cutoff2 <- function(mytable, p.value=.001, abslogFC=.5, up_or_down="all") {
  if (up_or_down == "all") {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)}
  else if (up_or_down == "up") {
    subset(mytable, mytable$logFC>abslogFC & mytable$adj.P.Val<p.value)}
  else if (up_or_down == "down") {
    subset(mytable, mytable$logFC< (0-abslogFC) & mytable$adj.P.Val<p.value)}
}

#for one gene
annoone <- function(geneid) {
genestable <- read.csv("./SGD retrieves/genelist2.csv", sep=",")
    geneidspace=paste(geneid, "")
    anno <-subset (genestable, genestable$gene_id==geneidspace)
    return(anno)
}

#for a list
annolist <-function (genelist) {
  anno<-c()
  for (gene in (1:length(genelist))) {
    geneid=paste(genelist[gene], "")
    anno <-rbind(anno, subset(genestable, genestable$gene_id==geneid))}
  return (anno)
}
    

numbers <-function (geneid, compar) {
comparison <-c()
  if (compar=="pm") comparison<-plus_minus_cut
  if (compar=="pd") comparison <-plus_del_cut
  if (compar=="md") comparison<-minus_del_cut
  if (compar=="dm") comparison<-del_minus_cut
nums <- subset(comparison, comparison$ID==geneid)
return(nums)
}


numberslist <-function (genelist, compar) {
  comparison <-c()
  if (compar=="pm") comparison<-plus_minus_cut
  if (compar=="pd") comparison <-plus_del_cut
  if (compar=="md") comparison<-minus_del_cut
  if (compar=="dm") comparison<-del_minus_cut
  nums<-c()
  for (gene in (1:length(genelist))) {
    geneid<-genelist[gene]
  nums <- rbind(nums, subset(comparison, comparison$ID==geneid))
    }
  return(nums)
}

annogenelist <- function (genelist, compar) {
comparison<-c()
mylist<-c()
  if (compar=="pm") comparison<-plus_minus_cut
  if (compar=="pd") comparison <-plus_del_cut
  if (compar=="md") comparison<-minus_del_cut
  else stop
   for (orf in (1:length(genelist))) {
listforone<-c(genelist[orf], numbers(genelist[orf], compar), 
              annoone(genelist[orf]))
mylist<-rbind(mylist, listforone)}
return(mylist)
}


#SGD intergenic spaces

#sfp1upstream<-scan(file="./SGD retrieves/list_aaaawtttt upstream.csv", what="character")


#a4wt4<-function (var) {
 # var2<-intersect(var[,1],sfp1upstream)
  #write.table(var2, file=paste("./SGD retrieves/Intersect/",
                               deparse(substitute(var)), "_has_a4wt4", ".csv"))
  #return (var2)
#}



library(seqinr)


threehund<-function (fasta_file) {
  fasta_file2<-lapply(fasta_file, function (elt) elt[701:1000])
  write.fasta(sequences=fasta_file2, names=names(fasta_file), 
              file.out=paste(deparse(substitute(fasta_file)), "_300", ".fasta", sep=""))
}


#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}


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

#function for tables (genelist, compar1, compar2)
exttable <- function (genelist, compar1, compar2) {
write.table(numberslist(genelist, compar1), 
             paste("./for LN/", deparse(substitute(genelist)), substring(deparse(substitute(compar1)),2,3), ".xls", sep=""), 
             col.names=NA, sep="\t") 
write.table(numberslist(genelist, compar2), 
             paste("./for LN/", deparse(substitute(genelist)), substring(deparse(substitute(compar2)),2,3), ".xls", sep=""), 
             col.names=NA, sep="\t") 
write.table(annolist(genelist), paste("./for LN/", deparse(substitute(genelist)), ".xls", sep=""), sep="\t", col.names=NA)
}



#install.packages(ggplot2)
library(ggplot2)



alldata <-read.table("../bioinfo/raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, del-minus, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
minus_del<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
del_minus<-topTable(fit2, coef=4, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
