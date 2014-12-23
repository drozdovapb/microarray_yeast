#annotation function
# for linux 
setwd("/media/sda6/Polina/Work/Ingenuity")

#source("http://bioconductor.org/biocLite.R")
#biocLite("yeast.db0")
#library(yeast.db0)
#couldn't find out how to use it

library (limma)

alldata <-read.table("./data/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, del-minus, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

#extracting variables
plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="p")
minus_del<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="p")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="p")
del_minus<-topTable(fit2, coef=4, number=nrow(fit2), adjust.method="BH", sort.by="logFC")

#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}

plus_minus_cut <- cutoff(plus_minus, 0.5, 0.001)
plus_del_cut <- cutoff(plus_del, 0.5, 0.001)
minus_del_cut <- cutoff(minus_del, 0.5, 0.001)
del_minus_cut<-cutoff(del_minus, 0.5, 0.001)


genestable <- read.csv("./SGD retrieves/genelist2.csv", sep=",")

#for one gene
annoone <- function(geneid) {
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