#gene list size challenge


alldata <-read.table("./data/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, del-minus, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
del_minus<-topTable(fit2, coef=4, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="logFC")



pval<-c(.05, .01, .001, .0001, 1e-5, 1e-6)
for (i in 1:length(pval)) {
  a<-nrow(subsetplus_minus
  
}