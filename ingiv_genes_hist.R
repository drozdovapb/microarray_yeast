# for linux 
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/")

library(limma)

alldata <-read.table("./raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="p")
minus_del<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="p")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="p")

#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}

#applying function!!!

plus_minus_cut <- cutoff(plus_minus, 0.5, 0.001)
plus_del_cut <- cutoff(plus_del, 0.5, 0.001)
minus_del_cut <- cutoff(minus_del, 0.5, 0.001)

# splitting into upregulated and downregulated
plus_minus_cut_up <- subset(plus_minus_cut, plus_minus_cut$logFC > 0)
plus_minus_cut_down <- subset(plus_minus_cut, plus_minus_cut$logFC < 0)
plus_del_cut_up <-subset(plus_del_cut, logFC > 0) 
plus_del_cut_down <- subset(plus_del_cut, logFC < 0)
minus_del_cut_up <- subset(minus_del_cut, logFC > 0)
minus_del_cut_down <- subset(minus_del_cut, logFC <0)


#sorting by FC, by descending abs(!) logFC
plus_minus_cut_up <- plus_minus_cut_up[order(-plus_minus_cut_up$logFC),]
plus_minus_cut_down <- plus_minus_cut_down[order(plus_minus_cut_down$logFC),]
plus_del_cut_up <-plus_del_cut_up[order(-plus_del_cut_up$logFC),] 
plus_del_cut_down <- plus_del_cut_down[order (plus_del_cut_down$logFC),]
minus_del_cut_up <- minus_del_cut_up[order(-minus_del_cut_up$logFC),]
minus_del_cut_down <- minus_del_cut_down [order(minus_del_cut_down$logFC),]


allgenesregulated <- union(plus_minus_cut$ID, union(plus_del_cut$ID,minus_del_cut$ID))

genestable <- read.csv("./SGD retrieves/S_cerevisiae_AROS_genelist.csv")
allyeastgenes <- genestable$systematic_name

ind_barplot <- function(ID) {
  datasubset<-alldata[ID,]
  datasubsetmeans<-c(mean(as.numeric(datasubset[1:4])),
                     mean(as.numeric(datasubset[5:8])),
                     mean(as.numeric(datasubset[9:12])))
  names(datasubsetmeans)<-c("[ISP+]", "[isp-]", "sfp1Δ")  
  sds<-c(sd(as.numeric(datasubset[1:4])),
              sd(as.numeric(datasubset[5:8])),
              sd(as.numeric(datasubset[9:12])))
  confints<-sds*1.96/2
  myceiling<-ceiling(max((datasubsetmeans))+2)
  myfilename=paste("./my/",ID,".png")
  #png(myfilename, bg="transparent")
  mybarplot<-barplot(datasubsetmeans, ylab="log2 интенсивности свечения",
                     col="chocolate3",
                     ylim=c(0,myceiling), las=1, yaxp=c(0,myceiling, myceiling), cex.axis=0.8, 
                     xlab=ID)
    for (i in (1:length(mybarplot))) {
      arrows(x0=mybarplot[i], y0=datasubsetmeans[i]+confints[i],
             y1=datasubsetmeans[i]-confints[i], angle=90, code=3)
  text(x=mybarplot[i], y=datasubsetmeans[i]+2, labels="***")
    }
#  dev.off()
  #cropped
  myfilenamecrop=paste("./my/",ID,"_crop.png")
  png(myfilenamecrop, bg="transparent")
  mybarplot<-barplot(datasubsetmeans, col="chocolate3",
                     ylim=c(0,myceiling), yaxp=c(0,myceiling,myceiling), las=1, xaxt='n',lwd=6, cex.axis=2)
  for (i in (1:length(mybarplot))) {
    arrows(x0=mybarplot[i], y0=datasubsetmeans[i]+confints[i],
           y1=datasubsetmeans[i]-confints[i], angle=90, code=3,lwd=6)
    }
  dev.off()
}

#RPS9A
ind_barplot("YPL081W")
#PRS9B
ind_barplot("YBR189W")
#CRF1

ind_barplot("YDR223W")

#SUP35
ind_barplot("YDR172W")

#SUP45
ind_barplot("YBR143C")

#RPI1
ind_barplot("YIL119C")

#FIT2
svg("FIT2.svg")
ind_barplot("YOR382W")
#FIG1
svg("FIG1.svg")
ind_barplot("YOR382W")
#AFT1
ind_barplot("YGL071W")



ind_barplot2 <- function(ID) {
  datasubset<-alldata[ID,]
  datasubsetmeans<-c(mean(as.numeric(datasubset[5:8])),
                     mean(as.numeric(datasubset[1:4])))
  names(datasubsetmeans)<-as.factor(c("[isp-]", "[ISP+]"))  
  sds<-c(sd(as.numeric(datasubset[5:8])),
         sd(as.numeric(datasubset[1:4])))
  confints<-sds*1.96/2
  myceiling<-ceiling(max((datasubsetmeans))+2)
  mybarplot<-barplot(datasubsetmeans, col="chocolate3",
                     ylim=c(0,myceiling), yaxp=c(0,myceiling,myceiling), las=1, xaxt='n',lwd=6, cex.axis=2)
  for (i in (1:length(mybarplot))) {
    arrows(x0=mybarplot[i], y0=datasubsetmeans[i]+confints[i],
           y1=datasubsetmeans[i]-confints[i], angle=90, code=3,lwd=6)}
}
  

setwd("/media/drozdovapb//441DB6652D7ED741//Work//Current work/Theses/Master's thesis/exhibitory")
#FIT2
png("2014-04-13_FIT2.png"); ind_barplot2("YOR382W"); dev.off()
#FIG1
png("2014-04-13_FIG1.png"); ind_barplot2("YOR382W"); dev.off()
#AFT1
png("2014-04-13_AFT1.png"); ind_barplot2("YGL071W"); dev.off()
#RPS9A
png("2014-04-13_RPS9A.png"); ind_barplot2("YPL081W"); dev.off()
#PRS9B
png("2014-04-13_RPS9B.png"); ind_barplot2("YBR189W"); dev.off()
#CRF1
png("2014-04-13_CRF1.png"); ind_barplot2("YDR223W"); dev.off()
#SUP35
png("2014-04-13_SUP35.png"); ind_barplot2("YDR172W"); dev.off()
#SUP45
png("2014-04-13_SUP45.png"); ind_barplot2("YBR143C"); dev.off()
