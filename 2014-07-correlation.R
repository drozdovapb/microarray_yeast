#is there any correlation?

#sfp1delta

setwd("/media/drozdovapb//441DB6652D7ED741/Work/Current work//Projects/bioinfo/processed_tables_R/")
sfp1delta<-read.table("./logdata.csv")

slowgrowth.signature <- "/media/drozdovapb/441DB6652D7ED741/Work/BigData/ODubhir/Supplementary_data4.txt"


read.tab <- function(file, ...)read.table(file=file, sep="\t", as.is=TRUE, quote="",
                                          header=TRUE, comment.char="", row.names=1, ...)
write.tab <- function(x, file="", ...)write.table(x, file, sep="\t", quote=FALSE, na="",
                                                  eol="\r\n",
                                                  row.names=TRUE, col.names=NA, ...)

sig.data <- read.tab(file=slowgrowth.signature)
numerical.cols <- names(which(unlist(lapply(sig.data,class)=="numeric")))
signature <- as.matrix(sig.data[,numerical.cols, drop=FALSE])

my<-merge(sfp1delta, signature, by="row.names")

library(graphics)
plot(my$plus_min, my$M)
smoothScatter(my$plus_min, my$M, colram=colorRampPalette(c("white", "black")))
plot(my$del_min, my$M)
cor(my$del_min, my$M)
smoothScatter(my$del_min, my$M, colram=colorRampPalette(c("white", "black")))
smoothScatter(my$plus_del, my$M, colram=colorRampPalette(c("white", "black")))
plot(my$plus_del, my$M)

#sup45
sup45mut <- read.table("/media/drozdovapb//441DB6652D7ED741//Work//Current work/Projects/SvM/2014-08-07-orderratio_sup45_wt_man_edits.csv", sep="\t", head=T)
#for (i in (1:nrow(sup45mut))) {
#    if (sup45mut$ORF[i] == "Unknown") sup45mut$ORF[i]<-sup45mut$Gene[i]
#}

sup45mut <- subset(sup45mut, !duplicated(ORF))
row.names(sup45mut)<-sup45mut$ORF
mytwo <- merge(sup45mut, signature, by="row.names")
plot(log(mytwo$FC), mytwo$M)


#Jorgensen

Jorg<- read.table("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/sfp1delta/manually processed_DEG and logFC_by authors/Jorgensen_04.csv", sep=",", head=T)
row.names(Jorg)<- Jorg$ORFs
mythree <- merge(Jorg, signature, by="row.names")

plot(mythree[,9], mythree$M)
plot(mythree[,9], mythree$M)
    cor(mythree[,9], mythree$M, use="complete")

par(mfrow=c(3,4))
for (i in 3:13) {
    plot(mythree[,i], mythree$M, xlab=colnames(mythree[i]), main=cor(mythree[,i], mythree$M, use="complete"))
}
         

#Reimand
Reim <- read.table("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/sfp1delta/manually processed_DEG and logFC_by authors/Reimand_10.txt", head=T)
row.names(Reim)<-Reim$ORF
ReimSig <- merge(Reim, signature, by="row.names")
plot(ReimSig$Reimand_sfp1, ReimSig$M)
cor(ReimSig$Reimand_sfp1, ReimSig$M)

