par(las=1)

#microarray
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/")
ma_data <- read.table("./raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", head=T)
RPS9A_ma <- t(subset(ma_data, locus=="YPL081W")[,2:9])
FIG1_ma <- t(subset(ma_data, locus=="YBR040W")[,2:9]) 
FIT2_ma <- t(subset(ma_data, locus=="YOR382W")[,2:9]) 
SUP45_ma <- t(subset(ma_data, locus=="YBR143C")[,2:9]) 

ma_sel<- cbind(RPS9A_ma, FIG1_ma, FIT2_ma, SUP45_ma)
colnames(ma_sel) <- c("RPS9A", "FIG1", "FIT2", "SUP45")
strains<-as.factor(c(rep("ISP+",4), rep("isp-", 4)))

#FIG1_FIT2_qPCR
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/Fe/")
pmd<-read.csv("./plots_tables/2014-03-06-qPCRs.csv")
pmonly<-droplevels(subset(pmd, pmd$strain=="[ISP+]" | pmd$strain=="[isp-]"))

#RPS9A_qPCR
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/ribosomes/")
alldata<-read.csv("./plots_tables/2014-03-26-values.csv", sep="\t", head=T)
pmdatasyntol<-droplevels(subset(alldata, alldata$kit=="Syntol"))
pmdatasyntol<-droplevels(subset(pmdatasyntol, pmdatasyntol$strain=="[ISP+]" | pmdatasyntol$strain=="[isp-]"))
pmdRPS9Asyntol<- -pmdatasyntol$meanRPS9A+pmdatasyntol$meanADH1
boxplot(pmdRPS9Asyntol ~ pmdatasyntol$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(RPS9A))), las=1, ylim=c(-12,-4))
points(pmdRPS9Asyntol ~ pmdatasyntol$strain, pch=19)

#SUP45_qPCR
setwd("../Sup35, Sup45, Sfp1/")
alldata<-read.table("./plots_tables/2014-03-28-SUP35, SUP45_ratios.csv", head=T, sep=",")
#library(coin); wilcox_test(alldata$dmeanSUP45~alldata$strain) #p = .04953
dmeanSUP45<-alldata$dmeanSUP45
S45strain<-alldata$strain


#для диплома
setwd("/media/drozdovapb//441DB6652D7ED741//Work//Current work/Theses/Master's thesis//exhibitory")
#RPS9A_both
svg("./2014-05-12_RPS9A2.svg", width=6, height=4)
par(mfrow=(c(1,2)), las=1)
boxplot(ma_sel[,1] ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, ylim=c(8,16)) #main="RPS9A (транскриптом)"
points(ma_sel[,1] ~ strains, pch=6)
boxplot(pmdRPS9Asyntol ~ pmdatasyntol$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, ylim=c(-12,-4)) #main="RPS9A (транскриптом)"
points(pmdRPS9Asyntol ~ pmdatasyntol$strain, pch=19)
dev.off()

#FIG1, FIT2_both
svg("2014-05-12_FIG1,FTI2.svg", width=6, height=8)
par(mfrow=(c(2,2)), las=1)
boxplot(ma_sel[,2] ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, ylim=c(4,12)) #main="FIG1 (транскриптом)"
points(ma_sel[,2] ~ strains, pch=6)
boxplot(pmonly$dFIG1 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, ylim=c(-18, -10)) #main="FIG1 (ОТ-ПЦРРВ)"
points(pmonly$dFIG1 ~ pmonly$strain, pch=19)
boxplot(ma_sel[,3] ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, ylim=c(6,14)) #main="FIT2 (транскриптом)"
points(ma_sel[,3] ~ strains, pch=6)
boxplot(pmonly$dFIT2 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, ylim=c(-8, 0)) #main="FIT2 (ОТ-ПЦРРВ)"
points(pmonly$dFIT2 ~ pmonly$strain, pch=19)
dev.off()

#SUP45_both
svg("2014-05-12_SUP45.svg", width=6, height=4)
par(mfrow=(c(1,2)), las=1)
boxplot(ma_sel[,4] ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, ylim=c(12,16)) #main="FIG1 (транскриптом)"
points(ma_sel[,4] ~ strains, pch=6)
boxplot(dmeanSUP45 ~ S45strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, ylim=c(-8, -4)) #main="FIG1 (ОТ-ПЦРРВ)"
points(dmeanSUP45 ~ S45strain, pch=19)
dev.off()


#ma_AFT1 2014-05-21
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/")
ma_data <- read.table("./raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", head=T)
AFT1_ma <- t(subset(ma_data, locus=="YGL071W")[,2:9])
boxplot(AFT1_ma ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, main="AFT1 (транскриптом)", ylim=c(8,16))
points(AFT1_ma ~ strains, ylab="log2 интенсивности свечения", main="AFT1 (транскриптом)", pch=19)

#qPCR_Aft1
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/Fe/")
pmd<-read.csv("./plots_tables/2014-03-06-qPCRs.csv")
pmonly<-droplevels(subset(pmd, pmd$strain=="[ISP+]" | pmd$strain=="[isp-]"))


setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Theses//Master's thesis/")
svg("2014-05-26_AFT1.svg", width=6, height=4)
par(mfrow=(c(1,2)), las=1)
boxplot(AFT1_ma ~ strains, ylab="log2 интенсивности свечения", boxlty=0, whisklty=0, staplelty=0, main="AFT1 (транскриптом)", ylim=c(8,16))
points(AFT1_ma ~ strains, ylab="log2 интенсивности свечения", main="AFT1 (транскриптом)", pch=19)
boxplot(pmonly$dAFT1 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(AFT1))))
points(pmonly$dAFT1 ~ pmonly$strain, pch=19)
dev.off()