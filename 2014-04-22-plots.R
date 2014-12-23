par(las=1)

#microarray
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/")
ma_data <- read.table("./raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", head=T)
RPS9A_ma <- t(subset(ma_data, locus=="YPL081W")[,2:9])
FIG1_ma <- t(subset(ma_data, locus=="YBR040W")[,2:9]) 
FIT2_ma <- t(subset(ma_data, locus=="YOR382W")[,2:9]) 
ma_sel<- cbind(RPS9A_ma, FIG1_ma, FIT2_ma)
colnames(ma_sel) <- c("RPS9A", "FIG1", "FIT2")
strains<-as.factor(c(rep("ISP+",4), rep("isp-", 4)))
#svg("./figures_R/2014-04-23_RPS9A.svg", width=4, height=3)
boxplot(ma_sel[,1] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(RPS9A))), ylim=c(8,16))
points(ma_sel[,1] ~ strains, pch=6)
#dev.off()
#svg("./figures_R/2014-04-23_FIG1.svg", width=4, height=3)
boxplot(ma_sel[,2] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(FIG1))), ylim=c(4,12))
points(ma_sel[,2] ~ strains, pch=6)
#dev.off()
#svg("./figures_R/2014-04-23_FIT2.svg", width=4, height=3)
boxplot(ma_sel[,3] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(FIT2))), ylim=c(6,14))
points(ma_sel[,3] ~ strains, pch=6)
#dev.off()



#FIG1_FIT2_qPCR
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/Fe/")
pmd<-read.csv("./plots_tables/2014-03-06-qPCRs.csv")
pmonly<-droplevels(subset(pmd, pmd$strain=="[ISP+]" | pmd$strain=="[isp-]"))
#svg( filename="./plots_tables/2014-03-09-FIG1 and FIT2_pm.svg", width=8, height=3)
par(mfrow=(c(1,2)))
boxplot(pmonly$dFIG1 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(FIG1))), ylim=c(-18, -10))
points(pmonly$dFIG1 ~ pmonly$strain, pch=19)
boxplot(pmonly$dFIT2 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(FIT2))), ylim=c(-8, 0))
points(pmonly$dFIT2 ~ pmonly$strain, pch=19)
dev.off()

#RPS9A_qPCR
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/ribosomes/")
alldata<-read.csv("./plots_tables/2014-03-26-values.csv", sep="\t", head=T)
pmdatasyntol<-droplevels(subset(alldata, alldata$kit=="Syntol"))
pmdatasyntol<-droplevels(subset(pmdatasyntol, pmdatasyntol$strain=="[ISP+]" | pmdatasyntol$strain=="[isp-]"))
pmdRPS9Asyntol<- -pmdatasyntol$meanRPS9A+pmdatasyntol$meanADH1
boxplot(pmdRPS9Asyntol ~ pmdatasyntol$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main=substitute(expression(italic(RPS9A))), las=1, ylim=c(-12,-4))
points(pmdRPS9Asyntol ~ pmdatasyntol$strain, pch=19)


#for the paper
setwd("/media/drozdovapb//441DB6652D7ED741//Work//Current work/Theses/Our papers/something-never-to-be-published/2014-04-the second round/")
#RPS9A_both
svg("2014-04-23_RPS9A.svg", width=6, height=3)
par(mfrow=(c(1,2)), las=1)
boxplot(ma_sel[,1] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main="RPS9A (microarray)", ylim=c(8,16))
points(ma_sel[,1] ~ strains, pch=6)
boxplot(pmdRPS9Asyntol ~ pmdatasyntol$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main="RPS9A (qPCR)", ylim=c(-12,-4))
points(pmdRPS9Asyntol ~ pmdatasyntol$strain, pch=19)
dev.off()

#FIG1, FIT2_both
svg("2014-04-23_FIG1,FTI2.svg", width=6, height=6)
par(mfrow=(c(2,2)), las=1)
boxplot(ma_sel[,2] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main="FIG1 (microarray)", ylim=c(4,12))
points(ma_sel[,2] ~ strains, pch=6)
boxplot(pmonly$dFIG1 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main="FIG1 (qPCR)", ylim=c(-18, -10))
points(pmonly$dFIG1 ~ pmonly$strain, pch=19)
boxplot(ma_sel[,3] ~ strains, ylab="log2 fluorescence intensity", boxlty=0, whisklty=0, staplelty=0, main="FIT2 (microarray)", ylim=c(6,14))
points(ma_sel[,3] ~ strains, pch=6)
boxplot(pmonly$dFIT2 ~ pmonly$strain, ylab="ΔCq", boxlty=0, whisklty=0, staplelty=0, main="FIT2 (qPCR)", ylim=c(-8, 0))
points(pmonly$dFIT2 ~ pmonly$strain, pch=19)
dev.off()