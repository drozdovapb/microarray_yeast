# for linux 
setwd("/media/sda6/Polina/Work/Ingenuity")

source("./2. programs and scripts/functions.R")

SCRAndRibi<-read.table("./3. exploratory/new ideas/GO RP Ribi/RP+Ribi.csv", sep="\t", head=TRUE)
RPLs <-read.table("./3. exploratory/new ideas/GO RP Ribi/RPLgenes.txt", sep="\t", head=TRUE)
RPSs <-read.table("./3. exploratory/new ideas/GO RP Ribi/RPSgenes.txt", sep="\t", head=TRUE)

SCR<-union(SCRAndRibi[,1], SCRAndRibi[,2])
Ribi<-union(SCRAndRibi[,3], SCRAndRibi[,4])
RibiMan<-as.vector(SCRAndRibi[1:21,3])
RibiComp<-as.vector(SCRAndRibi[1:179,4])

comparisons<-c("plus_minus", "minus_del", "plus_del")
pvalues<-c(.05, .01, .001)
abslogFCs<-c(0,.5,1,2,3,4)

#SCR for pm
resultsSCR_pm<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsSCR_pm[j,k]<-length(intersect(SCR, cutoff2(
  mytable=plus_minus, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsSCR_pm)<-abslogFCs
#write.csv(resultsSCR_pm,"./new ideas/GO RP Ribi/Results/resultsSCR_pm")

#SCR for pd
resultsSCR_pd<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsSCR_pd[j,k]<-length(intersect(SCR, cutoff2(
  mytable=plus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsSCR_pd)<-abslogFCs
#write.csv(resultsSCR_pd,"./new ideas/GO RP Ribi/Results/resultsSCR_pd")

#SCR for md
resultsSCR_md<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsSCR_md[j,k]<-length(intersect(SCR, cutoff2(
  mytable=minus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsSCR_md)<-abslogFCs
#write.csv(resultsSCR_md,"./new ideas/GO RP Ribi/Results/resultsSCR_md")

#Ribi for pm
resultsRibi_pm<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibi_pm[j,k]<-length(intersect(Ribi, cutoff2(
  mytable=plus_minus, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibi_pm)<-abslogFCs
#write.csv(resultsRibi_pm,"./new ideas/GO RP Ribi/Results/resultsRibi_pm")

intersect(Ribi, cutoff2(mytable=plus_minus, p.value=.05, abslogFC=1)[,1])

#Ribi for pd
resultsRibi_pd<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibi_pd[j,k]<-length(intersect(Ribi, cutoff2(
  mytable=plus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibi_pd)<-abslogFCs
#write.csv(resultsRibi_pd,"./new ideas/GO RP Ribi/Results/resultsRibi_pd")

#Ribi for md
resultsRibi_md<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibi_md[j,k]<-length(intersect(Ribi, cutoff2(
  mytable=minus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibi_md)<-abslogFCs
#write.csv(resultsRibi_md,"./new ideas/GO RP Ribi/Results/resultsRibi_md")



intersect(SCR, cutoff2(plus_minus, abslogFC=1)[,1])


length(Ribi)
length(SCR)
RPs<-union(RPLs[,1], RPSs[,1])
RPsIntersect<-intersect(RPLs[,1], RPSs[,1])
length(RPs)

#RPs for pm
resultsRPs_pm<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRPs_pm[j,k]<-length(intersect(RPs, cutoff2(
  mytable=plus_minus, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRPs_pm)<-abslogFCs
write.csv(resultsRPs_pm,"./new ideas/GO RP Ribi/Results/resultsRPs_pm")

#RPs for pd
resultsRPs_pd<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRPs_pd[j,k]<-length(intersect(RPs, cutoff2(
  mytable=plus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRPs_pd)<-abslogFCs
#write.csv(resultsRPs_pd,"./new ideas/GO RP Ribi/Results/resultsRPs_pd")

#RPs for md
resultsRPs_md<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRPs_md[j,k]<-length(intersect(RPs, cutoff2(
  mytable=minus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRPs_md)<-abslogFCs
#write.csv(resultsRPs_md,"./new ideas/GO RP Ribi/Results/resultsRPs_md")


#RibiMan for pm
resultsRibiMan_pm<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiMan_pm[j,k]<-length(intersect(RibiMan, cutoff2(
  mytable=plus_minus, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiMan_pm)<-abslogFCs
#write.csv(resultsRibiMan_pm,"./new ideas/GO RP Ribi/Results/resultsRibiMan_pm")

#RibiMan for pd
resultsRibiMan_pd<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiMan_pd[j,k]<-length(intersect(RibiMan, cutoff2(
  mytable=plus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiMan_pd)<-abslogFCs
#write.csv(resultsRibiMan_pd,"./new ideas/GO RP Ribi/Results/resultsRibiMan_pd")

#RibiMan for md
resultsRibiMan_md<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiMan_md[j,k]<-length(intersect(RibiMan, cutoff2(
  mytable=minus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiMan_md)<-abslogFCs
#write.csv(resultsRibiMan_md,"./new ideas/GO RP Ribi/Results/resultsRibiMan_md")


#RibiComp for pm
resultsRibiComp_pm<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiComp_pm[j,k]<-length(intersect(RibiComp, cutoff2(
  mytable=plus_minus, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiComp_pm)<-abslogFCs
#write.csv(resultsRibiComp_pm,"./new ideas/GO RP Ribi/Results/resultsRibiComp_pm")

#RibiComp for pd
resultsRibiComp_pd<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiComp_pd[j,k]<-length(intersect(RibiComp, cutoff2(
  mytable=plus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiComp_pd)<-abslogFCs
#write.csv(resultsRibiComp_pd,"./new ideas/GO RP Ribi/Results/resultsRibiComp_pd")

#RibiComp for md
resultsRibiComp_md<-data.frame(row.names=pvalues)
{for (j in 1:length(pvalues))
{for (k in 1:length(abslogFCs))
{resultsRibiComp_md[j,k]<-length(intersect(RibiComp, cutoff2(
  mytable=minus_del, p.value=pvalues[j], abslogFC=abslogFCs[k])[,1]))}}}
names(resultsRibiComp_md)<-abslogFCs
#write.csv(resultsRibiComp_md,"./new ideas/GO RP Ribi/Results/resultsRibiComp_md")
