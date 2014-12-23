setwd("/media/sda6/Polina/Work/Ingenuity/raw_data_web/GEO")

#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")

library(GEOquery)
chipollina2007 <-getGEO("GSE5238")
chua2006<-getGEO("GSE5499")

#show(chua2006)
#show(pData(phenoData(chua2006)))