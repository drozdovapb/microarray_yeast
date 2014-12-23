setwd("/media/sda6/Polina/Work/Ingenuity/")

##http://www.bioconductor.org/help/workflows/gene-regulation-tfbs/

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("MotifDb",  "GenomicFeatures", 
#           "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
#           "org.Sc.sgd.db", "BSgenome.Scerevisiae.UCSC.sacCer3",
#           "motifStack", "seqLogo"))

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Sc.sgd.db)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

##1000 bp upstream, 0 bp downstream for the sake of simplicity, as in the original example

pfm.sfp1 <- query(MotifDb, "SFP1")

