#control promoters

setwd("/media/4D9308A6770B38AB/Polina/Work/Ingenuity")

allgenesregulated <- union(plus_minus_cut$ID, union(plus_del_cut$ID,minus_del_cut$ID))
                
genestable <- read.csv("./SGD retrieves/S_cerevisiae_AROS_genelist.csv")
allyeastgenes <- genestable$systematic_name

nonregulatedgenes <- setdiff(allyeastgenes, allgenesregulated)
                              
length(allyeastgenes)
length(nonregulatedgenes)
length(allgenesregulated)                              
                              
write.table(nonregulatedgenes, "./SGD retrieves/nonregulatedgenes.csv")
write.table(allgenesregulated, "./SGD retrieves/regulatedgenes.csv")

# get frequencies
#install.packages("seqinr")
library(seqinr)
                              
promoters_nonreg <- read.fasta("./SGD retrieves/nonreg_promoters.fasta", 
                               seqtype="DNA")
                              
                              
nonreg_seq<-getSequence(promoters_nonreg)
nonreg_freq_table<-sapply(nonreg_seq, count, wordsize=2, freq=TRUE)
nonreg_freq<-rowMeans(nonreg_freq_table)
nonreg_freq<-nonreg_freq*100
write.table(nonreg_freq,"./SGD retrieves/nonreg_freq.csv")


promoters_reg <- read.fasta("./SGD retrieves/reg_promoters.fasta", 
                               seqtype="DNA")


reg_seq<-getSequence(promoters_reg)
reg_freq_table<-sapply(reg_seq, count, wordsize=2, freq=TRUE)
reg_freq<-rowMeans(reg_freq_table, na.rm=T)
reg_freq<-reg_freq*100
write.table(reg_freq,"./SGD retrieves/reg_freq.csv")