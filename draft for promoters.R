setwd("/media/sda6/Polina/Work/Ingenuity/TFBSs")

#install.packages("seqinr")
library(seqinr)

promoters_plus_del <- read.fasta("promoters_plus_del.fasta", seqtype="DNA")
promoters_plus_minus<-read.fasta("promoters_plus_minus.fasta", seqtype="DNA")
promoters_minus_del <- read.fasta("promoters_minus_del.fasta", seqtype="DNA")

library(Biostrings)

#monont
pro_plus_minus_seq<-getSequence(promoters_plus_minus)
plus_minus_freq_table<-sapply(pro_plus_minus_seq, count, wordsize=1, freq=TRUE)
plus_minus_freq<-rowMeans(plus_minus_freq_table)
write.table(plus_minus_freq,"plus_minus_freq.csv")


pro_plus_del_seq<-getSequence(promoters_plus_del)
plus_del_freq_table<-sapply(pro_plus_del_seq, count, wordsize=1, freq=TRUE)
plus_del_freq<-rowMeans(plus_del_freq_table)
write.table(plus_del_freq,"plus_del_freq.csv")

pro_minus_del_seq<-getSequence(promoters_minus_del)
minus_del_freq_table<-sapply(pro_minus_del_seq, count, wordsize=1, freq=TRUE)
minus_del_freq<-rowMeans(minus_del_freq_table)
write.table(minus_del_freq,"minus_del_freq.csv")

#di
pro_plus_minus_seq<-getSequence(promoters_plus_minus)
plus_minus_freq_table2<-sapply(pro_plus_minus_seq, count, wordsize=2, freq=T)
plus_minus_freq2<-rowMeans(plus_minus_freq_table2, na.rm=T)
plus_minus_freq2<-plus_minus_freq2*100
write.table(plus_minus_freq2,"plus_minus_freq2.csv")


pro_plus_del_seq<-getSequence(promoters_plus_del)
plus_del_freq_table2<-sapply(pro_plus_del_seq, count, wordsize=2, freq=T)
plus_del_freq2<-rowMeans(plus_del_freq_table2, na.rm=T)
plus_del_freq2<-plus_del_freq2*100
write.table(plus_del_freq2,"plus_del_freq2.csv")

pro_minus_del_seq<-getSequence(promoters_minus_del)
minus_del_freq_table2<-sapply(pro_minus_del_seq, count, wordsize=2, freq=T)
minus_del_freq2<-rowMeans(minus_del_freq_table2, na.rm=T)
minus_del_freq2<-minus_del_freq2*100
write.table(minus_del_freq2,"minus_del_freq2.csv")




#group 1, coexpressed in plus and gel

group1up <- intersect(plus_minus_cut_up$ID, del_minus_cut_up$ID)
group1down <- intersect(plus_minus_cut_down$ID, del_minus_cut_down$ID)

#group 2, differentially expressed in plus and del
group2up<-intersect(plus_minus_cut_up$ID, minus_del_cut_up$ID)
group2down<-intersect(plus_minus_cut_down$ID, minus_del_cut_down$ID)

#group 3, prion specific
group3up <- setdiff(plus_del_cut_up$ID, minus_del_cut_up$ID)
group3down <- setdiff(plus_del_cut_down$ID, minus_del_cut_down$ID)


setwd("/media/sda6/Polina/Work/Ingenuity/TFBSs/promoters_new")

threehund<-function (fasta_file) {
  fasta_file2<-lapply(fasta_file, function (elt) elt[701:1000])
  write.fasta(sequences=fasta_file2, names=names(fasta_file), 
              file.out=paste(deparse(substitute(fasta_file)), "_300", ".fasta", sep=""))
}


promoters_plus_minus_up<-read.fasta("promoters_plus_minus_up.fasta", as.string=FALSE)
threehund(promoters_plus_minus_up)

promoters_plus_minus_down<-read.fasta("promoters_plus_minus_down.fasta", as.string=FALSE)
threehund(promoters_plus_minus_down)

promoters_minus_del_up<-read.fasta("promoters_minus_del_up.fasta", as.string=FALSE)
threehund(promoters_minus_del_up)

promoters_minus_del_down<-read.fasta("promoters_minus_del_down.fasta", as.string=FALSE)
threehund(promoters_minus_del_down)

promoters_plus_del_up<-read.fasta("promoters_plus_del_up.fasta", as.string=FALSE)
threehund(promoters_plus_del_up)

promoters_plus_del_down<-read.fasta("promoters_plus_del_down.fasta", as.string=FALSE)
threehund(promoters_plus_del_down)

promoters_group1_up<-read.fasta("promoters_group1_up.fasta", as.string=FALSE)
threehund(promoters_group1_up)

promoters_group1_down<-read.fasta("promoters_group1_down.fasta", as.string=FALSE)
threehund(promoters_group1_down)

promoters_group2_up<-read.fasta("promoters_group2_up.fasta", as.string=FALSE)
threehund(promoters_group2_up)

promoters_group2_down<-read.fasta("promoters_group2_down.fasta", as.string=FALSE)
threehund(promoters_group2_down)

promoters_group3_up<-read.fasta("promoters_group3_up.fasta", as.string=FALSE)
threehund(promoters_group3_up)

promoters_group3_down<-read.fasta("promoters_group3_down.fasta", as.string=FALSE)
threehund(promoters_group3_down)