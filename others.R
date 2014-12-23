#We and others

chowdhury<-read.table("./for LN/others/Chowdhury 2010.csv", header=TRUE, sep="\t")
chowdhury<-chowdhury[order(chowdhury[2]),]
for (i in (1:nrow(chowdhury))) {
  chowdhury[i,2]=-chowdhury[i,2]
}
chowdhury_inv<-chowdhury[order(chowdhury[2]),]
chowdhury_inv<-chowdhury_inv[order(chowdhury_inv[2]),]
write.table(chowdhury_inv, "./for LN/others/chowdhury_inv.xls", sep="\t", col.names=NA)



jorgensen1<-read.table("./for LN/others/Jorgensen 2002-1.csv", header=TRUE, sep="\t")
jorgensen1<-jorgensen1[order(jorgensen1[2]),]
write.table(jorgensen1, "./for LN/others/jorgensen1.xls", sep="\t", col.names=NA)

jorgensen2<-read.table("./for LN/others/Jorgensen 2002-2.csv", header=TRUE, sep="\t")
jorgensen2<-jorgensen2[order(jorgensen2[2]),]
write.table(jorgensen2, "./for LN/others/jorgensen2.xls", sep="\t", col.names=NA)
