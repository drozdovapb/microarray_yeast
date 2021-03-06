# for linux 
setwd("/media/drozdovapb//441DB6652D7ED741//Work//Current work//Projects//bioinfo")

library(limma)

genestable <- read.csv("./raw_data_web/AROS_genelist.csv")
allyeastgenes <- genestable$systematic_name


alldata <-read.table("./raw_data_our/Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, minus-del, plus-del, del-minus, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)

plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
minus_del<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
del_minus<-topTable(fit2, coef=4, number=nrow(fit2), adjust.method="BH", sort.by="logFC")
                    
#subsetting function
cutoff <- function(mytable, abslogFC, p.value) {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)
}

plus_minus_cut <- cutoff(plus_minus, 0.5, 0.001)
plus_del_cut <- cutoff(plus_del, 0.5, 0.001)
minus_del_cut <- cutoff(minus_del, 0.5, 0.001)
del_minus_cut <- cutoff(del_minus, 0.5, 0.001)


# splitting into upregulated and downregulated
plus_minus_cut_up <- subset(plus_minus_cut, plus_minus_cut$logFC > 0)
write.table(plus_minus_cut_up, "plus_minus_cut_up.xls", col.names=NA)
plus_minus_cut_down <- subset(plus_minus_cut, plus_minus_cut$logFC < 0)
plus_del_cut_up <-subset(plus_del_cut, logFC > 0) 
plus_del_cut_down <- subset(plus_del_cut, logFC < 0)
minus_del_cut_up <- subset(minus_del_cut, logFC > 0)
minus_del_cut_down <- subset(minus_del_cut, logFC <0)
del_minus_cut_up <- subset(del_minus_cut, logFC > 0)
del_minus_cut_down <- subset(del_minus_cut, logFC <0)

#group 1, coexpressed in plus and gel

group1up <- intersect(plus_minus_cut_up$ID, del_minus_cut_up$ID)
group1down <- intersect(plus_minus_cut_down$ID, del_minus_cut_down$ID)

#group 2, differentially expressed in plus and del
group2up<-intersect(plus_minus_cut_up$ID, minus_del_cut_up$ID)
group2down<-intersect(plus_minus_cut_down$ID, minus_del_cut_down$ID)

#group 3, prion specific
group3up <- setdiff(plus_del_cut_up$ID, minus_del_cut_up$ID)
group3down <- setdiff(plus_del_cut_down$ID, minus_del_cut_down$ID)

group3all<-setdiff(plus_del_cut$ID, minus_del_cut$ID)

#function for tables (genelist, compar1, compar2)
exttable <- function (genelist, compar1, compar2) {
write.table(numberslist(genelist, compar1), 
             paste("./for LN/", deparse(substitute(genelist)), substring(deparse(substitute(compar1)),2,3), ".xls", sep=""), 
             col.names=NA, sep="\t") 
write.table(numberslist(genelist, compar2), 
             paste("./for LN/", deparse(substitute(genelist)), substring(deparse(substitute(compar2)),2,3), ".xls", sep=""), 
             col.names=NA, sep="\t") 
write.table(annolist(genelist), paste("./for LN/", deparse(substitute(genelist)), ".xls", sep=""), sep="\t", col.names=NA)
}

#exttable(group1up, "pm", "dm")
#exttable(group1down, "pm", "dm")
#exttable(group2up, "pm", "md")
#exttable(group2down, "pm", "md")
#exttable(group3up, "pd", "md")
#exttable(group3down, "pd", "md")
#exttable(group3all, "pd", "md")
#in the case of group3all there should be nothing in "md"...


#merging

#group1
group1downanno<-read.table("./for LN/group1down.xls")
group1downdm<-read.table("./for LN/group1downdm.xls")
group1downpm<-read.table("./for LN/group1downpm.xls")
group1downmerged <- merge(group1downanno, group1downpm, group1downdm)

#group 3

group3allup <- intersect(plus_del_cut_up$ID, group3all)
#exttable(group3allup, "pd", "md")
group3alldown <- intersect(plus_del_cut_down$ID, group3all)
#exttable(group3alldown, "pd", "md")



# new lists (see thoughts 18.01.13)

p_to_m_.001_.5 <- plus_minus_cut
p_to_d_.001_.5 <-plus_del_cut
m_to_d_.001_.5 <-minus_del_cut
d_to_m_.001_.5 <-del_minus_cut
p_to_m_.001_.5_up <- subset(p_to_m_.001_.5, logFC > 0)
p_to_m_.001_.5_down <- subset(p_to_m_.001_.5, logFC < 0)
p_to_d_.001_.5_up <- subset(p_to_d_.001_.5, logFC > 0)
p_to_d_.001_.5_down <- subset(p_to_d_.001_.5, logFC < 0)
m_to_d_.001_.5_up <- subset(m_to_d_.001_.5, logFC > 0)
m_to_d_.001_.5_down <- subset(m_to_d_.001_.5, logFC < 0)
d_to_m_.001_.5_up <- subset(d_to_m_.001_.5, logFC > 0)
d_to_m_.001_.5_down <- subset(d_to_m_.001_.5, logFC < 0)


# addition
p_to_d_.001_1 <-cutoff(plus_del, 1, 0.001)
 p_to_d_.001_1_up<-subset(p_to_d_.001_1, logFC>0)
 p_to_d_.001_1_down<-subset(p_to_d_.001_1, logFC<0)
 m_to_d_.001_1 <-cutoff(minus_del, 1, 0.001)
 m_to_d_.001_1_up<-subset(m_to_d_.001_1, logFC>0)
 m_to_d_.001_1_down<-subset(m_to_d_.001_1, logFC<0)
d_to_m_.001_1<-cutoff(del_minus, 1, 0.001)
d_to_m_.001_1_up<-subset(d_to_m_.001_1, logFC>0)
d_to_m_.001_1_down<-subset(d_to_m_.001_1, logFC<0)

p_to_m_.001_1 <-cutoff(plus_minus, 1, 0.001)
p_to_m_.001_1_up<-subset(p_to_m_.001_1, logFC>0)
p_to_m_.001_1_down<-subset(p_to_m_.001_1, logFC<0)


p_to_d_.001_2 <-cutoff(plus_del, 2, 0.001)
m_to_d_.001_2 <-cutoff(minus_del, 2, 0.001)
p_to_d_.001_3 <-cutoff(plus_del, 3, 0.001)
m_to_d_.001_3 <-cutoff(minus_del, 3, 0.001)
p_to_d_.001_4 <-cutoff(plus_del, 4, 0.001)
m_to_d_.001_4 <-cutoff(minus_del, 4, 0.001)
p_to_d_.001_5 <-cutoff(plus_del, 5, 0.001)
m_to_d_.001_5 <-cutoff(minus_del, 5, 0.001)
p_to_d_.001_5.5 <-cutoff(plus_del, 5.5, 0.001)
m_to_d_.001_5.5 <-cutoff(minus_del, 5.5, 0.001)
p_to_d_.001_6 <-cutoff(plus_del, 6, 0.001)
m_to_d_.001_6 <-cutoff(minus_del, 6, 0.001)

p_to_d_e2_.5 <-cutoff(plus_del, 0.5, 1e-2)
  p_to_d_e2_.5_up <-subset(p_to_d_e2_.5, logFC>0)
  p_to_d_e2_.5_down <-subset(p_to_d_e2_.5, logFC<0)
m_to_d_e2_.5 <-cutoff(minus_del, 0.5, 1e-2)
  m_to_d_e2_.5_up <-subset(m_to_d_e2_.5, logFC>0)
  m_to_d_e2_.5_down <-subset(m_to_d_e2_.5, logFC<0)
p_to_d_e3_.5 <-cutoff(plus_del, 0.5, 1e-3)
  p_to_d_e3_.5_up <-subset(p_to_d_e3_.5, logFC>0)
  p_to_d_e3_.5_down <-subset(p_to_d_e3_.5, logFC<0)
m_to_d_e3_.5 <-cutoff(minus_del, 0.5, 1e-3)
  m_to_d_e3_.5_up <-subset(m_to_d_e3_.5, logFC>0)
  m_to_d_e3_.5_down <-subset(m_to_d_e3_.5, logFC<0)
p_to_d_e4_.5 <-cutoff(plus_del, 0.5, 1e-4)
  p_to_d_e4_.5_up <-subset(p_to_d_e4_.5, logFC>0)
  p_to_d_e4_.5_down <-subset(p_to_d_e4_.5, logFC<0)
m_to_d_e4_.5 <-cutoff(minus_del, 0.5, 1e-4)
  m_to_d_e4_.5_up <-subset(m_to_d_e4_.5, logFC>0)
  m_to_d_e4_.5_down <-subset(m_to_d_e4_.5, logFC<0)
p_to_d_e5_.5 <-cutoff(plus_del, 0.5, 1e-5)
  p_to_d_e5_.5_up <-subset(p_to_d_e5_.5, logFC>0)
  p_to_d_e5_.5_down <-subset(p_to_d_e5_.5, logFC<0)
m_to_d_e5_.5 <-cutoff(minus_del, 0.5, 1e-5)
  m_to_d_e5_.5_up <-subset(m_to_d_e5_.5, logFC>0)
  m_to_d_e5_.5_down <-subset(m_to_d_e5_.5, logFC<0)
p_to_d_e6_.5 <-cutoff(plus_del, 0.5, 1e-6)
  p_to_d_e6_.5_up <-subset(p_to_d_e6_.5, logFC>0)
  p_to_d_e6_.5_down <-subset(p_to_d_e6_.5, logFC<0)
m_to_d_e6_.5 <-cutoff(minus_del, 0.5, 1e-6)
  m_to_d_e6_.5_up <-subset(m_to_d_e6_.5, logFC>0)
  m_to_d_e6_.5_down <-subset(m_to_d_e6_.5, logFC<0)
p_to_d_e7_.5 <-cutoff(plus_del, 0.5, 1e-7)
  p_to_d_e7_.5_up <-subset(p_to_d_e7_.5, logFC>0)
  p_to_d_e7_.5_down <-subset(p_to_d_e7_.5, logFC<0)
m_to_d_e7_.5 <-cutoff(minus_del, 0.5, 1e-7)
  m_to_d_e7_.5_up <-subset(m_to_d_e7_.5, logFC>0)
  m_to_d_e7_.5_down <-subset(m_to_d_e7_.5, logFC<0)
p_to_d_e8_.5 <-cutoff(plus_del, 0.5, 1e-8)
  p_to_d_e8_.5_up <-subset(p_to_d_e8_.5, logFC>0)
  p_to_d_e8_.5_down <-subset(p_to_d_e8_.5, logFC<0)
m_to_d_e8_.5 <-cutoff(minus_del, 0.5, 1e-8)
  m_to_d_e8_.5_up <-subset(m_to_d_e8_.5, logFC>0)
  m_to_d_e8_.5_down <-subset(m_to_d_e8_.5, logFC<0)
p_to_d_e9_.5 <-cutoff(plus_del, 0.5, 1e-9)
  p_to_d_e9_.5_up <-subset(p_to_d_e9_.5, logFC>0)
  p_to_d_e9_.5_down <-subset(p_to_d_e9_.5, logFC<0)
m_to_d_e9_.5 <-cutoff(minus_del, 0.5, 1e-9)
  m_to_d_e9_.5_up <-subset(m_to_d_e9_.5, logFC>0)
  m_to_d_e9_.5_down <-subset(m_to_d_e9_.5, logFC<0)




#intersections
#coexpressed in p and d
p_to_m_up_AND_d_to_m_up <-intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_up[,1])
p_to_m_down_AND_d_to_m_down <-intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_down[,1])

#differentially expressed in p and d
p_to_m_down_AND_d_to_m_up <-intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_up[,1])
p_to_m_up_AND_d_to_m_down <-intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_down[,1])

#"prion specific" genes
p_to_d_.001_.5_AND_NOT_m_to_d_.001_.5 <-setdiff(p_to_d_.001_.5[,1], m_to_d_.001_.5[,1])
p_to_d_.001_1_AND_NOT_m_to_d_.001_1 <- setdiff(p_to_d_.001_1[,1], m_to_d_.001_1[,1])
p_to_d_.001_2_AND_NOT_m_to_d_.001_2 <- setdiff(p_to_d_.001_2[,1], m_to_d_.001_2[,1])
p_to_d_.001_3_AND_NOT_m_to_d_.001_3 <- setdiff(p_to_d_.001_3[,1], m_to_d_.001_3[,1])
p_to_d_.001_4_AND_NOT_m_to_d_.001_4 <- setdiff(p_to_d_.001_4[,1], m_to_d_.001_4[,1])
p_to_d_.001_5_AND_NOT_m_to_d_.001_5 <- setdiff(p_to_d_.001_5[,1], m_to_d_.001_5[,1])
p_to_d_.001_5.5_AND_NOT_m_to_d_.001_5.5 <- setdiff(p_to_d_.001_5.5[,1], m_to_d_.001_5.5[,1])
p_to_d_.001_6_AND_NOT_m_to_d_.001_6 <- setdiff(p_to_d_.001_6[,1], m_to_d_.001_6[,1])

p_to_d_.001_.5_up_AND_NOT_m_to_d_.001_.5_up <- 
                setdiff(p_to_d_.001_.5_up[,1],m_to_d_.001_.5_up[,1])
p_to_d_.001_.5_down_AND_NOT_m_to_d_.001_.5_down <- 
  setdiff(p_to_d_.001_.5_down[,1],m_to_d_.001_.5_down[,1])

p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5 <- 
  setdiff(p_to_d_e2_.5[,1],m_to_d_e2_.5[,1])
p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5 <- 
  setdiff(p_to_d_e3_.5[,1],m_to_d_e3_.5[,1])
p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5 <- 
  setdiff(p_to_d_e4_.5[,1],m_to_d_e4_.5[,1])
p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5 <- 
  setdiff(p_to_d_e5_.5[,1],m_to_d_e5_.5[,1])
p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5 <- 
  setdiff(p_to_d_e6_.5[,1],m_to_d_e6_.5[,1])
p_to_d_e7_.5_AND_NOT_m_to_d_e7_.5 <- 
  setdiff(p_to_d_e7_.5[,1],m_to_d_e7_.5[,1])
p_to_d_e8_.5_AND_NOT_m_to_d_e8_.5 <- 
  setdiff(p_to_d_e8_.5[,1],m_to_d_e8_.5[,1])
p_to_d_e9_.5_AND_NOT_m_to_d_e9_.5 <- 
  setdiff(p_to_d_e9_.5[,1],m_to_d_e9_.5[,1])



p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_up <- 
  setdiff(p_to_d_e2_.5_up[,1],m_to_d_e2_.5_up[,1])
p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_up <- 
  setdiff(p_to_d_e3_.5_up[,1],m_to_d_e3_.5_up[,1])
p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_up <- 
  setdiff(p_to_d_e4_.5_up[,1],m_to_d_e4_.5_up[,1])
p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_up <- 
  setdiff(p_to_d_e5_.5_up[,1],m_to_d_e5_.5_up[,1])
p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_up <- 
  setdiff(p_to_d_e6_.5_up[,1],m_to_d_e6_.5_up[,1])



p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_down <- 
  setdiff(p_to_d_e2_.5_down[,1],m_to_d_e2_.5_down[,1])
p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_down <- 
  setdiff(p_to_d_e3_.5_down[,1],m_to_d_e3_.5_down[,1])
p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_down <- 
  setdiff(p_to_d_e4_.5_down[,1],m_to_d_e4_.5_down[,1])
p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_down <- 
  setdiff(p_to_d_e5_.5_down[,1],m_to_d_e5_.5_down[,1])
p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_down <- 
  setdiff(p_to_d_e6_.5_down[,1],m_to_d_e6_.5_down[,1])
#complementary for prion specific genes
m_to_d_.001_.5_AND_NOT_p_to_d_.001_.5 <-setdiff(m_to_d_.001_.5[,1], p_to_d_.001_.5[,1])
m_to_d_.001_1_AND_NOT_p_to_d_.001_1 <- setdiff(m_to_d_.001_1[,1], p_to_d_.001_1[,1])
m_to_d_.001_2_AND_NOT_p_to_d_.001_2 <- setdiff(m_to_d_.001_2[,1], p_to_d_.001_2[,1])
m_to_d_.001_3_AND_NOT_p_to_d_.001_3 <- setdiff(m_to_d_.001_3[,1], p_to_d_.001_3[,1])
m_to_d_.001_4_AND_NOT_p_to_d_.001_4 <- setdiff(m_to_d_.001_4[,1], p_to_d_.001_4[,1])
m_to_d_.001_5_AND_NOT_p_to_d_.001_5 <- setdiff(m_to_d_.001_5[,1], p_to_d_.001_5[,1])
m_to_d_.001_5.5_AND_NOT_p_to_d_.001_5.5 <- setdiff(m_to_d_.001_5.5[,1], p_to_d_.001_5.5[,1])
m_to_d_.001_6_AND_NOT_p_to_d_.001_6 <- setdiff(m_to_d_.001_6[,1], p_to_d_.001_6[,1])

m_to_d_.001_.5_up_AND_NOT_p_to_d_.001_.5_up <- 
  setdiff(m_to_d_.001_.5_up[,1], p_to_d_.001_.5_up[,1])

m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5 <- 
  setdiff(m_to_d_e5_.5[,1], p_to_d_e5_.5[,1])
m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5 <- 
  setdiff(m_to_d_e5_.5[,1],p_to_d_e5_.5[,1])
m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5 <- 
  setdiff(m_to_d_e2_.5[,1],p_to_d_e2_.5[,1])
m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5 <- 
  setdiff(m_to_d_e3_.5[,1],p_to_d_e3_.5[,1])
m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5 <- 
  setdiff(m_to_d_e4_.5[,1],p_to_d_e4_.5[,1])
m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5 <- 
  setdiff(m_to_d_e6_.5[,1],p_to_d_e6_.5[,1])
m_to_d_e7_.5_AND_NOT_p_to_d_e7_.5 <- 
  setdiff(m_to_d_e7_.5[,1],p_to_d_e7_.5[,1])
m_to_d_e8_.5_AND_NOT_p_to_d_e8_.5 <- 
  setdiff(m_to_d_e8_.5[,1],p_to_d_e8_.5[,1])
m_to_d_e9_.5_AND_NOT_p_to_d_e9_.5 <- 
  setdiff(m_to_d_e9_.5[,1],p_to_d_e9_.5[,1])


m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_up <- 
  setdiff(m_to_d_e2_.5_up[,1],p_to_d_e2_.5_up[,1])
m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_up <- 
  setdiff(m_to_d_e3_.5_up[,1],p_to_d_e3_.5_up[,1])
m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_up <- 
  setdiff(m_to_d_e4_.5_up[,1],p_to_d_e4_.5_up[,1])
m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_up <- 
  setdiff(m_to_d_e5_.5_up[,1],p_to_d_e5_.5_up[,1])
m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_up <- 
  setdiff(m_to_d_e6_.5_up[,1],p_to_d_e6_.5_up[,1])



m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_down <- 
  setdiff(m_to_d_e2_.5_down[,1],p_to_d_e2_.5_down[,1])
m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_down <- 
  setdiff(m_to_d_e3_.5_down[,1],p_to_d_e3_.5_down[,1])
m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_down <- 
  setdiff(m_to_d_e4_.5_down[,1],p_to_d_e4_.5_down[,1])
m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_down <- 
  setdiff(m_to_d_e5_.5_down[,1],p_to_d_e5_.5_down[,1])
m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_down <- 
  setdiff(m_to_d_e6_.5_down[,1],p_to_d_e6_.5_down[,1])

# plot
#png("./new ideas/Prion-specific genes_adjusting FC.png") 
abslogFC<-c(.5, 1, 2, 3, 4, 5, 5.5, 6)
plot (NA, NA, xlim=c(0,6.5), ylim=c(0, 600), xlab="abslogFC threshold", 
      ylab="number of specific genes")
lines (abslogFC, c(length(p_to_d_.001_.5_AND_NOT_m_to_d_.001_.5),
                   length(p_to_d_.001_1_AND_NOT_m_to_d_.001_1),
                   length(p_to_d_.001_2_AND_NOT_m_to_d_.001_2),
                   length(p_to_d_.001_3_AND_NOT_m_to_d_.001_3),
                   length(p_to_d_.001_4_AND_NOT_m_to_d_.001_4),
                   length(p_to_d_.001_5_AND_NOT_m_to_d_.001_5),
                   length(p_to_d_.001_5.5_AND_NOT_m_to_d_.001_5.5),
                   length(p_to_d_.001_6_AND_NOT_m_to_d_.001_6)), col="violetred4")

points (abslogFC, c(length(p_to_d_.001_.5_AND_NOT_m_to_d_.001_.5),
                   length(p_to_d_.001_1_AND_NOT_m_to_d_.001_1),
                   length(p_to_d_.001_2_AND_NOT_m_to_d_.001_2),
                   length(p_to_d_.001_3_AND_NOT_m_to_d_.001_3),
                   length(p_to_d_.001_4_AND_NOT_m_to_d_.001_4),
                   length(p_to_d_.001_5_AND_NOT_m_to_d_.001_5),
                   length(p_to_d_.001_5.5_AND_NOT_m_to_d_.001_5.5),
                   length(p_to_d_.001_6_AND_NOT_m_to_d_.001_6)), col="violetred4")

lines (abslogFC, c(length(m_to_d_.001_.5_AND_NOT_p_to_d_.001_.5),
                   length(m_to_d_.001_1_AND_NOT_p_to_d_.001_1),
                   length(m_to_d_.001_2_AND_NOT_p_to_d_.001_2),
                   length(m_to_d_.001_3_AND_NOT_p_to_d_.001_3),
                   length(m_to_d_.001_4_AND_NOT_p_to_d_.001_4),
                   length(m_to_d_.001_5_AND_NOT_p_to_d_.001_5),
                   length(m_to_d_.001_5.5_AND_NOT_p_to_d_.001_5.5),
                   length(m_to_d_.001_6_AND_NOT_p_to_d_.001_6)), col="blue")

points (abslogFC, c(length(m_to_d_.001_.5_AND_NOT_p_to_d_.001_.5),
                   length(m_to_d_.001_1_AND_NOT_p_to_d_.001_1),
                   length(m_to_d_.001_2_AND_NOT_p_to_d_.001_2),
                   length(m_to_d_.001_3_AND_NOT_p_to_d_.001_3),
                   length(m_to_d_.001_4_AND_NOT_p_to_d_.001_4),
                   length(m_to_d_.001_5_AND_NOT_p_to_d_.001_5),
                   length(m_to_d_.001_5.5_AND_NOT_p_to_d_.001_5.5),
                   length(m_to_d_.001_6_AND_NOT_p_to_d_.001_6)), col="blue")
legend(x="topright", legend=c("plus to del and not minus to del, p=.001",
                            "minus to del and not plus to del, p=.001"),
       fill=c("violetred4", "blue"))
#dev.off()


par(col=1)
#plot2

# plot

#png("./new ideas/Prion-specific genes: Adjusting p-value.png") 
ppval<-c(2, 3, 4, 5, 6, 7, 8, 9)
plot (NA, NA, xlim=c(1,9.5), ylim=c(0, 800), xlab="-lg(p-value threshold)", 
      ylab="number of specific genes, abslogFC>0.5" )
lines (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5),
                   length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5),
                   length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5),
                   length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5),
                   length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5),
                   length(p_to_d_e7_.5_AND_NOT_m_to_d_e7_.5),
                   length(p_to_d_e8_.5_AND_NOT_m_to_d_e8_.5),
                   length(p_to_d_e9_.5_AND_NOT_m_to_d_e9_.5)),
                   col="violetred4")

       points (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5),
                          length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5),
                          length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5),
                          length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5),
                          length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5),
                        length(p_to_d_e7_.5_AND_NOT_m_to_d_e7_.5),
                        length(p_to_d_e8_.5_AND_NOT_m_to_d_e8_.5),
                        length(p_to_d_e9_.5_AND_NOT_m_to_d_e9_.5)),
                          col="violetred4")

               
               
               
               lines (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5),
                                  length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5),
                                  length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5),
                                  length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5),
                                  length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5),
                      length(m_to_d_e7_.5_AND_NOT_p_to_d_e7_.5),
                      length(m_to_d_e8_.5_AND_NOT_p_to_d_e8_.5),
                      length(m_to_d_e9_.5_AND_NOT_p_to_d_e9_.5)),
                                  col="blue")
                      
              
                      
                      points (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5),
                                      length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5),
                                      length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5),
                                      length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5),
                                      length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5),
                              length(m_to_d_e7_.5_AND_NOT_p_to_d_e7_.5),
                              length(m_to_d_e8_.5_AND_NOT_p_to_d_e8_.5),
                                       length(m_to_d_e9_.5_AND_NOT_p_to_d_e9_.5)),
                                      col="blue")
                             
                             
               
lines(0:9, rep.int(0,10))

legend(x="topright", legend=c("plus to del and not minus to del",
                              "minus to del and not plus to del"),
       fill=c("violetred4", "blue"))
#dev.off()


a<-setdiff(p_to_d_.001_.5[,1], m_to_d_.001_.5[,1])
b<-setdiff(m_to_d_.001_.5[,1], p_to_d_.001_.5[,1])
c<-intersect(a, p_to_m_.001_.5[,1])
d<-intersect(b, p_to_m_.001_.5[,1])


           
           
#png("./new ideas/Prion-specific genes: Adjusting p-value_upregs.png") 
           ppval<-c(2, 3, 4, 5, 6)
           plot (NA, NA, xlim=c(1,6.5), ylim=c(0, 500), xlab="-lg(p-value threshold)", 
                 ylab="number of specific genes, abslogFC>0.5" )
           lines (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_up),
                           length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_up),
                           length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_up),
                           length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_up),
                           length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_up)),
                  col="violetred4")
           
           points (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_up),
                            length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_up),
                            length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_up),
                            length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_up),
                            length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_up)),
                   col="violetred4")
           
           
           
           
           lines (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_up),
                           length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_up),
                           length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_up),
                           length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_up),
                           length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_up)),
                  col="blue")
           
           
           
           points (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_up),
                            length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_up),
                            length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_up),
                            length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_up),
                            length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_up)),
                   col="blue")
           
           
           
           
           legend(x="topright", legend=c("plus to del and not minus to del",
                                         "minus to del and not plus to del"),
                  fill=c("violetred4", "blue"))
#           dev.off()
           
           
           
#upregs           
#png("./new ideas/Prion-specific genes: upregs.png")
plot(NA, NA, xlim=c(1,6), ylim=c(1,2000), ylab="number of genes", xlab="-lg(p-value)")
           lines (ppval, c(nrow(p_to_d_e2_.5_up),
                           nrow(p_to_d_e3_.5_up),
                           nrow(p_to_d_e4_.5_up),
                           nrow(p_to_d_e5_.5_up),
                           nrow(p_to_d_e6_.5_up)),
                  col="violetred4")
           
           points (ppval, c(nrow(p_to_d_e2_.5_up),
                            nrow(p_to_d_e3_.5_up),
                            nrow(p_to_d_e4_.5_up),
                            nrow(p_to_d_e5_.5_up),
                            nrow(p_to_d_e6_.5_up)),
                   col="violetred4")
           
           lines (ppval, c(nrow(m_to_d_e2_.5_up),
                           nrow(m_to_d_e3_.5_up),
                           nrow(m_to_d_e4_.5_up),
                           nrow(m_to_d_e5_.5_up),
                           nrow(m_to_d_e6_.5_up)),
                  col="blue")
           
           points (ppval, c(nrow(m_to_d_e2_.5_up),
                            nrow(m_to_d_e3_.5_up),
                            nrow(m_to_d_e4_.5_up),
                            nrow(m_to_d_e5_.5_up),
                            nrow(m_to_d_e6_.5_up)),
                   col="blue")

legend(x="topright", legend=c("plus to del_up",
                                         "minus to del_up"),
                  fill=c("violetred4", "blue"))
#dev.off()           
           
#downregs           
png("./new ideas/Prion-specific genes: downregs.png")
plot(NA, NA, xlim=c(1,6), ylim=c(1,2000), ylab="number of genes", xlab="-lg(p-value)")
           lines (ppval, c(nrow(p_to_d_e2_.5_down),
                           nrow(p_to_d_e3_.5_down),
                           nrow(p_to_d_e4_.5_down),
                           nrow(p_to_d_e5_.5_down),
                           nrow(p_to_d_e6_.5_down)),
                  col="violetred4")
           
           points (ppval, c(nrow(p_to_d_e2_.5_down),
                            nrow(p_to_d_e3_.5_down),
                            nrow(p_to_d_e4_.5_down),
                            nrow(p_to_d_e5_.5_down),
                            nrow(p_to_d_e6_.5_down)),
                   col="violetred4")
           
           lines (ppval, c(nrow(m_to_d_e2_.5_down),
                           nrow(m_to_d_e3_.5_down),
                           nrow(m_to_d_e4_.5_down),
                           nrow(m_to_d_e5_.5_down),
                           nrow(m_to_d_e6_.5_down)),
                  col="blue")
           
           points (ppval, c(nrow(m_to_d_e2_.5_down),
                            nrow(m_to_d_e3_.5_down),
                            nrow(m_to_d_e4_.5_down),
                            nrow(m_to_d_e5_.5_down),
                            nrow(m_to_d_e6_.5_down)),
                   col="blue")
           
           legend(x="topright", legend=c("plus to del_down",
                                         "minus to del_down"),
                  fill=c("violetred4", "blue"))
#dev.off()           
           
#setdiffs for downregs
           
#png("./new ideas/Prion-specific genes: Adjusting p-value_downregs.png") 
           ppval<-c(2, 3, 4, 5, 6)
           plot (NA, NA, xlim=c(1,6.5), ylim=c(0, 500), xlab="-lg(p-value threshold)", 
                 ylab="number of specific genes, abslogFC>0.5" )
           lines (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_down),
                           length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_down),
                           length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_down),
                           length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_down),
                           length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_down)),
                  col="violetred4")
           
           points (ppval, c(length(p_to_d_e2_.5_AND_NOT_m_to_d_e2_.5_down),
                            length(p_to_d_e3_.5_AND_NOT_m_to_d_e3_.5_down),
                            length(p_to_d_e4_.5_AND_NOT_m_to_d_e4_.5_down),
                            length(p_to_d_e5_.5_AND_NOT_m_to_d_e5_.5_down),
                            length(p_to_d_e6_.5_AND_NOT_m_to_d_e6_.5_down)),
                   col="violetred4")
           
           
           
           
           lines (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_down),
                           length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_down),
                           length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_down),
                           length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_down),
                           length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_down)),
                  col="blue")
           
           
           
           points (ppval, c(length(m_to_d_e2_.5_AND_NOT_p_to_d_e2_.5_down),
                            length(m_to_d_e3_.5_AND_NOT_p_to_d_e3_.5_down),
                            length(m_to_d_e4_.5_AND_NOT_p_to_d_e4_.5_down),
                            length(m_to_d_e5_.5_AND_NOT_p_to_d_e5_.5_down),
                            length(m_to_d_e6_.5_AND_NOT_p_to_d_e6_.5_down)),
                   col="blue")
           
           
           
           
           legend(x="topright", legend=c("plus to del and not minus to del",
                                         "minus to del and not plus to del"),
                  fill=c("violetred4", "blue"))
#           dev.off()



#sampling
p_to_d_.001_.5_up_sample<-sample(p_to_d_.001_.5_up[,1], 300)
write.table(p_to_d_.001_.5_up_sample, "./new ideas/p_to_d_.001_.5_sample_300_2.csv")

m_to_d_.001_.5_up_sample<-sample(m_to_d_.001_.5_up[,1], 300)
write.table(m_to_d_.001_.5_up_sample, "./new ideas/m_to_d_.001_.5_sample_300.csv")

genestable <- read.csv("./SGD retrieves/genelist2.csv", sep=",")
allyeastgenes <- genestable$systematic_name
allyeastgenes2<-allyeastgenes[grep('Y.*', allyeastgenes)]
genesSample<-sample(allyeastgenes2, 300)
#write.table(genesSample, "./new ideas/genesSample.csv")


#for LN
length(intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_up[,1])
length(intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_up[,1]))
length(intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_down[,1]))
length(intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_down[,1]))
       
length(setdiff(p_to_m_.001_.5[,1], d_to_m_.001_.5[,1]))
#write.table(setdiff(p_to_m_.001_.5[,1], d_to_m_.001_.5[,1]), "df.csv")
       
write.table(setdiff(p_to_m_.001_.5_up[,1], d_to_m_.001_.5[,1]), "df.csv")
write.table(setdiff(p_to_m_.001_.5_down[,1], d_to_m_.001_.5[,1]), "df.csv")
       
p_m_AND_d_m_same <- c(intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_up[,1]),
                      intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_down[,1]))
length(p_m_AND_d_m_same)
length(intersect(p_m_AND_d_m_same, m_to_d_.001_.5_AND_NOT_p_to_d_.001_.5))

       
p_m_AND_d_m_diff<-c(intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_down[,1]),
                   intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_up[,1]))
length(p_m_AND_d_m_diff)
length(intersect(p_m_AND_d_m_diff, p_to_d_.001_.5_AND_NOT_m_to_d_.001_.5))
       
p_m_AND_d_m_same_up <- c(intersect(p_to_m_.001_.5_up[,1], d_to_m_.001_.5_up[,1]))
p_m_AND_d_m_same_down <- c(intersect(p_to_m_.001_.5_down[,1], d_to_m_.001_.5_down[,1]))
       

p_m_AND_d_m_same <- c(intersect(p_to_m_.001_1_up[,1], d_to_m_.001_1_up[,1]),
                             intersect(p_to_m_.001_1_down[,1], d_to_m_.001_1_down[,1]))

       
p_m_AND_d_m_same <- c(intersect(subset([,1], d_to_m_.001_1_up[,1]),
                             intersect(p_to_m_.001_1_down[,1], d_to_m_.001_1_down[,1]))
       
                      
                      
p_m_AND_d_m_same_2<-c(intersect(cutoff2(plus_minus, .001, 2, "up")[,1],
                      cutoff2(del_minus, .001, 2, "up")[,1]), 
                      intersect(cutoff2(plus_minus, .001, 2, "down")[,1],
                      cutoff2(del_minus, .001, 2, "down")[,1]))
                      
                      
                      
                      
logFC<-seq(0, 4, 0.1)             
plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(plus_del_cut)), xlab="abslogFC", ylab="number of genes")  
p_m_num <-c()
for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
lines (logFC, p_m_num, col="violetred")
#points(logFC, p_m_num, col="violetred")
p_d_num <-c()
for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
lines (logFC, p_d_num, col="dark violet")
#points(logFC, p_d_num, col="dark violet")
m_d_num <-c()
for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
lines (logFC, m_d_num, col="blue4")
#points(logFC, m_d_num, col="blue4")
legend(x="topright", legend=c("p/d", "m/d", "p/m"), fill=c("dark violet", "blue4", "violetred"))


plot(NA, NA, xlim=c(2,max(pval)), ylim=c(1,nrow(cutoff2(plus_del, p.value=.01))),
     xlab="-lg(p-value)", ylab="number of genes")

pval<-c(2, 3, 4, 5, 6)
p_m_num2<-c()
for (i in 1:length(pval)) {p_m_num2[i]<-nrow(cutoff2(plus_minus, p.value=10^(0-pval[i]))) }
lines (pval, p_m_num2, col="violetred")
points(pval, p_m_num2, col="violetred")
                      
p_d_num2<-c()
for (i in 1:length(pval)) {p_d_num2[i]<-nrow(cutoff2(plus_del, p.value=10^(0-pval[i]))) }
lines (pval, p_d_num2, col="dark violet")
points(pval, p_d_num2, col="dark violet")
                      
                      
m_d_num2<-c()
for (i in 1:length(pval)) {m_d_num2[i]<-nrow(cutoff2(minus_del, p.value=10^(0-pval[i]))) }
lines (pval, m_d_num2, col="blue4")
points(pval, m_d_num2, col="blue4")
                      
                      
legend(x="topright", legend=c("p/d", "m/d", "p/m"), fill=c("dark violet", "blue4", "violetred"))
                      

                      
#logFC for Tanya
logFC<-seq(0, 4, 0.1)             
plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(plus_del_cut)), xlab="log2(кратность изменения)", ylab="Число генов", main="Зависимость числа генов с изменённой экспрессией \nот выбранного порогового значения кратного изменения экспрессии", cex.main=0.9)  
p_m_num <-c()
for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
lines (logFC, p_m_num, col="violetred", lwd=3)
#points(logFC, p_m_num, col="violetred")
p_d_num <-c()
for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
lines (logFC, p_d_num, col="dark violet", lwd=3)
#points(logFC, p_d_num, col="dark violet")
m_d_num <-c()
for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
lines (logFC, m_d_num, col="blue4", lwd=3)
#points(logFC, m_d_num, col="blue4")
legend(x="topright", legend=c("[ISP+]/sfp1Δ", "[isp-]/sfp1Δ", "[ISP+]/[isp-]"), fill=c("dark violet", "blue4", "violetred"))

                      
#for the poster

intersect(cutoff2(plus_del, p.value=.001, abslogFC=.5, up_or_down="up")[,1],
                      cutoff2(minus_del, p.value=.001, abslogFC=.5, up_or_down="all")[,1])
                      
                      
write.table(intersect(cutoff2(minus_del)[,1], 
          cutoff2(plus_minus)[,1]), "md_pm.csv")