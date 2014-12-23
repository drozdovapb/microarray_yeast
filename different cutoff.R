  setwd("/media/sda6/Polina/Work/Ingenuity")
  
  volcanoplot(fit2, coef=1, col="blue4")
  
  
  #plotting logFC
  logFC<-seq(0, 4, 0.1)             
  plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(cutoff2(plus_del))), xlab="abslogFC", ylab="number of genes")
  p_m_num <-c()
  for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_m_num, col="violetred")
  #points(logFC, p_m_num, col="violetred", cex=.5, pch=19)
  p_d_num <-c()
  for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_d_num, col="dark violet")
  #points(logFC, p_d_num, col="dark violet", cex=.5, pch=19)
  m_d_num <-c()
  for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
  lines (lwd=2, logFC, m_d_num, col="blue4")
  #points(logFC, m_d_num, col="blue4", cex=.5, pch=19)
  legend(x="topright", legend=c("p/d", "m/d", "p/m"), fill=c("dark violet", "blue4", "violetred"))
  
  pval<-seq(2, 6, 1)
  plot(NA, NA, xlim=c(2,max(pval)), ylim=c(1,nrow(cutoff2(plus_del, p.value=.01))),
       xlab="-lg(p-value)", ylab="number of genes")
  
  p_m_num2<-c()
  for (i in 1:length(pval)) {p_m_num2[i]<-nrow(cutoff2(plus_minus, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_m_num2, col="violetred")
  #points(pval, p_m_num2, col="violetred")
  p_d_num2<-c()
  for (i in 1:length(pval)) {p_d_num2[i]<-nrow(cutoff2(plus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_d_num2, col="dark violet")
  #points(pval, p_d_num2, col="dark violet")
  m_d_num2<-c()
  for (i in 1:length(pval)) {m_d_num2[i]<-nrow(cutoff2(minus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, m_d_num2, col="blue4")
  #points(pval, m_d_num2, col="blue4")
  legend(x="topright", legend=c("p/d", "m/d", "p/m"), fill=c("dark violet", "blue4", "violetred"))
  
  pvalSel<-c(.05, .01, .001)
  abslogFCSel<-c(.5, 1, 2)
  for (i in 1:length(pvalSel)) {
    for (j in 1:length(abslogFCSel)) {
  #    for (k in 1:3) {
      write.csv(cutoff2(mytable=plus_minus, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="up"),
                  file=paste("./processed_tables_R/different cutoffs/plus_minus", 
                             abslogFCSel[j], pvalSel[i], "up", sep="_")) 
      write.csv(cutoff2(mytable=plus_minus, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="down"),
                file=paste("./processed_tables_R/different cutoffs/plus_minus", 
                           abslogFCSel[j], pvalSel[i], "down", sep="_"))     }
  }
  
  
  
  pvalSel<-c(.05, .01, .001)
  abslogFCSel<-c(.5, 1, 2)
  for (i in 1:length(pvalSel)) {
    for (j in 1:length(abslogFCSel)) {
      #    for (k in 1:3) {
      write.csv(cutoff2(mytable=minus_del, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="up"),
                file=paste("./processed_tables_R/different cutoffs/minus_del", 
                           abslogFCSel[j], pvalSel[i], "up", sep="_")) 
      write.csv(cutoff2(mytable=minus_del, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="down"),
                file=paste("./processed_tables_R/different cutoffs/minus_del", 
                           abslogFCSel[j], pvalSel[i], "down", sep="_"))     }
  }
  
  
  
  pvalSel<-c(.05, .01, .001)
  abslogFCSel<-c(.5, 1, 2)
  for (i in 1:length(pvalSel)) {
    for (j in 1:length(abslogFCSel)) {
      #    for (k in 1:3) {
      write.csv(cutoff2(mytable=plus_del, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="up"),
                file=paste("./processed_tables_R/different cutoffs/plus_del", 
                           abslogFCSel[j], pvalSel[i], "up", sep="_")) 
      write.csv(cutoff2(mytable=plus_del, abslogFC=abslogFCSel[j], p.value=pvalSel[i], up_or_down="down"),
                file=paste("./processed_tables_R/different cutoffs/plus_del", 
                           abslogFCSel[j], pvalSel[i], "down", sep="_"))     }
  }
  
  nrow(cutoff2(minus_del, abslogFC=.5, p.value=.001, up_or_down="up"))
  
  
  
  #plots in Russian
  #logFC in Russian
  
  #plotting logFC
  logFC<-seq(0, 4, 0.1)             
  plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(cutoff2(plus_del))), xlab="abslogFC", ylab="Число генов")
  p_m_num <-c()
  for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_m_num, col="violetred")
  #points(logFC, p_m_num, col="violetred", cex=.5, pch=19)
  p_d_num <-c()
  for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_d_num, col="dark violet")
  #points(logFC, p_d_num, col="dark violet", cex=.5, pch=19)
  m_d_num <-c()
  for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
  lines (lwd=2, logFC, m_d_num, col="blue4")
  #points(logFC, m_d_num, col="blue4", cex=.5, pch=19)
  legend(x="topright", legend=c("p/d", "m/d", "p/m"), fill=c("dark violet", "blue4", "violetred"))
  
#p-value Russian
png("../figures_R/different cutoffs/p-value Rus.png", width=700, height=400, res=96)
#svg("../figures_R/different cutoffs/p-value Rus.svg")
  pval<-seq(0, 6, .01)
 # pval2<-seq(0,6,.5)
#for presentation: change cex. Cex 1.5 would be okay.
  #par(ps=12, cex=1, cex.main=1)
  plot(NA, NA, xlim=c(min(pval),max(pval)), ylim=c(1,nrow(cutoff2(plus_del, p.value=10^(-min(pval))))),
       xlab="-lg(p-value)", ylab="Число генов с изменённой экспрессией")
    p_m_num2<-c()
  for (i in 1:length(pval)) {p_m_num2[i]<-nrow(cutoff2(plus_minus, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_m_num2, col="violetred")
#points(pval, p_m_num2, col="violetred")
  p_d_num2<-c()
  for (i in 1:length(pval)) {p_d_num2[i]<-nrow(cutoff2(plus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_d_num2, col="dark violet")
#points(pval, p_d_num2, col="dark violet")
  m_d_num2<-c()
  for (i in 1:length(pval)) {m_d_num2[i]<-nrow(cutoff2(minus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, m_d_num2, col="blue4")
#points(pval, m_d_num2, col="blue4")
  legend(x="topright", legend=c("[ISP+]/sfp1Δ", "[isp-]/sfp1Δ", "[ISP+]/[isp-]"), fill=c("dark violet", "blue4", "violetred"))
abline (v=c(-log(.05, 10), -log(.01, 10), -log(.001, 10)))
  dev.off()
  
  #plotting logFC — Russian
png("../figures_R/different cutoffs/abslogFC Rus.png", width=700, height=400, res=96)
#svg("../figures_R/different cutoffs/abslogFC Rus.svg")
  logFC<-seq(0, 4, 0.1)             
  plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(cutoff2(plus_del))), xlab="abslogFC", ylab="Число генов с изменённой экспрессией")
  p_m_num <-c()
  for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_m_num, col="violetred")
  p_d_num <-c()
  for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_d_num, col="dark violet")
  m_d_num <-c()
  for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
  lines (lwd=2, logFC, m_d_num, col="blue4")
  legend(x="topright", legend=c("[ISP+]/sfp1Δ", "[isp-]/sfp1Δ", "[ISP+]/[isp-]"), fill=c("dark violet", "blue4", "violetred"))
  abline (v=c(0.5, 1))
  dev.off()
  
  
  
  
  #plotting logFC — English
  png("../figures_R/different cutoffs/abslogFC Eng.png", width=700, height=400, res=96)
  logFC<-seq(0, 4, 0.1)             
  plot(NA, NA, xlim=c(0, max(logFC)), ylim=c(0, nrow(cutoff2(plus_del))), xlab="abslogFC", ylab="# of DE genes")
  p_m_num <-c()
  for (i in 1:length(logFC)) {p_m_num[i] <- nrow(cutoff2(mytable=plus_minus, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_m_num, col="violetred")
  p_d_num <-c()
  for (i in 1:length(logFC)) {p_d_num[i] <- nrow(cutoff2(mytable=plus_del, abslogFC=logFC[i]))}
  lines (lwd=3, logFC, p_d_num, col="dark violet")
  m_d_num <-c()
  for (i in 1:length(logFC)) {m_d_num[i] <- nrow(cutoff2(mytable=minus_del, abslogFC=logFC[i]))}
  lines (lwd=2, logFC, m_d_num, col="blue4")
  legend(x="topright", legend=c("[ISP+]/sfp1Δ", "[isp-]/sfp1Δ", "[ISP+]/[isp-]"), fill=c("dark violet", "blue4", "violetred"))
  abline (v=c(0.5, 1))
  dev.off()
  
  #p-value English
  png("../figures_R/different cutoffs/p-value Eng.png", width=700, height=400, res=96)
  pval<-seq(0, 6, .01)
  plot(NA, NA, xlim=c(min(pval),max(pval)), ylim=c(1,nrow(cutoff2(plus_del, p.value=10^(-min(pval))))),
       xlab="-lg(p-value)", ylab="# of DE genes")
  p_m_num2<-c()
  for (i in 1:length(pval)) {p_m_num2[i]<-nrow(cutoff2(plus_minus, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_m_num2, col="violetred")
  #points(pval, p_m_num2, col="violetred")
  p_d_num2<-c()
  for (i in 1:length(pval)) {p_d_num2[i]<-nrow(cutoff2(plus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, p_d_num2, col="dark violet")
  #points(pval, p_d_num2, col="dark violet")
  m_d_num2<-c()
  for (i in 1:length(pval)) {m_d_num2[i]<-nrow(cutoff2(minus_del, p.value=10^(0-pval[i]))) }
  lines (lwd=3, pval, m_d_num2, col="blue4")
  #points(pval, m_d_num2, col="blue4")
  legend(x="topright", legend=c("[ISP+]/sfp1Δ", "[isp-]/sfp1Δ", "[ISP+]/[isp-]"), fill=c("dark violet", "blue4", "violetred"))
  abline (v=c(-log(.05, 10), -log(.01, 10), -log(.001, 10)))
  dev.off()
  