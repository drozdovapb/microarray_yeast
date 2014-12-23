#extracting lists


write.table(setdiff(p_to_d_.001_.5[,1], m_to_d_.001_.5[,1]), "plus to del AND NOT minus to del_.001_.5.xls")
 write.table(setdiff(p_to_d_.001_.5[,1], m_to_d_.001_.5[,1]), "plus to del AND NOT minus to del_.001_.5.csv")
 write.table(setdiff(m_to_d_.001_.5[,1], p_to_d_.001_.5[,1]), "minus to del AND NOT plus to del_.001_.5.csv")
 write.table(setdiff(subset(m_to_d_.001_.5, logFC>0)[,1], subset(p_to_d_.001_.5, logFC>0)[,1]), "minus to del AND NOT plus to del_.001_.5_upregs.csv")
 setwd("./new ideas")
 write.table(setdiff(subset(m_to_d_.001_.5, logFC<0)[,1], subset(p_to_d_.001_.5, logFC<0)[,1]), "minus to del AND NOT plus to del_.001_.5_downregs.csv")
 write.table(setdiff(subset(p_to_d_.001_.5, logFC<0)[,1], subset(m_to_d_.001_.5, logFC<0)[,1]), "plus to del AND NOT minus to del_.001_.5_downregs.csv")
 write.table(setdiff(subset(p_to_d_.001_.5, logFC>0)[,1], subset(m_to_d_.001_.5, logFC>0)[,1]), "plus to del AND NOT minus to del_.001_.5_upregs.csv")
 write.table(intersect(subset(p_to_m_.001_.5, logFC>0)[,1], subset(d_to_m_.001_.5, logFC<0)[,1]), "p_to_m_up_and_d_to_m_down_.001_.5.csv")
 write.table(intersect(subset(p_to_m_.001_.5, logFC<0)[,1], subset(d_to_m_.001_.5, logFC>0)[,1]), "p_to_m_down_and_d_to_m_up_.001_.5.csv")
 write.table(intersect(subset(p_to_m_.001_.5, logFC>0)[,1], subset(d_to_m_.001_.5, logFC>0)[,1]), "p_to_m_up_and_d_to_m_up_.001_.5.csv")
 write.table(intersect(subset(p_to_m_.001_.5, logFC<0)[,1], subset(d_to_m_.001_.5, logFC<0)[,1]), "p_to_m_down_and_d_to_m_down_.001_.5.csv")