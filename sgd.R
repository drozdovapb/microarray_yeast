# for linux 
setwd("/media/sda6/Polina/Work/Ingenuity")

#for windows
#setwd("D:/Polina/Work/Ingenuity")

sfp1upstream<-scan(file="./SGD retrieves/list_aaaawtttt upstream.csv", what="character")


a4wt4<-function (var) {
  var2<-intersect(var[,1],sfp1upstream)
  write.table(var2, file=paste("./SGD retrieves/Intersect/",
                               deparse(substitute(var)), "_has_a4wt4", ".csv"))
  return (var2)
}

m_d_a4wt4<-a4wt4(minus_del_cut)
p_d_a4wt4<-a4wt4(plus_del_cut)
p_m_a4wt4<-a4wt4(plus_minus_cut)

a4wt4(minus_del_cut_up)
a4wt4(minus_del_cut_down)
a4wt4(plus_del_cut_up)
a4wt4(plus_del_cut_down)
a4wt4(plus_minus_cut_up)
a4wt4(plus_minus_cut_down)

numbers<-rbind(a4wt4(minus_del_cut),
                      a4wt4(plus_del_cut),
                      a4wt4(plus_minus_cut),
                      a4wt4(minus_del_cut_up),
                      a4wt4(minus_del_cut_down),
                      a4wt4(plus_del_cut_up),
                      a4wt4(plus_del_cut_down),
                      a4wt4(plus_minus_cut_up),
                      a4wt4(plus_minus_cut_down))                

names<-rbind(deparse(substitute(minus_del)),
                  deparse(substitute(plus_del)),
                  deparse(substitute(plus_minus)),
                  deparse(substitute(minus_del_up)),
                  deparse(substitute(minus_del_down)),
                  deparse(substitute(plus_del_up)),
                 deparse(substitute(plus_del_down)),
                 deparse(substitute(plus_minus_up)),
                deparse(substitute(plus_minus_down)))                

mytable<-cbind(names, numbers)

write.table(mytable, "./SGD retrieves/Intersect/genes with a4wt4.csv",
         col.names=c("comparison", "number of genes with a4wt4"))

write.table(intersect(plus_minus$ID, sfp1upstream), 
            file="plus_minus_sfp1siteupstream.csv")

source("http://bioconductor.org/biocLite.R")
biocLite("yeast.db0")

getannotation <- function (comparison) {

  output<-data.frame()

  
  
for (i in (1:length(sfp1upstream))) {
  output<-rbind(output,subset (comparison, comparison$ID==sfp1upstream[i]))
     
} 
  write.table(output, file=paste(deparse(substitute
  return (output)
}
