# for linux 
setwd("/media/drozdovapb/441DB6652D7ED741/Work//Current work/Projects/bioinfo/")

#for windows
#setwd("D:/Polina/Work/Ingenuity")

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

#don't forget to execute functions.R

contrast.matrix8<-makeContrasts(plus-minus, del-minus, plus-del, levels=design)
fit8<-contrasts.fit(fit, contrast.matrix8)
fit8<-eBayes(fit8)
res8<-decideTests(fit8, p.value=.001, lfc=.5)
venn8<-vennCounts(res8)
vennDiagram(venn8)
venn8up<-vennCounts(res8, include="up")
par(col="green4")
vennDiagram(venn8up)

venn8down<-vennCounts(res8, include="down")
par(col="green4")
vennDiagram(venn8down)

venn8all<-vennCounts(res8, include="both")
par(col="blue")
vennDiagram(venn8all)

#venn diagrams for 3 usual comparisons

#png("3_triple_diagrams.png")
res<-decideTests(fit2, p.value=0.001, lfc=0.5)
allvenn <- vennCounts(res)

#par(mfrow=c(3,3), height=1600, width=5800)

png("all.png")
par(col="blue4")
vennDiagram(allvenn)
dev.off()

png("up.png")
par(col="red3")
upvenn <-vennCounts(res, include="up")
vennDiagram(upvenn)
dev.off()

png("down.png")
par(col="green4")
downvenn <-vennCounts(res, include="down")
vennDiagram(downvenn)
dev.off()

#everything above taken from limma user's guide


#venn diagrams for inverted plus-minus

contrast.matrix2 <-makeContrasts(minus-plus, minus-del, plus-del, levels=design)
fit3 <- contrasts.fit(fit, contrast.matrix2)
fit3 <-eBayes(fit3)

res2<-decideTests(fit3, p.value=0.001, lfc=0.5)
allvenninv <- vennCounts(res2)
png("inv_all.png")
par(col="blue4")
vennDiagram(allvenninv)
dev.off()

png("inv_up.png")
upvenninv <-vennCounts(res2, include="up")
par(col="red3")
vennDiagram(upvenninv)
dev.off()

png("inv_down.png")
downvenninv <-vennCounts(res2, include="down")
par(col="green4")
vennDiagram(downvenninv)
dev.off()

#install.packages("VennDiagram")
library(VennDiagram)

#venn.diagram(allvenn) returns error: Incorrect number of elements

dev.off()
draw.triple.venn(area1=332, area2=1773, area3=1942,
                 n12=213, n23=1385, n13=175, n123=92, 
                 category=c("[ISP+]/[isp-]","[isp-]/sfp1Δ", "[ISP+]/sfp1Δ"),
                 scaled=TRUE,
                 col=c("blue4", "blue4", "blue4"),
                 fill=c("blue4", "blue4", "blue4"),
                 alpha=c(332/1942, 1773/1942, 1),
                 cat.fontface=3,
                 cat.cex=2,
                 cex=1.5,
                 cat.just=list(c(0.5,-1), c(1,-1), c(0,1)),
                 mar=c(0.1,0.1,0.1,0.1))
              

#venn diagrams for two comparisons only

contrast.matrix <-makeContrasts(plus-minus, del-minus, levels=design)
fit4 <- contrasts.fit(fit, contrast.matrix)
fit4 <-eBayes(fit4)

res4<-decideTests(fit4, p.value=0.001, lfc=0.5)
allvenn4 <- vennCounts(res4)
png("2_comp_all.png")
par(col="blue4")
vennDiagram(allvenn4)
dev.off()

png("2_comp_up.png")
upvenn4 <-vennCounts(res4, include="up")
par(col="red3")
vennDiagram(upvenn4)
dev.off()

png("2_comp_down.png")
par(col="green4")
downvenn4 <-vennCounts(res4, include="down")
vennDiagram(downvenn4)
dev.off()


#venn diagrams for two comparisons only

contrast.matrix <-makeContrasts(plus-minus, minus-del, levels=design)
fit5 <- contrasts.fit(fit, contrast.matrix)
fit5 <-eBayes(fit5)

res5<-decideTests(fit5, p.value=0.001, lfc=0.5)
allvenn5 <- vennCounts(res5)
png("2comp_inv_all.png")
par(col="blue4")
vennDiagram(allvenn5)
dev.off()

png("2comp_inv_up.png")
upvenn5 <-vennCounts(res5, include="up")
par(col="red3")
vennDiagram(upvenn5)
dev.off()

png("2comp_inv_down.png")
downvenn5 <-vennCounts(res5, include="down")
par(col="green4")
vennDiagram(downvenn5)
dev.off()


#two comparisons for the sake of my mind

contrast.matrix <-makeContrasts(plus-minus, del-minus, levels=design)
fit5 <- contrasts.fit(fit, contrast.matrix)
fit5 <-eBayes(fit5)
res5<-decideTests(fit5, p.value=0.001, lfc=0.5)
venninv <- vennCounts(res5)
#png("2venninv.png")
par(col="blue4")
vennDiagram(venninv)
#dev.off()

#ups

venninvup <- vennCounts(res5, include="up")
#png("2venninvup.png")
par(col="red3")
vennDiagram(venninvup)
#dev.off()

#downs
venninvdown <- vennCounts(res5, include="down")
#png("2venninvdown.png")
par(col="green4")
vennDiagram(venninvdown)
#dev.off()

#to deletant

contrast.matrix <-makeContrasts(plus-del, minus-del, levels=design)
fit8 <- contrasts.fit(fit, contrast.matrix)
fit8 <-eBayes(fit8)
res8<-decideTests(fit8, p.value=0.001, lfc=0.5)
todel <- vennCounts(res8)


png("todelall.png", bg="transparent")
par(col="blue4")
vennDiagram(todel)
dev.off()

todel <- vennCounts(res8, include="down")

png("todeldown.png")
par(col="green4")
vennDiagram(todel)
dev.off()

todel <- vennCounts(res8, include="up")

png("todelup.png")
par(col="red3")
vennDiagram(todel)
dev.off()


#Venn Diagram for the poster

png("p_d_m_d.png", bg="transparent")
draw.pairwise.venn(area1=1942, area2=1773, cross.area=1385, 
                 scaled=TRUE,
                 col=c("blue4", "blue4"),
                 fill=c("dark violet","cyan"),
                 alpha=c(.8, .5),
                 cat.fontface=3,
                 cat.cex=2,
                 cex=1.5,
                 mar=rep(.025,4),
                   category=c("ISP+/Δsfp1", "isp-/Δsfp1"),
                   cat.just=list(c(0.37,-1.5), c(0.57,-1)))
dev.off()



png("p_d_up_m_d_up.png", bg="transparent")
draw.pairwise.venn(area1=1043, area2=986, cross.area=721, 
                   scaled=TRUE,
                   col=c("red4", "red4"),
                   fill=c("dark violet","cyan"),
                   alpha=c(.8, .5),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("ISP+/Δsfp1", "isp-/Δsfp1"),
                   cat.just=list(c(0.37,-1.5), c(0.57,-1)),
                   cat.col="red3")
dev.off()

png("p_d_down_m_d_down.png", bg="transparent")
draw.pairwise.venn(area1=899, area2=787, cross.area=659, 
                   scaled=TRUE,
                   col=c("green4", "green4"),
                   fill=c("dark violet","cyan"),
                   alpha=c(.8, .5),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("ISP+/Δsfp1", "isp-/Δsfp1"),
                   cat.just=list(c(0.37,-1.5), c(0.57,-1)),
                   cat.col="green4")
dev.off()


#m_to_d and p_to_m

contrast.matrix <-makeContrasts(minus-del, plus-minus, levels=design)
fit11 <- contrasts.fit(fit, contrast.matrix)
fit11 <-eBayes(fit11)
res11<-decideTests(fit11, p.value=0.001, lfc=0.5)
m_to_d_p_to_m <- vennCounts(res11)
#png("todelall.png", bg="transparent")
par(col="blue4")
vennDiagram(m_to_d_p_to_m)
#dev.off()

png("m_d_p_m.png", bg="transparent")
draw.pairwise.venn(area1=1773, area2=332, cross.area=213, 
                   scaled=TRUE,
                   col=c("blue4", "blue4"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("isp-/Δsfp1", "ISP+/isp-"),
                   cat.just=list(c(0.37,-1.5), c(0.8,-1)),
                   cat.col="blue4")
dev.off()


#down

m_to_d_p_to_m_down <- vennCounts(res11, include="down")
#png("m_to_d_p_to_m_down.png", bg="transparent")
par(col="green4")
vennDiagram(m_to_d_p_to_m_down)
#dev.off()

#png("p_m_down_m_d_down.png", bg="transparent")
svg("pm_down_dm_up.svg", bg="transparent")
draw.pairwise.venn(area1=787, area2=258, cross.area=33, 
                   scaled=TRUE,
                   col=c("green4", "green4"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("isp-/Δsfp1", "ISP+/isp-"),
                   cat.just=list(c(0.37,-1.5), c(0.8,-1)),
                   cat.col="green4")
dev.off()


#up

m_to_d_p_to_m_up <- vennCounts(res11, include="up")
#png("m_to_d_p_to_m_up.png", bg="transparent")
par(col="red4")
vennDiagram(m_to_d_p_to_m_up)
#dev.off()


#png("m_d_up_p_m_up.png", bg="transparent")
svg("pm_up_dm_down.svg", bg="transparent")
draw.pairwise.venn(area1=986, area2=74, cross.area=19, 
                   scaled=TRUE,
                   col=c("red3", "red3"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("isp-/Δsfp1", "ISP+/isp-"),
                   cat.just=list(c(0.37,-1.5), c(0.8,-1)),
                   cat.col="red3")
dev.off()


contrast.matrix <-makeContrasts(del-minus, plus-minus, levels=design)
fit11 <- contrasts.fit(fit, contrast.matrix)
fit11 <-eBayes(fit11)
res11<-decideTests(fit11, p.value=0.001, lfc=0.5)
m_to_d_p_to_m <- vennCounts(res11)
#png("todelall.png", bg="transparent")
par(col="blue4")
vennDiagram(m_to_d_p_to_m)
#dev.off()

svg("./figures_R/d_m_p_m.svg")
draw.pairwise.venn(area1=332, area2=1773, cross.area=213, 
                   scaled=TRUE,
                   col=c("blue4", "blue4"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=3,
                   cex=3,
                   mar=rep(.025,4),
                   category=c("Δsfp1/isp-", "ISP+/isp-"),
                   cat.just=list(c(0.2,-1.5), c(0.8,-2)),
                   cat.col="blue4")
dev.off()

svg("./figures_R/d_m_p_m_down.svg")
draw.pairwise.venn(area1=258, area2=986, cross.area=141, 
                   scaled=TRUE,
                   col=c("green4", "green4"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=3,
                   cex=3,
                   mar=rep(.025,4),
                   category=c("Δsfp1/isp-", "ISP+/isp-"),
                   cat.just=list(c(0.2,-1.5), c(0.8,-2)),
                   cat.col="green4")
dev.off()

svg("./figures_R/d_m_p_m_up.svg")
draw.pairwise.venn(area1=74, area2=787, cross.area=20, 
                   scaled=TRUE,
                   col=c("red4", "red4"),
                   fill=c("cyan","dark violet"),
                   alpha=c(.5, .8),
                   cat.fontface=3,
                   cat.cex=3,
                   cex=3,
                   mar=rep(.025,4),
                   category=c("Δsfp1/isp-", "ISP+/isp-"),
                   cat.just=list(c(0.2,-1.5), c(0.8,-2)),
                   cat.col="red4")
dev.off()



c(m_to_d_.001_.5[,1], p_to_m_.001_.5[,1])
length(intersect(c(m_to_d_.001_.5[,1], p_to_m_.001_.5[,1]), p_to_d_.001_.5[,1]))

png("sizes of gene lists.png")
draw.pairwise.venn(area1=2105, area2=1942, cross.area=1468, 
                   scaled=TRUE,
                   col=c("blue4", "blue4"),
                   fill=c("violet","dark violet"),
                   alpha=c(.5, .5),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c("(isp-/Δspf1) + (ISP+/isp-)", "ISP+/Δsfp1"),
                   cat.just=list(c(.1,-3.5), c(.8,-2.5)),
                   cat.col="blue4")
dev.off()



png("p_m_AND_d_m_same_m_d_AND_NOT_p_to_d_1.png")
draw.pairwise.venn(area1=length(p_m_AND_d_m_same), area2=length(m_to_d_.001_1_AND_NOT_p_to_d_.001_1), 
                   cross.area=length(intersect(p_m_AND_d_m_same, m_to_d_.001_1_AND_NOT_p_to_d_.001_1)), 
                   scaled=TRUE,
                   col=c("blue4", "blue4"),
                   fill=c("violet","dark violet"),
                   alpha=c(.5, .5),
                   cat.fontface=3,
                   cat.cex=2,
                   cex=1.5,
                   mar=rep(.025,4),
                   category=c( "m_d_AND_NOT_p_d", "p_m_AND_d_m_same"),
                   cat.just=list(c(-.1, 0),c(1, 0)),
                   cat.col="blue4")
dev.off()
dev.new()





contrast.matrix <-makeContrasts(del-minus, plus-minus, levels=design)
fit <- eBayes(contrasts.fit(fit, contrast.matrix))
res<-decideTests(fit, p.value=0.001, lfc=0.5)
par(col="blue4")
vennDiagram(vennCounts(res))

res<-decideTests(fit, p.value=0.001, lfc=1)
par(col="blue4")
vennDiagram(vennCounts(res))

res<-decideTests(fit, p.value=0.001, lfc=2)
par(col="blue4")
vennDiagram(vennCounts(res))

res<-decideTests(fit, p.value=0.001, lfc=3)
par(col="blue4")
vennDiagram(vennCounts(res))

#png("todelall.png", bg="transparent")


contrast.matrix <-makeContrasts(del-minus, plus-minus, levels=design)
fit <- eBayes(contrasts.fit(fit, contrast.matrix))
res<-decideTests(fit, p.value=0.001, lfc=0.5)
par(col="red4")
vennDiagram(vennCounts(res, include="up"))

res<-decideTests(fit, p.value=0.001, lfc=1)
par(col="red4")
vennDiagram(vennCounts(res, include="up"))

res<-decideTests(fit, p.value=0.001, lfc=2)
par(col="red4")
vennDiagram(vennCounts(res, include="up"))

res<-decideTests(fit, p.value=0.001, lfc=3)
par(col="red4")
vennDiagram(vennCounts(res, include="up"))


contrast.matrix <-makeContrasts(del-minus, plus-minus, levels=design)
fit <- eBayes(contrasts.fit(fit, contrast.matrix))
res<-decideTests(fit, p.value=0.001, lfc=0.5)
par(col="green4")
vennDiagram(vennCounts(res, include="down"))

res<-decideTests(fit, p.value=0.001, lfc=1)
par(col="green4")
vennDiagram(vennCounts(res, include="down"))

res<-decideTests(fit, p.value=0.001, lfc=2)
par(col="green4")
vennDiagram(vennCounts(res, include="down"))

res<-decideTests(fit, p.value=0.001, lfc=3)
par(col="green4")
vennDiagram(vennCounts(res, include="down"))


install.packages("colorfulVennPlot")
