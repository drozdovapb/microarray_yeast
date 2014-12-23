phyper(q=213-1, m=332, n=6225-332, k=1773, lower.tail=F, log.p=F); phyper(q=213-1, m=1773, n=6225-1773, k=332, lower.tail=F)
phyper(q=20-1, m=74, n=6225-74, k=787, lower.tail=F, log.p=F)
phyper(q=141-1, m=258, n=6225-258, k=986, lower.tail=F, log.p=F)
phyper(q=33-1, m=258, n=6225-258, k=787, lower.tail=F, log.p=F)
phyper(q=19-1, m=74, n=6225-74, k=979, lower.tail=F, log.p=F)


phyper(10, 1000, 20000, 1000, log.p=F, lower.tail=F)

library(Vennerable)

pm_and_dm<-list(c(row.names(pm_cut_up),row.names(pm_cut_down)), c(row.names(dm_cut_down), row.names(dm_cut_up)))
pm_and_dm2<-Venn(pm_and_dm, SetNames=c("ISP+", "sfp1dela"))
plot(pm_and_dm2)

#the problem of colors!!!!!!!!!!! It chooses colors itself!

pm_dm_up <- Venn(list(row.names(pm_cut_up), row.names(dm_cut_up)), SetNames=c("prion_up","sfp1_up"))
plot(pm_dm_up)

pm_dm_down <- Venn(list(row.names(pm_cut_down), row.names(dm_cut_down)), SetNames=c("prion_down","sfp1_down"))
plot(pm_dm_down)

pm_dm_ud <- Venn(list(row.names(pm_cut_up), row.names(dm_cut_down)), SetNames=c("prion_up","sfp1_up"))
plot(pm_dm_ud)

pm_dm_du <- Venn(list(row.names(pm_cut_down), row.names(dm_cut_up)), SetNames=c("prion_up","sfp1_up"))
plot(pm_dm_du)
