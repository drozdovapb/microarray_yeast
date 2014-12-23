

plus_minus_cut <- cutoff(plus_minus, .5, 0.001)
plus_del_cut <- cutoff(plus_del, .5, 0.001)
minus_del_cut <- cutoff(minus_del, .5, 0.001)
del_minus_cut <- cutoff(del_minus, .5, 0.001)


plus_minus_cut_up <- subset(plus_minus_cut, plus_minus_cut$logFC > 0)
plus_minus_cut_down <- subset(plus_minus_cut, plus_minus_cut$logFC < 0)
plus_del_cut_up <-subset(plus_del_cut, logFC > 0) 
plus_del_cut_down <- subset(plus_del_cut, logFC < 0)
minus_del_cut_up <- subset(minus_del_cut, logFC > 0)
minus_del_cut_down <- subset(minus_del_cut, logFC <0)
del_minus_cut_up <- subset(del_minus_cut, logFC > 0)
del_minus_cut_down <- subset(del_minus_cut, logFC <0)

intersect(row.names(plus_minus_cut_up), row.names(plus_minus_cut_up))
intersect(row.names(plus_minus_cut), row.names(del_minus_cut))
write.table(intersect(row.names(plus_minus_cut_up), row.names(del_minus_cut_down)), "2014-06-07-pm_up_dm_down.csv")
write.table(intersect(row.names(plus_minus_cut_down), row.names(del_minus_cut_up)), "2014-06-07-pm_down_dm_up.csv")