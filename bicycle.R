allDataPlus<-alldata[,1:4]
allDataMinus<-alldata[,5:8]
allDataDel<-alldata[,9:12]


allDataPlusZ<-apply(allDataPlus, 1, sd)/apply(allDataPlus, 1, mean) 
allDataMinusZ<-apply(allDataMinus, 1, sd)/apply(allDataPlus, 1, mean) 
allDataDelZ<-apply(allDataDel, 1, sd)/apply(allDataDel, 1, mean) 

mean(allDataPlusZ)
mean(allDataMinusZ)
mean(allDataDelZ)