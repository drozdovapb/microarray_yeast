#todos
my function for phyper!
my function for intersection??? Not easy
Finctional classification


Intersections
========================================================

This analysis is intended to find out whether the genes affected by [ISP+] or deletion of SFP1 are the same.

## Introduction and abbreviations
p = [ISP+] (plus)
m = [isp-] (minus)
d = sfp1delta (delta)
s = subtracted (corrected for slow growth)

## Part Zero: Processing my data


```r
setwd("/media/drozdovapb//441DB6652D7ED741/Work/Current work/Projects/bioinfo/raw_data_our/")
#source("http://bioconductor.org/biocLite.R"); biocLite("limma")
library(limma)

alldata <-read.table("./Statistical data analysis (Sophie Lamarre)/normdata.csv", dec=",", sep="\t", head=TRUE, row.names="locus")
alldata.matrix <-as.matrix(alldata)
design <- model.matrix(~0 + factor(c(1,1,1,1,2,2,2,2,3,3,3,3)))
colnames(design) <- c("plus","minus","del")
fit <-lmFit(alldata, design)
contrast.matrix <-makeContrasts(plus-minus, del-minus, plus-del, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <-eBayes(fit2)
#everything above taken from limma user's guide
#extracting variables and writing tables
plus_minus<-topTable(fit2, coef=1, number=nrow(fit2), adjust.method="BH", sort.by="p")
del_minus<-topTable(fit2, coef=2, number=nrow(fit2), adjust.method="BH", sort.by="p")
plus_del<-topTable(fit2, coef=3, number=nrow(fit2), adjust.method="BH", sort.by="p")
#write.table(plus_minus, "../processed_tables_R/plus-minus.csv", sep="\t", dec=".", col.names=NA)
#write.table(del_minus, "../processed_tables_R/del-minus.csv", sep="\t", dec=".", col.names=NA)
#write.table(plus_del, "../processed_tables_R/minus-del.csv", sep="\t", dec=".", col.names=NA)

logdata<-cbind(plus_minus$logFC, del_minus$logFC, plus_del$logFC)
colnames(logdata)<-c("plus_min", "del_min", "plus_del")
row.names(logdata)<-row.names(plus_minus)
#write.table(logdata, "../processed_tables_R/logdata.csv", sep="\t", dec=".", col.names=NA, quote=FALSE)
```

Now we have a dataframe 'logdata' containing log2 fold changes for 3 comparisons. 

## Part 1: Extracting DEG lists from initial data


```r
cutoff2 <- function(mytable, p.value=.001, abslogFC=.5, up_or_down="all") {
  if (up_or_down == "all") {
  subset(mytable, abs(mytable$logFC)>abslogFC & mytable$adj.P.Val<p.value)}
  else if (up_or_down == "up") {
    subset(mytable, mytable$logFC>abslogFC & mytable$adj.P.Val<p.value)}
  else if (up_or_down == "down") {
    subset(mytable, mytable$logFC< (0-abslogFC) & mytable$adj.P.Val<p.value)}
}

pm_cut_up<-cutoff2(plus_minus, up_or_down="up")
pm_cut_down<-cutoff2(plus_minus, up_or_down="down")
dm_cut_up<-cutoff2(del_minus, up_or_down="up")
dm_cut_down<-cutoff2(del_minus, up_or_down="down")
```

## Intersections

Now we intersect the lists of genes affected by [ISP+] or SFP1 deletion.


```r
#install.packages("bvenn")
library(bvenn)
```

```
## Loading required package: grid
```

```r
library(Vennerable)
```

```
## Loading required package: graph
## Loading required package: RBGL
## Loading required package: lattice
## Loading required package: RColorBrewer
## Loading required package: reshape
## Loading required package: gtools
## Loading required package: xtable
```

```r
#install.packages("colorfulVennPlot")
library(colorfulVennPlot)
pm_and_dm<-list(c(row.names(pm_cut_up),row.names(pm_cut_down)), c(row.names(dm_cut_down), row.names(dm_cut_up)))
pm_and_dm2<-Venn(pm_and_dm, SetNames=c("pm", "dm"))
plot(pm_and_dm2)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-31.png) 

```r
#the problem of colors!!!!!!!!!!! It chooses colors itself!

pm_dm_up <- Venn(list(row.names(pm_cut_up), row.names(dm_cut_up)), SetNames=c("pm_up","dm_up"))
plot(pm_dm_up)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-32.png) 

```r
pm_dm_down <- Venn(list(row.names(pm_cut_down), row.names(dm_cut_down)), SetNames=c("pm_down","dm_down"))
plot(pm_dm_down)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-33.png) 

```r
pm_dm_ud <- Venn(list(row.names(pm_cut_up), row.names(dm_cut_down)), SetNames=c("pm_up","dm_down"))
plot(pm_dm_ud)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-34.png) 

```r
pm_dm_du <- Venn(list(row.names(pm_cut_down), row.names(dm_cut_up)), SetNames=c("pm_down","dm_up"))
plot(pm_dm_du)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-35.png) 

All intersections are > 0.

## Hypergeometrical distribution

Here we apply hypergeometrical distribution to answer the question whether these intersections could be due to chance.


```r
pm_all_dm_all <- phyper(q=213-1, m=332, n=6225-332, k=1773, lower.tail=F, log.p=F) #phyper(q=213-1, m=1773, n=6225-1773, k=332, lower.tail=F)
pm_up_dm_up <- phyper(q=20-1, m=74, n=6225-74, k=787, lower.tail=F, log.p=F)
pm_down_dm_down <- phyper(q=141-1, m=258, n=6225-258, k=986, lower.tail=F, log.p=F)
pm_down_dm_up <- phyper(q=33-1, m=258, n=6225-258, k=787, lower.tail=F, log.p=F)
pm_up_dm_down <- phyper(q=19-1, m=74, n=6225-74, k=986, lower.tail=F, log.p=F)
```

We see that intersection of all genes in both comparions is not likely to be a coincidence: p = 6.0811 &times; 10<sup>-44</sup>. The same applies to genes upregulated in both comparisons (p = 6.249 &times; 10<sup>-4</sup>) and for genes downregulated in both comparisons (p = 3.5537 &times; 10<sup>-49</sup>). 

In addition we applied the same test to genes diversely regulated by [ISP+] and deletion of SFP1. Interestingly, while for the genes downregulated by [ISP+] and upregulated by deletion seem to intersect due to chance p = 0.5002; for the other group (genes upregulated by [ISP+] and downregulated when SFP1 is deleted) we see p = 0.0193. 

This may mean that while the [ISP+] prion does not affect the direct targets of Sfp1p some of the functions are the same...
I don't know what the second part means.

## Correcting for slow growth

Our data (Rogoza et al., 2010) as well as data by other authors (e.g. Jorgensen et al., 2002) indicate that deletion of SFP1 significantly decreases growth rate.
A paper recently published by O'Duibhir et al. suggests a methods for removing slow grower's signature from microarray data. 


Intro by Philip Lijnzaad: 

> This script removes the slow growth signature from a dataset as
> described in O'Duibhir et al., Mol. Sys. Biol. 2014, and was used
> e.g. for the amino acid starvation data shown in Fig. 6A.

> For this, the slow growth signature is needed; in this script we use
> the one described in the publication, called
> 'Supplementary_data4.txt'. (You can easily calculate the signature from the
> original data; see script 'svd-transform.R' for that).

> The current script applies the signature removal to
> the original data, called 'Supplementary_data3.txt'.  (This is just
> by way of example, since that data had already its slow growth
> signature removed, with the script 'svd-transform.R' ...)

> It has been tested on Linux and Windows.
> On Unix, it can also work as a stand-alone script:

> chmod a+x remove.slowgrowth-profile.R ## needed once to make the script executable

> ./remove.slowgrowth-profile.R

> We use 'Supplementary_data4.txt' as the slow growth signature to be
> removed, and 'Supplementary_data3.txt' as the data set to remove it from.
> For simplicity, the script expects these files in the
> directory from which it is run. The correspondence between the
> transcripts in the slow growth signature and the data set is
> established based on gene names. The output will be only for
> transcripts present in both the input files.

> The script produces, again in the current directory, the files
> data-signatureremoved.txt (the input data that has now been made
> independent of the signature) and 'data-signaturestrength.txt' (the
> degree to which the signature was present in each of the columns of
> the the input data).

> Written by Philip Lijnzaad <p.lijnzaad@umcutrecht.nl>, december 2013



```r
#my working directory
setwd("/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/processed_tables_R/slow_growth/")
#my data
data.with.slowgrowth.signature <- "/media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/processed_tables_R/logdata.csv"
slowgrowth.signature <- "/media/drozdovapb/441DB6652D7ED741/Work/BigData/ODubhir/Supplementary_data4.txt"
```

Philip's code

> definitions of proper inner product and Euclidean norm


```r
dot.product <- function(x,y=x)drop(crossprod(x,y))
eucl.norm <- function(x)sqrt(dot.product(x))
```

> the main function:
> general function to correct the columns of 'dataset' by getting rid
> of a component indicated by 'signature'. Typical usage would be 
> e.g.:

>   removal <- remove.signature(timeseries, signature=svd$u[,1])
>   timeseries.new <- removal$data
>   signature.presence.per.timepoint <- removal$strength

> After this operation, all the dataset columns are orthogonal to the
> signature, and their correlation to the signature is zero.
> (note: zero-correlation equals orthogonality *only* for centered
> data, which in this case is true due to the normalization)
> The function returns a list with the transformed data, and a vector
> of the degree to which the signature was present in the input data



```r
remove.signature <- function(dataset, signature) {

  if(is.null(dim(dataset)))
    dataset <- as.matrix(dataset,ncol=1)

  if(is.null(rownames(dataset)))
    stop("dataset must have rownames")

  if(is.matrix(signature)) {
    if(ncol(signature)!=1)
      stop("signature must be a simple vector, or a single-column matrix")
    signature <- drop(signature)
  }
  if(is.null(names(signature)))
    stop("signature must be a named vector, or a single-column matrix with rownames")
  
  use <- intersect(rownames(dataset), names(signature))
  
  if ( length(use) == 0)
    stop("Data set and signature have no names in common")

  signature <- signature / eucl.norm(signature)

  column.factors <- dot.product(dataset[use,,drop=FALSE], signature[use])
  correction <- signature[use,drop=FALSE] %*% t(column.factors)
  return(list(data=dataset[use,,drop=FALSE] - correction, strength=column.factors))
}                                       #remove.signature
```



> convenience functions to read/write tab-delimited data with row+column names; topleft cell is empty.



```r
read.tab <- function(file, ...)read.table(file=file, sep="\t", as.is=TRUE, quote="",
                                          header=TRUE, comment.char="", row.names=1, ...)
write.tab <- function(x, file="", ...)write.table(x, file, sep="\t", quote=FALSE, na="",
                                                  eol="\r\n",
                                                  row.names=TRUE, col.names=NA, ...)

message("Reading data from ", data.with.slowgrowth.signature," ...\n") 
```

```
## Reading data from /media/drozdovapb/441DB6652D7ED741/Work/Current work/Projects/bioinfo/processed_tables_R/logdata.csv ...
```

```r
data <- read.tab(file=data.with.slowgrowth.signature)
numerical.cols <- names(which(unlist(lapply(data,class)=="numeric")))
M <- as.matrix(data[,numerical.cols])
### Matrix of log2(ratio)'s. Rows: transcripts; columns: deletion strains
annot.cols <- setdiff(names(data), numerical.cols) ## columns with annotation such as commonName

sig.data <- read.tab(file=slowgrowth.signature)
numerical.cols <- names(which(unlist(lapply(sig.data,class)=="numeric")))
signature <- as.matrix(sig.data[,numerical.cols, drop=FALSE])

message("Transforming the data ...\n")
```

```
## Transforming the data ...
```

```r
removal <- remove.signature(M, signature)

file <- "data-signatureremoved.txt"
message("Writing transformed data to ",file, " ...\n")
```

```
## Writing transformed data to data-signatureremoved.txt ...
```

```r
if(length(annot.cols)==0)  {
  write.tab(file=file, removal$data)
} else {
  ## put back the original annotations:
  transformed <- data.frame(data[rownames(removal$data),annot.cols],
                            removal$data)
  names(transformed)[1:length(annot.cols)] <- annot.cols
  write.tab(file=file, transformed)
}

file <- "data-signaturestrength.txt"
message("Writing signature presence to ",file, " ...\n")
```

```
## Writing signature presence to data-signaturestrength.txt ...
```

```r
write.tab(file=file,data.frame(signature.strength=removal$strength))

message("Done")
```

```
## Done
```


Philip's script ends here

Since now we have new fold change values => re-extraction of DE genes. 
pvalues are taken from the original comparisons that are in the 2014-07-07-short_limma.R script. Well, they could be but the number of rows don't match now. To overcome that, we merge data sets by row names. 


```r
plus_minus_subtr<-merge(as.data.frame(removal$data), plus_minus, by="row.names")
    row.names(plus_minus_subtr)<-plus_minus_subtr$Row.names
del_minus_subtr<-merge(as.data.frame(removal$data), del_minus, by="row.names")
    row.names(del_minus_subtr)<-del_minus_subtr$Row.names
plus_del_subtr<-merge(as.data.frame(removal$data), plus_del, by="row.names")
    row.names(plus_del_subtr)<-plus_del_subtr$Row.names


pms_cut_up<-subset(plus_minus_subtr, plus_minus_subtr$adj.P.Val<.001 & plus_minus_subtr$plus_min>.5)
#write.table(pms_cut_up, "2014-07-07-pm_subtr_up.csv", row.names=F)
pms_cut_down<-subset(plus_minus_subtr, plus_minus_subtr$adj.P.Val<.001 & plus_minus_subtr$plus_min< (-.5))
#write.table(pms_cut_down, "2014-07-07-pm_subtr_down.csv", row.names=F)
dms_cut_up<-subset(del_minus_subtr, del_minus_subtr$adj.P.Val<.001 & del_minus_subtr$del_min>.5)
#write.table(dms_cut_up, "2014-07-07-dm_subtr_up.csv")
dms_cut_down<-subset(del_minus_subtr, del_minus_subtr$adj.P.Val<.001 & del_minus_subtr$del_min < (-.5))
#write.table(dms_cut_down, "2014-07-07-dm_subtr_down.csv", row.names=F)
pds_cut_up<-subset(plus_del_subtr, plus_del_subtr$adj.P.Val<.001 & plus_del_subtr$plus_del>.5)
#write.table(pds_cut_up, "2014-07-07-pd_subtr_up.csv", row.names=F)
pds_cut_down<-subset(plus_del_subtr, plus_del_subtr$adj.P.Val<.001 & plus_del_subtr$abslogFC<(-.5), up_or_down="up")
#write.table(pds_cut_up, "2014-07-07-pd_subtr_down.csv")
```


## Intersections II

The first question is whether this transformation affects pm (it shouldn't) and dm (it should).


```r
setdiff(row.names(pm_cut_up), row.names(pms_cut_up)) #one gene "YHR213W" #it's a pseudogene #http://www.yeastgenome.org/cgi-bin/locus.fpl?dbid=S000001256
```

```
## [1] "YHR213W"
```

```r
length(setdiff(row.names(pm_cut_down), row.names(pms_cut_down))) #11 genes  #"YBR093C" "YJL153C" "YIL117C" "YOL082W" "YOL083W" "YIR017C" "YGR043C" "YBR046C" "YOL048C" "YKL103C" "YBR196C"
```

```
## [1] 11
```

```r
length(setdiff(row.names(dm_cut_up), row.names(dms_cut_up))) #578! genes
```

```
## [1] 578
```

```r
length(setdiff(row.names(dm_cut_down), row.names(dms_cut_down))) #607! genes
```

```
## [1] 607
```

```r
#what's happening?
#Venn diagrams
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("graph", "RBGL", "gtools", "xtable"))
#install.packages("Vennerable", repos="http://R-Forge.R-project.org")

vennd<-list(row.names(dm_cut_up), row.names(dms_cut_up))
venndd<-Venn(vennd, SetNames=c("dm_up", "dms_up"))
plot(venndd)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) 

```r
pms_dms<-list(c(row.names(pms_cut_up), row.names(pms_cut_down)), c(row.names(dms_cut_up), row.names(dms_cut_down)))
pms_dms_Venn<-Venn(pms_dms, SetNames=c("pms", "dms"))
plot(pms_dms_Venn)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) 

```r
pms_dms_uu<-list(row.names(pms_cut_up), row.names(dms_cut_up))
pms_dms_uu_Venn<-Venn(pms_dms_uu, SetNames=c("pms_up", "dms_up"))
plot(pms_dms_uu_Venn)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-103.png) 

```r
pms_dms_dd<-list(row.names(pms_cut_down), row.names(dms_cut_down))
pms_dms_dd_Venn<-Venn(pms_dms_dd, SetNames=c("pms_down", "dms_down"))
plot(pms_dms_dd_Venn)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-104.png) 

```r
pms_dms_ud<-list(row.names(pms_cut_up), row.names(dms_cut_down))
pms_dms_ud_Venn<-Venn(pms_dms_ud, SetNames=c("pms_up", "dms_down"))
plot(pms_dms_ud_Venn)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-105.png) 

```r
pms_dms_du<-list(row.names(pms_cut_down), row.names(dms_cut_up))
pms_dms_du_Venn<-Venn(pms_dms_du, SetNames=c("pms_down", "dms_up"))
plot(pms_dms_du_Venn)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-106.png) 

## Hypergeometrical distribution II


```r
pms_all_dms_all <- phyper(q=206-1, m=320, n=6225-320, k=1148, lower.tail=F, log.p=F) #phyper(q=213-1, m=1773, n=6225-1773, k=332, lower.tail=F)
pms_up_dms_up <- phyper(q=19-1, m=73, n=6225-73, k=552, lower.tail=F, log.p=F)
pms_down_dms_down <- phyper(q=94-1, m=249, n=6225-249, k=596, lower.tail=F, log.p=F)
pms_down_dms_up <- phyper(q=74-1, m=249, n=6225-249, k=552, lower.tail=F, log.p=F)
pms_up_dms_down <- phyper(q=19-1, m=73, n=6225-73, k=596, lower.tail=F, log.p=F)
```


We see that intersection of all genes in both comparions is not likely to be a coincidence: p = 2.7433 &times; 10<sup>-78</sup>. The same applies to genes upregulated in both comparisons (p = 1.247 &times; 10<sup>-5</sup>) and for genes downregulated in both comparisons (p = 6.3438 &times; 10<sup>-35</sup>) and to the genes downregulated by [ISP+] and upregulated by deletion seem to intersect due to chance p = 3.1515 &times; 10<sup>-22</sup>, as well as for the other group (genes upregulated by [ISP+] and downregulated when SFP1 is deleted) we see p = 3.7019 &times; 10<sup>-5</sup>. 

So???
Some potentially misinterpreted results.
Now we skip transforming pm.

Intersections III


```r
pm_dms<-list(c(row.names(pm_cut_up), row.names(pm_cut_down)), c(row.names(dms_cut_up), row.names(dms_cut_down)))
pm_dms_Venn<-Venn(pm_dms, SetNames=c("pm", "dms"))
plot(pm_dms_Venn)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) 

```r
pm_dms_uu<-list(row.names(pm_cut_up), row.names(dms_cut_up))
pm_dms_uu_Venn<-Venn(pm_dms_uu, SetNames=c("pms_up", "dms_up"))
plot(pm_dms_uu_Venn)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) 

```r
pm_dms_dd<-list(row.names(pm_cut_down), row.names(dms_cut_down))
pm_dms_dd_Venn<-Venn(pm_dms_dd, SetNames=c("pms_down", "dms_down"))
plot(pm_dms_dd_Venn)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) 

```r
pm_dms_ud<-list(row.names(pm_cut_up), row.names(dms_cut_down))
pm_dms_ud_Venn<-Venn(pm_dms_ud, SetNames=c("pm_up", "dms_down"))
plot(pm_dms_ud_Venn)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-124.png) 

```r
pm_dms_du<-list(row.names(pm_cut_down), row.names(dms_cut_up))
pm_dms_du_Venn<-Venn(pm_dms_du, SetNames=c("pm_down", "dms_up"))
plot(pm_dms_du_Venn)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-125.png) 

```r
write.table(intersect(row.names(pm_cut_up), row.names(dms_cut_up)), "pm_dms_uu.csv")
write.table(intersect(row.names(pm_cut_down), row.names(dms_cut_down)), "pm_dms_dd.csv")
```

## Hypergeometrical distribution III


```r
pm_all_dms_all <- phyper(q=212-1, m=332, n=6225-332, k=1148, lower.tail=F, log.p=F) 
pm_up_dms_up <- phyper(q=19-1, m=74, n=6225-74, k=552, lower.tail=F, log.p=F)
pm_down_dms_down <- phyper(q=94-1, m=258, n=6225-248, k=596, lower.tail=F, log.p=F)
pm_down_dms_up <- phyper(q=75-1, m=258, n=6225-258, k=552, lower.tail=F, log.p=F)
pm_up_dms_down <- phyper(q=19-1, m=74, n=6225-74, k=596, lower.tail=F, log.p=F)
```

All genes: p = 1.0453 &times; 10<sup>-79</sup>. 
Genes upregulated in both comparisons (p = 1.5418 &times; 10<sup>-5</sup>) and genes downregulated in both comparisons (p = 1.6403 &times; 10<sup>-33</sup>) 
Genes downregulated by [ISP+] and upregulated by deletion: p = 7.1079 &times; 10<sup>-22</sup>.
Genes upregulated by [ISP+] and downregulated when SFP1 is deleted): p = 4.5446 &times; 10<sup>-5</sup>. 
All intersections are significant. 
