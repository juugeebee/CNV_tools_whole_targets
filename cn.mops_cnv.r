# Author: Julie BOGOIN

library(cn.mops)

segments <- read.table(file="/media/Data1/jbogoin/ref/gencode/v34_hg38/gencode.v34.basic.annotation.CDS.merged.bed",
                    header=FALSE, sep="\t", as.is=TRUE) 

BAMFiles <- list.files(path=".", pattern=".dedup.bam$")

gr <- GRanges(segments[,1], IRanges(segments[,2],segments[,3]))

X <- getSegmentReadCountsFromBAM(BAMFiles, GR=gr)

resCNMOPS <- exomecn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

#plot(resCNMOPS, which=5)

setwd(dir="./cn.mops_output")

segm <- as.data.frame(segmentation(resCNMOPS))
write.csv(segm, file="segmentation.csv")

CNVs <- as.data.frame(cnvs(resCNMOPS))
write.csv(CNVs, file="cnvs.csv")

CNVRegions <- as.data.frame(cnvr(resCNMOPS))
write.csv(CNVRegions, file="cnvr.csv")
