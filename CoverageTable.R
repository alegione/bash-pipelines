#!/bin/R

#Allow command line arguments
args = commandArgs(trailingOnly=TRUE)

varSWD<-args[1]	#working directory for depth files
setwd(varSWD)
CoverageTable<-read.table(file = "Results_Summary/CoverageTable.tsv", header = TRUE)


varfilename<-args[2]	#depth file in .txt format

varDepthFile<-paste0("Coverage/",varfilename,".depth.txt")

dat <- read.table(varDepthFile)

dat$v4 <- log2(dat$V3 + 1) #will take a few seconds
colnames(dat) <- c("ref", "pos", "depth", "log2_depth")

len<-max(dat$pos)
aveDepth<-round(mean(dat$depth), digits = 2)
aveLog2Depth<-round(mean(dat$log2_depth), digits = 2)
coverage<-sum(dat$depth > 0)
coverageAbove10<-sum(dat$depth > 10)
missingBases<-(len-coverage)
percentCoverage<-round(coverage/len*100, 2)
percentCoverageAbove10<-round(coverageAbove10/len*100, 2)
medianDepth<-median(dat$depth)
minDepth<-as.numeric(min(dat$depth))
maxDepth<-as.numeric(max(dat$depth))
q1depth<-as.numeric(quantile(dat$depth,0.25))
q3depth<-as.numeric(quantile(dat$depth,0.75))

#output results
varNewtable<-rbind(CoverageTable, data.frame(Name=varfilename,ReferenceLength=len,AverageDepth=aveDepth,MinDepth=minDepth,Q1=q1depth,MedianDepth=medianDepth,Q3=q3depth,MaxDepth=maxDepth,Coverage=percentCoverage,GoodCoverage=percentCoverageAbove10))
write.table(x = varNewtable, quote = FALSE, file = "Results_Summary/CoverageTable.tsv", sep = "\t", row.names = FALSE)

