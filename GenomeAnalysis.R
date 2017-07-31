#!/usr/local/bin/r
#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

#Allow command line arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input file).n", call.=FALSE)
#} else if (length(args)==1) {
#  # default output file
#  args[2] = "out.txt"
#}

varSWD<-args[1]	#working directory for depth files
setwd(varSWD)
CoverageTable<-read.table(file = "coverage/CoverageTable.txt", header = TRUE)
#View(CoverageTable)


#if function here to get see if a table exists?

varfilename<-args[2]	#depth file in .txt format
varDepthFile<-paste0("coverage/",varfilename,".depth.txt")

dat <- read.table(varDepthFile)
#View(dat)
dat$v4 <- log2(dat$V3 + 1) #will take a few seconds
colnames(dat) <- c("ref", "pos", "depth", "log2_depth")

len<-max(dat$pos)
#len
aveDepth<-round(mean(dat$depth), digits = 2)
#aveDepth
aveLog2Depth<-round(mean(dat$log2_depth), digits = 2)
#aveLog2Depth
coverage<-sum(dat$depth > 0)
#coverage
coverageAbove10<-sum(dat$depth > 10)
#coverageAbove10
missingBases<-(len-coverage)
#missingBases
percentCoverage<-round(coverage/len*100, 2)
#percentCoverage
percentCoverageAbove10<-round(coverageAbove10/len*100, 2)
#percentCoverageAbove10
medianDepth<-median(dat$depth)
#medianDepth
minDepth<-as.numeric(min(dat$depth))
maxDepth<-as.numeric(max(dat$depth))
q1depth<-as.numeric(quantile(dat$depth,0.25))
q3depth<-as.numeric(quantile(dat$depth,0.75))

#output results
varNewtable<-rbind(CoverageTable, data.frame(Name=varfilename,ReferenceLength=len,AverageDepth=aveDepth,MinDepth=minDepth,Q1=q1depth,MedianDepth=medianDepth,Q3=q3depth,MaxDepth=maxDepth,Coverage=percentCoverage,GoodCoverage=percentCoverageAbove10))
write.table(x = varNewtable, quote = FALSE, file = "coverage/CoverageTable.txt", sep = "\t", row.names = FALSE)
#View(varNewtable)
#varPDFname<-paste0(varfilename,".pdf") #concatentates a string, here the file name with the file type
varTIFFname<-paste0("coverage/",varfilename,".tiff") #concatentates a string, here the file name with the file type


#pdf(varPDFname, width = 600, height = 600, useDingbats = FALSE)

if(!file.exists(varTIFFname)) {
  options(scipen = 10)
  tiff(filename = varTIFFname, width = 1800, height = 1200, res = 300)
  #tiff(filename = "19234_1_nametestsmall.tiff", res = 600)
  #plot(dat$pos, dat$log2_depth, type = "l", xlab = varfilename, ylab = "Normalised depth of coverage")
  par(bty = "l",xaxs="i",yaxs="i")
  plot(dat$pos, dat$log2_depth, type = "l", xlab = "", ylab = "Normalised depth of coverage")
#  plot(dat$pos, dat$log2_depth, type = "l", xlab = "", ylim = c(0,8), ylab = "Normalised depth of coverage")
  lines(lowess(dat$pos,dat$log2_depth, f = 0.01), col = "red")
  dev.off()
}


