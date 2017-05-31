#!/home/genomics/.linuxbrew/bin/R

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
#library("seqinr")

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
#setwd("PhD/Sequences/CpecHope/")
setwd(varSWD)
#getwd()
CoverageTable<-read.table(file = "Results_Summary/CoverageTable.tsv", header = TRUE)
#View(CoverageTable)

#varFasta<-args[2]

#varFasta <- read.fasta(file = "PhD/Sequences/Cpec_whole_genome/03RayIs_3F3_B_UGT.fasta")
#varFasta
#varFastaseq <- varFasta[[1]]  # Put the sequence in a vector.
#varFastaseq
#varGC<-round(GC(varFastaseq)*100,2)
#varGC<-00
#if function here to get see if a table exists?
#varDataOld<-data.frame(Length=numeric(),'Percentage Coverage'=numeric(),'PercentCoverage where depth above 10'=numeric(),'GC'=numeric(),'Average Depth'=numeric(),'Median Depth'=numeric())

varfilename<-args[2]	#depth file in .txt format
#varfilename<-"19238_1_02"
varDepthFile<-paste0("Coverage/",varfilename,".depth.txt")

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
write.table(x = varNewtable, quote = FALSE, file = "Results_Summary/CoverageTable.tsv", sep = "\t", row.names = FALSE)
#View(varNewtable)
#varPDFname<-paste0(varfilename,".pdf") #concatentates a string, here the file name with the file type
