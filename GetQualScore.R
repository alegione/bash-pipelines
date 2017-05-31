#Takes in quality scores table and outputs and average quality score for a fastqc file
#Allow command line arguments
args = commandArgs(trailingOnly=TRUE)

varSWD<-args[1]		#current working directory
#varfilename<-args[2]	#depth file in .tsv format

setwd(varSWD)
QualScores<-read.table(file = "Qual.txt", header = TRUE)
QualScores$Product<-(QualScores$Quality*QualScores$Count)
AveQual <- sum(QualScores$Product)/sum(QualScores$Count)
AveQual <- round(AveQual, digits = 1)
cat(AveQual)
