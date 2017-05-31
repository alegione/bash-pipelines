#!/bin/R

#Allow command line arguments
args = commandArgs(trailingOnly=TRUE)

varSWD<-args[1]		#current working directory
varfilename<-args[2]	#depth file in .txt format

setwd(varSWD)

CoverageTable<-read.table(file = "Results_Summary/CoverageTable.tsv", header = TRUE)
varYLim<-round(max(CoverageTable$MaxDepth),digits = 0)
varTIFFname<-paste0("Results_Summary/plots/",varfilename,".tiff") #concatentates a string, here the file name with the file type

if(!file.exists(varTIFFname)) {
  varDepthFile<-paste0("Coverage/",varfilename,".depth.txt")
  dat <- read.table(varDepthFile)
  dat$v4 <- log2(dat$V3 + 1) #will take a few seconds
  dat$v5 <- (dat$V3 + 1)
  
  colnames(dat) <- c("ref", "pos", "depth", "log2_depth", "depth_and_1")
  options(scipen = 10)
  tiff(filename = varTIFFname, width = 1800, height = 1200, res = 300)
  par(bty = "o",xaxs="i",yaxs="i", cex = 1)
  plot(dat$pos, dat$depth_and_1, type = "l", xlab = "E58 reference genome (Mb)", ylim = c(1,varYLim), ylab = expression("Depth of coverage (log "[10]*")"), log = "y", yaxt = "n", xaxt = "n")
  lines(lowess(dat$pos,dat$depth_and_1, f = 0.01), col = "red")
  title(main = varfilename, line = 2.5)
  # The below axes are specific to a Chlamydia pecorum project and will not be applicable to others
  axis(side = 2, at = c(1, 10, 100, 1000, 10000, 100000), labels = c("", "1","2", "3", "4", "5"), tick = TRUE)
  axis(side = 3, at = c(159117, 163832, 362801, 363976, 590994, 627950,900088,939121), labels = FALSE, tick = TRUE)
  axis(side = 3, at = c(((362801+363976)/2), ((590994+627950)/2)),labels = c("ompA","pmp"), font = 3, tick = FALSE)
  axis(side = 3, at = c(((159117+163832)/2), ((900088+939121)/2)), labels = c("23S/16S","PZ"), tick = FALSE)
  axis(side = 1, at = c(seq(0,1100000, 100000)), labels = c(seq(0,1.1,0.1)))
  dev.off()
}
