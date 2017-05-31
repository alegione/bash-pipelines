#Allow command line arguments
library(ape, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(methods, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
library(ade4, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)
suppressMessages(library(adegenet, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
library(pegas, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE)


args = commandArgs(trailingOnly=TRUE)

varSWD<-args[1]		#current working directory
varfilename<-args[2]	#depth file in .txt format
setwd(varSWD)

fastaFile<-paste0("core_genome_alignments/",varfilename,".aln")

fastaFile<-read.dna(file = fastaFile, format = "fasta")
# get nucleotide diversity
Div<-nuc.div(x = fastaFile)
# get number of haplotypes
h<-haplotype(x = fastaFile)
h<-attr(x = h,"dim")[[1]]
SegSites<-length(seg.sites(fastaFile))
#Test of the Neutral Mutation Hypothesishw
TajD<-tajima.test(fastaFile)[[1]]
TajDPval<-tajima.test(fastaFile)[[2]]

cat(varfilename,h,SegSites,Div,TajD,TajDPval, sep = "\t")

#The below can be used on the final output table if you run the above in a loop
#align_metric<-read.table(file = "Alignment_metrics.tsv", header = TRUE, sep = "\t")
#align_metric$PAdjust_BH<-p.adjust(p=align_metric$P.value, method = "BH")
#align_metric$PAdjust_Bonf<-p.adjust(p=align_metric$P.value, method = "bonferroni")
