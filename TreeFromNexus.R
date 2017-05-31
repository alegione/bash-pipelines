### To convert MrBayes consensus tree output to Newick format
### Inspired/adapted from code from Dr Neil Young
#install.packages("ape")
library(ape)
args = commandArgs(trailingOnly=TRUE)
options(digits=4)
varSWD<-args[1]		#Current file directory
varFileName<-args[2] #Current tree being converted (should end in *.con.tre)
#varNexusInput<-paste0(varFileName,".con.tre")
#varTreeOutput<-paste0(varFileName,".newick.tre")
setwd(varSWD)

varStatInput<-read.table(file = paste0(varFileName,".pstat"), header = TRUE, skip = 1)
statMin<-min(varStatInput$PSRF)
statMax<-max(varStatInput$PSRF)
statMedian<-median(varStatInput$PSRF)

if((statMedian < 0.96) || (statMedian > 1.04)) {
  varControl<-"No"
} else {
  varControl<-"Yes"
}

cat(varFileName,statMedian, statMin, statMax, varControl, sep = "\t")

tree <- read.nexus(file=paste0(varFileName,".con.tre"))

write.tree(tree,file=paste0(varFileName,".newick.tre"))
