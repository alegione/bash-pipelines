### To convert to Newick format
## Made a tree conversion program in R, saved this as convert_tree.R
#install.packages("ape")
library(ape)
args = commandArgs(trailingOnly=TRUE)
options(digits=4)
varSWD<-args[1]		#Current file directory
varFileName<-args[2] #Current tree being converted
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
# Loop through the Mr Bayes output trees and save with a new name. For this to work, you just have to make sure your trees end in ".nexus.mb.con.tre".
#ls *nexus.mb.con.tre | sed 's/.nexus.mb.con.tre//1' | while read line; do cp $line.nexus.mb.con.tre in.tre; R --vanilla < convert_tree.R; mv temp.R.tre $line.mb.tre; done
