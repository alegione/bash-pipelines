#!/bin/bash

# Use the below to get a table of the quality scores from fastqc output, will then have to send it to an R script
# awk '/>>Per sequence quality scores/,/>>END_MODULE/' Desktop/Shared_folder/CpecFinal/fastQC_output/fastqc_data.txt  | grep -v ">>"
#echo $(Rscript Desktop/Shared_folder/Scripts/R/GetQualScore.R Desktop/Shared_folder/CpecFinal)
echo -e "Genome\tOriginal Average Quality\tTrimmed Average Quality" > Qualtable.tsv
rm -r tmp
while read i; do
	mkdir tmp
	unzip -q Desktop/Shared_folder/CpecFinal/fastQC_output/Original/${i}_1_fastqc.zip -d tmp
	awk '/>>Per sequence quality scores/,/>>END_MODULE/' tmp/${i}_1_fastqc/fastqc_data.txt  | grep -v ">>" | sed 's/#//g' > tmp/Qual.txt
	OriginalQual=$(Rscript Desktop/Shared_folder/Scripts/R/GetQualScore.R tmp/)
	rm -r tmp
	mkdir tmp
	unzip -q Desktop/Shared_folder/CpecFinal/fastQC_output/Trimmed/${i}_1_fastqc.zip -d tmp
	awk '/>>Per sequence quality scores/,/>>END_MODULE/' tmp/${i}_1_fastqc/fastqc_data.txt  | grep -v ">>" | sed 's/#//g' > tmp/Qual.txt
	TrimmedQual=$(Rscript Desktop/Shared_folder/Scripts/R/GetQualScore.R tmp)
	rm -r tmp
	echo -e "$i\t$OriginalQual\t$TrimmedQual" >> Qualtable.tsv
done < Desktop/Shared_folder/CpecFinal/MetaData/GenomeList.txt

