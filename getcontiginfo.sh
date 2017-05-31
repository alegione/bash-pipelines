#!/bin/bash


#echo -e "Name\tcontigs >= 1000 bp\tlength range\tsum of lengths\t%identity average\t%identity range" > Shared_Folder/CpecHope/spadesoutput/contiginfo.txt
#1=Shared_Folder/CpecAssemble/spadesoutput/
while read i; do
#contigs >= 1000 bp
	contignumber=$(grep "number of" Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt | cut -f13 -d " ")

	contigmin=$(grep "length range" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f4 -d " " | cut -f1 -d "-")

	contigmax=$(grep "length range" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f4 -d " " | cut -f2 -d "-")
	
	contigrange=$(echo "$contigmin - $contigmax")
	contigsum=$(grep "sum of contigs" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f6 -d " ")

	averageid=$(grep "average id" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f4 -d " ")

#id range low
	idlow=$(grep "id% range" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f4 -d " " | cut -f1 -d "-")
#id range high
	idhigh=$(grep "id% range" "Shared_Folder/CpecAssemble/spadesoutput/$i/contig_summary/$i.blastn.summary.txt" | cut -f4 -d " " | cut -f2 -d "-")
	idrange=$(echo "$idlow - $idhigh")
	
	if ! grep -q $i "Shared_Folder/CpecHope/spadesoutput/contiginfo.txt"; then
		echo -e "$i\t$contignumber\t$contigrange\t$contigsum\t$averageid\t$idrange" >> Shared_Folder/CpecHope/spadesoutput/contiginfo.txt
	fi
done < Shared_Folder/CpecAssemble/spadescontents.txt
