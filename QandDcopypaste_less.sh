#!/bin/bash

ls "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/" > meta2.txt

sed -i 's/\.fasta//g' meta2.txt
mkdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/genomes_noquestion"

while read fileline; do
	cp -v "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/$fileline/$fileline.gbk" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/genomes_noquestion/$fileline.gbk" 
done < meta2.txt

