#!/bin/bash

ls "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/" > meta2.txt

sed -i 's/\.fasta//g' meta2.txt

while read fileline; do
	cp -v "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes/$fileline/${fileline}.gbk" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes/genomes/${fileline}.gbk" 
	cp -v "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/$fileline/${fileline}_noquestion.gbk" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/genomes_noquestion/${fileline}_noquestion.gbk" 
done < meta.txt

