#!/bin/bash

ls "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/" > meta2.txt

sed -i 's/\.fasta//g' meta2.txt
mkdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/Annotated_genomes/"

while read fileline; do
	AlBin/prokka_bin/prokka --outdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/Annotated_genomes/$fileline" --prefix "$fileline" --locustag "$fileline" --kingdom Bacteria --addgenes "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/$fileline.fasta"

done < meta2.txt

