#!/bin/bash

ls "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/Annotated_genomes/" > meta2.txt

sed -i 's/\.fasta//g' meta2.txt
mkdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/genomes_gff"

while read fileline; do
	cp -vf "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/Annotated_genomes/$fileline/$fileline.gff" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/genomes_gff/$fileline.gff" 
done < meta2.txt

