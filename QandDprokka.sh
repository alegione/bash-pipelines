#!/bin/bash

ls "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/" > meta2.txt

sed -i 's/\.fasta//g' meta2.txt

while read fileline; do
	AlBin/prokka_bin/prokka --outdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes/$fileline" --prefix $fileline --genus Chlamydia --addgenes "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes/$fileline.fasta" --force --usegenus

	AlBin/prokka/bin/prokka --outdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/prokka_genomes_noquestion/$fileline" --prefix "${fileline}_noquestion" --genus Chlamydia --addgenes "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes_noquestion/$fileline.fasta" --force --usegenus

done < meta.txt

