#!/bin/bash


mkdir "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/renamed_draft_genomes"
while read old new; do
	cp "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/draft_genomes/$old.fasta" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/renamed_draft_genomes/$new.fasta"
done < renamefile.txt
