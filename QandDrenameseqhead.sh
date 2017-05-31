#!/bin/bash

while read old new; do
	sed -i "s/>.*/>$new/" "Shared_Folder/Cpec_whole_genome/ReferenceGenomes/renamed_draft_genomes/$new.fasta"
done < renamefile.txt
