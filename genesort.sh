#!/bin/sh
########################################
#
# run pipeline 
#
########################################

#things to add: how does grep deal with duplicates?? assuming it takes EVERYTHING
#		need to tweak my prompt, not that it matters too much


if [ ! -d "/home/qiime/Desktop/Shared_Folder/Genes/" ]; then
	mkdir /home/qiime/Desktop/Shared_Folder/Genes/
fi


echo "What project are you running"
read Project
echo "You entered: $Project"


Dir="/home/qiime/Desktop/Shared_Folder/Genes/$Project" #the location to store the project files



if [ ! -d "$Dir" ]; then #if the project directory doesn't exist, create all appropriate folders and sub-folders
	mkdir $Dir
	mkdir $Dir/Metadata
fi

Meta="$Dir/Metadata"
Mapping="$Meta/Genelist"
AllGenes="$Meta/AllGenes"
echo "" > "$AllGenes.txt"
FileLocation="/home/qiime/Desktop/Shared_Folder/Cpec_whole_genome/Koalagenes.txt"
echo "I'm here"

#read -p "Is your target $FileLocation? (Y/N)" Switch

#if [[ $Switch != "y"* ]]; then
#	echo "target file location:"
#	#/home/qiime/Desktop/Shared_Folder/Cpec_whole_genome/Koalagenes.txt
#	read FileLocation
#fi

echo "$FileLocation"

if [ ! -a "$AllGenes.txt" ]; then
	sed 's/[] [\/,:"@-]/_/g' $FileLocation | sed '1d' | sort -k2 > "$AllGenes.txt" #remove special characters, important to have the hyphen at the start or end, otherwise it will treat it as signifying a range
	cut -f 1 "$AllGenes.txt" | sort -u > "$Meta/GenomeList.txt"
fi
GenomeCount=$(wc -l < "$Meta/GenomeList.txt")

echo "Genomes detect = $GenomeCount"


#if all ORFS have novel names, can sort as per genome list above

if [ ! -a "$Mapping.uniq.txt" ]; then
	cut -f 2 "$AllGenes.txt" > "$Mapping.txt"
	uniq -i "$Mapping.txt"  > "$Mapping.uniq.txt"
fi


Count="1"
TotalORFs=$(wc -l < "$Mapping.uniq.txt")
echo "Total unique ORFs = $TotalORFs"

if [ ! -d "$Dir/Genes" ]; then #if the genes directory doesn't exist, create all appropriate folders and sub-folders
	mkdir $Dir/Genes
	mkdir $Dir/Genes/Individual
	mkdir $Dir/Proteins
	mkdir $Dir/Proteins/Individual
	while read FileLine; do
		echo "processing ORF $Count/$TotalORFs => $FileLine"
		#look for the gene name (FileLine) from the mapping file, within the list of unique genes, add > to start for fasta and merge sample and gene name, then remove either gene or protein and add a new line between the fasta heading and the sequence
		grep -h -w "$FileLine" "$AllGenes.txt" | sed 's/^/>/' | cut -f 1,5 | sed -e 's/\t/\n/' > "$Dir/Genes/Individual/$Count.$FileLine.fa" #-h stops the file name being printed
		grep -h -w "$FileLine" "$AllGenes.txt" | sed 's/^/>/' | cut -f 1,6 | sed -e 's/\t/\n/' > "$Dir/Proteins/Individual/$Count.$FileLine.fa" #-h stops the file name being printed
		Count=$((Count+1))
	done < "$Mapping.uniq.txt"
fi
# -e 's/>/\n>/'


if [ ! -e "$Dir/Metadata/NumberedGeneList.txt" ]; then
	Count="1"
	ls "$Dir/Genes/Individual/" > "$Dir/Metadata/NumberedGeneList.txt"
	ls "$Dir/Proteins/Individual/" > "$Dir/Metadata/NumberedProteinList.txt"
	Low=$GenomeCount
	High=$GenomeCount
	FileArray=($GenomeCount)
	while read FileLine; do
		Count=$(grep -c ">" "$Dir/Genes/Individual/$FileLine")
		FileArray=("${FileArray[@]}" "$Count")
		if [ "$Count" -gt "$High" ]; then
			High=$Count
		elif [ "$Count" -lt "$Low" ]; then
			Low=$Count
		fi
	done < "$Dir/Metadata/NumberedGeneList.txt"
	
	sorted_unique_ids=($(echo "${FileArray[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ' | sed 's/://g'))
	echo "${sorted_unique_ids[@]}"
	
	echo "The most lines are $High and the fewest are $Low"

	Count="1"
	High=${#sorted_unique_ids[@]}

	echo "Processing $High new folders"
fi


if [ ! -d "$Dir/Genes/Counted/" ]; then
	mkdir "$Dir/Genes/Counted/"
	for i in "${sorted_unique_ids[@]}"; do
		mkdir "$Dir/Genes/Counted/$i.genes"
		echo "$Count/$High"
		while read FileLine; do
			Genes=$(grep -c ">" "$Dir/Genes/Individual/$FileLine")
			if [ "$i" -eq "$Genes" ]; then
				cp "$Dir/Genes/Individual/$FileLine" "$Dir/Genes/Counted/$i.genes"
			fi
		done < "$Dir/Metadata/NumberedGeneList.txt"
		Count=$((Count+1))
	done
fi

mkdir "$Dir/Aligned/"

i="0"
for i in "${sorted_unique_ids[@]}"; do
	WorkDir="$Dir/Genes/Counted/$i.genes"
	echo "$WorkDir"
	ls $WorkDir > "$Meta/TargetFiles.txt"
	High=$(wc -l < "$Meta/TargetFiles.txt")
	Count="1"
	mkdir "$Dir/Aligned/$i.aligned"
	while read FileLine; do
		echo "Aligning $FileLine with Muscle (File $Count/$High)"
		sed -i '/^$/d' "$WorkDir/$FileLine"
		muscle -quiet -in "$WorkDir/$FileLine" -out "$Dir/Aligned/$i.aligned/$FileLine"
		Count=$((Count+1))
	done < "$Meta/TargetFiles.txt"
done

#Add function to remove reference genome and copy to new location of haplotype diversity
#Add function to add metadata manually for sample?? Would be useful for microbiome work too

#ls "$Dir/Aligned/" > "$Meta/Aligned.txt"
#while read FileLine; do
#	sed -i 's/|.*/\n/' "$Dir/Aligned/$FileLine"
#done < "$Meta/Aligned.txt"

# ADD FUNCTION TO COUNT HOW MANY FILES HAVE BEEN ADDED TO CHECK WITH GENES IN TOTAL




#while [ "genes in file" -lt "some arbitrary number OR find out which file has most > in it and use that number"]
	#make a folder for least number of lines file
	#put all gene files in it that have X genes
	#theoretically, the folder with 42 genes should be the one to focus on
	#	Count=$((Count+1))
	#	ls "$Dir/Genes/" > "$Dir/Metadata/NumberedGeneList.txt"
	#	ls "$Dir/Proteins/" > "$Dir/Metadata/NumberedProteinList.txt"
	


	#CDP-diacylglycerol--serine O-phosphatidyltransferase (EC 2.7.8.8) CDS
	# - trunkation?






