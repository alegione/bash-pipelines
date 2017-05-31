#!/bin/bash
########################################
#
# Basic genome assembly bash pipeline written by Alistair Legione (alistair.legione@gmail.com) and Alyce Taylor-Brown
# - unzip reads
# - iterate reads through SPAdes
# - record quality of assembly through QUAST
#
########################################

#Set colour variables
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
PURPLE='\033[1;35m'
NOCOLOUR='\033[0m'

#ask user the name of the project for file name/directory purposes
#set this path to the target directory you want all files kept, from multiple pipeline runs. Be careful not to put a / at the end
Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please enter your home directory ${YELLOW}(ie: Documents/YOURNAME):${NOCOLOUR}"
	read -e Dir
	Dir=$(printf $Dir | sed 's/\/$//')
	echo $Dir
	if [ ! -d $Dir ]; then
		echo -e "${RED}No directory of that name exists${NOCOLOUR}"
		echo -e "Would you like to create it? (Y/N)${NOCOLOUR}"
		read -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ $yesno = "y" ]; then
			mkdir $Dir
			Switch=1
		fi
	else	
		echo -e "${BLUE}You entered: ${GREEN}$Dir${NOCOLOUR}"		
		Switch=1
	fi
done

# ask user the name of the project for file name/directory purposes
# echo prints text to the terminal, read pauses the terminal and awaits input from the user
# in this case the input will be saved to a variable called 'Project' and repeats it to the user
# Variables in bash are invoked with a $ at the beginning. Case matters, 'Project' and 'project' would be two different variables. 
# I typically use a capital for each word for clarity.

echo -e "${BLUE}Please enter a project title:${NOCOLOUR}"
read Project
echo -e "${BLUE}You entered: ${GREEN}$Project${NOCOLOUR}"


Dir="$Dir/$Project" #the location to store the project files

#File variables
Meta="$Dir/Metadata" #variable for metadata folder, used for shorter downstream coding
#Reads="$Dir/rawreads"
TrimOut="$Dir/trimmedoutput"
BowTieOut="$Dir/BAMfiles"
#file location for trimmomatic adapters
AdapterLocations="/home/genomics/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/"



#if the project directory doesn't exist, create project folders and initial sub-folders (metadata folder in this case)
if [ ! -d "$Dir" ]; then 
	mkdir $Dir
	mkdir $Meta
fi


#Asks for path to the target directory where the read files are (eg zipped read files)
Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please enter the file location of your reads ${RED}(your files MUST be paired end reads ending in a #.fastq):${NOCOLOUR}"
	read -e ReadDir
	ReadDir=$(printf $ReadDir | sed 's/\/$//')
	if [ -d $ReadDir ]; then
		Switch=1
		echo -e "${BLUE}You entered: ${GREEN}$ReadDir${NOCOLOUR}"
	else
		echo -e "${RED}Directory does not exist: ${GREEN}$ReadDir${NOCOLOUR}"
	fi
done

Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please enter number corresponding to your adapter sequences (1-3)${NOCOLOUR}"
	echo -e "${YELLOW}	1 - TruSeq3 - Single ends"
	echo -e "${YELLOW}	2 - TruSeq3 - Paired ends"
	echo -e "${YELLOW}	3 - TruSeq3 - Paired ends version 2"
	read -e -N Adapter
	if [ $Adapter -eq "1" ]; then
		Adapter="$AdapterLocation/TruSeq3-SE.fa"
		echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Single ends${NOCOLOUR}"
		Switch=1
	else if [ $Adapter -eq "2" ]; then
		Adapter="$AdapterLocation/TruSeq3-PE.fa"
		echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Paired ends${NOCOLOUR}"
		Switch=1
	else if [ $Adapter -eq "3" ]; then
		Adapter="$AdapterLocation/TruSeq3-PE-2.fa"
		echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Paired ends version 2${NOCOLOUR}"
		Switch=1
	else
		echo -e "${RED}Option does not exist: ${GREEN}$Adapter${NOCOLOUR}"
		echo -e "${RED}Please select one of the available options${NOCOLOUR}"
	fi
done

#Asks for path a reference genome(s)
Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please enter the file location of your reference genomes (they must be in .fasta format):${NOCOLOUR}"
	read -e Refdir
	if [ -d $Refdir ]; then
		Switch=1
		echo -e "${BLUE}You entered: ${GREEN}$RefGenome${NOCOLOUR}"
		echo -e "${BLUE}List of reference genomes created at ${GREEN}$Meta${NOCOLOUR}"
		ls $RefDir/*.fasta | tr '\n' '\0' | xargs -0 -n 1 basename | sed 's/\.fasta//g' > "$Meta/ReferenceList.txt"
	else
		echo -e "${RED}Directory does not exist: ${GREEN}$RefGenome${NOCOLOUR}"
	fi
done

while read Ref; do
	if [ ! -e "$Ref.1.bt2" ]; then
		echo -e "${BLUE}Building bowtie index from reference genome:${GREEN} $Ref ${NOCOLOUR}"
		bwa index $RefGenome
		bowtie2-build --threads 4 $Ref.fa
	fi	
done < "$Meta/ReferenceList.txt"


###########################################################################3

echo -e "${BLUE}Preparing metadata file (your files MUST be paired end reads ending in a #.fastq.gz)${NOCOLOUR}"
ls $ReadDir/*1.fastq.gz | tr '\n' '\0' | xargs -0 -n 1 basename | sed 's/_1\.fastq\.gz//' > "$Meta/GenomeList.txt" #lists all files ending in 1.fastq.gz (* is a wildcard) in ReadDir, removes the extension from the list and saves the list to a text file in the metadata folder

GenomeCount=$(wc -l < "$Meta/GenomeList.txt") #counts (wc) the number of lines (-l) in an input file (<). Putting the command within $() allows you to save the value to a variable, which we can then print to the terminal
echo -e "${BLUE}Genomes detected = ${GREEN}$GenomeCount${NOCOLOUR}"

if [ ! -e "$Meta/ReadCount.txt" ]; then
	> "$Meta/ReadCount.txt"
	> "$Meta/subsample.txt"
fi

TooLarge=0
if [ $GenomeCount -ne $((wc -l "$Meta/ReadCount.txt")) ]; then
	while read FileLine; do
		echo -e "${BLUE}Counting reads within ${GREEN}$FileLine ${BLUE}${NOCOLOUR}"
		ReadCount=$(zcat "$ReadDir/${FileLine}_1.fastq.gz" | echo $((`wc -l`/4)))
		echo "$FileLine = $ReadCount" >> "$Meta/ReadCount.txt"
		echo -e "		${YELLOW}$FileLine = $ReadCount${NOCOLOUR}"
		if [ $ReadCount -gt "20000000" ]; then
			echo $FileLine >> "$Meta/subsample.txt"
			TooLarge=1
		fi
	done <  "$Meta/GenomeList.txt"
fi


if [ TooLarge -eq "1" ]; then
	echo -e "${BLUE}What depth would you like to subsample reads at?${NOCOLOUR}"
	read Depth
	Depth=$((Depth*4))

	while read FileLine; do
		if [ ! -e "$ReadDir/${FileLine}_1.sub.fastq.gz" ]; then
			echo -e "${BLUE}Sampling forward reads within ${GREEN}$FileLine${BLUE} to depth of ${YELLOW}$Depth${BLUE} reads${NOCOLOUR}"
			zcat "$ReadDir/${FileLine}_1.fastq.gz" | head -$Depth | gzip > "$ReadDir/${FileLine}_1.sub.fastq.gz"
			echo -e "${BLUE}Sampling reverse reads within ${GREEN}$FileLine${BLUE} to depth of ${YELLOW}$Depth${BLUE} reads${NOCOLOUR}"
			zcat "$ReadDir/${FileLine}_2.fastq.gz" | head -$Depth | gzip > "$ReadDir/${FileLine}_2.sub.fastq.gz"
		fi
	done < "$Meta/subsample.txt"
fi



if [ ! -d $TrimOut ]; then
	mkdir $TrimOut
fi


Count="1"
while read FileLine; do
	if [ ! -e "$TrimOut/${FileLine}_1.fastq.gz" ]; then
		if [ -e "$ReadDir/${FileLine}_1.sub.fastq.gz" ]; then
			echo -e "${BLUE}Trimming ${GREEN}$Count ${BLUE}of ${GREEN}$GenomeCount ${BLUE}paired end reads with Trimmomatic${NOCOLOUR}"
			echo -e "${BLUE}Read1 = ${GREEN}"$ReadDir/$FileLine"_1.sub.fastq.gz${NOCOLOUR}"
			echo -e "${BLUE}Read2 = ${GREEN}"$ReadDir/$FileLine"_2.sub.fastq.gz${NOCOLOUR}"
			trimmomatic PE -phred33 $ReadDir/$FileLine"_1.sub.fastq.gz" $ReadDir/$FileLine"_2.sub.fastq.gz" $TrimOut/$FileLine"_1.fastq.gz" $TrimOut/$FileLine"_1_unpaired.fastq.gz" $TrimOut/$FileLine"_2.fastq.gz" $TrimOut/$FileLine"_2_unpaired.fastq.gz" ILLUMINACLIP:$Adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
			
		else
			echo -e "${BLUE}Trimming ${GREEN}$Count ${BLUE}of ${GREEN}$GenomeCount ${BLUE}paired end reads with Trimmomatic${NOCOLOUR}"
			echo -e "${BLUE}Read1 = ${GREEN}"$ReadDir/$FileLine"_1.fastq.gz${NOCOLOUR}"
			echo -e "${BLUE}Read2 = ${GREEN}"$ReadDir/$FileLine"_2.fastq.gz${NOCOLOUR}"
			trimmomatic PE -phred33 $ReadDir/$FileLine"_1.fastq.gz" $ReadDir/$FileLine"_2.fastq.gz" $TrimOut/$FileLine"_1.fastq.gz" $TrimOut/$FileLine"_1_unpaired.fastq.gz" $TrimOut/$FileLine"_2.fastq.gz" $TrimOut/$FileLine"_2_unpaired.fastq.gz" ILLUMINACLIP:$Adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
		fi
	fi
	Count=$((Count+1))
done < "$Meta/GenomeList.txt"

while read genome; do
	while read FileLine; do
		#bowtie2 on genome/reads using reference
		#samtools to bam? remove unmapped reads and de-interlace based on name		
		#point read variable to new reads and send it through the loop again
		
	done < "$Meta/ReferenceList.txt"
done < "$Meta/GenomeList.txt"


