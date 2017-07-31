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

#get screen size and maximise
LINE=`xrandr -q | grep Screen`
WIDTH=`echo ${LINE} | awk '{ print $8 }'`
HEIGHT=`echo ${LINE} | awk '{ print $10 }' | awk -F"," '{ print $1 }'`
echo -e "\e[4;$HEIGHT;${WIDTH}t"


#program location variables (if not in bin)
Spades="spades.py"
Quast="quast.py"
bwa="bwa"
trimmomatic="trimmomatic"
srst2="srst2"
#srst2="srst2.py"
coverageTableR="bin/CoverageTable.R"
coveragePlotsR="bin/CoveragePlots.R"


#Requires my fork of getmlst.py that produces an output folder containing all downloaded data
#getmlst="/home/genomics/Shared_folder/Alistair_Legione/Algetmlst.py"
getmlst="Algetmlst.py"


#ask user the name of the project for file name/directory purposes
#set this path to the target directory you want all files kept, from multiple pipeline runs.Be careful not to put a / at the end

if [ -e "$1/Metadata/Parameters.txt" ]; then
	echo -e "${BLUE}Parameter file detected...obtaining previously entered options${NOCOLOUR}"
	ParFile="$1/Metadata/Parameters.txt"
	Meta="$1/Metadata"
	if grep -i -q "Project" $ParFile; then Dir=$(grep -i "Project" $ParFile | cut -f2); echo -e "${GREEN}Project directory: $Dir${NOCOLOUR}";else Dir="nil"; fi
	if grep -i -q "Original reads" $ParFile; then ReadDir=$(grep -i "Original reads" $ParFile | cut -f2); echo -e "${GREEN}Reads directory: $ReadDir${NOCOLOUR}";else ReadDir="nil";fi
	if grep -i -q "Adapters" $ParFile; then Adapter=$(grep -i "Adapters" $ParFile | cut -f2); echo -e "${GREEN}Adapters to trim: $Adapter${NOCOLOUR}";else Adapters="nil";fi
	if grep -i -q "Reference genome" $ParFile; then RefGenome=$(grep -i "Reference genome" $ParFile | cut -f2); echo -e "${GREEN}Reference genome: $RefGenome${NOCOLOUR}";else RefGenome="nil";fi
	if grep -i -q "Prior contigs" $ParFile; then PriorContig=$(grep -i "Prior contigs" $ParFile | cut -f2); echo -e "${GREEN}Prior contigs: $PriorContig${NOCOLOUR}";else PriorContig="nil";fi
	if grep -i -q "Prior sanger" $ParFile; then PriorSanger=$(grep -i "Prior sanger" $ParFile | cut -f2); echo -e "${GREEN}Prior sanger sequence: $PriorSanger${NOCOLOUR}";else PriorSanger="nil";fi
	if grep -i -q "Blast database" $ParFile; then BlastDB=$(grep -i "Blast database" $ParFile | cut -f2); echo -e "${GREEN}Blast database: $BlastDB${NOCOLOUR}";else BlastDB="nil";fi
	if grep -i -q "Depth" $ParFile; then Depth=$(grep -i "Depth" $ParFile | cut -f2); echo -e "${GREEN}Subsampling depth: $Depth${NOCOLOUR}"; else Depth="nil";fi
	if grep -i -q "MLST database" $ParFile && grep -i -q "MLST definitions" $ParFile && grep -i -q "MLST delimiter" $ParFile && grep -i -q "MLST fasta" $ParFile; then
		MLSTDB=$(grep -i "MLST database" $ParFile | cut -f2)
		MLSTDelimiter=$(grep -i "MLST delimiter" $ParFile | cut -f2)
		MLSTfasta=$(grep -i "MLST fasta" $ParFile | cut -f2)
		MLSTdefinitions=$(grep -i "MLST definitions" $ParFile | cut -f2)
		echo -e "${GREEN}MLST database: $MLSTDB${NOCOLOUR}"
	else
		MLSTDB="nil"
		MLSTfasta="nil"
		MLSTdefinitions="nil"
		MLSTDelimiter="nil"
	fi
	sleep 1
else
	Dir="nil"
	ReadDir="nil"
	Adapters="nil"
	RefGenome="nil"
	PriorSanger="nil"
	PriorContig="nil"
	BlastDB="nil"
	Depth="nil"
	MLSTDB="nil"
	MLSTfasta="nil"
	MLSTdefinitions="nil"
	MLSTDelimiter="nil"
fi


if [ "$Dir" == "nil" ]; then
	
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter your home directory ${YELLOW}(ie: Documents/YOURNAME):${NOCOLOUR}"
		read -e HomeDir
		HomeDir=$(printf $HomeDir | sed 's/\/$//')
		echo $HomeDir
		if [ ! -d $HomeDir ]; then
			echo -e "${RED}No directory of that name exists${NOCOLOUR}"
			echo -e "Would you like to create it? (Y/N)${NOCOLOUR}"
			read -N 1 yesno
			yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
			if [ $yesno = "y" ]; then
				mkdir $HomeDir
				Switch=1
			fi
		else	
			echo -e "${BLUE}You entered: ${GREEN}$HomeDir${NOCOLOUR}"		
			Switch=1
		fi
	done

	echo -e "${BLUE}Please enter a project title:${NOCOLOUR}"
	read Project
	echo -e "${BLUE}You entered: ${GREEN}$Project${NOCOLOUR}"
	Dir="$HomeDir/$Project" #the location to store the project files
	Meta="$Dir/Metadata"
	ParFile="$Meta/Parameters.txt"
	if [ ! -d $Dir ]; then
		mkdir $Dir
	fi
	if [ ! -d $Meta ]; then	
		mkdir $Meta
	fi
	echo "Home	$HomeDir" >> $ParFile
	echo "Project	$Dir" >> $ParFile
fi



#File variables
if [ ! -z $Meta ]; then Meta="$Dir/Metadata"; fi
if [ ! -z $ParFile ]; then ParFile="$Meta/Parameters.txt"; fi
RawReads="$Dir/Original_reads"
nameconvert="$Dir/Metadata/rename_file.txt"
SpadeOut="$Dir/spades_output"
TrimOut="$Dir/trimmed_output"
bwaOut="$Dir/BWA_output"
coverageOut="$Dir/Coverage"
Pileup="$Dir/Pileup"
fastqcOut="$Dir/fastQC_output"
MLSTout="$Dir/MLST_output"
ResultsSum="$Dir/Results_Summary"
#file location for trimmomatic adapters
	AdapterLocation="/home/genomics/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters"
	#file location for SRST2 MLST databases
	SRST2MLST="/home/genomics/MLST_databases"
	#SRST2MLST="/home/qiime/MLST_databases"
#file location for SRST2 gene databases
	SRST2gene="/home/genomics/SRST2_gene_databases"
	#SRST2gene="/home/qiime/SRST2_gene_databases"


#if the project directory doesn't exist, create project folders and initial sub-folders (metadata folder in this case)
if [ ! -d "$Dir" ]; then 
	mkdir $Dir
fi

if [ ! -d "$Meta" ]; then 
	mkdir $Meta
fi

if [ ! -d "$SRST2MLST" ]; then 
	mkdir $SRST2MLST
fi

if [ ! -d "$ResultsSum" ]; then 
	mkdir $ResultsSum
fi


#Asks for path to the target directory where the read files are (eg zipped read files)
if [ "$ReadDir" == "nil" ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter the file location of your reads ${RED}(your files MUST be paired end reads ending in a #.fastq.gz):${NOCOLOUR}"
		read -e ReadDir
		ReadDir=$(printf $ReadDir | sed 's/\/$//')
		if [ -d $ReadDir ]; then
			Switch=1
			echo -e "${BLUE}You entered: ${GREEN}$ReadDir${NOCOLOUR}"
			echo -e "Original reads	$ReadDir" >> $ParFile
		else
			echo -e "${RED}Directory does not exist: ${GREEN}$ReadDir${NOCOLOUR}"
		fi
	done
fi

if [ "$Adapters" == "nil" ]; then
	echo "in adapter section"
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter number corresponding to your adapter sequences (1-3)${NOCOLOUR}"
		echo -e "${YELLOW}	1 - TruSeq3 - Single ends"
		echo -e "${YELLOW}	2 - TruSeq3 - Paired ends"
		echo -e "${YELLOW}	3 - TruSeq3 - Paired ends version 2"
		read -e -N 1 Adapter
		if [ $Adapter -eq "1" ]; then
			Adapter="$AdapterLocation/TruSeq3-SE.fa"
			echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Single ends${NOCOLOUR}"
			echo -e "Adapters	$Adapter" >> $ParFile
			Switch=1
		elif [ $Adapter -eq "2" ]; then
			Adapter="$AdapterLocation/TruSeq3-PE.fa"
			echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Paired ends${NOCOLOUR}"
			echo -e "Adapters	$Adapter" >> $ParFile
			Switch=1
		elif [ $Adapter -eq "3" ]; then
			Adapter="$AdapterLocation/TruSeq3-PE-2.fa"
			echo -e "${BLUE}You entered: ${GREEN}TruSeq3 - Paired ends version 2${NOCOLOUR}"
			echo -e "Adapters	$Adapter" >> $ParFile
			Switch=1
		else
			echo -e "${RED}Option does not exist: ${GREEN}$Adapter${NOCOLOUR}"
			echo -e "${RED}Please select one of the available options${NOCOLOUR}"
		fi
	done
fi


if [ "$RefGenome" == "nil" ]; then
	#Asks for path a reference genome
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to use a reference genome in downstream analysis (Y/N)?	${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			echo -e "${BLUE}Please enter the file location of your reference genome (in fasta format):${NOCOLOUR}"
			read -e RefGenome
			if [ -e $RefGenome ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$RefGenome${NOCOLOUR}"
				echo -e "Reference genome	$RefGenome" >> $ParFile
				if [ ! -e "$RefGenome.bwt" ]; then
					echo -e "${BLUE}Building bwa index from reference genome:${GREEN} $RefGenome ${NOCOLOUR}"
					bwa index $RefGenome
				fi				
			else
				echo -e "${RED}Directory does not exist: ${GREEN}$RefGenome${NOCOLOUR}"
			fi
		else	
			Switch=1
			RefGene="nil"
		fi
	done
fi


#Asks for path to previous contigs
Switch=0
if [ "$PriorContig" == "nil" ]; then
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to use previous contigs in downstream analysis (Y/N)?${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			echo -e "${BLUE}Please enter the directory location of your prior contigs:${NOCOLOUR}"
			read -e PriorContig
			PriorContig=$(printf $PriorContig | sed 's/\/$//')
			if [ -d $PriorContig ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$PriorContig${NOCOLOUR}"
				echo -e "Prior contigs	$PriorContig" >> $ParFile
				
			else
				echo -e "${RED}Directory does not exist: ${GREEN}$PriorContig${NOCOLOUR}"
			fi
		else	
			Switch=1
			PriorContig="none"
			echo -e "Prior contigs	$PriorContig" >> $ParFile
		fi
	done
fi
#Asks for path to previous Sanger reads
if [ "$PriorSanger" == "nil" ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to use previous sanger sequencing in downstream analysis (Y/N)?${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			echo -e "${BLUE}Please enter the directory location of your prior sanger sequencing:${NOCOLOUR}"
			read -e PriorSanger
			PriorSanger=$(printf $PriorSanger | sed 's/\/$//')
			if [ -d $PriorSanger ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$PriorSanger${NOCOLOUR}"
				echo -e "Prior sanger	$PriorSanger" >> $ParFile
			else
				echo -e "${RED}Directory does not exist: ${GREEN}$PriorSanger${NOCOLOUR}"
			fi
		else	
			Switch=1
			PriorSanger="none"
			echo -e "Prior sanger	$PriorSanger" >> $ParFile
		fi
	done
fi


#Asks for path to the blast database
if [ "$BlastDB" == "nil" ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to use a custom BLAST database in downstream analysis (Y/N)?${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			echo -e "${BLUE}Please enter the file location of your custom blast database (in fasta format):${NOCOLOUR}"
			read -e BlastDB
			if [ -e $BlastDB ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$BlastDB${NOCOLOUR}"
				echo -e "Blast database	$BlastDB" >> $ParFile
			else
				echo -e "${RED}Database does not exist: ${GREEN}$BlastDB${NOCOLOUR}"
			fi
		else	
			Switch=1
			BlastDB="none"
			echo -e "Blast database	$BlastDB" >> $ParFile
		fi
	done
fi


#file location for SRST2 MLST databases
#SRST2MLST="/home/genomics/MLST_databases"
#file location for SRST2 gene databases
#SRST2gene="/home/genomics/SRST2_gene_databases"
#MLSTDelimiter
#Asks for path to the MLST database
if [ "$MLSTDB" == "nil" ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to use a MLST database for short read analysis${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			if [ "$(ls -A "$SRST2MLST/")" ]; then			
				echo -e "${BLUE}Please type in the full path for one of the following databases, located at ${GREEN}$SRST2MLST${NOCOLOUR}"
				echo -e "${BLUE}If the database you enter isn't available we will attempt to download it${NOCOLOUR}"
				ls -d "$SRST2MLST/"*"/"
			else
				echo echo -e "${BLUE}Please enter an MLST database to download${NOCOLOUR}"
			fi
			read -e MLSTDB
			MLSTDB=$(printf $MLSTDB | sed 's/\/$//')
			if [ -d $MLSTDB ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$MLSTDB${NOCOLOUR}"
				echo -e "MLST database	$MLSTDB" >> $ParFile
				MLSTfasta=$(ls "$MLSTDB/"*".fasta")
				MLSTdefinitions=$(ls "$MLSTDB/"*".txt")
				MLSTDelimiter=$(grep -i "Suggested mlst delimiter" "$MLSTDB/"*".delimiter" | cut -f2)
				
				echo -e	"MLST fasta file	$MLSTfasta" >> $ParFile
				echo -e "MLST definitions file	$MLSTdefinitions" >> $ParFile
				echo -e "Suggested mlst delimiter	$MLSTDelimiter" >> $ParFile
			else
				echo -e "${RED}Database does not exist: ${GREEN}$MLSTDB${NOCOLOUR}"
				echo -e "${RED}Attempting download of ${GREEN}$MLSTDB${RED} with srst2 script getmlst.py${NOCOLOUR}"
				if [ -d "$SRST2MLST/tmp" ]; then 
					rm -r "$SRST2MLST/tmp"
					mkdir "$SRST2MLST/tmp"
				fi
				MLSTDB=$(printf $MLSTDB | sed 's!.*/!!')
				$getmlst --species "$MLSTDB" --output "$SRST2MLST/tmp/"
				if [ "$(ls -A "$SRST2MLST/tmp/")" ]; then
					newMLST=$(ls "$SRST2MLST/tmp/"*".fasta" | sed 's!.*/!!' | sed 's/.fasta//g')
					MLSTDB="$SRST2MLST/$newMLST"
					mv "$SRST2MLST/tmp" $MLSTDB
					MLSTfasta=$(ls "$MLSTDB/"*".fasta")
					MLSTdefinitions=$(ls "$MLSTDB/"*".txt")
					MLSTDelimiter=$(grep -i "Suggested mlst delimiter" "$MLSTDB/"*".delimiter" | cut -f2)
					echo -e "${BLUE}Successfully downloaded: ${GREEN}$MLSTfasta${NOCOLOUR}"
					
					echo -e "MLST database	$MLSTDB" >> $ParFile
					echo -e	"MLST fasta file	$MLSTfasta" >> $ParFile
					echo -e "MLST definitions file	$MLSTdefinitions" >> $ParFile
					echo -e "Suggested mlst delimiter	$MLSTDelimiter" >> $ParFile
					Switch=1
				else
					echo -e "${BLUE}No file of that name available for downloaded, please try one from this list${NOCOLOUR}"
					echo -e "${YELLOW}"
					cat "$SRST2MLST/MLST_schemes_list.txt"
					echo -e "${NOCOLOUR}"
				fi
			fi
		else	
			Switch=1
			MLSTDB="none"
			MLSTdefinitions="none"
			MLSTDelimiter="none"
			MLSTfasta="none"
			echo -e "MLST database	$MLSTDB" >> $ParFile
			echo -e "MLST definitions	$MLSTdefinitions" >> $ParFile
			echo -e "MLST delimiter	$MLSTDelimiter" >> $ParFile
			echo -e "MLST fasta	$MLSTfasta" >> $ParFile
		fi
	done
fi


if [ ! -e "$Meta/GenomeList.txt" ]; then
	echo -e "${BLUE}Preparing metadata file (your files MUST be paired end reads ending in a #.fastq.gz)${NOCOLOUR}"

	#lists all files ending in 1.fastq.gz (* is a wildcard) in ReadDir, removes the extension from the list and saves the list to a text file in the metadata folder
	ls $ReadDir/*1.fastq.gz | tr '\n' '\0' | xargs -0 -n 1 basename | sed 's/_1\.fastq\.gz//' > "$Meta/GenomeList.txt"
fi


spadeslist="$Meta/GenomeList.txt"


#count (wc) the number of lines (-l) in an input file (<). Putting the command within $() allows you to save the value to a variable, which we can then print to the terminal
GenomeCount=$(wc -l < "$Meta/GenomeList.txt") 




usingGenomeCount=$GenomeCount
echo -e "${BLUE}Genomes detected = ${GREEN}$GenomeCount${NOCOLOUR}"

if [ ! -e $nameconvert ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Would you like to rename your files before beginning${NOCOLOUR}"
		echo -e "${RED}PLEASE NOTE - THIS WILL CREATE NEW READ FILES, NOT OVERWRITE YOUR OLD FILES${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" == "y" ]; then
			echo -e "${BLUE}Please enter a location of a rename file (Must have one column with old file basenames (no extensions) and new file basenames (no extensions)${NOCOLOUR}"
			read -e renamefile
			if [ -e $renamefile ]; then
				Switch=1
				echo -e "${BLUE}You entered: ${GREEN}$renamefile${NOCOLOUR}"
				cp $renamefile $nameconvert
				mkdir $RawReads
				while read old new; do
					if [ -e "$ReadDir/${old}_1.fastq.gz" ]; then
						echo -e ${PURPLE}$(date)${NOCOLOUR}			
						echo -e "${BLUE}Moving ${GREEN}$old ${BLUE}to ${GREEN}$RawReads ${BLUE}and renaming to ${GREEN}$new${NOCOLOUR}"
						sed -i "s/$old/$new/g" "$Meta/GenomeList.txt"
						cp "$ReadDir/${old}_1.fastq.gz" "$RawReads/${new}_1.fastq.gz"
						cp "$ReadDir/${old}_2.fastq.gz" "$RawReads/${new}_2.fastq.gz"
					fi
				done < $nameconvert
				if grep -i -q "Original reads" $ParFile; then
					sed -i "s|$ReadDir|$RawReads|" $ParFile
				fi
				ReadDir=$RawReads
			else
				echo -e "${RED}File does not exist: ${GREEN}$renamefile${NOCOLOUR}"
			fi
		else	
			Switch=1
			echo -e "${BLUE}Files will retain current names throughout all analysis${NOCOLOUR}"
			echo -e "The original files have not been renamed before analyis" >> $nameconvert
		fi
	done
fi

if [ ! -d $TrimOut ]; then
	mkdir $TrimOut
	mkdir $TrimOut/paired
	mkdir $TrimOut/unpaired
	mkdir $fastqcOut
	mkdir $fastqcOut/Original
	mkdir $fastqcOut/Trimmed
fi


#fastqc if raw reads don't match fastqc results (this effectively skips it if equal, but does it all again if different, not ideal if only one new edition, perhaps a list of new genomes?
#	Could remove those already in fastqc folder from genome list with sed? and output a txt, remove \n and replace with ', ' and plug that into fastqc?
qcCount="0"
readCount=$(find "$ReadDir" -type f | wc -l)
qcCount=$(find "$fastqcOut/Original" -type f | wc -l)
readCount=$(echo "$readCount*2" | bc)

if [ $qcCount -eq "0" ]; then
	echo -e ${PURPLE}$(date)${NOCOLOUR}
	echo -e "${BLUE}Producing quality metrics for all raw reads ${NOCOLOUR}"
	fastqc -o "$fastqcOut/Original" --quiet -t 4 "$ReadDir/"*".fastq.gz"
elif [ $readCount -ne $qcCount ]; then
	while read i; do
		if [ ! -e "$fastqcOut/Original/${i}_1_fastqc.html" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Producing quality metrics for $i ${NOCOLOUR}"
			fastqc -o "$fastqcOut/Original" --quiet -t 4 "$ReadDir/${i}"*".fastq.gz"
		fi
	done < "$Meta/GenomeList.txt"
fi	
Count="1"
while read FileLine; do
	if [ ! -e "$TrimOut/paired/${FileLine}_1.fastq.gz" ]; then
		echo -e ${PURPLE}$(date)${NOCOLOUR}
		echo -e "${BLUE}Trimming ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}paired end reads with Trimmomatic${NOCOLOUR}"
		echo -e "${BLUE}Read1 = ${GREEN}"$ReadDir/$FileLine"_1.fastq.gz${NOCOLOUR}"
		echo -e "${BLUE}Read2 = ${GREEN}"$ReadDir/$FileLine"_2.fastq.gz${NOCOLOUR}"
		$trimmomatic PE -phred33 $ReadDir/$FileLine"_1.fastq.gz" $ReadDir/$FileLine"_2.fastq.gz" $TrimOut/paired/$FileLine"_1.fastq.gz" $TrimOut/unpaired/$FileLine"_1_unpaired.fastq.gz" $TrimOut/paired/$FileLine"_2.fastq.gz" $TrimOut/unpaired/$FileLine"_2_unpaired.fastq.gz" ILLUMINACLIP:$Adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi
	Count=$((Count+1))
done < "$Meta/GenomeList.txt"	

#fastqc if trimmed reads don't match fastqc results (this effectively skips it if equal, but does it all again if different, not ideal if only one new edition, perhaps a list of new genomes?
#	Could remove those already in fastqc folder from genome list with sed? and output a txt, remove \n and replace with ', ' and plug that into fastqc?
qcCount="0"
readCount=$(find "$TrimOut/paired/" -type f | wc -l)
qcCount=$(find "$fastqcOut/Trimmed" -type f | wc -l)
readCount=$(echo "$readCount*2" | bc)


if [ $qcCount -eq "0" ]; then
	echo -e ${PURPLE}$(date)${NOCOLOUR}
	echo -e "${BLUE}Producing quality metrics for all trimmed reads ${NOCOLOUR}"
	fastqc -o "$fastqcOut/Trimmed" --quiet -t 4 "$TrimOut/paired/"*".fastq.gz"
elif [ $readCount -ne $qcCount ]; then
	while read i; do
		if [ ! -e "$fastqcOut/Trimmed/${i}_1_fastqc.html" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Producing quality metrics for trimmed $i ${NOCOLOUR}"
			fastqc -o "$fastqcOut/Trimmed" --quiet -t 4 "$TrimOut/paired/${i}"*".fastq.gz"
		fi
	done < "$Meta/GenomeList.txt"
fi

if [ ! -e "$ResultsSum/trimmedreadtable.tsv" ]; then
	echo -e "Sample	Raw read count	Raw read GC%	Trimmed read count	Trimmed read GC%" > "$ResultsSum/trimmedreadtable.tsv"
fi

if [ ! -e "$Meta/genomesubsamplelist.txt" ]; then
	TooLarge="0"
	Count=0
else
	TooLarge="1"
	Count=$(wc -l < "$Meta/genomesubsamplelist.txt")
fi


echo -e ${PURPLE}$(date)${NOCOLOUR}
echo -e "${BLUE}Producing table of raw and trimmed read data${NOCOLOUR}"

#probably don't need this
echo -e "${YELLOW}Sample	Raw read count	Raw read GC%	Trimmed read count	Trimmed read GC% ${NOCOLOUR}"


tmp=$fastqcOut/tmp		
if [ ! -d $tmp ]; then
	mkdir $tmp
	mkdir $tmp/original
	mkdir $tmp/trimmed
fi

while read i; do 
	if ! grep -i -q $i "$ResultsSum/trimmedreadtable.tsv"; then

		unzip -q "$fastqcOut/Original/${i}_1_fastqc.zip" -d $tmp/original
		unzip -q "$fastqcOut/Trimmed/${i}_1_fastqc.zip" -d $tmp/trimmed
		originalreadcount=$(grep "Total Sequences" "$tmp/original/${i}_1_fastqc/fastqc_data.txt" | cut -f2)
		trimmedreadcount=$(grep "Total Sequences" "$tmp/trimmed/${i}_1_fastqc/fastqc_data.txt" | cut -f2)
		originalgc=$(grep "%GC" "$tmp/original/${i}_1_fastqc/fastqc_data.txt" | cut -f2)
		trimmedgc=$(grep "%GC" "$tmp/trimmed/${i}_1_fastqc/fastqc_data.txt" | cut -f2)
		if [ $trimmedreadcount -gt "20000000" ]; then
			echo $i >> "$Meta/genomesubsamplelist.txt"
			TooLarge=1
			Count=$((Count+1))
		fi
		echo -e "$i\t$originalreadcount\t$originalgc\t$trimmedreadcount\t$trimmedgc" >> "$ResultsSum/trimmedreadtable.tsv"
		echo -e "${YELLOW}$i\t$originalreadcount\t$originalgc\t$trimmedreadcount\t$trimmedgc${NOCOLOUR}"
	fi
done < "$Meta/GenomeList.txt"

if [ -d $tmp ]; then 
	rm -r $tmp
fi

if [ "$TooLarge" -eq "1" ]; then
	if [ ! -d "$Dir/LargeReads" ]; then
		mkdir "$Dir/LargeReads"
	fi
	if [ "$Depth" == "nil" ]; then
		echo -e "${BLUE}There were ${GREEN}$Count ${BLUE}read files found with >20,000,000 reads. It is recommended you subsample to at least this value${NOCOLOUR}"
		echo -e "${BLUE}What depth would you like to subsample reads to?${NOCOLOUR}"
		read Depth
		DepthOld=$Depth
		Depth=$((Depth*4))
		echo -e "Depth	$Depth" >> $ParFile
	else
		DepthOld=$((Depth/4))
	fi
	while read FileLine; do
		if [ ! -e "$Dir/LargeReads/${FileLine}_1.fastq.gz" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Large read files will be stored in ${GREEN}$Dir/LargeReads${NOCOLOUR}"			
			mv -f "$TrimOut/paired/${FileLine}_1.fastq.gz" "$Dir/LargeReads/${FileLine}_1.fastq.gz"
			mv -f "$TrimOut/paired/${FileLine}_2.fastq.gz" "$Dir/LargeReads/${FileLine}_2.fastq.gz"
		fi
		if [ ! -e "$TrimOut/paired/${FileLine}_2.fastq.gz" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Subsampled reads will move to ${GREEN}$TrimOut/paired/${NOCOLOUR}"			
			echo -e "${BLUE}Sampling forward reads within ${GREEN}$FileLine${BLUE} to depth of ${YELLOW}$DepthOld${BLUE} reads${NOCOLOUR}"
			zcat "$Dir/LargeReads/${FileLine}_1.fastq.gz" | head -$Depth | gzip > "$TrimOut/paired/${FileLine}_1.fastq.gz"
			echo -e "${BLUE}Sampling reverse reads within ${GREEN}$FileLine${BLUE} to depth of ${YELLOW}$DepthOld${BLUE} reads${NOCOLOUR}"
			zcat "$Dir/LargeReads/${FileLine}_2.fastq.gz" | head -$Depth | gzip > "$TrimOut/paired/${FileLine}_2.fastq.gz"
		fi
	done < "$Meta/genomesubsamplelist.txt"
fi

#	
#	Need to ask for what database at some point and if it is needed then download it
#		Probably list those available, with number selection? or download as secondary option?
#

if [ "$MLSTDB" != "none" ]; then
	Count=1
	if [ ! -d $MLSTout ]; then
		mkdir $MLSTout
		mkdir "$MLSTout/consensus_sequences"
	fi
	MLSTtype=$(echo "$MLSTfasta" | sed 's!.*/!!' | sed 's/.fasta//g')
	echo -e "${BLUE}MLST database = ${GREEN}$MLSTtype${NOCOLOUR}"
	while read i; do
		if [ ! -e "$MLSTout/$i/${i}__mlst__${MLSTtype}__results.txt" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Mapping MLST for ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}paired end reads with srst2${NOCOLOUR}"			
			mkdir "$MLSTout/$i"
			echo -e "${BLUE}MLST output will be ${GREEN}$MLSTout/$i/$i${NOCOLOUR}"
			echo -e "${BLUE}read 1 = ${GREEN}$TrimOut/paired/${i}_1.fastq.gz${NOCOLOUR}"
			echo -e "${BLUE}read 2 = ${GREEN}$TrimOut/paired/${i}_2.fastq.gz${NOCOLOUR}"

			echo $(date) > "$Meta/Progress.txt"
			echo -e "$srst2 --input_pe $TrimOut/paired/${i}_1.fastq.gz $TrimOut/paired/${i}_2.fastq.gz --mlst_db $MLSTfasta --mlst_definitions $MLSTdefinitions --mlst_delimiter $MLSTDelimiter --log --report_all_consensus --threads 4 --output $MLSTout/$i/$i"	> "$Meta/Progress.txt"
			$srst2 --input_pe "$TrimOut/paired/${i}_1.fastq.gz" "$TrimOut/paired/${i}_2.fastq.gz" --mlst_db $MLSTfasta --mlst_definitions "$MLSTdefinitions" --mlst_delimiter $MLSTDelimiter --log --report_all_consensus --threads 4 --output "$MLSTout/$i/$i"
			if [ -e "$MLSTout/$i/"*"all_consensus_alleles.fasta" ]; then
				cp "$MLSTout/$i/"*"all_consensus_alleles.fasta" "$MLSTout/consensus_sequences"
			fi
			if [ ! -e "$ResultsSum/MLSTcollated.tsv" ]; then
				cat "$MLSTout/$i/${i}__mlst__${MLSTtype}__results.txt" >> "$ResultsSum/MLSTcollated.tsv"
			else
				grep "$i" "$MLSTout/$i/${i}__mlst__${MLSTtype}__results.txt" >> "$ResultsSum/MLSTcollated.tsv"
			fi
			rm $MLSTDB/*.fai
			rm $MLSTDB/*.bt2
			
			#grep -v "^>" "$MLSTout/$i/$i.all_consensus_alleles.fasta" | awk -v var=$i 'BEGIN { ORS=""; print '>var\n' } { print }' > "$MLSTout/consensus_sequences/$i"
		fi
	Count=$((Count+1))
	done < "$Meta/GenomeList.txt"
fi

if [ ! -d $bwaOut ] && [ $RefGenome != "nil" ]; then
	mkdir $bwaOut
	mkdir $coverageOut
	mkdir $Pileup
fi

#Burrows wheeler alignment of all reads against a reference (only if a reference is provided)
Count="1"
if [ "$RefGenome" != "nil" ]; then
	while read i; do
		if [ ! -e "$Pileup/$i.vcf.gz" ]; then
			if [ ! -e "$bwaOut/$i.sort.bam" ]; then
				echo -e ${PURPLE}$(date)${NOCOLOUR}				
				echo -e "${BLUE}Aligning ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}paired end reads with burrows wheeler alignment${NOCOLOUR}"
				echo -e "${BLUE}Total paired end reads = ${GREEN}$(grep "$i" $ResultsSum/trimmedreadtable.tsv | cut -f4)${NOCOLOUR}"
				echo -e "${BLUE}Read1 = ${GREEN}"$TrimOut/paired/$i"_1.fastq.gz${NOCOLOUR}"
				echo -e "${BLUE}Read2 = ${GREEN}"$TrimOut/paired/$i"_2.fastq.gz${NOCOLOUR}"
				echo -e "${BLUE}Reference = ${GREEN}$RefGenome${NOCOLOUR}"
				bwa mem -v 0 -t 4 $RefGenome $TrimOut/paired/$i"_1.fastq.gz" $TrimOut/paired/$i"_2.fastq.gz" > "$bwaOut/$i.bwa.sam"
				echo -e ${PURPLE}$(date)${NOCOLOUR}			
				echo -e "${BLUE}Converting ${GREEN}$i.bwa.sam ${BLUE}to bam format${NOCOLOUR}"
				samtools view -bS "$bwaOut/$i.bwa.sam" > "$bwaOut/$i.bwa.bam"
				echo -e "${RED}Deleting ${GREEN}$i.bwa.sam ${RED}to save space${NOCOLOUR}"		
				rm "$bwaOut/$i.bwa.sam"
				
				echo -e ${PURPLE}$(date)${NOCOLOUR}	
				echo -e "${BLUE}Sorting ${GREEN}$i.bwa.bam ${NOCOLOUR}"
				samtools sort -@ 4 -o "$bwaOut/$i.sort.bam" --reference $RefGenome "$bwaOut/$i.bwa.bam"
				rm "$bwaOut/$i.bwa.bam"
			fi
			if [ ! -e "$coverageOut/$i.depth.txt" ]; then
				echo -e ${PURPLE}$(date)${NOCOLOUR}
				echo -e "${BLUE}Producing depth of coverage statistics of ${GREEN}$i ${BLUE}against ${GREEN}$RefGenome ${NOCOLOUR}"				
				genomeCoverageBed -ibam "$bwaOut/$i.sort.bam" -d -g $RefGenome > "$coverageOut/$i.depth.txt"
			fi
			if [ ! -e "$Pileup/$i.vcf.gz" ]; then
				echo -e ${PURPLE}$(date)${NOCOLOUR}
				echo -e "${BLUE}Variant calling (mPileUp) of ${GREEN}$i ${BLUE}against ${GREEN}$RefGenome ${NOCOLOUR}"						
				samtools mpileup -v -o "$Pileup/$i.vcf.gz" -f $RefGenome "$bwaOut/$i.sort.bam"
			fi
		fi
		Count=$((Count+1))
	done < "$Meta/GenomeList.txt"
	if [ ! -e "$ResultsSum/CoverageTable.tsv" ]; then
		echo -e "Name\tReferenceLength\tAverageDepth\tMinDepth\tQ1\tMedianDepth\tQ3\tMaxDepth\tCoverage\tGoodCoverage\tGC" > "$ResultsSum/CoverageTable.tsv"
	fi

	while read i; do
		if ! grep -i -q $i "$ResultsSum/CoverageTable.tsv"; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Producing coverage statistics for ${GREEN}$i${NOCOLOUR}"
			Rscript $coverageTableR $Dir $i
			goodcov=$(grep $i "$ResultsSum/CoverageTable.tsv" | cut -f10)
			echo "$goodcov"
			goodcov=$(echo "($goodcov+0.5)/1" | bc)
			echo "goodcov rounded? = $goodcov"
			if [ "$goodcov" -ge "90" ]; then #should make this threshold an option
				echo $i >> "$Meta/denovoGenomeList.txt"
			fi
		fi
	done < "$Meta/GenomeList.txt"
	if [ ! -d $ResultsSum/plots ]; then	
		mkdir $ResultsSum/plots	
	fi
	while read i; do
		if [ ! -e "$ResultsSum/plots/$i.tiff" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Producing coverage plot for ${GREEN}$i${NOCOLOUR}"
			Rscript "$coveragePlotsR" $Dir $i
		fi
	done < "$Meta/GenomeList.txt"
	denovoGenomeCount=0
	if [ -e "$Meta/denovoGenomeList.txt" ]; then
		echo -e "${BLUE}In total there are ${GREEN}$denovoGenomeCount ${BLUE}genomes with >=90% coverage${NOCOLOUR}"
		denovoGenomeCount=$(wc -l < "$Meta/denovoGenomeList.txt")
		usingGenomeCount=$denovoGenomeCount
		spadeslist="$Meta/denovoGenomeList.txt"
		echo -e "Original genomes\t$GenomeCount" >> "$ResultsSum/goodcoverage.txt"
		echo -e "Genomes used in spades\t$denovoGenomeCount" >> "$ResultsSum/goodcoverage.txt"
	else
		echo -e "${RED}In total there are ${GREEN}$denovoGenomeCount ${RED}genomes with >=90% coverage${NOCOLOUR}"
		echo -e "${RED}Would you like to continue with SPAdes assembly ${GREEN}(Y/N)${RED}?\nPOOR COVERAGE WILL RESULT IN POOR CONTIGS${NOCOLOUR}"
		read -e -N 1 yesno
		yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
		if [ "$yesno" != "y" ]; then
			exit
		fi
	fi

		

	sleep 2
fi


if [ ! -d $SpadeOut ]; then
	mkdir $SpadeOut
	mkdir "$Dir/contigs"
	mkdir "$Dir/graphs"
	mkdir "$Dir/scaffolds"
fi

Count="1"
while read FileLine; do
	if [ ! -d "$SpadeOut/$FileLine" ]; then
		echo -e ${PURPLE}$(date)${NOCOLOUR}
		echo -e "${BLUE}Assembling ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}paired end reads with SPAdes${NOCOLOUR}"
		echo -e "${BLUE}Total paied end reads = ${GREEN}$(grep "$FileLine" $ResultsSum/trimmedreadtable.tsv | cut -f4)${NOCOLOUR}"
		echo -e "${BLUE}Read1 = ${GREEN}"$TrimOut/paired/$FileLine"_1.fastq.gz${NOCOLOUR}"
		echo -e "${BLUE}Read2 = ${GREEN}"$TrimOut/paired/$FileLine"_2.fastq.gz${NOCOLOUR}"
		mkdir "$SpadeOut/$FileLine"

		if [ $PriorSanger != "none" ] && [ $PriorContig  != "none" ]; then			
			if [ -e "$PriorSanger/$FileLine.fasta" ] && [  -e "$PriorContig/$FileLine.CLC.fa" ]; then	
				echo -e "${BLUE}Using prior sanger sequencing and untrusted contigs${NOCOLOUR}"				
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful --sanger "$PriorSanger/$FileLine.fasta" --untrusted-contigs "$PriorContig/$FileLine.CLC.fa" -o "$SpadeOut/$FileLine"
			elif [ -e "$PriorSanger/$FileLine.fasta" ]; then
				echo -e "${BLUE}Using prior sanger sequencing${NOCOLOUR}"			
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful --sanger "$PriorSanger/$FileLine.fasta" -o "$SpadeOut/$FileLine"
			elif [ -e "$PriorContig/$FileLine.CLC.fa" ]; then
				echo -e "${BLUE}Using prior untrusted contigs${NOCOLOUR}"			
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful --untrusted-contigs "$PriorContig/$FileLine.CLC.fa" -o "$SpadeOut/$FileLine"
			else
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful -o "$SpadeOut/$FileLine"
			fi
		elif [ $PriorSanger  != "none" ]; then
			if [ -e "$PriorSanger/$FileLine.fasta" ]; then
				echo -e "${BLUE}Using prior Sanger sequencing${NOCOLOUR}"
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful --sanger "$PriorSanger/$FileLine.fasta" -o "$SpadeOut/$FileLine"
			else
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful -o "$SpadeOut/$FileLine"
			fi
		elif [ $PriorContig  != "none" ]; then
			if [ -e "$PriorContig/$FileLine.CLC.fa" ]; then
				echo -e "${BLUE}Using prior untrusted contigs${NOCOLOUR}"
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful --untrusted-contigs "$PriorContig/$FileLine.CLC.fa" -o "$SpadeOut/$FileLine"
			else
				$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful -o "$SpadeOut/$FileLine"
			fi
		else
			$Spades -t 4 -1 $TrimOut/paired/$FileLine"_1.fastq.gz" -2 $TrimOut/paired/$FileLine"_2.fastq.gz" -k 21,33,55,67 --careful -o "$SpadeOut/$FileLine"
		fi


		echo -e ${PURPLE}$(date)${NOCOLOUR}
		cp "$SpadeOut/$FileLine/assembly_graph.fastg" "$Dir/graphs/$FileLine.assembly_graph.fastg"
		mkdir "$SpadeOut/$FileLine/contig_summary"
# Moved the blast step to within the spades assembly loop so genomes are completed as the loop continues
		echo -e "${BLUE}comparing contigs of $SpadeOut/$FileLine to custom database:${NOCOLOUR}"
		echo -e "${GREEN}$BlastDB${NOCOLOUR}"
		
		blastn -query "$SpadeOut/$FileLine/contigs.fasta" -db $BlastDB -outfmt 6 -out "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.blastn.tab"
		get_top_hit_from_blast.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab"
		blast_summary.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.blastn.summary.txt"
		get_contig_name_only_from_blast.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.names.tab"
		#can probably just cut output from cpec blast results rather than running script?
		for i in `less $SpadeOut/$FileLine/contig_summary/$FileLine.names.tab`; do
			fasta_get_keyword_matches.pl "$SpadeOut/$FileLine/contigs.fasta" $i >> "$SpadeOut/$FileLine/contig_summary/$FileLine.contigs.fasta"
		echo -e "${YELLOW}"	
		head -5 "$SpadeOut/$FileLine/contig_summary/$FileLine.blastn.summary.txt"
		echo -e "${NOCOLOUR}"
		done
		cp "$SpadeOut/$FileLine/contig_summary/$FileLine.contigs.fasta" "$Dir/contigs/$FileLine.contigs.fasta"
	fi
	Count=$((Count+1))
done < $spadeslist


while read FileLine; do
	if [ -d "$SpadeOut/$FileLine" ]; then
		if [ -e "$SpadeOut/$FileLine/contigs.fasta" ] && [ ! -e "$Dir/contigs/$FileLine.contigs.fasta" ]; then
			echo -e "${BLUE}comparing contigs of $SpadeOut/$FileLine to custom database:${NOCOLOUR}"
			echo -e "${GREEN}$BlastDB${NOCOLOUR}"
			
			blastn -query "$SpadeOut/$FileLine/contigs.fasta" -db $BlastDB -outfmt 6 -out "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.blastn.tab"
			get_top_hit_from_blast.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab"
			blast_summary.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.blastn.summary.txt"
			get_contig_name_only_from_blast.pl "$SpadeOut/$FileLine/contig_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/contig_summary/$FileLine.names.tab"
			#can probably just cut output from cpec blast results rather than running script?
			for i in `less $SpadeOut/$FileLine/contig_summary/$FileLine.names.tab`; do
				fasta_get_keyword_matches.pl "$SpadeOut/$FileLine/contigs.fasta" $i >> "$SpadeOut/$FileLine/contig_summary/$FileLine.contigs.fasta"
			echo -e "${YELLOW}"	
			head -5 "$SpadeOut/$FileLine/contig_summary/$FileLine.blastn.summary.txt"
			echo -e "${NOCOLOUR}"
			done
			cp "$SpadeOut/$FileLine/contig_summary/$FileLine.contigs.fasta" "$Dir/contigs/$FileLine.contigs.fasta"
		fi
		if [ -e "$SpadeOut/$FileLine/scaffolds.fasta" ] && [ ! -d "$SpadeOut/$FileLine/scaffold_summary" ]; then
			echo -e ${PURPLE}$(date)${NOCOLOUR}			
			mkdir "$SpadeOut/$FileLine/scaffold_summary"
			echo -e "${BLUE}comparing scaffolds of $SpadeOut/$FileLine to custom database:${NOCOLOUR}"
			echo -e "${GREEN}$BlastDB${NOCOLOUR}"
			blastn -query "$SpadeOut/$FileLine/scaffolds.fasta" -db $BlastDB -outfmt 6 -out "$SpadeOut/$FileLine/scaffold_summary/$FileLine.ref.blastn.tab"
			get_top_hit_from_blast.pl "$SpadeOut/$FileLine/scaffold_summary/$FileLine.ref.blastn.tab" > "$SpadeOut/$FileLine/scaffold_summary/$FileLine.ref.TOP.blastn.tab"
			blast_summary.pl "$SpadeOut/$FileLine/scaffold_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/scaffold_summary/$FileLine.blastn.summary.txt"
			echo -e "${YELLOW}"	
			head -5 "$SpadeOut/$FileLine/scaffold_summary/$FileLine.blastn.summary.txt"
			echo -e "${NOCOLOUR}"			

			get_contig_name_only_from_blast.pl "$SpadeOut/$FileLine/scaffold_summary/$FileLine.ref.TOP.blastn.tab" > "$SpadeOut/$FileLine/scaffold_summary/$FileLine.names.tab"
			#can probably just cut output from cpec blast results rather than running script?
			for i in `less $SpadeOut/$FileLine/scaffold_summary/$FileLine.names.tab`; do
				fasta_get_keyword_matches.pl "$SpadeOut/$FileLine/scaffolds.fasta" $i >> "$SpadeOut/$FileLine/scaffold_summary/$FileLine.scaffolds.fasta"
			done
			cp "$SpadeOut/$FileLine/scaffold_summary/$FileLine.scaffolds.fasta" "$Dir/scaffolds/$FileLine.scaffolds.fasta"
		fi
	fi
done < $spadeslist





Count="1"

if [ ! -d "$Dir/quastoutput" ]; then
	mkdir "$Dir/quastoutput"
	if [ "$RefGenome" == "nil" ]; then	
		while read FileLine; do
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Analysing assembly ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}with Quast${NOCOLOUR}"
			$Quast -f -o "$Dir/quastoutput/$FileLine" "$Dir/contigs/$FileLine.contigs.fasta"
			$Quast -f -o "$Dir/quastoutput/$FileLine.finalQuast" "$SpadeOut/$FileLine/contigs.fasta"
			Count=$((Count+1))
		done < $spadeslist
	else
		while read FileLine; do
			echo -e ${PURPLE}$(date)${NOCOLOUR}
			echo -e "${BLUE}Analysing assembly ${GREEN}$Count ${BLUE}of ${GREEN}$usingGenomeCount ${BLUE}with Quast${NOCOLOUR}"
			echo -e "${BLUE}Reference = ${GREEN}$RefGenome${NOCOLOUR}"
			$Quast -f -o "$Dir/quastoutput/$FileLine" -R $RefGenome "$Dir/contigs/$FileLine.contigs.fasta"
			$Quast -f -o "$Dir/quastoutput/$FileLine.finalQuast" -R $RefGenome "$SpadeOut/$FileLine/contigs.fasta"
			Count=$((Count+1))
		done < $spadeslist
	fi
	echo -e ${PURPLE}$(date)${NOCOLOUR}
	echo -e "${BLUE}Producing quast summary reports${NOCOLOUR}"
	while read i; do
		
		if [ -e "$Dir/quastoutput/$i/contigs_reports/contigs_report_$i.contigs.mis_contigs.info" ]; then
			echo "$i" >> "$ResultsSum/Misassemblies_details.txt"
			cat "$Dir/quastoutput/$i/contigs_reports/contigs_report_$i.contigs.mis_contigs.info" >> "$ResultsSum/Misassemblies_details.txt"
			echo -e "\n===================================\n" >> "$ResultsSum/Misassemblies_details.txt"
		fi
		if [ ! -e "$ResultsSum/Misassemblies_summary.tsv" ]; then
			cat "$Dir/quastoutput/$i/contigs_reports/transposed_report_misassemblies.tsv" > "$ResultsSum/Misassemblies_summary.tsv"
		else
			grep $i "$Dir/quastoutput/$i/contigs_reports/transposed_report_misassemblies.tsv" >> "$ResultsSum/Misassemblies_summary.tsv"
		fi
		if [ ! -e "$ResultsSum/quast_summary.tsv" ]; then	
		cat "$Dir/quastoutput/$i/transposed_report.tsv" > "$ResultsSum/quast_summary.tsv"
		else
			grep $i "$Dir/quastoutput/$i/transposed_report.tsv" >> "$ResultsSum/quast_summary.tsv"
		fi
		columns=$(awk '{print NF}' "$Dir/quastoutput/$i/contigs_reports/alignments_$i.contigs.tsv" | sort -nu | tail -n 1)
		columns=$((columns-1))
		gaps=$(wc -l < "$Dir/quastoutput/$i/genome_stats/$i.contigs_gaps.txt")
		gaps=$((gaps-1))
		if [ ! -e "$ResultsSum/contig_coverage.tsv" ]; then
			echo -e "Sample\tcontigs aligned to reference\tGaps detected" > "$ResultsSum/contig_coverage.tsv"
			echo -e "$i\t$columns\t$gaps" >> "$ResultsSum/contig_coverage.tsv"
		else
			echo -e "$i\t$columns\t$gaps" >> "$ResultsSum/contig_coverage.tsv"
		fi
		if [ ! -d "$Dir/quastoutput/gapfiles" ]; then
			mkdir "$Dir/quastoutput/gapfiles"
		fi
		cp "$Dir/quastoutput/$i/genome_stats/$i.contigs_gaps.txt" "$Dir/quastoutput/gapfiles/$i.contigs_gaps.txt"
	done < $spadeslist
	#this needs editing to remove contig name from each file
	while read i; do
		if [ -e "$Dir/quastoutput/$i.finalQuast/contigs_reports/contigs_report_contigs.mis_contigs.info" ]; then
			echo "$i" >> "$ResultsSum/Misassemblies_details_notfiltered.txt"
			cat "$Dir/quastoutput/$i.finalQuast/contigs_reports/contigs_report_contigs.mis_contigs.info" >> "$ResultsSum/Misassemblies_details_notfiltered.txt"
			echo -e "\n===================================\n" >> "$ResultsSum/Misassemblies_details_notfiltered.txt"
		fi
		if [ ! -e "$ResultsSum/Misassemblies_summary_notfiltered.tsv" ]; then
			cat "$Dir/quastoutput/$i.finalQuast/contigs_reports/transposed_report_misassemblies.tsv" > "$ResultsSum/Misassemblies_summary_notfiltered.tsv"
		else
			grep $i "$Dir/quastoutput/$i.finalQuast/contigs_reports/transposed_report_misassemblies.tsv" >> "$ResultsSum/Misassemblies_summary_notfiltered.tsv"
		fi
		if [ ! -e "$ResultsSum/quast_summary_notfiltered.tsv" ]; then	
		cat "$Dir/quastoutput/$i.finalQuast/transposed_report.tsv" > "$ResultsSum/quast_summary_notfiltered.tsv"
		else
			grep $i "$Dir/quastoutput/$i.finalQuast/transposed_report.tsv" >> "$ResultsSum/quast_summary_notfiltered.tsv"
		fi
		columns=$(awk '{print NF}' "$Dir/quastoutput/$i.finalQuast/contigs_reports/alignments_contigs.tsv" | sort -nu | tail -n 1)
		columns=$((columns-1))
		gaps=$(wc -l < "$Dir/quastoutput/$i.finalQuast/genome_stats/contigs_gaps.txt")
		gaps=$((gaps-1))
		if [ ! -e "$ResultsSum/contig_coverage_notfiltered.tsv" ]; then
			echo -e "Sample\tcontigs aligned to reference\tGaps detected" > "$ResultsSum/contig_coverage_notfiltered.tsv"
			echo -e "$i\t$columns\t$gaps" >> "$ResultsSum/contig_coverage_notfiltered.tsv"
		else
			echo -e "$i\t$columns\t$gaps" >> "$ResultsSum/contig_coverage_notfiltered.tsv"
		fi
		if [ ! -d "$Dir/quastoutput/gapfiles_unfilteredcontigs" ]; then
			mkdir "$Dir/quastoutput/gapfiles_unfilteredcontigs"
		fi
		cp "$Dir/quastoutput/$i.finalQuast/genome_stats/contigs_gaps.txt" "$Dir/quastoutput/gapfiles_unfilteredcontigs/$i.contigs_gaps.txt"
	done < $spadeslist
fi
echo -e "${BLUE}Assembly pipeline complete${NOCOLOUR}"

#prokka to annotate genes
#if reference provided, map to reference and concatenate? (perhaps do this before prokka?)
#move gff files to stand alone directory
#run roary without deleting final files?
#build trees for individual genes from roary?
#run codeml on individual genes from roary using alignment and trees (can use fast tree or some other max likilihood?)
	#need to convert gene files from fasta to paml compatible
