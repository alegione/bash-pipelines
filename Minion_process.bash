#!/bin/bash

#get CPU and RAM metrics
CORES=$(nproc)
RAM=$(free -g | tr -s "[:space:]" "\t" | cut -f11)

#Set colour variables
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
PURPLE='\033[1;35m'
NOCOLOUR='\033[0m'
#USE THIS TO STRIP COLOUR FROM FILES: sed 's/\x1b\[[0-9;]*m//g'

#get screen size and maximise
LINE=`xrandr -q | grep Screen`
WIDTH=`echo ${LINE} | awk '{ print $8 }'`
HEIGHT=`echo ${LINE} | awk '{ print $10 }' | awk -F"," '{ print $1 }'`
echo -e "\e[4;$HEIGHT;${WIDTH}t"

#fresh Screen
clear

if [ -e "$1/Metadata/Parameters.txt" ]; then
	echo -e "${BLUE}Parameter file detected...obtaining previously entered options${NOCOLOUR}"
	ParFile="$1/Metadata/Parameters.txt"
	Meta="$1/Metadata"

	if grep -i -q "User" $ParFile; then HomeDir=$(grep -i "User" $ParFile | cut -f2); echo -e "${GREEN}User: $HomeDir${NOCOLOUR}";else HomeDir="nil"; fi
	if grep -i -q "Project" $ParFile; then Project=$(grep -i "Project" $ParFile | cut -f2); echo -e "${GREEN}Project directory: $Project${NOCOLOUR}";else Project="nil"; fi
	if grep -i -q "Original reads" $ParFile; then DIR_RawReads=$(grep -i "Original reads" $ParFile | cut -f2); echo -e "${GREEN}Raw reads directory: $DIR_RawReads${NOCOLOUR}";else DIR_RawReads="nil";fi
	if grep -i -q "Called reads" $ParFile; then DIR_calledreads=$(grep -i "Called reads" $ParFile | cut -f2); echo -e "${GREEN}Called reads directory: $DIR_calledreads${NOCOLOUR}";else DIR_calledreads="nil";fi

	if grep -i -q "Flow cell" $ParFile; then Minion_Flow=$(grep -i "Flow cell" $ParFile | cut -f2); echo -e "${GREEN}Flow cell: $Minion_Flow${NOCOLOUR}";else Minion_Flow="nil";fi
	if grep -i -q "Kit" $ParFile; then Minion_Kit=$(grep -i "Kit" $ParFile | cut -f2); echo -e "${GREEN}Kit used: $Minion_Kit${NOCOLOUR}";else Minion_Kit="nil";fi
	if grep -i -q "Direction" $ParFile; then Direction=$(grep -i "Direction" $ParFile | cut -f2); echo -e "${GREEN}Direction: $Direction${NOCOLOUR}";else Direction="nil";fi
	if grep -i -q "Product" $ParFile; then Product=$(grep -i "Product" $ParFile | cut -f2); echo -e "${GREEN}Product: $Product${NOCOLOUR}";else Product="nil";fi

#	if grep -i -q "Sample names" $ParFile; then Species=$(grep -i "Sample names" $ParFile | cut -f2);echo -e "${GREEN}Sample prefix: $Species${NOCOLOUR}";else Species="nil";fi
#	if grep -i -q "Diversity profile target" $ParFile; then divprotarget=$(grep -i "Diversity profile target" $ParFile | cut -f2);echo -e "${GREEN}Diversity profile target: $divprotarget${NOCOLOUR}";else divprotarget="nil";fi
#	if grep -i -q "Adapters" $ParFile; then Adapter=$(grep -i "Adapters" $ParFile | cut -f2); echo -e "${GREEN}Adapters to trim: $Adapter${NOCOLOUR}";else Adapters="nil";fi
#	if grep -i -q "Taxonomy database" $ParFile; then taxa=$(grep -i "Taxonomy database" $ParFile | cut -f2); echo -e "${GREEN}Taxonomy database: $taxa${NOCOLOUR}";else taxa="nil";fi
#	if grep -i -q "Rarefraction Low" $ParFile; then taxa=$(grep -i "Rarefraction Low" $ParFile | cut -f2); echo -e "${GREEN}Rarefraction Low: $rarefractionlow${NOCOLOUR}";else rarefractionlow="nil";fi
#	if grep -i -q "Rarefraction High" $ParFile; then taxa=$(grep -i "Rarefraction High" $ParFile | cut -f2); echo -e "${GREEN}Rarefraction High: $rarefractionhigh${NOCOLOUR}";else rarefractionhigh="nil";fi
#	if grep -i -q "Iterations" $ParFile; then iterations=$(grep -i "Iterations" $ParFile | cut -f2); echo -e "${GREEN}Iterations: $iterations${NOCOLOUR}";else iterations="nil";fi
#	if grep -i -q "Depth" $ParFile; then Depth=$(grep -i "Depth" $ParFile | cut -f2); echo -e "${GREEN}Subsampling depth: $Depth${NOCOLOUR}"; else Depth="nil";fi
#	if grep -i -q "Keep Singletons? (Y/N=1/0)" $ParFile; then singles=$(grep -i "Keep Singletons? (Y/N=1/0)" $ParFile | cut -f2); echo -e "${GREEN}Keep Singletons? (Y/N=1/0): $singles${NOCOLOUR}"; else singles="nil";fi
#	if grep -i -q "Minimum length for clustering" $ParFile; then Truncate=$(grep -i "Minimum length for clustering" $ParFile | cut -f2); echo -e "${GREEN}Minimum length for clustering: $Depth${NOCOLOUR}"; else Truncate="nil";fi
# forward primer
# reverse primer (could also be in mapping file)
# diversity profiling target (ITS or 16S)

	sleep 1
else
	HomeDir="nil"
	Project="nil"
	DIR_RawReads="nil"
	DIR_calledreads="nil"
	Minion_Flow="nil"
	Minion_Kit="nil"
#	Species="nil"
#	divprotarget="nil"
#	Adapters="nil"
#	taxa="nil"
#	rarefractionlow="nil"
#	rarefractionhigh="nil"
#	iterations="nil"
#	Depth="nil"
#	singles="nil"
#	Truncate="nil"
fi

if [ "$HomeDir" == "nil" ]; then

	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter your home directory ${YELLOW}(eg: Student/YOURNAME):${NOCOLOUR}"
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
fi


if [ "$Project" == "nil" ]; then
	echo -e "${BLUE}Please enter a project title:${NOCOLOUR}"
	read Project
	echo -e "${BLUE}You entered: ${GREEN}$Project${NOCOLOUR}"
fi



Dir="$HomeDir/$Project"
Progress="$Dir/Progress.txt"

if [ ! -d $Dir ]; then
	mkdir $Dir
fi

if [ -z $Meta ]; then Meta="$Dir/Metadata"; fi

if [ ! -d $Meta ]; then
	mkdir $Meta
fi

if [ -z $ParFile ]; then ParFile="$Meta/Parameters.txt"; fi

if [ ! -e $ParFile ]; then touch $ParFile; fi

if ! grep -i -q "User" $ParFile; then echo -e "User	$HomeDir" >> $ParFile; fi

if ! grep -i -q "Project" $ParFile; then echo -e "Project	$Project" >> $ParFile; fi

#File variables

#DIR_RawReads="$Dir/Original_reads"
DIR_calledreads="$Dir/Called_reads"
DIR_TrimmedReads="$Dir/PoreChop"
DIR_FilteredReads="$Dir/FiltLong"
DIR_ResultsSum="$Dir/Results_Summary"
DIR_canu="$Dir/canu"
#Minion_Kit="SQK-DCS108"
WorkerThreads="16"

if [ ! -d $DIR_ResultsSum ]; then
	mkdir $DIR_ResultsSum
fi

Switch=0
if [ $Minion_Flow == "nil" ]; then
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter the number corresponding to the flow cell in use:${NOCOLOUR}"
		echo -e "${YELLOW}1 - FLO-MIN106 (R9.4)${NOCOLOUR}"
		echo -e "${YELLOW}2 - FLO-MIN107 (R9.5)${NOCOLOUR}"
		read -e Minion_Flow
		case $Minion_Flow in
			1)
				echo -e "${GREEN}Flow cell FLO-MIN106 selected${NOCOLOUR}" | tee -a $Progress
				Minion_Flow="FLO-MIN106"
				Switch=1
				;;
			2)
				echo -e "${GREEN}Flow cell FLO-MIN107 selected${NOCOLOUR}" | tee -a $Progress
				Minion_Flow="FLO-MIN107"
				Switch=1
				;;
			*)
				echo -e "${RED}ERROR: Please only type 1 or 2${NOCOLOUR}"
				;;
		esac
	done
fi
if ! grep -i -q "Flow cell" $ParFile; then echo -e "Flow cell	$Minion_Flow" >> $ParFile; fi

# ADD IN KIT SELECTION CASES
Switch=0
if [ $Minion_Kit == "nil" ]; then
	echo -e "${BLUE}Please enter the number corresponding to the kit in use:${NOCOLOUR}"
	echo -e "${BLUE}DNA${NOCOLOUR}"
	echo -e "${YELLOW}1 - Ligation sequencing kit 1D${GREEN}\tSQK-LSK108${NOCOLOUR}"
	echo -e "${YELLOW}2 - Ligation sequencing kit${GREEN}\tSQK-LSK109${NOCOLOUR}"
	echo -e "${YELLOW}3 - Rapid sequencing kit${GREEN}\tSQK-RAD004${NOCOLOUR}"
	echo -e "${YELLOW}4 - PCR sequencing kit${GREEN}\tSQK-PSK004${NOCOLOUR}"
	echo -e "${YELLOW}5 - 1D^2 sequencing kit${GREEN}\tSQK-LSK308${NOCOLOUR}"
	echo -e "${YELLOW}6 - PCR Barcoding kit${GREEN}\tSQK-PBK004${NOCOLOUR}"
	echo -e "${YELLOW}7 - Rapid PCR Barcoding kit${GREEN}\tSQK-RPB004${NOCOLOUR}"
	echo -e "${YELLOW}8 - 16S Barcoding kit${GREEN}\tSQK-RAB204${NOCOLOUR}"
	echo -e "${BLUE}RNA${NOCOLOUR}"
	echo -e "${YELLOW}9 - Direct cDNA sequencing kit${GREEN}\tSQK-DCS108${NOCOLOUR}"
	echo -e "${YELLOW}10 - Direct RNA sequencing kit${GREEN}\tSQK-RNA001${NOCOLOUR}"
	echo -e "${YELLOW}11 - cDNA-PCR Sequencing Kit${GREEN}\tSQK-PCS108${NOCOLOUR}"
	while [ "$Switch" -eq "0" ]; do
		read -e Minion_Kit
		case $Minion_Kit in
			1)
				echo -e "${GREEN}Selected Ligation sequencing kit 1D${GREEN}\tSQK-LSK108${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK108"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			2)
				echo -e "${GREEN}Selected Ligation sequencing kit${GREEN}\tSQK-LSK109${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK109"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			3)
				echo -e "${GREEN}Selected Rapid sequencing kit${GREEN}\tSQK-RAD004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RAD004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			4)
				echo -e "${GREEN}Selected PCR sequencing kit${GREEN}\tSQK-PSK004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-PSK004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			5)
				echo -e "${GREEN}Selected 1D^2 sequencing kit${GREEN}\tSQK-LSK308${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK308"
				Direction="1D2"
				Product="DNA"
				Switch=1
			;;
			6)
				echo -e "${GREEN}Selected PCR Barcoding kit${GREEN}\tSQK-PBK004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-PBK004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			7)
				echo -e "${GREEN}Selected Rapid PCR Barcoding kit${GREEN}\tSQK-RPB004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RPB004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			8)
				echo -e "${GREEN}Selected 16S Barcoding kit${GREEN}\tSQK-RAB204${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RAB204"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			9)
				echo -e "${GREEN}Selected Direct cDNA sequencing kit${GREEN}\tSQK-DCS108${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-DCS108"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			10)
				echo -e "${GREEN}Selected Direct RNA sequencing kit${GREEN}\tSQK-RNA001${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RNA001"
				Direction="1D"
				Product="RNA"
				Switch=1
			;;
			11)
				echo -e "${GREEN}Selected cDNA-PCR Sequencing Kit${GREEN}\tSQK-PCS108${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-PCS108"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			*)
				echo -e "${RED}ERROR: Please only type a number from 1 to 11${NOCOLOUR}"
			;;
		esac
	done
fi
if ! grep -i -q "Kit" $ParFile; then echo -e "Kit	$Minion_Kit" >> $ParFile; fi
if ! grep -i -q "Direction" $ParFile; then echo -e "Direction	$Minion_Kit" >> $ParFile; fi
if ! grep -i -q "Product" $ParFile; then echo -e "Product	$Minion_Kit" >> $ParFile; fi

if [ $DIR_RawReads == "nil" ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter the file location of your fast5 reads:${NOCOLOUR}"
		read -e DIR_RawReads
		if [ -d $DIR_RawReads ]; then
			Switch=1
			echo -e "${BLUE}You entered: ${GREEN}$DIR_RawReads${NOCOLOUR}"
		else
			echo -e "${RED}Directory does not exist: ${GREEN}$DIR_RawReads${NOCOLOUR}"
		fi
	done
fi

if ! grep -i -q "Original reads" $ParFile; then echo -e "Original reads	$DIR_RawReads" >> $ParFile; fi


if [ ! -d $DIR_calledreads ]; then
	mkdir $DIR_calledreads
fi

if [ ! -d $DIR_calledreads/workspace ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running Albacore on fast5 files${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Minion flow cell:${GREEN} $Minion_Flow${NOCOLOUR}"
	echo -e "${BLUE}Minion kit:${GREEN} $Minion_Kit${NOCOLOUR}"
	echo -e "${BLUE}Worker threads:${GREEN} $WorkerThreads${NOCOLOUR}"
	read_fast5_basecaller.py -v | tee -a $Progress
	read_fast5_basecaller.py -i $DIR_RawReads -s $DIR_calledreads -r -t $WorkerThreads -f $Minion_Flow -k $Minion_Kit -o fastq --disable_filtering
fi

if [ ! -e "$DIR_ResultsSum/Stats-01-Called_reads.txt" ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running NanoStat on called reads${NOCOLOUR}" | tee -a $Progress
	NanoStat --version | tee -a $Progress
	NanoStat --fastq "$DIR_calledreads/workspace/"* --readtype $Direction -t $WorkerThreads -n "$DIR_ResultsSum/Stats-01-Called_reads.txt"
fi

if [ ! -d $DIR_TrimmedReads ]; then
	mkdir $DIR_TrimmedReads
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running PoreChop on called reads to remove adapters${NOCOLOUR}" | tee -a $Progress
	echo -e "Porechop version: $(porechop --version)" | tee -a $Progress

# NEED TO CHECK STDERR FOR OUTPUT ALSO OR PLACE IN TEMP AND REMOVE ADAPTERS FOR PROGRESS PAGE
	if [ $Barcoded == "TRUE" ]; then
		porechop -i "$DIR_calledreads" -b "$DIR_TrimmedReads" -v 1
	else
		porechop -i "$DIR_calledreads" -o "$DIR_TrimmedReads/$Project.porechop.fastq"  -v 1
	fi
fi

	basename -a $DIR_TrimmedReads/* > "$Meta/ReadFileNames.txt"

if [ ! -e "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt" ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running NanoStat on trimmed reads${NOCOLOUR}" | tee -a $Progress
	NanoStat --version | tee -a $Progress
	while read i; do
		echo -e $i >> "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt"
		NanoStat --fastq "$DIR_TrimmedReads/$i" --readtype $Direction -t $WorkerThreads >> "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt"
	done < "$Meta/ReadFileNames.txt"
fi

if [ ! -d $DIR_FilteredReads ]; then
	mkdir $DIR_FilteredReads
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running FiltLong on called reads to remove low quality reads${NOCOLOUR}" | tee -a $Progress
	filtlong --version | tee -a $Progress
	while read i; do
		echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
		echo -e "${GREEN}$i${NOCOLOUR}" | tee -a $Progress
		filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 "$DIR_TrimmedReads/$i" > "$DIR_FilteredReads/$i"
	done < "$Meta/ReadFileNames.txt"
fi

if [ ! -e "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt" ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running NanoStat on trimmed reads${NOCOLOUR}" | tee -a $Progress
	NanoStat --version | tee -a $Progress
	while read i; do
		echo -e $i >> "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt"
		NanoStat --fastq "$DIR_FilteredReads/$i" --readtype $Direction -t $WorkerThreads >> "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt"
	done < "$Meta/ReadFileNames.txt"
fi

if [ ! -d $DIR_canu ]; then
	mkdir $DIR_canu
fi
while read i; do
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running canu on filtered $i${NOCOLOUR}" | tee -a $Progress
	canu --version | tee -a $Progress
	canu -d $DIR_canu.$Project -p $i genomeSize=500k -nanopore-raw "$DIR_FilteredReads/$i"
done < "$Meta/ReadFileNames.txt"