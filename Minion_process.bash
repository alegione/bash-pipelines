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
  DelayNeeded="FALSE"

  if grep -i -q "User" $ParFile; then HomeDir=$(grep -i "User" $ParFile | cut -f2); echo -e "${GREEN}User:\t$HomeDir${NOCOLOUR}";else HomeDir="nil"; fi
  if grep -i -q "Project" $ParFile; then Project=$(grep -i "Project" $ParFile | cut -f2); echo -e "${GREEN}Project directory:\t$Project${NOCOLOUR}";else Project="nil"; fi
  if grep -i -q "Original reads" $ParFile; then DIR_RawReads=$(grep -i "Original reads" $ParFile | cut -f2); echo -e "${GREEN}Raw reads directory:\t$DIR_RawReads${NOCOLOUR}";else DIR_RawReads="nil";fi
  if grep -i -q "Called reads" $ParFile; then DIR_calledreads=$(grep -i "Called reads" $ParFile | cut -f2); echo -e "${GREEN}Called reads directory:\t$DIR_calledreads${NOCOLOUR}";else DIR_calledreads="nil";fi

  if grep -i -q "Flow cell" $ParFile; then Minion_Flow=$(grep -i "Flow cell" $ParFile | cut -f2); echo -e "${GREEN}Flow cell:\t$Minion_Flow${NOCOLOUR}";else Minion_Flow="nil";fi
  if grep -i -q "Kit" $ParFile; then Minion_Kit=$(grep -i "Kit" $ParFile | cut -f2); echo -e "${GREEN}Kit:\t$Minion_Kit${NOCOLOUR}";else Minion_Kit="nil";fi
  if grep -i -q "Barcoded" $ParFile; then Barcoded=$(grep -i "Barcoded" $ParFile | cut -f2); echo -e "${GREEN}Barcoded:\t$Barcoded${NOCOLOUR}";else Barcoded="nil";fi
  if grep -i -q "Direction" $ParFile; then Direction=$(grep -i "Direction" $ParFile | cut -f2); echo -e "${GREEN}Direction:\t$Direction${NOCOLOUR}";else Direction="nil";fi
  if grep -i -q "Product" $ParFile; then Product=$(grep -i "Product" $ParFile | cut -f2); echo -e "${GREEN}Product:\t$Product${NOCOLOUR}";else Product="nil";fi

  if grep -i -q "Reference genome" $ParFile; then RefGenome=$(grep -i "Reference genome" $ParFile | cut -f2); echo -e "${GREEN}Reference genome:\t$RefGenome${NOCOLOUR}";else RefGenome="nil";fi
  if grep -i -q "Reference length" $ParFile; then RefLength=$(grep -i "Reference length" $ParFile | cut -f2); echo -e "${GREEN}Reference length:\t$RefLength${NOCOLOUR}";else RefLength="nil";fi
  sleep 1

else
  DelayNeeded="TRUE"
  HomeDir="nil"
  Project="nil"
  DIR_RawReads="nil"
  DIR_calledreads="nil"
  Minion_Flow="nil"
  Minion_Kit="nil"
  Barcoded="nil"
  Direction="nil"
  Product="nil"
  RefGenome="nil"
  RefLength="nil"
fi

if [ "$HomeDir" == "nil" ]; then

  Switch=0
  while [ "$Switch" -eq "0" ]; do
    echo -e "${BLUE}Please enter your home directory ${YELLOW}(eg: Student/YOURNAME):${NOCOLOUR}"
    read -e HomeDir
    HomeDir=$(printf $HomeDir | sed 's/\/$//')
    if [ ! -d $HomeDir ]; then
      echo -e "${RED}No directory of that name exists${NOCOLOUR}"
      echo -e "Would you like to create it? (Y/N)${NOCOLOUR}"
      read -N 1 yesno
      echo -e ""
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
DIR_calledreads="$Dir/Called_reads"
DIR_TrimmedReads="$Dir/PoreChop"
DIR_FilteredReads="$Dir/FiltLong"
DIR_ResultsSum="$Dir/Results_Summary"
DIR_Alignment="$Dir/Alignments"
DIR_canu="$Dir/canu"
WorkerThreads=$(($CORES * 2))

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
				echo -e "${GREEN}Selected Ligation sequencing kit 1D\tSQK-LSK108${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK108"
				Direction="1D"
				Product="DNA"
        Switch=1
			;;
			2)
				echo -e "${GREEN}Selected Ligation sequencing kit\tSQK-LSK109${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK109"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			3)
				echo -e "${GREEN}Selected Rapid sequencing kit\tSQK-RAD004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RAD004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			4)
				echo -e "${GREEN}Selected PCR sequencing kit\tSQK-PSK004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-PSK004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			5)
				echo -e "${GREEN}Selected 1D^2 sequencing kit\tSQK-LSK308${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-LSK308"
				Direction="1D2"
				Product="DNA"
				Switch=1
			;;
			6)
				echo -e "${GREEN}Selected PCR Barcoding kit\tSQK-PBK004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-PBK004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			7)
				echo -e "${GREEN}Selected Rapid PCR Barcoding kit\tSQK-RPB004${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RPB004"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			8)
				echo -e "${GREEN}Selected 16S Barcoding kit\tSQK-RAB204${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RAB204"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			9)
				echo -e "${GREEN}Selected Direct cDNA sequencing kit\tSQK-DCS108${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-DCS108"
				Direction="1D"
				Product="DNA"
				Switch=1
			;;
			10)
				echo -e "${GREEN}Selected Direct RNA sequencing kit\tSQK-RNA001${NOCOLOUR}" | tee -a $Progress
				Minion_Kit="SQK-RNA001"
				Direction="1D"
				Product="RNA"
				Switch=1
			;;
			11)
				echo -e "${GREEN}Selected cDNA-PCR Sequencing Kit\tSQK-PCS108${NOCOLOUR}" | tee -a $Progress
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

if [ $Barcoded == "nil" ]; then
  echo -e "${GREEN}Did you include barcoding (Yes/No)${NOCOLOUR}"
  read -N 1 yesno
  yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
  if [ $yesno = "y" ]; then
    Barcoded="TRUE"
    echo -e "${BLUE}\nReads will be binned by barcode in PoreChop${NOCOLOUR}" | tee -a $Progress
  else
    Barcoded="FALSE"
    echo -e "${BLUE}\nReads will NOT be binned by barcode in PoreChop${NOCOLOUR}" | tee -a $Progress
  fi
fi

if ! grep -i -q "Kit" $ParFile; then echo -e "Kit\t$Minion_Kit" >> $ParFile; fi
if ! grep -i -q "Direction" $ParFile; then echo -e "Direction\t$Direction" >> $ParFile; fi
if ! grep -i -q "Product" $ParFile; then echo -e "Product\t$Product" >> $ParFile; fi
if ! grep -i -q "Barcoded" $ParFile; then echo -e "Barcoded\t$Barcoded" >> $ParFile; fi


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

if ! grep -i -q "Original reads" $ParFile; then echo -e "Original reads\t$DIR_RawReads" >> $ParFile; fi

if [ "$RefLength" == "nil" ]; then
  #Asks for path a reference genome
  Switch=0
  while [ "$Switch" -eq "0" ]; do
    echo -e "${BLUE}Would you like to use a reference genome in downstream analysis (Y/N)?${NOCOLOUR}"
    read -N 1 yesno
    echo "" #adds a new line
    yesno=$(echo -e "$yesno" | tr '[:upper:]' '[:lower:]')
    if [ "$yesno" == "y" ]; then
      echo -e "${BLUE}Please enter the file location of your reference genome (in fasta format):${NOCOLOUR}"
      read -e RefGenome
      if [ -e $RefGenome ]; then
        Switch=1
        echo -e "${BLUE}You entered: ${GREEN}$RefGenome${NOCOLOUR}"
        echo -e "Reference genome\t$RefGenome" >> $ParFile
        echo -e "${BLUE}Building reference genome index with minimap2${NOCOLOUR}"
#       cp $RefGenome $Meta/ref.fa
        minimap2 -d $Meta/ref.mmi $RefGenome > /dev/null 2>&1
#       if [ ! -e "$RefGenome.bwt" ]; then
#         echo -e "${BLUE}Building bwa index from reference genome:${GREEN} $RefGenome ${NOCOLOUR}"
#         bwa index $RefGenome
#       fi
        RefLength=$(cat $RefGenome | awk '$0 !~ ">" {c+=length($0);} END { print c; }')
	echo -e "Reference length\t$RefLength" >> $ParFile
      else
        echo -e "${RED}File does not exist: ${GREEN}$RefGenome${NOCOLOUR}"
      fi
    else
      Switch=1
      RefGenome="None"
      echo -e "${BLUE}Please enter an approximate length of your target (in bp):${NOCOLOUR}"
      read -e RefLength
      echo -e "Reference length\t$RefLength" >> $ParFile
    fi
  done
fi

#### SLEEP FOR HOW LONG BEFORE STARTING??
#if [ $DelayNeeded == "TRUE" ]; then
  echo -e "${BLUE}How long to sleep for before starting (eg. for 2 hours type: 2h): ${NOCOLOUR}"
  read sleeptime
  sleep "$sleeptime"
#fi

# START WORKFLOW HERE

if [ ! -d $DIR_calledreads ]; then
	mkdir $DIR_calledreads
fi

if [ ! -d $DIR_calledreads/workspace ]; then
	#DETERMINE SCRIPT FOR 1D or 1D^2
	if [ "$Direction" == "1D2" ]; then
		Albacore_Script="full_1dsq_basecaller.py"
	else
		Albacore_Script="read_fast5_basecaller.py"
	fi
	#PRINT INFO TO TERMINAL AND PROGRESS FILE
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running Albacore on fast5 files${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Minion flow cell:${GREEN} $Minion_Flow${NOCOLOUR}"
	echo -e "${BLUE}Minion kit:${GREEN} $Minion_Kit${NOCOLOUR}"
	echo -e "${BLUE}Worker threads:${GREEN} $WorkerThreads${NOCOLOUR}"
	$Albacore_Script -v | tee -a $Progress
 	#RUN ALBACORE
	$Albacore_Script -i $DIR_RawReads -s $DIR_calledreads -r -t $WorkerThreads -f $Minion_Flow -k $Minion_Kit -o fastq --disable_filtering
fi

if [ ! -e "$DIR_ResultsSum/Stats-01-Called_reads.txt" ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running NanoStat on called reads${NOCOLOUR}" | tee -a $Progress
	NanoStat --version | tee -a $Progress
  if [ $Direction != "1D2" ]; then
	   NanoStat --fastq "$DIR_calledreads/workspace/"*.fastq --readtype $Direction -t $WorkerThreads -n "$DIR_ResultsSum/Stats-01-Called_reads.txt"
   else
     NanoStat --fastq "$DIR_calledreads/workspace/"*.fastq --readtype $Direction -t $WorkerThreads -n "$DIR_ResultsSum/Stats-01-Called_reads.txt"
     NanoStat --fastq "$DIR_calledreads/1dsq_analysis/workspace/"*.fastq --readtype $Direction -t $WorkerThreads -n "$DIR_ResultsSum/Stats-01-Called_reads_1D2.txt"
   fi
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
    if [ $Direction != "1D2" ]; then
      porechop -i "$DIR_calledreads/workspace/" -o "$DIR_TrimmedReads/$Project.porechop.fastq"  -v 1
    else
      porechop -i "$DIR_calledreads/workspace/" -o "$DIR_TrimmedReads/$Project.1D.porechop.fastq"  -v 1
      porechop -i "$DIR_calledreads/1dsq_analysis/workspace/" -o "$DIR_TrimmedReads/$Project.1D2.porechop.fastq"  -v 1
    fi
  fi
fi

# COLLECT NAMES OF FILES FROM PORECHOP OUTPUT
  basename -a -s ".porechop.fastq" $DIR_TrimmedReads/* > "$Meta/ReadFileNames.txt"

# PRODUCE STATS FOR EACH PORECHOP RESULT
if [ ! -e "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt" ]; then
  echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
  echo -e "${BLUE}Running NanoStat on trimmed reads${NOCOLOUR}" | tee -a $Progress
  NanoStat --version | tee -a $Progress
  while read i; do
#      if [ $Direction != "1D2" ]; then
    echo -e "$i.porechop.fastq" >> "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt"
    NanoStat --fastq "$DIR_TrimmedReads/$i.porechop.fastq" --readtype $Direction -t $WorkerThreads >> "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt"
#      else
#          echo -e "$i.fastq" >> "$DIR_ResultsSum/Stats-02-PoreChop_reads.txt"
#      fi

  done < "$Meta/ReadFileNames.txt"
fi

# FILTER READS WITH FILTLONG
if [ ! -d $DIR_FilteredReads ]; then
	mkdir $DIR_FilteredReads
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running FiltLong on called reads to remove low quality reads${NOCOLOUR}" | tee -a $Progress
	filtlong --version | tee -a $Progress
	while read i; do
		echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
		echo -e "${GREEN}Filtering $i.porechop.fastq and saving as $i.filtered.fastq${NOCOLOUR}" | tee -a $Progress
		filtlong --min_length 500 --keep_percent 90 "$DIR_TrimmedReads/$i.porechop.fastq" > "$DIR_FilteredReads/$i.filtered.fastq"
	done < "$Meta/ReadFileNames.txt"
fi

# PRODUCE STATS FOR EACH FILTLONG RESULT
if [ ! -e "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt" ]; then
	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
	echo -e "${BLUE}Running NanoStat on filtered reads${NOCOLOUR}" | tee -a $Progress
	NanoStat --version | tee -a $Progress
	while read i; do
		echo -e "$i.filtered.fastq" >> "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt"
		NanoStat --fastq "$DIR_FilteredReads/$i.filtered.fastq" --readtype $Direction -t $WorkerThreads >> "$DIR_ResultsSum/Stats-03-FiltLong_reads.txt"
	done < "$Meta/ReadFileNames.txt"
fi

# PRODUCE ALIGNMENT FOR EACH FILTLONG RESULT TO THE REFERENCE GENOME
if [ "$RefGenome" != "nil" ] || [ "$RefGenome" != "None" ]; then
  if [ ! -d $DIR_Alignment ]; then
    mkdir $DIR_Alignment
  fi
  echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
  echo -e "${BLUE}Running Minimap2 on filtered reads against the provided reference: $RefGenome${NOCOLOUR}" | tee -a $Progress
  while read i; do
    if [ ! -e "$DIR_Alignment/$i.sam" ]; then
      echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
      echo -e "${BLUE}Running Minimap2 on $i${NOCOLOUR}" | tee -a $Progress
      minimap2 -L -ax map-ont "$Meta/ref.mmi" "$DIR_FilteredReads/$i.filtered.fastq" > "$DIR_Alignment/$i.sam"
            samtools view -h -S -F 4 "$DIR_Alignment/$i.sam" > "$DIR_Alignment/$i.aligned.sam"
            samtools sort --output-fmt SAM --threads $CORES "$DIR_Alignment/$i.aligned.sam" > "$DIR_Alignment/$i.sort.sam"
            cut -f1 "$DIR_Alignment/$i.sort.sam" | sort | uniq -c | grep "2 " | cut -f8 -d " " > "$DIR_Alignment/$i.dups.txt"
            grep -v --file="$DIR_Alignment/$i.dups.txt" "$DIR_Alignment/$i.sort.sam" > "$DIR_Alignment/$i.nodups.sam"
    fi
  done < "$Meta/ReadFileNames.txt"
fi

# MAKE THE DIRECTORY FOR EACH ASSEMBLY
if [ ! -d $DIR_canu ]; then
  mkdir $DIR_canu
fi

if [ "$Product" == "RNA" ]; then
  if [ ! -d $DIR_canu/RNAtoDNA ]; then
    mkdir $DIR_canu/RNAtoDNA
  fi
  while read i; do
    if [ ! -e "$DIR_canu/RNAtoDNA/$i.DNA.fastq" ]; then
      echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
      echo -e "${BLUE}Converting U to T in $i filtered fastq file${NOCOLOUR}" | tee -a $Progress
      awk 'NR % 4 == 2 {$1;for (i=1;i<=NF;i++) {gsub(/U/,"T",$i); printf "%s\n",$i}} NR % 4 != 2 { print }' "$DIR_FilteredReads/$i.filtered.fastq" > "$DIR_canu/RNAtoDNA/$i.DNA.fastq"
    fi
  done < "$Meta/ReadFileNames.txt"
fi

# ASSEMBLE EACH FILE IN CANU, NEED TO UPDATE genomeSize BASED ON REFERENCE, OR ASK AT THE SAME INITIAL PROMPT, THIS WON'T WORK FOR MULTIPLE GENOMES OR REFERENCES THOUGH
while read i; do
  if [ ! -e "$DIR_canu/$i/$i.contigs.fasta" ]; then
  	echo -e "${PURPLE}$(date)${NOCOLOUR}" | tee -a $Progress
  	echo -e "${BLUE}Running canu on filtered $i.fastq${NOCOLOUR}" | tee -a $Progress
  	canu --version | tee -a $Progress
	if [ "$Product" == "DNA" ]; then
#          canu -d $DIR_canu/$i -p $i genomeSize=$RefLength minMemory-nanopore-raw "$DIR_FilteredReads/$i.filtered.fastq"
          canu -d $DIR_canu/$i -p $i batMemory=32 genomeSize=$RefLength minMemory-nanopore-raw "$DIR_FilteredReads/$i.filtered.fastq"
	else
	  canu -d $DIR_canu/$i -p $i genomeSize=$RefLength -nanopore-raw "$DIR_canu/RNAtoDNA/$i.DNA.fastq"
	fi
  fi
done < "$Meta/ReadFileNames.txt"

# CLEAN UP PROGRESS TEXT FILE TO REMOVE SPECIAL CHARACTERS THAT DENOTE COLOUR
sed -i 's/\x1b\[[0-9;]*m//g' $Progress
