#!/bin/bash

#get CPU metrics
CORES=$(nproc)
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

#Get start time
SCRIPTSTART=`date +%s`

walltime () {
	RED='\033[1;31m'
	GREEN='\033[1;32m'
	YELLOW='\033[1;33m'
	BLUE='\033[1;34m'
	PURPLE='\033[1;35m'
	NOCOLOUR='\033[0m'	
	TIMENOW=`date +%s`
	STARTTIME=$1
	TOTAL=$((TIMENOW-STARTTIME))
	
	if [ "$TOTAL" -lt "120" ]; then
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}seconds${NOCOLOUR}"
	elif [ "$TOTAL" -lt "7200" ]; then
		TOTAL=$(echo "$TOTAL/60" | bc)
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}minutes${NOCOLOUR}"
	else
		TOTAL=$(echo "$TOTAL/3600" | bc)
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}hours${NOCOLOUR}"
	fi
}




#Determine location and set environmental variables
if [ -d /home/alistair ]; then
	LOCATION="laptop"
else
	LOCATION="communal"
fi
echo "$LOCATION"


if [ "$LOCATION" = "laptop" ]; then
	if [ ! -e "Desktop/Shared_folder/mount_test.txt" ]; then
		sudo mount -t vboxsf Sequences Desktop/Shared_folder/
	fi
	DirReads="Desktop/Shared_folder/CpecFinal/Original_Reads"
	DirOutput="haptest"
	DirBWAOutput="$DirOutput/Sorted_BAM"
	DirFreebayes="$DirOutput/VCF"
	DirFlow="$DirOutput/FlowFiles"
	Reference="$DirOutput/E58_complete.fa"
	BWA="bwa mem"
	FREEBAYES="freebayes -F 0.05 -p 10"
elif [ "$LOCATION" = "communal" ]; then
	Dir="DataDrive/Alistair_Legione/CpecFinal"
	DirReads="DataDrive/Alistair_Legione/CpecFinal/Original_Reads"
	DirOutput="HapFlow"
	DirBWAOutput="$DirOutput/Sorted_BAM"
	DirFreebayes="$DirOutput/VCF"
	DirFlow="$DirOutput/FlowFiles"
	Reference="$DirOutput/E58_complete.fa"
	BWA="bwa mem"
	FREEBAYES="freebayes -F 0.05 -p 10"
fi

if [ ! -d $DirOutput ]; then
	mkdir $DirOutput
fi
if [ ! -d $DirBWAOutput ]; then
	mkdir $DirBWAOutput
fi
if [ ! -d $DirFreebayes ]; then
	mkdir $DirFreebayes
fi
if [ ! -d $DirFlow ]; then
	mkdir $DirFlow
fi

basename -a -s "_1.fastq.gz" $DirReads/*_1.fastq.gz > $DirOutput/metadata.txt

GenomeCount=$(wc -l < $DirOutput/metadata.txt)

echo -e "Processing $GenomeCount reads"

count="1"
while read i; do
	if [ ! -e "$DirBWAOutput/$i.bam" ] && [ ! -e "$DirBWAOutput/$i.sorted.bam" ]; then
		if [ ! -e $DirBWAOutput/$i.sam ]; then
			STARTTIME=`date +%s`
			echo -e "${PURPLE}$(date)${NOCOLOUR}"
			echo -e "${BLUE}Undertaking bwa mem alignment to reference for genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"
			$BWA -v 0 -t $CORES "$Reference" "$DirReads/${i}_1.fastq.gz" "$DirReads/${i}_2.fastq.gz" > "$DirBWAOutput/${i}.sam"
			walltime $STARTTIME
		fi
		STARTTIME=`date +%s`
		echo -e "${PURPLE}$(date)${NOCOLOUR}"
		echo -e "${BLUE}Converting SAM file to BAM FILE for genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"
		samtools import "$Reference.fai" "$DirBWAOutput/${i}.sam" "$DirBWAOutput/${i}.bam"
		walltime $STARTTIME		
		if [ -e "$DirBWAOutput/$i.sam" ]; then
			echo -e "${RED}DELETING REDUNDANT SAM FILE TO SAVE SPACE"
			rm "$DirBWAOutput/$i.sam"
		fi
	fi
	if [ ! -e "$DirBWAOutput/${i}.sorted.bam" ]; then
		STARTTIME=`date +%s`
		echo -e "${PURPLE}$(date)${NOCOLOUR}"
		echo -e "${BLUE}Sorting bam file for genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"	
		samtools sort -o "$DirBWAOutput/${i}.sorted.bam" "$DirBWAOutput/${i}.bam"
		walltime $STARTTIME
		echo -e "${RED}DELETING REDUNDANT UNSORTED BAM FILE TO SAVE SPACE"
		rm "$DirBWAOutput/${i}.bam"
	fi
	if [ ! -e "$DirFreebayes/${i}.vcf" ] && [ ! -e "$DirFreebayes/${i}.filt.vcf" ]; then
			STARTTIME=`date +%s`
			echo -e "${PURPLE}$(date)${NOCOLOUR}"
			echo -e "${BLUE}Producing vcf file from genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"
			samtools index "$DirBWAOutput/${i}.sorted.bam"
			$FREEBAYES -f "$Reference" -b "$DirBWAOutput/${i}.sorted.bam" > "$DirFreebayes/${i}.vcf"
			walltime $STARTTIME
	fi
	if [ ! -e "$DirFreebayes/${i}.filt.vcf" ]; then
		STARTTIME=`date +%s`
		echo -e "${PURPLE}$(date)${NOCOLOUR}"
		echo -e "${BLUE}Filtering vcf file for genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"		
		vcffilter -f "QUAL > 20" "$DirFreebayes/${i}.vcf" > "$DirFreebayes/${i}.filt.vcf"
		walltime $STARTTIME
		echo -e "${RED}DELETING REDUNDANT UNFILTERED VCF FILE TO SAVE SPACE${NOCOLOUR}"
		rm "$DirFreebayes/${i}.vcf" 
	fi
	if [ ! -e "$DirFlow/${i}.0.flw" ]; then
		STARTTIME=`date +%s`
		echo -e "${PURPLE}$(date)${NOCOLOUR}"
		echo -e "${BLUE}Producing flow file for genome ${YELLOW}$count/$GenomeCount\t${GREEN}${i}${NOCOLOUR}"
		HapFlow -b "$DirBWAOutput/${i}.sorted.bam" -v "$DirFreebayes/${i}.filt.vcf" -o "$DirFlow/${i}"
		walltime $STARTTIME	
	fi
	count=$((count+1))
done < $DirOutput/metadata.txt


echo -e "${BLUE}Hapflow loop finished${NOCOLOUR}"
walltime $SCRIPTSTART

#Get end time


