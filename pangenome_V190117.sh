#!/bin/sh
########################################
#
# run pipeline 
#
########################################

##
# BELOW ARE RANDOM COMMANDS TO BE ADDED TO EITHER THIS FILE OR TO Assemble.sh
##
# while read i; do echo $i; genomeCoverageBed -ibam Shared_folder/Alistair_Legione/Cpec_Assemble/BWAbamsorted/$i.sort.bam -d > Shared_folder/Alistair_Legione/Cpec_Assemble/$i.depth.txt; done < Shared_folder/Alistair_Legione/Cpec_Assemble/BWAbamsorted/contents.txt 
#echo -e "Genome\tDepth\tCoverage\tCoverage>10" > "Shared_folder/Alistair_Legione/Cpec_Assemble/coverage.txt"; while read i; do tot=$(awk '{s+=$3}END{print s}' Shared_folder/Alistair_Legione/Cpec_Assemble/$i.depth.txt); covered=$(awk 'BEGIN{count=0}{if ($3!=0) count++;}END{print count}' Shared_folder/Alistair_Legione/Cpec_Assemble/$i.depth.txt); coveredwell=$(awk 'BEGIN{count=0}{if ($3>9) count++;}END{print count}' Shared_folder/Alistair_Legione/Cpec_Assemble/$i.depth.txt); len=$(wc -l < Shared_folder/Alistair_Legione/Cpec_Assemble/$i.depth.txt); depth=$(bc <<< "scale=2; $tot/$len"); coverage=$(bc <<< "scale=4; $covered/$len*100");  goodcoverage=$(bc <<< "scale=4; $coveredwell/$len*100"); echo -e "$i\t$depth\t$coverage\t$goodcoverage" >> Shared_folder/Alistair_Legione/Cpec_Assemble/coverage.txt; done < Shared_folder/Alistair_Legione/Cpec_Assemble/BWAbamsorted/contents.txt 



#Set colour variables
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
PURPLE='\033[1;35m'
NOCOLOUR='\033[0m'

#Edit for path to prokka (if prokka is not in your path)
Prokka="prokka"
Baseml="/home/qiime/Documents/paml/src/baseml"
Codeml="/home/qiime/Documents/paml/src/codeml"
PamlConverter="/home/qiime/Documents/phylip2paml.pl"
PhylipConverter="/home/qiime/Documents/fasta2phylip.pl"


Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please enter your project directory ${YELLOW}(ie: Documents/YOURNAME/Project):${NOCOLOUR}"
	read -e Dir
	Dir=$(printf $Dir | sed 's/\/$//')
	echo $Dir
	if [ ! -d $Dir ]; then
		echo -e "${RED}No project of that name exists${NOCOLOUR}"
	else	
		echo -e "${BLUE}You entered: ${GREEN}$Dir${NOCOLOUR}"		
		Switch=1
	fi
done


Dir="$Dir/pangenome"
Meta="$Dir/Metadata" #variable for metadata folder, used for shorter downstream coding
FileLocation="$Dir/genomes_fasta"

if [ ! -d "$Dir" ]; then
	mkdir "$Dir"
	mkdir "$Dir/Metadata"
fi



#Asks for path to the target directory where the read files are (eg zipped read files)
if [ ! -d $FileLocation ]; then
	Switch=0
	while [ "$Switch" -eq "0" ]; do
		echo -e "${BLUE}Please enter the file location of your genomes ${RED}(your files MUST be a #.fasta format):${NOCOLOUR}"
		read -e GenomeDir
		GenomeDir=$(printf $GenomeDir | sed 's/\/$//')
		if [ -d $GenomeDir ]; then
			Switch=1
			echo -e "${BLUE}You entered: ${GREEN}$GenomeDir${NOCOLOUR}"
		else
			echo -e "${RED}Directory does not exist: ${GREEN}$GenomeDir${NOCOLOUR}"
		fi
	echo -e "${BLUE}Preparing genome metadata file${NOCOLOUR}"
	ls $GenomeDir/*.fasta | tr '\n' '\0' | xargs -0 -n 1 basename | sed 's/\.fasta\//' > "$Meta/PanGenomeList.txt" #lists all files ending in 1.fastq.gz (* is a wildcard) in GenomeDir, removes the extension from the list and saves the list to a text file in the metadata folder
	mkdir $FileLocation
	cp $GenomeDir/*.fasta $FileLocation
	done
fi

GenomeCount=$(wc -l < "$Meta/PanGenomeList.txt")

Switch=0
while [ "$Switch" -eq "0" ]; do
	echo -e "${BLUE}Please indicate the location of your Metadata/Mapping file${NOCOLOUR}"
	read -e Mapping
	if [ -e $Mapping ]; then
		Switch=1
		echo -e "${BLUE}You entered ${GREEN}$Mapping${NOCOLOUR}"
		echo -e "${BLUE}Moving mapping file to ${GREEN}$Dir/Metadata/metadata.txt${NOCOLOUR}"
		cat $Mapping > "$Dir/Metadata/$Project.txt"
		MapHead="$Dir/Metadata/$Project.head.txt"				
		head -1 $Mapping > $MapHead
	else
		echo -e "${RED}File does not exist: ${GREEN}$Mapping${NOCOLOUR}"
	fi
done

MappingCount=$(wc -l < "$Meta/PanGenomeList.txt")
MappingCount=$((MappingCount-1))

if [ $GenomeCount -ne $MappingCount ]; then
	echo -e "${RED}The number of genomes in the mapping file (${GREEN}$MappingCount${RED}) does not match the number detected in the genome directory (${GREEN}$GenomeCount${RED}), please rectify to continue${NOCOLOUR}"
	echo -e "${RED}Exiting program${NOCOLOUR}"
	exit 1
else
	echo -e "${BLUE}The number of genomes in the mapping file (${GREEN}$MappingCount${BLUE}) matches the number detected in the genome directory (${GREEN}$GenomeCount${BLUE}), good for you!${NOCOLOUR}"
fi


#Mapping="$Meta/Genelist"
#AllGenes="$Meta/AllGenes"
#echo "" > "$AllGenes.txt"
#echo "I'm here"

#read -p "Is your target $FileLocation? (Y/N)" Switch

#if [[ $Switch != "y"* ]]; then
#	echo "target file location:"
#	#/home/qiime/Desktop/Shared_Folder/Cpec_whole_genome/Koalagenes.txt
#	read FileLocation
#fi

#echo "$FileLocation"

#if [ ! -a "$AllGenes.txt" ]; then
#	sed 's/[] [\/,:"@-]/_/g' $FileLocation | sed '1d' | sort -k2 > "$AllGenes.txt" #remove special characters, important to have the hyphen at the start or end, otherwise it will treat it as signifying a range
#	cut -f 1 "$AllGenes.txt" | sort -u > "$Meta/GenomeList.txt"
#fi

#if all ORFS have novel names, can sort as per genome list above










##### add protein database of trusted pecorum proteins?!?!?! E58! where to get? what format?

#Making a Core Databases

#If you want to modify these core databases, the included script prokka-uniprot_to_fasta_db, along with the official uniprot_sprot.dat, can be used to generate a new database to put in /path/to/prokka/db/kingdom/. If you add new ones, the command prokka --listdb will show you whether it has been detected properly.

#The Genus Databases

#If you enable --usegenus and also provide a Genus via --genus then it will first use a BLAST database which is Genus specific. Prokka comes with a set of databases for the most common Bacterial genera; type prokka --listdb to see what they are.


#Adding a Genus Databases

#If you have a set of Genbank files and want to create a new Genus database, Prokka comes with a tool called prokka-genbank_to_fasta_db to help. For example, if you had four annotated "Coccus" genomes, you could do the following:

#% quokka-genbank_to_fasta_db Coccus1.gbk Coccus2.gbk Coccus3.gbk Coccus4.gbk > Coccus.faa
#% cd-hit -i Coccus.faa -o Coccus -T 0 -M 0 -g 1 -s 0.8 -c 0.9
#% rm -fv Coccus.faa Coccus.bak.clstr Coccus.clstr
#% makeblastdb -dbtype prot -in Coccus
#% mv Coccus.p* /path/to/prokka/db/genus/










#if the annotated directory doesn't exist, create all appropriate folders and sub-folders
if [ ! -d "$Dir/Annotated_Genomes" ]; then 
	mkdir "$Dir/Annotated_Genomes"
	mkdir "$Dir/gff_genomes"
fi

#run through the list of genomes run prokka if it hasn't been done before
Count=1
while read FileLine; do
	
	if [ ! -d "$Dir/Annotated_Genomes/$FileLine" ]; then
		echo -e "${BLUE}Annotating ${YELLOW}$FileLine ${BLUE}(${GREEN}$Count ${BLUE}of ${GREEN}$GenomeCount ${BLUE}genomes) with prokka${NOCOLOUR}"
		perl $Prokka --outdir "$Dir/Annotated_Genomes/$FileLine" --addgenes --force --kingdom Bacteria --genus Chlamydia --species pecorum --usegenus --proteins [X] --locustag "$FileLine" "$FileLocation/$FileLine"
		DateExtension=$(date +%d%m%Y)
		echo -e "${BLUE}Moving .gff file to group directory${NOCOLOUR}"
		
		cp "$Dir/Annotated_Genomes/$FileLine/"$FileLine"_${DateExtension}.gff" "$Dir/gff_genomes"
		
	else
		echo -e "${YELLOW}$FileLine ${BLUE}(${GREEN}$Count ${BLUE}of ${GREEN}$GenomeCount ${BLUE}genomes) has already been annotated and will be skipped${NOCOLOUR}"
	fi
	Count=$((Count+1))
done < "$Meta/PanGenomeList.txt"

#produce a pan genome with roary
#Usage:   roary [options] *.gff
#Options: -p INT    number of threads [1]
#         -o STR    clusters output filename [clustered_proteins]
#         -f STR    output directory [.]
#         -e        create a multiFASTA alignment of core genes using PRANK
#         -n        fast core gene alignment with MAFFT, use with -e
#         -i        minimum percentage identity for blastp [95]
#         -cd FLOAT percentage of isolates a gene must be in to be core [99]
#         -qc       generate QC report with Kraken
#         -k STR    path to Kraken database for QC, use with -qc
#         -a        check dependancies and exit
#         -b STR    blastp executable [blastp]
#         -c STR    mcl executable [mcl]
#         -d STR    mcxdeblast executable [mcxdeblast]
#         -g INT    maximum number of clusters [50000]
#         -m STR    makeblastdb executable [makeblastdb]
#         -r        create R plots, requires R and ggplot2
#         -s        dont split paralogs
#         -t INT    translation table [11]
#         -z        dont delete intermediate files
#         -v        verbose output to STDOUT
#         -w        print version and exit
#         -y        add gene inference information to spreadsheet, doesnt work with -e
#         -h        this help message

if [ ! -d "$Dir/roary_output" ]; then 
	roary -f "$Dir/roary_output" -e -n -i 90 -cd 90 -z -p 4 -s -r -v $Dir/gff_genomes/*.gff
	
fi

if [ ! -d "$Dir/scoary_output" ]; then 
	mkdir "$Dir/scoary_output"
	scoary -t "$Dir/roary_output/gene_presence_absence.csv" -i "$Dir/Metadata/$Project.txt" -o "$Dir/scoary_output" -c BH -p 0.05 -u --threads 4 --no-time -e 10 --collapse
fi

#scoary -t Shared_folder/Alistair_Legione/traitdata3.csv -g DataDrive/Alistair_Legione/roary_koalas/gene_presence_absence.csv --no-time -e 10 -o DataDrive/Alistair_Legione/Scoary

	
if [ ! -d "$Dir/analysis" ]; then 
	mkdir "$Dir/analysis"
	FastTree -nt -gtr "$Dir/roary_output/core_gene_alignment.aln" > "$Dir/analysis/core_gene_alignment.newick"


	python "/home/qiime/Roary/contrib/roary_plots/roary_plots.py" "$Dir/analysis/core_gene_alignment.newick" "$Dir/roary/gene_presence_absence.csv"



fi

#srst2 --output testchlam --input_pe Shared_folder/Alistair_Legione/Cpec_reads/19238_1_16_1.fastq.gz Shared_folder/Alistair_Legione/Cpec_reads/19238_1_16_2.fastq.gz --mlst_db Chlamydiales_spp..fasta --mlst_definitions chlamydiales.txt --mlst_delimiter '_'


#while read i; do mkdir "Shared_folder/Alistair_Legione/shortreadMLST/$i"; srst2 --output "Shared_folder/Alistair_Legione/shortreadMLST/$i/$i" --input_pe "Shared_folder/Alistair_Legione/Cpec_reads/${i}_1.fastq.gz" "Shared_folder/Alistair_Legione/Cpec_reads/${i}_2.fastq.gz" --mlst_db Chlamydiales_spp..fasta --mlst_definitions chlamydiales.txt --mlst_delimiter '_' --log --gene_db ARGannot.r1.fasta --report_new_consensus --threads 4; done < Shared_folder/Alistair_Legione/cpecreadsforsrst2.txt 



#roary -p 4 -f DataDrive/Alistair_Legione/roary_koalas -e -n -r -z -v -i 90 -cd 90 Shared_folder/Alistair_Legione/genomes_gff/*.gff
#ls DataDrive/Alistair_Legione/roary_koalas/pan_genome_sequences/ > DataDrive/Alistair_Legione/roary_koalas/GeneList.txt
#sed -i 's/\.fa\.aln//g' DataDrive/Alistair_Legione/roary_koalas/GeneList.txt 
#while read i; do FastTree -nt -gtr -gamma -log "DataDrive/Alistair_Legione/fasttree_logs/$i.log" "DataDrive/Alistair_Legione/roary_koalas/pan_genome_sequences/$i.fa.aln" > "DataDrive/Alistair_Legione/fasttrees/$i.tre";done < DataDrive/Alistair_Legione/roary_koalas/GeneList.txt 
#scoary -t Shared_folder/Alistair_Legione/traitdata3.csv -g DataDrive/Alistair_Legione/roary_koalas/gene_presence_absence.csv --no-time -e 10 -o DataDrive/Alistair_Legione/Scoary

#scoary -t Shared_folder/Alistair_Legione/traitdata3.csv -g DataDrive/Alistair_Legione/roary_koalas/gene_presence_absence.csv --no-time -e 10 --collapse -o DataDrive/Alistair_Legione/Scoary


#Rscript Documents/Alistair_Legione/plottreeR_args.r Shared_folder/Alistair_Legione/traitdata3.csv Shared_folder/Alistair_Legione/fasttrees/tyrS.tre UGT0ocular1 Desktop/plot_test.pdf

#while read i; do count=$(grep -c ">" "Shared_folder/Alistair_Legione/Pangenome/roary_koalas/pan_genome_sequences/$i.fa.aln"); Rscript Documents/Alistair_Legione/plottreeR_args.r "Shared_folder/Alistair_Legione/traitdata3.csv" "Shared_folder/Alistair_Legione/Pangenome/fasttrees/$i.tre" "UGT1ocular0" "Shared_folder/Alistair_Legione/Pangenome/plotTrees/$count/$i.pdf"; done < Shared_folder/Alistair_Legione/Pangenome/roary_koalas/GeneList.txt 



Count="1"
if [ ! -d "$Dir/paml" ]; then
	mkdir "$Dir/paml"
	ls "$Dir/roary/pan_genome_sequences" | sed 's/.fa.aln//' > "$Meta/GeneList.txt"
	GeneCount=$(wc -l < "$Meta/GeneList.txt")
	while read FileLine; do
		if [ ! -d "$Dir/paml/$FileLine" ]; then
			echo "I'm currently looking for $FileLine"
			read  blank
			mkdir "$Dir/paml/$FileLine"
			cp "$Dir/roary/pan_genome_sequences/$FileLine.fa.aln" "$Dir/paml/$FileLine"
			mkdir "$Dir/paml/$FileLine/baseml"		
			cp $Baseml "$Dir/paml/$FileLine/baseml"
			mkdir "$Dir/paml/$FileLine/codeml"		
			cp $Codeml "$Dir/paml/$FileLine/codeml"
			echo "Building tree $Count/$GeneCount, using $FileLine alignment"
			perl $PhylipConverter "$Dir/paml/$FileLine/$FileLine.fa.aln" > "$Dir/paml/$FileLine/$FileLine.ph"
			FastTree -quiet -nosupport -nt -gtr $Dir/paml/$FileLine/$FileLine.ph > $Dir/paml/$FileLine/$FileLine.tre
			#phyml --quiet -i $Dir/paml/$FileLine/$FileLine.ph > $Dir/paml/$FileLine/$FileLine.tre
		fi
		Count=$((Count+1))
	done < "$Meta/GeneList.txt"
fi
Count="1"
GeneCount=$(wc -l < "$Meta/GeneList.txt")
if [ ! -d "$Dir/NoReference" ]; then 
	mkdir "$Dir/NoReference"
	mkdir "$Dir/Reference"
	while read FileLine; do
		echo "Reading gene $Count/$GeneCount and removing E58 reference"
		cat "$Dir/roary/pan_genome_sequences/$FileLine.fa.aln" | "/home/qiime/bioawk/bioawk" -c fastx '$name !~ /E58/ {print ">" $name "\n" $seq}' >> "$Dir/NoReference/$FileLine.fa"
		cp "$Dir/NoReference/$FileLine.fa" "$Dir/Reference/$FileLine.fa"
		cat "$Dir/roary/pan_genome_sequences/$FileLine.fa.aln" | "/home/qiime/bioawk/bioawk" -c fastx '$name ~ /E58/ {print ">" $name "\n" $seq}' >> "$Dir/Reference/$FileLine.fa"
		echo "Aligning $FileLine with Muscle (File $Count/$GeneCount)"
		muscle -quiet -in "$Dir/NoReference/$FileLine.fa" -out "$Dir/NoReference/$FileLine.aln.fa"
		muscle -quiet -in "$Dir/Reference/$FileLine.fa" -out "$Dir/Reference/$FileLine.fa.tmp"
		python "/home/qiime/Documents/stable.py" "$Dir/Reference/$FileLine.fa" "$Dir/Reference/$FileLine.fa.tmp" > "$Dir/Reference/$FileLine.aln.fa"
		sed -i 's/fasta.*/fasta/g' "$Dir/Reference/$FileLine.aln.fa" "$Dir/NoReference/$FileLine.aln.fa"
		rm "$Dir/Reference/$FileLine.fa"
		rm "$Dir/Reference/$FileLine.fa.tmp"
		rm "$Dir/NoReference/$FileLine.fa"
		Count=$((Count+1))
	done < "$Meta/GeneList.txt"
fi

if [ ! -d "$Dir/Collapse" ]; then
	mkdir "$Dir/Collapse/"
	mkdir "$Dir/Collapse/Reference"
	mkdir "$Dir/Collapse/NoReference"
	while read FileLine; do
		fastx_collapser -v "$Dir/Reference/$FileLine.aln.fa" "$Dir/Collapse/Reference/$FileLine.aln.fa"
		fastx_collapser -v "$Dir/NoReference/$FileLine.aln.fa" "$Dir/Collapse/NoReference/$FileLine.aln.fa"
	done < "$Meta/GeneList.txt"
fi

ls "$Dir/roary/pan_genome_sequences" | sed 's/.fa.aln//' > "$Meta/GeneList.txt"
Count="1"
GeneCount=$(wc -l < "$Meta/GeneList.txt")
	while read FileLine; do
		if [ -e "$Dir/paml/$FileLine/$FileLine.tre" ]; then
			perl $PamlConverter "$Dir/paml/$FileLine/$FileLine.ph" > "$Dir/paml/$FileLine/$FileLine.nuc"
			echo "Building ctl files $Count/$GeneCount, for $FileLine"

			echo "      seqfile = $Dir/paml/$FileLine/$FileLine.nuc					" > "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "     treefile = $Dir/paml/$FileLine/$FileLine.tre					" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo " 											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "      outfile = $Dir/paml/$FileLine/$FileLine.mlb       * main result file	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        noisy = 2   * 0,1,2,3: how much rubbish on the screen			" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "      verbose = 0   * 1: detailed output, 0: concise output			" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "                    * 3: StepwiseAddition; (4,5):PerturbationNNI 			" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        model = 4   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85			" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "                    * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "*        ndata = 100								" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        kappa = 5  * initial or fixed kappa					" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "       Malpha = 0   * 1: different alpha's for genes, 0: one alpha		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo " RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "   Small_Diff = 7e-6								" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "*        icode = 0  * (with RateAncestor=1. try GC in data,model=4,Mgene=4)	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed		" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			echo "       method = 0  * Optimization method 0: simultaneous; 1: one branch a time	" >> "$Dir/paml/$FileLine/baseml/baseml.ctl"
			cat "$Dir/paml/$FileLine/baseml/baseml.ctl" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
#		fi
			cd $Dir/paml/$FileLine/baseml/
			./baseml
			echo "	    seqfile = $Dir/paml/$FileLine/$FileLine.nuc					" > "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "     treefile = $Dir/paml/$FileLine/$FileLine.tre					" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "      outfile = $Dir/paml/$FileLine/$FileLine.mlc   * main result file name	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "      verbose = 1  * 0: concise; 1: detailed, 2: too much				" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs				" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "*        ndata = 10								" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        model = 2									" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * models for codons:						" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * models for AAs or codon-translated AAs:			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "      NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * 13:3normal>0							" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        icode = 10 * 0:universal code; 1:mammalian mt; 2-10:see below		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        Mgene = 0									" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * codon: 0:rates,1:separate;2:diff pi,3:diff kapa,4:all diff	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "                   * AA: 0:rates, 1:separate					" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        kappa = 2  * initial or fixed kappa					" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        omega = .4 * initial or fixed omega, for codons or codon-based AAs	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "       Malpha = 0  * different alphas for genes					" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        ncatG = 8  * # of categories in dG of NSsites models			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo " RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "   Small_Diff = .5e-6								" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "       method = 0  * Optimization method 0: simultaneous; 1: one branch a time	" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "											" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,		" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "* 10: blepharisma nu.								" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			echo "* These codes correspond to transl_table 1 to 11 of GENEBANK.			" >> "$Dir/paml/$FileLine/codeml/codeml.ctl"
			cd $Dir/paml/$FileLine/codeml/
			./codeml
		fi # remove me later?
		Count=$((Count+1))
	done < "$Meta/GeneList.txt"
#fi # remove me later?

#	query_pan_genome -a union "$Dir/gff/*.gff"
#	query_pan_genome -a complement "$Dir/gff/*.gff"
