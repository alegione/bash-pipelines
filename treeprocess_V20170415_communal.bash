#!/bin/bash
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



#need $DirRoary directory, modeltest.tsv and bayesblock.txt in home directory, SNPcore in Snippy directory, in turn in home directory

# Program, script and dependency location variables
PROKKA="prokka"
ROARY="roary"
MRBAYES="mb"
ETE="ete3"
#MODELTEST="java -jar AlPrograms/jmodeltest2/dist/jModelTest.jar"
MODELTEST="java -jar /home/genomics/jmodeltest2-master/dist/jModelTest.jar"
#NEWICKCONVERTER="Desktop/Shared_folder/Scripts/R/TreeFromNexus.R"
NEWICKCONVERTER="/home/genomics/TreeFromNexus.R"
ALIGNMENTMETRICS="/home/genomics/Shared_folder/Alistair_Legione/R/Alignment_metrics.R"
#REFERENCETREE="/home/alistair/Snippy-core_noref/Bayes/SNPcore.newick.tre"
REFERENCETREE="/home/genomics/Shared_folder/Alistair_Legione/CpecFinal/Genome_Pipeline/SNPcore/core.3.bayes.newick.tre"
GCODE="11"
GENUS="Chlamydia"
SPECIES="pecorum"




####################################################
#
#
#	Order genelist by haplotype diversity (or equivalent) and run ALL processes in that order, to obtain output
#
#
####################################################







# Program input and output directories
PROJECT=""
Dir="/home/genomics/Documents/Alistair_Legione/koala_only"
#Dir="/home/alistair/Genome_Pipeline/koala_only"
#DirRoary="/home/genomics/Documents/Alistair_Legione/roary_80p_newProkka"
#DirProkka="/home/genomics/Documents/Alistair_Legione/NewProkka_Annotated_noref"
DirProkka="$Dir/Prokka_annotated"
DirRoary="$Dir/roary_koalas"
DirFasta="Desktop/Shared_folder/CpecFinal/GenomesFasta"
DirResults="$Dir/Results_Summary"
DirMetadata="$Dir/Metadata"

#if [ -e "$1/Metadata/Parameters.txt" ]; then
#	echo -e "${BLUE}Parameter file detected...obtaining previously entered options${NOCOLOUR}"
#	ParFile="$1/Metadata/Parameters.txt"
#	Meta="$1/Metadata"
#	if grep -i -q "Project" $ParFile; then Dir=$(grep -i "Project" $ParFile | cut -f2); echo -e "${GREEN}Project directory: $Dir${NOCOLOUR}";else Dir="nil"; fi
#	if grep -i -q "Original reads" $ParFile; then ReadDir=$(grep -i "Original reads" $ParFile | cut -f2); echo -e "${GREEN}Reads directory: $ReadDir${NOCOLOUR}";else ReadDir="nil";fi
#	if grep -i -q "Adapters" $ParFile; then Adapter=$(grep -i "Adapters" $ParFile | cut -f2); echo -e "${GREEN}Adapters to trim: $Adapter${NOCOLOUR}";else Adapters="nil";fi
#	if grep -i -q "Reference genome" $ParFile; then RefGenome=$(grep -i "Reference genome" $ParFile | cut -f2); echo -e "${GREEN}Reference genome: $RefGenome${NOCOLOUR}";else RefGenome="nil";fi
#	if grep -i -q "Prior contigs" $ParFile; then PriorContig=$(grep -i "Prior contigs" $ParFile | cut -f2); echo -e "${GREEN}Prior contigs: $PriorContig${NOCOLOUR}";else PriorContig="nil";fi
#	if grep -i -q "Prior sanger" $ParFile; then PriorSanger=$(grep -i "Prior sanger" $ParFile | cut -f2); echo -e "${GREEN}Prior sanger sequence: $PriorSanger${NOCOLOUR}";else PriorSanger="nil";fi
#	if grep -i -q "Blast database" $ParFile; then BlastDB=$(grep -i "Blast database" $ParFile | cut -f2); echo -e "${GREEN}Blast database: $BlastDB${NOCOLOUR}";else BlastDB="nil";fi
#	if grep -i -q "Depth" $ParFile; then Depth=$(grep -i "Depth" $ParFile | cut -f2); echo -e "${GREEN}Subsampling depth: $Depth${NOCOLOUR}"; else Depth="nil";fi
#	if grep -i -q "MLST database" $ParFile && grep -i -q "MLST definitions" $ParFile && grep -i -q "MLST delimiter" $ParFile && grep -i -q "MLST fasta" $ParFile; then
#		MLSTDB=$(grep -i "MLST database" $ParFile | cut -f2)
#		MLSTDelimiter=$(grep -i "MLST delimiter" $ParFile | cut -f2)
#		MLSTfasta=$(grep -i "MLST fasta" $ParFile | cut -f2)
#		MLSTdefinitions=$(grep -i "MLST definitions" $ParFile | cut -f2)
#		echo -e "${GREEN}MLST database: $MLSTDB${NOCOLOUR}"
#	else
#		MLSTDB="nil"
#		MLSTfasta="nil"
#		MLSTdefinitions="nil"
#		MLSTDelimiter="nil"
#	fi




#sleep 1
#else
#	Dir="nil"
#	ReadDir="nil"
#	Adapters="nil"
#	RefGenome="nil"
#	PriorSanger="nil"
#	PriorContig="nil"
#	BlastDB="nil"
#	Depth="nil"
#	MLSTDB="nil"
#	MLSTfasta="nil"
#	MLSTdefinitions="nil"
#	MLSTDelimiter="nil"
#fi

if [ ! -d $DirMetadata ]; then
	mkdir $DirMetadata
fi

if [ ! -d $DirResults ]; then
	mkdir $DirResults
fi

if [ ! -d $DirProkka ]; then
	mkdir $DirProkka
fi
if [ ! -d $DirProkka/gff ]; then
	mkdir $DirProkka/gff
fi

if [ ! -d $DirProkka/gbk ]; then
	mkdir $DirProkka/gbk
fi

if [ ! -d $DirProkka/txt ]; then
	mkdir $DirProkka/txt
fi

#NewProkka=$(printf $(which prokka | sed 's/prokka//g')"Alprokka")
#sed 's/my $MAXCONTIGIDLEN = 20;  # Genbank rule/my $MAXCONTIGIDLEN = 200;  # altered by ALscript from: my $MAXCONTIGIDLEN = 20;  # Genbank rule/g' $(which prokka) | sed 's/contig%06d/contig%04d/g' > #$NewProkka
#chmod 777 $NewProkka

# Get list of genomes to annotate
if [ ! -e $DirMetadata/GenomeList.txt ]; then
	basename -a -s ".fasta" $DirFasta/*.fasta > $DirMetadata/GenomeList.txt
fi
GenomeCount=$(wc -l < "$DirMetadata/GenomeList.txt")

# Annotate genomes with prokka, copy the gff, gbk and txt files to new folders for easier future accessibility
if [ ! -e "$DirResults/Annotation_Summary.tsv" ]; then
	echo -e "Genome\tLength\tCDS\tgene\tmisc_RNA\tmRNA\trRNA\tsig_peptide\ttmRNA\ttRNA" > "$DirResults/Annotation_Summary.tsv"
fi

count="1"
while read i; do
	if [ ! -d $DirProkka/$i ]; then
		date
		echo -e "${BLUE}Annotating ${YELLOW}$i ${BLUE}(${GREEN}$count ${BLUE}of ${GREEN}$GenomeCount ${BLUE}genomes) with prokka${NOCOLOUR}"
		$PROKKA --outdir $DirProkka/$i --quiet --prefix $i --locustag $i --addgenes --genus $GENUS --species $SPECIES --strain $i --kingdom Bacteria --gcode $GCODE --rfam $DirFasta/$i.fasta
	fi
	count=$((count+1))
done < $DirMetadata/GenomeList.txt

while read i; do
	if [ ! -e $DirProkka/gbk/$i.gbk ]; then
		cp $DirProkka/$i/*.gbk $DirProkka/gbk/
	fi
	if [ ! -e $DirProkka/txt/$i.txt ]; then
		cp $DirProkka/$i/*.txt $DirProkka/txt/
	fi
	if [ ! -e $DirProkka/gff/$i.gff ]; then
		cp $DirProkka/$i/*.gff $DirProkka/gff/
	fi
done < $DirMetadata/GenomeList.txt

while read i; do
	if ! grep -q $i "$DirResults/Annotation_Summary.tsv"; then
		echo -e "$i\t$(tail -9 $DirProkka/txt/$i.txt | sort | cut -f2 -d " " | sed '$!{:a;N;s/\n/\t/;ta}')" >> "$DirResults/Annotation_Summary.tsv"
	fi
done < $DirMetadata/GenomeList.txt

if [ ! -d $DirRoary ]; then
	echo -e "${PURPLE}$(date)\n${BLUE}Running Roary pan-genome pipeline${NOCOLOUR}"
	$ROARY -p 4 -f $DirRoary -i 75 -s -e -n -r -v -z $DirProkka/gff/*
	# quiet version
	# roary -p 4 -f $DirRoary -s -e -n -r -z $DirProkka/gff/*
fi


#	basename -a -s ".txt" $DirProkka/txt/*.txt > $DirMetadata/GenomeList.txt

if [ ! -d $Dir/core_genome_alignments ]; then
	mkdir $Dir/core_genome_alignments
fi

if [ ! -e $DirMetadata/GeneList.txt ]; then
	basename -a -s ".fa.aln" $DirRoary/pan_genome_sequences/*.aln > $DirMetadata/GeneList.txt
fi

if [ ! -e $DirMetadata/presence_count.txt ]; then
	sed 's/\",\"/\t/g' "$DirRoary/gene_presence_absence.csv" | cut -f1,4 | sed 's/"//g' | sed s'|/|_|'g > "$DirMetadata/presence_count.txt"
fi

count="0"
echo -e "transferring core gene alignments"
while read i; do
	if [ ! -e $Dir/core_genome_alignments/$i.aln ]; then
		genomes=$(grep -w $i $DirMetadata/presence_count.txt | cut -f2)
		CoreCount=$(grep -c ">" $DirRoary/pan_genome_sequences/$i.fa.aln)
		if [ "$genomes" -eq "$GenomeCount" ] && [ "$CoreCount" -eq "$GenomeCount" ]; then
			printf $i
			printf ", "
			sed 's/\_[^_]*$//g' $DirRoary/pan_genome_sequences/$i.fa.aln > $Dir/core_genome_alignments/$i.aln
			#removes everything after the last underscore (in this case the locus number that prokka adds)
			count=$((count+1))
		fi
	fi
done < $DirMetadata/GeneList.txt
echo -e "\n$count gene files transferred"

#while read i; do genomes=$(grep -w $i $DirRoary/presencecount.txt | cut -f2); if [ "$genomes" -eq "$GenomeCount" ]; then sed 's/\_[^_]*$//g' $DirRoary/pan_genome_sequences/$i.fa.aln > $Dir/core_genome_alignments/$i.aln; fi; done < $DirRoary/GeneList.txt

if [ ! -e $DirMetadata/CoreGeneList.txt ]; then
	basename -a -s ".aln" $Dir/core_genome_alignments/*.aln > $DirMetadata/CoreGeneList.txt
fi

if [ ! -e $DirResults/Alignment_metrics.tsv ]; then
	echo -e "Gene\tHaplotypes\tSegregation Sites\tNucleotide Diversity\tTajima's D\tP-value" > $DirResults/Alignment_metrics.tsv
fi

echo -e "${PURPLE}$(date)\n${BLUE}Extracting alignment metrics from fasta files${NOCOLOUR}"
count="1"
while read i; do
	printf $i
	printf ", "
	if ! grep -w -q $i $DirResults/Alignment_metrics.tsv; then
		echo -e "$(Rscript $ALIGNMENTMETRICS "$Dir" "$i")" >> $DirResults/Alignment_metrics.tsv
	fi
done < $DirMetadata/CoreGeneList.txt





#model test each gene to get an idea of best models?

if [ ! -e "$DirResults/ModelTest.tsv" ]; then
	example=$(head -1 $DirMetadata/CoreGeneList.txt)
	echo -e "Gene\tLeaves\tMethod$($MODELTEST -tr 4 -w -a -v -p -AIC -s 3 -g 8 -i -f -d $DirRoary/pan_genome_sequences/$example.fa.aln | tail -3 | head -1 | sed 's/\t\t/\t/g')" > "$DirResults/ModelTest.tsv"
fi

if [ ! -d $Dir/CoreGene_modeltest ]; then	
	mkdir $Dir/CoreGene_modeltest
fi

if [ ! -d $Dir/nexus ]; then	
	mkdir $Dir/nexus
fi

if [ ! -d $Dir/phylip ]; then
	mkdir $Dir/phylip
fi

if [ ! -d $Dir/modeltest_trees ]; then	
	mkdir $Dir/modeltest_trees
fi	

if [ ! -d $Dir/Bayes ]; then
	mkdir $Dir/Bayes
fi


count="1"
echo -e "${PURPLE}$(date)\n${BLUE}Running jModeltest2 pipeline on core genes${NOCOLOUR}"
while read i; do
	if [ ! -e "$Dir/CoreGene_modeltest/$i.txt" ]; then
		echo -e "${PURPLE}$(date)\n${BLUE}Running jModeltest2 on gene ${YELLOW}$count ${GREEN}($i)${NOCOLOUR}"
		$MODELTEST -tr 4 -w -v -p -a -AIC -s 3 -g 8 -i -f -d $Dir/core_genome_alignments/$i.aln > $Dir/CoreGene_modeltest/$i.txt
	fi

	if ! grep -w -q $i $DirResults/ModelTest.tsv; then
		echo -e "${PURPLE}$(date)\n${BLUE}Converting alignment for gene ${YELLOW}$count ${GREEN}($i) ${BLUE}to nexus & phylip 
files${NOCOLOUR}"
		CoreCount=$(grep -c ">" $Dir/core_genome_alignments/$i.aln)
#		model=$(tail -1 $Dir/CoreGene_modeltest/$i.txt | sed 's/\+/\t/g' | sed 's/\t\t/\t/g')
		model=$(tail -1 $Dir/CoreGene_modeltest/$i.txt | sed 's/ /\t/g' | sed 's/\t\t/\t/g')
		paupblock=$(awk '/BEGIN PAUP/,/END/' $Dir/CoreGene_modeltest/$i.txt)
		#bayesblock=$(awk '/BEGIN PAUP/,/END/' | sed 's/PAUP/mrbayes/g')
		echo -e "$i\t${CoreCount}\t${model}" >> "$DirResults/ModelTest.tsv"
		seqret fasta::$Dir/core_genome_alignments/$i.aln nexus::$Dir/nexus/$i.nexus
		#echo -e "\n$bayesblock\n" >> $Dir/nexus/$i.nexus
		echo -e "\n$paupblock\n" >> $Dir/nexus/$i.nexus
		#gets best tree from jmodeltest
		grep "Tree for the best" $Dir/CoreGene_modeltest/$i.txt | sed 's/ = /=/g' | cut -f2 -d "=" > $Dir/modeltest_trees/$i.modeltest.tre
	fi

	if [ ! -d $Dir/Bayes/$i ]; then
		echo -e "${PURPLE}$(date)\n${BLUE}Building BayesBlock for gene ${YELLOW}$count ${GREEN}($i)${NOCOLOUR}"
		mkdir $Dir/Bayes/$i
#		sed -i 's/mrbayes/PAUP2/g' $Dir/nexus/$i.nexus
		cat $Dir/nexus/$i.nexus BayesBlock.txt | sed "s/FILE/$i/g" | sed "s|DIR|$Dir\/Bayes\/$i|g" > $Dir/Bayes/$i/$i.mb.nexus
	
		model=$(grep -w $i $DirResults/ModelTest.tsv | cut -f4)
		if [ `echo $model | grep -c "GTR"` -eq "1" ] || [ `echo $model | grep -c "SYM"` -eq "1" ]; then
			nst="6"
		elif [ `echo $model | grep -c "HKY"` -eq "1" ] || [ `echo $model | grep -c "K80"` -eq "1" ]; then
			nst="2"
		else
			nst="1"
		fi
		sed -i "s/NST/$nst/g" $Dir/Bayes/$i/$i.mb.nexus
		if [ `echo $model | grep -c "+G"` -eq "1" ] && [ `echo $model | grep -c "+I"` -eq "1" ]; then
			rates="invgamma"
		elif [ `echo $model | grep -c "+G"` -eq "1" ] && [ `echo $model | grep -c "+I"` -ne "1" ]; then
			rates="gamma"
		elif [ `echo $model | grep -c "+G"` -ne "1" ] && [ `echo $model | grep -c "+I"` -eq "1" ]; then
			rates="propinv"
		else
			rates="Equal"
		fi
		sed -i "s/RATES/$rates/g" $Dir/Bayes/$i/$i.mb.nexus
		if [ `echo $model | grep -c "JC"` -eq "1" ] || [ `echo $model | grep -c "SYM"` -eq "1" ] || [ `echo $model | grep -c "K80"` -eq "1" ]; then
			statefreqpr="fixed(equal)"
		else
			statefreqpr="Dirichlet(1.0,1.0,1.0,1.0)"
		fi
		sed -i "s/STATEFREQPR/$statefreqpr/g" $Dir/Bayes/$i/$i.mb.nexus
		echo -e "model=$model\tSubs = $nst\trates = $rates\tstatefreqpr = $statefreqpr"
	fi
	count=$((count+1))
done < $DirMetadata/CoreGeneList.txt

#perhaps run all samples through then tally the most common result and just use that for all genes?





#take roary output and convert it to a mr bayes file, using best model approach
#if [ -d $Dir/Bayes ]; then
#	mkdir $Dir/Bayes
#fi
#count="1"
#while read i; do
#	if [ ! -d $Dir/Bayes/$i ]; then
#		echo $i
#		echo $count
#		mkdir $Dir/Bayes/$i
#		sed -i 's/mrbayes/PAUP2/g' $Dir/nexus/$i.nexus
#		cat $Dir/nexus/$i.nexus BayesBlock.txt | sed "s/FILE/$i/g" | sed "s/DIR/$Dir\/Bayes\/$i/g" > $Dir/Bayes/$i/$i.mb.nexus
#	
#		model=$(grep -w $i $DirResults/ModelTest.tsv | cut -f4)
#		if [ `echo $model | grep -c "GTR"` -eq "1" ] || [ `echo $model | grep -c "SYM"` -eq "1" ]; then
#			nst="6"
#		elif [ `echo $model | grep -c "HKY"` -eq "1" ] || [ `echo $model | grep -c "K80"` -eq "1" ]; then
#			nst="2"
#		else
#			nst="1"
#		fi
#		sed -i "s/NST/$nst/g" $Dir/Bayes/$i/$i.mb.nexus
#		if [ `echo $model | grep -c "+G"` -eq "1" ] && [ `echo $model | grep -c "+I"` -eq "1" ]; then
#			rates="invgamma"
#		elif [ `echo $model | grep -c "+G"` -eq "1" ] && [ `echo $model | grep -c "+I"` -ne "1" ]; then
#			rates="gamma"
#		elif [ `echo $model | grep -c "+G"` -ne "1" ] && [ `echo $model | grep -c "+I"` -eq "1" ]; then
#			rates="propinv"
#		else
#			rates="Equal"
#		fi
#		sed -i "s/RATES/$rates/g" $Dir/Bayes/$i/$i.mb.nexus
#		
#		if [ `echo $model | grep -c "JC"` -eq "1" ] || [ `echo $model | grep -c "SYM"` -eq "1" ] || [ `echo $model | grep -c "K80"` -eq "1" ]; then
#			statefreqpr="fixed(equal)"
#		else
#			statefreqpr="Dirichlet(1.0,1.0,1.0,1.0)"
#		fi
#		sed -i "s/STATEFREQPR/$statefreqpr/g" $Dir/Bayes/$i/$i.mb.nexus
#	fi
#	count=$((count+1))
#done < $DirMetadata/CoreGeneList.txt


if [ ! -d $Dir/BayesTrees ]; then
	mkdir $Dir/BayesTrees
fi 

if [ ! -e $DirResults/PRSF_stats.tsv ]; then
		echo -e "Gene\tMedian\tMin\tMax\tConvergence" > $DirResults/PRSF_stats.tsv
fi


#run nexus alignments though mr bayes
count="1"
echo -e "${PURPLE}$(date)\n${BLUE}Running MrBayes pipeline on core genes${NOCOLOUR}"
while read i; do
	if [ ! -e $Dir/Bayes/$i/$i.con.tre ] && [ `grep -c ">" $Dir/core_genome_alignments/$i.aln` -eq "$GenomeCount" ]; then
		echo -e "${PURPLE}$(date)\n${BLUE}Running MrBayes on gene ${YELLOW}$count ${GREEN}($i)${NOCOLOUR}"
		mb $Dir/Bayes/$i/$i.mb.nexus
	fi	
	count=$((count+1))
done < $DirMetadata/CoreGeneList.invert.txt

count="1"
echo -e "${PURPLE}$(date)\n${BLUE}Extracting final trees in newick format${NOCOLOUR}"
while read i; do
		if [ ! -e $Dir/BayesTrees/$i.newick.tre ]; then
			printf $i
			printf ", "
			echo -e "$(Rscript $NEWICKCONVERTER $Dir/Bayes/$i/ $i)" >> $DirResults/PRSF_stats.tsv
			mv $Dir/Bayes/$i/$i.newick.tre $Dir/BayesTrees/$i.newick.tre
		fi
	count=$((count+1))
done < $DirMetadata/CoreGeneList.txt



#take the bayes tree and convert it to newick format
### To convert to Newick format
## Made a tree conversion program in R, saved this as convert_tree.R
## install.packages("ape")
## library(ape)
## tree <- read.nexus(file="in.tre")
## write.tree(tree,file="temp.R.tre")
## Loop through the Mr Bayes output trees and save with a new name. For this to work, you just have to make sure your trees end in ".nexus.mb.con.tre".




#/usr/local/bin/ete3 compare -t $words -r /home/alistair/Snippy-core_noref/Bayes/SNPcore.newick.tre --unrooted --taboutput
if [ ! -e $DirMetadata/TreeList.txt ]; then
	basename -a $Dir/BayesTrees/*.tre > $DirMetadata/TreeList.txt
fi

words=""
while read i; do 
	words="$words $i,"
done < $DirMetadata/TreeList.txt

echo "$words" > $DirMetadata/CommaTreeList.txt

sed -i 's/\,[^,]*$//g' $DirMetadata/CommaTreeList.txt


#run the trees through ete3 and compare with snp tree?
if [ ! -e $DirResults/treecompare.tsv ]; then
	echo -e "source\tref\tE.size\tnRF\tRF\tmaxRF\tsrc-branches\tref-branches\tsubtrees\ttreekoD" > $DirResults/treecompare.tsv
fi

count="1"
echo -e "${PURPLE}$(date)\n${BLUE}Running ete3 tree comparison pipeline${NOCOLOUR}"
while read i; do
	if ! grep -w -q $i $DirResults/treecompare.tsv; then
		echo -e "${BLUE}Comparing gene ${YELLOW}$count ${GREEN}($i) ${BLUE}to reference tree${NOCOLOUR}"
		echo -e "$i\tSNPtree\t$($ETE compare -t $Dir/BayesTrees/$i.newick.tre -r $REFERENCETREE --unrooted --taboutput | tail -1 | cut -f 3-10)" >> $DirResults/treecompare.tsv
	fi
	count=$((count+1))
done < $DirMetadata/CoreGeneList.txt










#run the trees through codeml with the alignment file

if [ ! -d $Dir/paml ]; then
	mkdir $Dir/paml
fi

if [ ! -d $Dir/paml/trees ]; then
	mkdir $Dir/paml/trees
fi
if [ ! -d $Dir/paml/alignment ]; then
	mkdir $Dir/paml/alignment
fi
if [ ! -d $Dir/paml/CTL ]; then
	mkdir $Dir/paml/CTL
fi
if [ ! -d $Dir/paml/output ]; then
	mkdir $Dir/paml/output
fi



count="1"
while read i; do
	if [ -e $Dir/BayesTrees/$i.newick.tre ] && [ ! -e "$Dir/paml/CTL/$i.ctl" ]; then
		echo -e "${PURPLE}$(date)\n${BLUE}Processing paml input for gene ${GREEN}$i ${YELLOW}(#$count)${NOCOLOUR}"
		if [ ! -e $Dir/paml/trees/$i.newick.tre ]; then
			cp $Dir/BayesTrees/$i.newick.tre $Dir/paml/trees/$i.newick.tre
			while read new old; do
				if [ "$new" != "$old" ]; then
				sed -i "0,/${old}:/{s/${old}:/$new:/}" $Dir/paml/trees/$i.newick.tre
				fi
			done < treecorrect.txt
		fi
		if [ ! -e $Dir/paml/alignment/$i.phylip ]; then
			seqret fasta::$Dir/core_genome_alignments/$i.aln phylip::$Dir/phylip/$i.phylip 
			sed -i '1!b;s/$/\ I/g' $Dir/phylip/$i.phylip
			while read new old; do
				if [ "$GCODE" -eq "11" ]; then
					sed -i "0,/${old}AT/{s/${old}AT/$new\ \ AT/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}GTG/{s/${old}GTG/$new\ \ \GTG/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}TTG/{s/${old}TTG/$new\ \ \TTG/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}CTG/{s/${old}CTG/$new\ \ \CTG/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}at/{s/${old}at/$new\ \ at/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}gtg/{s/${old}gtg/$new\ \ \gtg/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}ttg/{s/${old}ttg/$new\ \ \ttg/}" $Dir/phylip/$i.phylip
					sed -i "0,/${old}ctg/{s/${old}ctg/$new\ \ \ctg/}" $Dir/phylip/$i.phylip
				fi
				sed -i "0,/${old}\-/{s/${old}\-/$new\ \ \-/}" $Dir/phylip/$i.phylip		
			done < nexus2phylip.txt
			cp $Dir/phylip/$i.phylip $Dir/paml/alignment/$i.phylip
		fi
		sed "s|INPUTFILE|$Dir\/paml\/alignment\/$i.phylip|g" $Dir/paml/example.ctl | sed "s|TREEFILE|$Dir\/paml\/trees\/$i.newick.tre|g" | sed "s|OUTPUTFILE|$Dir\/paml\/output\/$i.mcl|g" > "$Dir/paml/CTL/$i.ctl"
	fi
	count=$((count+1))
done < $DirMetadata/CoreGeneList.txt





echo -e "${PURPLE}$(date)\n${BLUE}Running PAML pipeline${NOCOLOUR}"

#count="1"
#while read i; do
#	if [ -e "$Dir/paml/CTL/$i.ctl" ]; then
#		echo -e "conducting codeml for gene $i (#$count)"
#		codeml "$Dir/paml/CTL/$i.ctl"
#		count=$((count+1))
#done < $DirMetadata/CoreGeneList.txt

#ls "$Dir/paml/CTL/"*".ctl" | while read i; do echo $i; codeml $i; done
if [ ! -e $DirResults/Selection_metrics.tsv ]; then
	echo -e "Gene\tomega\tkappa" > $DirResults/Selection_metrics.tsv
fi
count="1"
for i in `less $DirMetadata/CoreGeneList.invert.txt`; do
	if [ -e "$Dir/paml/CTL/$i.ctl" ] && [ ! -e "$Dir/paml/output/$i.mcl" ]; then
		echo -e "${PURPLE}$(date)\n${BLUE}Processing paml input for gene ${GREEN}$i ${YELLOW}(#$count)${NOCOLOUR}"
		codeml "$Dir/paml/CTL/$i.ctl" <<< "\n"
	fi
	if ! grep -w -q $i $DirResults/Selection_metrics.tsv; then
		echo -e "${BLUE}Extracting selection metrics for ${GREEN}$i${NOCOLOUR}"
		echo -e "$i\t$(grep "(dN/dS) =" "$Dir/paml/output/$i.mcl" | cut -f2 -d "=" | sed 's/\ //g')\t$(grep "tv) =" "$Dir/paml/output/$i.mcl" | cut -f2 -d "=" | sed 's/\ //g')" >> $DirResults/Selection_metrics.tsv
	fi
	count=$((count+1))
done



#START=$(date +%s.%N); sleep 2s; END=$(date +%s.%N); dateshow=$(echo "$END - $START" | bc); echo -e "$START\n$END\n$dateshow"




