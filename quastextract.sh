#!/bin/bash

Dir=$1

if [ ! -d "$Dir/Results_Summary" ]; then
	mkdir "$Dir/Results_Summary"
fi

while read i; do
	if [ -e "$Dir/quastoutput/$i/contigs_reports/contigs_report_$i.contigs.mis_contigs.info" ]; then
		echo "$i" >> "$Dir/Results_Summary/Misassemblies_details.txt"
		echo "$Dir/quastoutput/$i/contigs_reports/contigs_report_$i.contigs.mis_contigs.info" >> "$Dir/Results_Summary/Misassemblies_details.txt"
		echo -e "\n===================================\n" >> "$Dir/Results_Summary/Misassemblies_details.txt"
	fi
	if [ ! -e "$Dir/Results_Summary/Misassemblies_summary.tsv" ]; then
		cat "$Dir/quastoutput/$i/contigs_reports/transposed_report_misassemblies.tsv" > "$Dir/Results_Summary/Misassemblies_summary.tsv"
	else
		grep $i "$Dir/quastoutput/$i/contigs_reports/transposed_report_misassemblies.tsv" >> "$Dir/Results_Summary/Misassemblies_summary.tsv"
	fi
	if [ ! -e "$Dir/Results_Summary/quast_summary.tsv" ]; then	
	cat "$Dir/quastoutput/$i/transposed_report.tsv" > "$Dir/Results_Summary/quast_summary.tsv"
	else
		grep $i "$Dir/quastoutput/$i/transposed_report.tsv" >> "$Dir/Results_Summary/quast_summary.tsv"
	fi


	columns=$(awk '{print NF}' "$Dir/quastoutput/$i/contigs_reports/alignments_$i.contigs.tsv" | sort -nu | tail -n 1)
	columns=$((columns-1))
	gaps=$(wc -l < "$Dir/quastoutput/$i/genome_stats/$i.contigs_gaps.txt")
	gaps=$((gaps-1))
	if [ ! -e "$Dir/Results_Summary/contig_coverage.tsv" ]; then
		echo -e "Sample\tcontigs aligned to reference\tGaps detected" > "$Dir/Results_Summary/contig_coverage.tsv"
		echo -e "$i\t$columns\t$gaps" >> "$Dir/Results_Summary/contig_coverage.tsv"
	else
		echo -e "$i\t$columns\t$gaps" >> "$Dir/Results_Summary/contig_coverage.tsv"
	fi

done < "$Dir/Metadata/GenomeList.txt"
