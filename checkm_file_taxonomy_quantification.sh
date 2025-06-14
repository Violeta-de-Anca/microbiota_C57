#!/bin/bash

input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
folder_path=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

for suffix in last_feces cecum_samples; do
       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
                #echo $F
                echo $b
                b=${F##*/}
                c=${F%/*}
                #echo $c
		sed -nE "s/^(Refined_[0-9]+)[[:space:]]+\{'marker lineage':[[:space:]]*'([^']+)'/\1 \2/p" $folder_path/$b/refined_libraries_megahit/quantification/$b.metaspades.taxonomy.tsv > $folder_path/$b/refined_libraries_megahit/quantification/$b.megahit.taxonomy.txt
		cut -d , -f 1 $folder_path/$b/refined_libraries_megahit/quantification/$b.megahit.taxonomy.txt > $folder_path/$b/refined_libraries_megahit/quantification/$b.final.txt
		mv $folder_path/$b/refined_libraries_megahit/quantification/$b.final.txt $folder_path/$b/refined_libraries_megahit/quantification/$b.megahit.taxonomy.txt
		mv $folder_path/$b/refined_libraries_megahit/quantification/$b.megahit.taxonomy.txt /proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/megahit/.
	done
done

