#!/bin/bash

# Replace 'filename.ext' with the file you are looking for
file_to_search="transposed_report.txt"
folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes
output_file="megahit_contig_assembly_report.txt"
# Loop through all subfolders
for dir in $folder/* ; do
	# Check if the file exists in the current folder
	if [[ -f "$dir/QUAST_out/$file_to_search" ]]; then
	# Echo the folder name and "yes"
		echo "$dir" >> "$output_file"
		cat $dir/QUAST_out/$file_to_search >> $output_file
	else
		echo "$dir: no" >> "$output_file"
fi
done

