#!/bin/bash
# Define the source directory and destination directory
src=/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter
dest=/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter
# Create the destination directory if it does not exist
mkdir -p "$dest"
# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with: - 'F' followed by 0, 1, or 2, - then an underscore, - then any number of characters (which should be the numeric part), -
# then an underscore, - then either A or B, - then _male_M.
#for d in "$src"/*.fastqc.zip.folder; do
    # Make sure that $d is a directory.
    # Define the file to be copied.
#        file="$d/*/fastqc_data.txt"
#	echo $file
#        # Check if the kronagram.html file exists.
#        # Get the base folder name.
#       base=$(basename "$d")
#	folder="${base%.fastqc.zip.folder}"
#		awk -v folder="$folder" 'NR==4 {
#			if ($2=="final_pure_reads_1.fastq" || $2=="final_pure_reads_2.fastq")
#			$2 = folder
#		}
#		{ print }' $file > $src/tmp && mv $src/tmp $file
#done

for d in "$src"/*.fastqc.zip.folder; do
	file="$d/*/fastqc_data.txt"
	base=$(basename "$d")
	folder="${base%.fastqc.zip.folder}"
	sed -i '4s/ /\t/' $file
done
