#!/bin/bash
# Define the source directory and destination directory
src=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
dest=/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter

# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with:
# - 'F' followed by 0, 1, or 2, - then an underscore, - then any number of characters (which should be
#the numeric part), -
# then an underscore, - then either A or B, - then _male_M.
for d in "$src"/cecum_samples/trimmed_host_removed/F[0-2]_*_[AB]_*_M; do
#    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        #file="$d/final_pure_reads_1_fastqc.zip"
        file=$d/final_pure_reads_1_fastqc.html
	# Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus the .html extension.
            #mv "$file" "$dest/${base}.1.fastqc.zip"
            mv "$file" "$dest/${base}.1.fastqc.html"
	else
            echo "Warning: final_assembly.fasta not found in $d"
        fi
    fi
    if [ -d "$d" ]; then
        # Define the file to be copied.
        #file="$d/final_pure_reads_2_fastqc.zip"
	file=$d/final_pure_reads_2_fastqc.html
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus $
            #mv "$file" "$dest/${base}.2.fastqc.zip"
            mv "$file" "$dest/${base}_2_fastqc.html
	else
            echo "Warning: final_assembly.fasta not found in $d"

        fi
    fi
done

for d in "$src"/last_feces/trimmed_host_removed/LF_*; do
    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        #file="$d/final_pure_reads_1_fastqc.zip"
	file=$d/final_pure_reads_1_fastqc.html
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus the .html extension.
            #mv "$file" "$dest/${base}.1.fastqc.zip"
            mv "$file" "$dest/${base}_1_fastqc.html
	else
            echo "Warning: final_assembly.fasta not found in $d"
        fi
    fi
    if [ -d "$d" ]; then
        # Define the file to be copied.
        #file="$d/final_pure_reads_2_fastqc.zip"
        file=$d/final_pure_reads_2_fastqc.html
	# Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus $
            #mv "$file" "$dest/${base}.2.fastqc.zip"
            mv "$file" "$dest/${base}_2_fastqc.html
	else
            echo "Warning: final_assembly.fasta not found in $d"

        fi
    fi
done

