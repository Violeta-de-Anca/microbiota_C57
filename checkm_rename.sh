#!/bin/bash
# Define the source directory and destination directory
src=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter/checkm_refined_bins_megahit
dest=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with: - 'F' followed by 0, 1, or 2, - then an underscore, - then any number of characters (which should be the numeric part), -
# then an underscore, - then either A or B, - then _male_M.
for d in "$src"/F[0-2]_*_[AB]_*_M; do
    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        file=$d/storage/bin_stats_ext.tsv
	c=${d##*/}
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
	    base=${c%.tsv}
            cp "$file" "$dest/$c/refined_libraries_megahit/quantification/${base}.metaspades.taxonomy.tsv"
        else
            echo "Warning: final_assembly.fasta not found in $d"
        fi
    fi
done

src=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter/checkm_refined_bins_megahit
dest=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics
# Create the destination directory if it does not exist
mkdir -p "$dest"
# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with: - 'F' followed by 0, 1, $
# then an underscore, - then either A or B, - then _male_M.
for d in "$src"/LF_*; do
    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        file=$d/storage/bin_stats_ext.tsv
	c=${d##*/}
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
	    base=${c%.tsv}
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus the .html extension.
            cp "$file" "$dest/$c/refined_libraries_megahit/quantification/${base}.metaspades.taxonomy.tsv"
        else
            echo "Warning: final_assembly.fasta not found in $d"
        fi
    fi
done
