#!/bin/bash
# Define the source directory and destination directory
src=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy
dest=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces
# Create the destination directory if it does not exist
mkdir -p "$dest"
# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with: - 'F' followed by 0, 1, or 2, - then an underscore, - then any number of characters (which should be the numeric part), - 
# then an underscore, - then either A or B, - then _male_M.
for d in "$src"/F[0-2]_*_[AB]_*_M; do
    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        file="$d/final_pure_reads.krona"
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus the .html extension.
            cp "$file" "$dest/${base}.krona"
        else
            echo "Warning: kronagram.html not found in $d"
        fi
    fi
done

mv $dest/*A* $dest/small_litter/.
mv $dest/*B* $dest/control/.

#do per generation the reports and per group
#small litter, cecum
#module load bioinfo-tools python3
#python KrakenTools/combine_kreports.py -r $dest/F0_103_A_female_M.krona $dest/F0_12_A_female_M.krona $dest/F0_22_A_male_M.krona $dest/F0_23_A_male_M.krona $dest/F0_24_A_male_M.krona $dest/F0_31_A_male_M.krona $dest/F0_33_A_female_M.krona $dest/F0_34_A_female_M.krona -o $dest/F0_small_litter.krona
#python KrakenTools/combine_kreports.py -r $dest/F1_1033_A_male_M.krona $dest/F1_111_A_male_M.krona $dest/F1_112_A_female_M.krona $dest/F1_121_A_female_M.krona $dest/F1_122_A_female_M.krona $dest/F1_123_A_male_M.krona $dest/F1_341_A_female_M.krona -o $dest/F1_small_litter.krona

#control, cecum

dest=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/last_feces
# Create the destination directory if it does not exist
mkdir -p "$dest"
# Loop over the subfolders that match the pattern. The glob F[0-2]_*_[AB]_male_M matches folders with: - 'F' followed by 0, 1, or 2, - then an underscore, - then any number of characters (which should be the numeric part), -
# then an underscore, - then either A or B, - then _male_M.
for d in "$src"/LF_F[0-2]_*_[AB]_*; do
    # Make sure that $d is a directory.
    if [ -d "$d" ]; then
        # Define the file to be copied.
        file="$d/final_pure_reads.krona"
        # Check if the kronagram.html file exists.
        if [ -f "$file" ]; then
            # Get the base folder name.
            base=$(basename "$d")
            # Copy the kronagram.html to the destination folder, renaming it to the folder name plus the .html extension.
            cp "$file" "$dest/${base}.krona"
        else
            echo "Warning: kronagram.html not found in $d"
        fi
    fi
done

mv $dest/*A* $dest/small_litter/.
mv $dest/*B* $dest/control/.
