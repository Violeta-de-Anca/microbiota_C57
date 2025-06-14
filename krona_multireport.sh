#!/bin/bash
# This script assumes all input files are in the current directory.
# It will create for each file an output file “filename.A” containing:
#   columns1-8: the original file’s columns 2–9 (joined by “|”)
#   column9: the ratio = (sum of col1 for that key)/(total sum of col1)
#input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter
#input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control
#input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/last_feces/small_litter
input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/last_feces/control

# Step 1. Process each file to create fileA (one per file)
#for file in $input_folder/F*M.krona; do
for file in $input_folder/LF*e.krona; do
  b=${file##*/}
  c=${b%.krona}
  # Only process plain files that do not already have an .krona extension
  if [[ -f "$file"  ]]; then
    awk '{
        total += $1;
        # build composite key from columns 2-10 separated by a pipe
        key =$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10; #all the levels until subspecies level
        sum[key] += $1;
      }
      END {
        for (k in sum)
          # Print the key and the ratio (6 decimal digits)
          printf "%.12f\t%s\n", sum[k], k;
      }' $file > $input_folder/${c}.ratio.krona
  fi
done

# Step 2. Get the union of all keys from all fileA files.
# (The key here is the composite string built from original columns 5–10)
#cat $input_folder/F*M.krona | cut -f 2-10 | sort -u | uniq | sed 's/\t/ /g' > $input_folder/keys.txt
cat $input_folder/LF*e.krona | cut -f 2-10 | sort -u | uniq | sed 's/\t/ /g' > $input_folder/keys.txt

# Create file B.
# The header row: first column is “key” then one column per input file.
header="key"
for f in $input_folder/*ratio.krona; do
  # Remove the .A extension for the header
  base="${f%.ratio.krona}"
  header="${header}\t${base}"
done
#echo -e "$header" > $input_folder/small_megahit_litter_report_cecum_feces.txt
#echo -e "$header" > $input_folder/control_megahit_report_cecum_feces.txt
#echo -e "$header" > $input_folder/small_litter_megahit_report_last_feces.txt
echo -e "$header" > $input_folder/control_megahit_report_last_feces.txt

# For each unique key, go through each ratio file and extract the ratio if present.
while IFS= read -r key; do
  line="${key}"
  for f in $input_folder/*ratio.krona; do
    # Look up the key in the current file’s .A file.
    ratio=$(awk -F'\t' -v k="$key" '$2==k {print $1}' "$f")
    # If not found, set ratio to 0.
    if [ -z "$ratio" ]; then ratio="NA"; fi
    line="${line}\t${ratio}"
  done
  #echo -e "$line" >> $input_folder/small_megahit_litter_report_cecum_feces.txt
  #echo -e "$line" >> $input_folder/control_megahit_report_cecum_feces.txt
  #echo -e "$line" >> $input_folder/small_litter_megahit_report_last_feces.txt
  echo -e "$line" >> $input_folder/control_megahit_report_last_feces.txt
done < $input_folder/keys.txt

echo "Done"

