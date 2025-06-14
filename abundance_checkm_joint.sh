#!/bin/bash
for tabfile in /proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/megahit/*_megahit_abundance_bins.tab; do
    # Extract unique individual name from tab file name
    id="${tabfile%%_megahit_abundance_bins.tab}"
    txtfile="${id}.megahit.taxonomy.txt"
    outfile="${id}_megahit_checkm_taxonomy_abundance.txt"
    echo $txtfile
    # Check if the taxonomy file exists
    if [[ -f "$txtfile" ]]; then
        # Read the taxonomy file first (FNR==NR) to build an array: key = genomic bin, value = taxonomy.
        # Then process the tab file: skip its header and, for each line, print the genomic bin,
        # its taxonomy (from the .txt file), and the abundance (2nd column of the .tab file).
        awk 'FNR==NR {tax[$1]=$2; next} FNR==1 {next} ($1 in tax) {print $1, tax[$1], $2}' "$txtfile" "$tabfile" > "$outfile"
    else
        echo "Taxonomy file for $id not found" >&2
    fi
done

