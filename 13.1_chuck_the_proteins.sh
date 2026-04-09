#!/bin/bash
input=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/all_proteins_multigen.faa
outdir=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks

mkdir -p $outdir

N_CHUNKS=600

total=$(grep -c ">" $input)
echo "Total sequences: $total"

awk -v chunk=$(( (total + N_CHUNKS - 1) / N_CHUNKS )) -v outdir="$outdir" '
  /^>/ {
    count++
    if ((count-1) % chunk == 0) {
      outfile = outdir "/chunk_" sprintf("%04d", int((count-1)/chunk)) ".faa"
    }
    print > outfile; next
  }
  { print > outfile }
' $input

echo "Chunks created: $(ls $outdir/*.faa | wc -l)"
