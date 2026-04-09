#!/bin/bash
##SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
#SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_15.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_15.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks

log_files=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files

rm $log_files/call_py_15_*

for i in $(seq 0 599); do
        F=$(printf "%04d" $i)
        output=$input_path/hssp_chunk_${F}.txt
	echo $output
        if [ -f "$output" ] && [ -s "$output" ]; then
                echo "[SKIP] $output already exists and is non-empty."
                continue
        fi
        sbatch --export=ALL,a=$F call_15_py_script.sh $F
done
