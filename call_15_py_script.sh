#!/bin/bash
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 50GB
#SBATCH -t 1:00:00
#SBATCH -J call.py
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/call_py_15_%j.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/call_py_15_%j.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks

module load Python/3.13.5-GCCcore-14.3.0 Biopython IPython Python-bundle-PyPI bx-python

F=$1
input=$input_path/psi_results_chunck_${F}_multigen.txt
echo $input
output=$input_path/hssp_chunk_${F}.txt
#check that the output does not exist
if [ -f "$output" ] && [ -s "$output" ]; then
	echo "[SKIP] $output already exists and is non-empty."
	echo "=========================================="
	echo "Array task : ${SLURM_ARRAY_TASK_ID}  (chunk ${F})"
	echo "Node       : $(hostname)"
	echo "Input      : $input"
	echo "Output     : $output"
	echo "Started    : $(date)"
	echo "=========================================="
	exit 0
fi

echo "=========================================="
echo "Array task : ${SLURM_ARRAY_TASK_ID}  (chunk ${F})"
echo "Node       : $(hostname)"
echo "Input      : $input"
echo "Output     : $output"
echo "Started    : $(date)"
echo "=========================================="

#check that the input does exist
if [ ! -f "$input" ]; then
	echo "[ERROR] Input file not found: $input"
	exit 1
fi

#run the python script
python 15_hssp_matrix.py --input $input --output $output --threshold 10
