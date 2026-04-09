#!/bin/bash
##SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
#SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_13.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_13.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks

#for i in $(seq 0 599); do
#	F=$(printf "%04d" $i)
#	echo $F
#	sbatch --export=ALL,a=$F 13_psi_blast.sh $F
#done

for i in $(seq 64 99); do
        F=$(printf "%04d" $i)
        echo $F
        sbatch --export=ALL,a=$F 13_psi_blast.sh $F
done
