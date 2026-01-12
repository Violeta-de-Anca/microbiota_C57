#!/bin/bash
#SBATCH -A uppmax2025-2-150
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_010.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_010.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

#for doing per generations/per sample type, but to get individual information
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

for suffix in cecum_samples last_feces; do
	for F in F0 F2 F1 F3 F4 F5; do
		generation=$input_path/${F}_${suffix}
		echo $generation
		sbatch --export=ALL,sample=$generation <> $generation
	done
done
