#!/bin/bash
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_12.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_12.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

for F in F0 F1 F2; do
	for i in cecum_samples last_feces; do
		[[ -d "$input_path/${F}_${i}" ]] || continue
        	echo $input_path/${F}_${i}
        	sbatch --export=ALL,a=$input_path/${F}_${i} 12_functional_table_clean.sh $input_path/${F}_${i}
	done
done

