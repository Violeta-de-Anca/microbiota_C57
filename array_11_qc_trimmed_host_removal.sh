#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_11.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_11.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

#for last_feces
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed

for F in $input_path/*; do
	[[ -d "$F" ]] || continue
        echo $F
	sbatch --export=ALL,a=$F qc_trimmed_host_removal.sh $F
done

#for cecum samples
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed

for F in $input_path/*; do
        [[ -d "$F" ]] || continue
        echo $F
        sbatch --export=ALL,a=$F qc_trimmed_host_removal.sh $F
done
