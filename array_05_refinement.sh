#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_05_refinement.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_05_refinement.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

#for doing per generations/per sample type
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics
for suffix in last_feces; do
	for F in F0 F2 F1; do
		generation=$input_path/${F}_${suffix}
		echo $generation
		sbatch --export=ALL,sample=$generation refinement_per_gen.sh $generation
	done
done

# for doing it individually
#input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

#for F in $(cat $input_path/samples_multigenerational.txt); do
#	echo $F
#	sbatch --export=ALL,a=$F refinement_bins.sh $F
#done
