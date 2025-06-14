#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -N 1
#SBATCH -t 10-00:00:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_03_kraken2.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_03_kraken2.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes
input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

for suffix in last_feces cecum_samples; do
       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
		echo $F
		echo $b
		b=${F##*/}
		c=${F%/*}
		echo $c
		sbatch --export=ALL,assembly=$F \
		/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/kraken2_taxonomic_composition.sh $F
	done
done
