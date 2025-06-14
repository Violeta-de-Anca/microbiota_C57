#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10-00:00:00
#SBATCH -J 02
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/02_assembling_metagenomes.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/02_assembling_metagenomes.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes

for suffix in cecum_samples last_feces; do
       for F in $(cat $input_main_path/$suffix/trimmed_host_removed/trimmed.files); do
		b=${F##*/}
		c=${F%/*}
		echo $b
		echo $c
		mkdir -p $output_path/${b}
		metawrap assembly -1 $c/${b}/final_pure_reads_1.fastq  -2 $c/${b}/final_pure_reads_2.fastq \
		-o $output_path/${b} -m 24
	done
done
