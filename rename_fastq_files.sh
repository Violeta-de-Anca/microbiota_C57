#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J rename
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/rename_fastq_files.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/rename_fastq_files.out

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

for suffix in last_feces; do
	for i in $input_main_path/$suffix/trimmed_host_removed/LF_*; do
		b=${i##*/} # file name
		c=${i%/*} #the path
		echo $b
		echo $i
		echo $c
		echo $input_main_path/$suffix/${b}_1.fastq.gz
		cp $i/final_pure_reads_1.fastq.gz $input_main_path/$suffix/${b}_1.fastq.gz
		#cp $i/final_pure_reads_1.fastq $input_main_path/$suffix/${b}_1.fastq
		cp $i/final_pure_reads_2.fastq.gz $input_main_path/$suffix/${b}_2.fastq.gz
		#cp $i/final_pure_reads_2.fastq $input_main_path/$suffix/${b}_2.fastq
	done
done
