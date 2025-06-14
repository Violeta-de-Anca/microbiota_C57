#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00
#SBATCH -J pure
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/host_zip.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/host_zip.out


input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

for suffix in last_feces cecum_samples; do
       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
                echo $F
                echo $b
                b=${F##*/}
                c=${F%/*}
                echo $c
		gzip $F/final_pure_reads_1.fastq
		gzip $F/final_pure_reads_2.fastq
	done
done
