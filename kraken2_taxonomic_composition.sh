#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J kraken2_%j
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/kraken2_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/kraken2_%j.out


sample=$1 #path till the fasta files
echo $sample
c=${sample%/*} #this is with the path to the samples
b=${sample##*/} #this is the last folder, which is the sample name


input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes

fastq_1=$c/$b/final_pure_reads_1.fastq

echo $fastq_1

fastq_2=$c/$b/final_pure_reads_2.fastq

echo $fastq_2

assembly=$input_main_path/$b/final_assembly.fasta

module load bioinfo-tools metaWRAP/1.3.2

output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy

mkdir -p $output_path/$b
metawrap kraken2 -o $output_path/$b $assembly $fastq_1 $fastq_2
