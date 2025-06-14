#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH --ntasks-per-core 8
#SBATCH -J quant
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_%j.out

module load bioinfo-tools metaWRAP/1.3.2

a=${1##*/} #without path
b=${1%/*} #the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly="$temp_dir/assembly.fasta"
temp_assembly1="$temp_dir/assembly_1.fastq"
temp_assembly2="$temp_dir/assembly_2.fastq"

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}.megahit.assembly.fasta
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/quantification
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/binsABC

echo $input_megahit
echo $a

# Decompress the assembly into the temporary file
zcat "$input_megahit" > "$temp_assembly"
zcat $1/final_pure_reads_1.fastq.gz > "$temp_assembly1"
zcat $1/final_pure_reads_2.fastq.gz > "$temp_assembly2"


mkdir -p $output_folder


metawrap quant_bins -b $bins_folder -o $output_folder -t 8 -a $temp_assembly $temp_assembly1 $temp_assembly2

