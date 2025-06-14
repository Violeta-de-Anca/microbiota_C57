#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00:00
##SBATCH --ntasks-per-core 8
#SBATCH -J reasssembly
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/reassembly/reassembly_megahit_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/reassembly/reassembly_megahit_%j.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path
b=${1%/*} #the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly="$temp_dir/assembly.fasta"
temp_assembly1="$temp_dir/assembly_1.fastq"
temp_assembly2="$temp_dir/assembly_2.fastq"

pre_output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}/reassembled_genomes
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/binsABC
output_folder=$pre_output_folder/megahit

echo $a

# Decompress the assembly into the temporary file
zcat $1/final_pure_reads_1.fastq.gz > "$temp_assembly1"
zcat $1/final_pure_reads_2.fastq.gz > "$temp_assembly2"

rm -r $output_folder
mkdir -p $pre_output_folder
mkdir -p $output_folder

#metawrap reassemble_bins -o $output_folder -b $bins_folder -1 $temp_assembly1 -2 $temp_assembly2

