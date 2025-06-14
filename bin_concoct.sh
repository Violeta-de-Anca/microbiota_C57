#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH --ntasks-per-core 8
#SBATCH -J b_BAT_%j
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bins_concoct_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bins_concoct_%j.out

module load bioinfo-tools metaWRAP/1.3.2

a=${1##*/} #without path
b=${1%/*} #the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly=$(mktemp "$temp_dir/assembly.XXXXXX.fasta")

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}.megahit.assembly.fasta
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

echo $input_megahit
echo $a

# Decompress the assembly into the temporary file
zcat "$input_megahit" > "$temp_assembly"

mkdir -p $output_folder/$a/concoct

metawrap binning -o $output_folder/$a/concoct -a $temp_assembly --concoct -m 64 -t 8 $1/final_pure_reads_*

rm "$temp_assembly"
