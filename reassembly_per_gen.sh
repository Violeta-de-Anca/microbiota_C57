#!/bin/bash
#SBATCH -A uppmax2025-2-150
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J reassem
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/reassembly_megahit_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/reassembly_megahit_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path (F 0/1/2)
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed/)
c=${1%/*/*} #sample type with path
type=${c##*/} #only sample type

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

#create a temporary file for the fastq files
temp_fastq1=$temp_dir/fastq_1.fastq
temp_fastq2=$temp_dir/fastq_2.fastq

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}_${type}/final_assembly.fasta
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}_${type}/reassembly
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}_${type}/refined_libraries_megahit/metawrap_70_10_bins
folder_fastq=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/${type}/trimmed_host_removed

echo $input_megahit
echo $a
echo $output_folder
echo $bins_folder

# Decompress the assembly into the temporary file
zcat $folder_fastq/${a}_${type}_1.fastq.gz > $temp_fastq1
zcat $folder_fastq/${a}_${type}_2.fastq.gz > $temp_fastq2

mkdir -p $output_folder

metawrap reassemble_bins -o $output_folder -b $bins_folder -1 $temp_fastq1 -2 $temp_fastq2
