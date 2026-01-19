#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 64GB
#SBATCH -t 2-00:00:00
#SBATCH -J 01plate7
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_plate7_host_removal.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_plate7_host_removal.out


module load metaWRAP/1.4-20230728-foss-2024a-Python-2.7.18

echo $1

#input files:
a=${1##*/} #without path
b=${1%/*} #the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

echo $temp_dir

# Create a temporary file on the scratch storage for the decompressed assembly
temp_fastq_1=$temp_dir/fastq_1.fastq
temp_fastq_2=$temp_dir/fastq_2.fastq

input_1=$1/${a}_1.fastq.gz
input_2=$1/${a}_2.fastq.gz

output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed
mkdir -p $output_folder/${a}

echo $input_1
echo $input_2
echo $a

# Decompress the fastq files into the temporary file
zcat $input_1 > $temp_fastq_1 || { echo "failed decompression";exit 1; }
zcat $input_2 > $temp_fastq_2 || { echo "failed decompression";exit 1; }

metaWRAP read_qc -1 $temp_fastq_1 -2 $temp_fastq_2 -o $output_folder/${a} -x mm39
