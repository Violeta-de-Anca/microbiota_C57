#!/bin/bash -l -x
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J 01plate6
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_plate6_host_removal_%j.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_plate6_host_removal_%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2

#input files:
a=${1##*/} #without path
c=${a#*_}
b=${1%/*} #the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

# Create a temporary file on the scratch storage for the decompressed assembly
temp_fastq_1=$temp_dir/fastq_1.fastq
temp_fastq_2=$temp_dir/fastq_2.fastq


input_1=$1/${c}_1.merged.fastq.gz
input_2=$1/${c}_2.merged.fastq.gz
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed

echo $input_1
echo $input_2
echo $a

# Decompress the assembly into the temporary file
zcat $input_1 > $temp_fastq_1 || { echo "failed decompression";exit 1; }
zcat $input_2 > $temp_fastq_2 || { echo "failed decompression";exit 1; }

metawrap read_qc -1 $temp_fastq_1 -2 $temp_fastq_2 -o $output_folder/${a} -x mm39
