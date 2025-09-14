#!/bin/bash
#SBATCH -A uppmax2025-2-150
#SBATCH -p node
##SBATCH -p core
#SBATCH -n 1
##SBATCH -t 1:00:00
#SBATCH -t 10-00:00:00
#SBATCH --ntasks-per-core 8
#SBATCH -J b_all_%j
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bins_all_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bins_all_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2

#module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

#we start with the direction to the fastq files
#_last_feces_1.fastq.gz
sample=$1
echo $sample
a=${sample##*/} #without path (F0)
b=${sample%/*} #the path
g=${sample%/*/*} #type of sample with the path
c=${g##*/} #type of sample without the path

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

# Create a temporary file on the scratch storage for the decompressed fastq files
temp_fastq1=$temp_dir/fastq_1.fastq
temp_fastq2=$temp_dir/fastq_2.fastq
temp_assembly=$temp_dir/assembly.fasta

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}_${c}/final_assembly.fasta.gz
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

echo $input_megahit
echo $a
echo $c
echo ${1}_${c}_1.fastq.gz

# Decompress the assembly and the fastq files into the temporary files
zcat ${1}_${c}_1.fastq.gz > $temp_fastq1 || { echo "failed decompression";exit 1; }
zcat ${1}_${c}_2.fastq.gz > $temp_fastq2 || { echo "failed decompression";exit 1; }
zcat $input_megahit > $temp_assembly || { echo "failed decompression";exit 1; }

mkdir -p $output_folder/${a}_$c

metawrap binning -o $output_folder/${a}_$c \
	-a $temp_assembly --concoct --metabat2 \
	--maxbin2 -m 64 -t 8 $temp_fastq1 $temp_fastq2

rm $temp_fastq1
rm $temp_fastq2

