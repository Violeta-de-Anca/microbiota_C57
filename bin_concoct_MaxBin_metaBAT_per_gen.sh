#!/bin/bash
#SBATCH -A uppmax2025-2-151
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
temp_assembly1=$temp_dir/assembly_1.fastq
temp_assembly2=$temp_dir/assembly_2.fastq

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}_${c}/final_assembly.fasta
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

echo $input_megahit
echo $a
echo $c
echo ${1}_last_feces_1.fastq.gz

# Decompress the assembly into the temporary file
zcat ${1}_${c}_1.fastq.gz > $temp_assembly1
zcat ${1}_${c}_2.fastq.gz > $temp_assembly2

mkdir -p $output_folder/${a}_$c

metawrap binning -o $output_folder/${a}_$c \
	-a $input_megahit --concoct --metabat2 \
	--maxbin2 -m 64 -t 8 $temp_assembly1 $temp_assembly2

rm $temp_assembly1
rm $temp_assembly2

