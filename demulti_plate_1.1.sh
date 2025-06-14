#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J dem1.1
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/00_demul_plate1.1.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/00_demul_plate1.1.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2

#input folder
input_fol=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/data_F_09/01.RawData/Undetermined
barcode_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

# Create a temporary file on the scratch storage for the decompressed assembly
temp_fastq_1=$temp_dir/fastq_1.fasta
temp_fastq_2=$temp_dir/fastq_2.fasta


input_1=$input_fol/*_*_L8_1.fq.gz
input_2=$input_fol/*_*_L8_2.fq.gz
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/whole_families_plate1.1


echo $input_1
echo $input_2
echo $a

# Decompress the assembly into the temporary file
zcat $input_1 > $temp_fastq_1
zcat $input_2 > $temp_fastq_2


./demultiplex demux -p $output_folder $barcode_folder/plate_1_1_barcodes.txt $temp_fastq_1 $temp_fastq_2
