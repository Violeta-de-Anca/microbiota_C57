#!/bin/bash -l
#SBATCH -A naiss2023-22-162
#SBATCH -p core -n 1
##SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J QC
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/QC_raw_fastq_files.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/QC_raw_fastq_files.out
#SBATCH --mail-type=FAIL,COMPLETED
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
#module load Stacks/2.62
#module load cutadapt/4.0
module load MultiQC
module load samtools/1.14
#module load bowtie2/2.3.5.1
module load FastQC/0.11.9
#module load picard
module load bcftools/1.17
module load QualiMap/2.2.1
#module unload java

# Define working directory
working_dir=/proj/naiss2024-23-57/C57_female_lineage_microbiota
raw_samples=$working_dir/samples/cecum_samples/raw_files_F0_till_f2/X204SC23116322-Z01-F001/undetermined

# Define output directory
qc_out=$working_dir/quality_control
fastqc_out=$qc_out/fastqc
multiqc_out=$qc_out/multiqc
mkdir -p $qc_out
mkdir -p $fastqc_out
mkdir -p $multiqc_out

# First QC to see how the sequecing looks like #
cd $raw_samples
fastqc -o $fastqc_out $raw_samples/*
multiqc -o $multiqc_out $raw_samples/*

file_list=($trimmed_input/*.fq.gz)
