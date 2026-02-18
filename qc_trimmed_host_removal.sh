#!/bin/bash -l
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 24:00:00
#SBATCH -J QC_trim
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/qc_trimmed_host_removed.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/qc_trimmed_host_removed.out

module load FastQC/0.12.1-Java-17

fastqc_out=/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter/trimmed_host_removed
a=${1##*/} #without path

mkdir -p $fastqc_out/${a}

fastqc -o $fastqc_out/${a} $1/*
