#!/bin/bash
#SBATCH -A uppmax2025-2-150
#SBATCH -p node
##SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J quant
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path, generation and type of sample
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics)
gen=${a%%_*}
sample=${a#*_}

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

cp 

checkm lineage_wf
