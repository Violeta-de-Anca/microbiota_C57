#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J functanno
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/functional_annotation_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/functional_annotation_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

mkdir -p $output_path

a=${1##*/} #without path (F 0/1/2_last_feces)
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes)

metawrap annotate_bins -o $output_path/${a} -b $1/reassembly/reassemblies
