#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH --ntasks-per-core 8
#SBATCH -J refin_gen
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/refinement/bins_refinement_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/refinement/bins_refinement_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

#for doing it per generations/sample type
a=${1##*/} #without path
b=${1%/*} #the path

#folder where there are the different bins
folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

#output folder
mkdir -p $1/refined_libraries_megahit

output_folder=$1/refined_libraries_megahit

metawrap bin_refinement -o $output_folder -t 8 -m 64 -A $1/metabat2_bins/ -B $1/maxbin2_bins/ -C $1/concoct_bins/
