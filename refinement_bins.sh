#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH --ntasks-per-core 8
#SBATCH -J refin_%j
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/refinement/bins_refinement_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/refinement/bins_refinement_%j.out

module load bioinfo-tools metaWRAP/1.3.2

#for doing it per generations/sample type



#for doing the refinement individually

a=${1##*/} #without path
b=${1%/*} #the path

#folder where there are the different bins
folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

#output folder
mkdir -p $folder/$a/refined_libraries_megahit

output_folder=$folder/$a/refined_libraries_megahit

metawrap bin_refinement -o $output_folder -t 8 -m 64 -A $folder/$a/metaBAT2/metabat2_bins/ -B $folder/$a/maxbin2_bins/ -C $folder/$a/concoct/concoct_bins/
