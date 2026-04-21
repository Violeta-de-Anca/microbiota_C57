#!/bin/bash
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 20GB
#SBATCH -t 1:00:00
#SBATCH -J table_mcl
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/format_16_mcl_py.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/format_16_mcl_py.out

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/mcl
output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

module load Python/3.13.5-GCCcore-14.3.0 Biopython IPython Python-bundle-PyPI bx-python

#run the python script
python 16.1_formatting_mcl.py --input $input_path/mcl_clusters.txt --output $output_path/relational_table_functional_bins_multigenerational_microbiota.txt
