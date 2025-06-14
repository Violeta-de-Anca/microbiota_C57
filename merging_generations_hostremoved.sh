#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J rezip
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_generations_hostremoved.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_generations_hostremoved.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2

