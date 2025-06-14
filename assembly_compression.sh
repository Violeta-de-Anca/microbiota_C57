#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J zip_assembly

gzip /proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/*fasta
