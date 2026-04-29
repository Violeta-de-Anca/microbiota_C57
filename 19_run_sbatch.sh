#!/bin/bash
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 3:00:00
#SBATCH -J weight
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/relational_table_abundances_functionals_multigen_C57.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/relational_table_abundances_functionals_multigen_C57.out


BASE=/proj/naiss2024-23-57/C57_female_lineage_microbiota
FUNC=${BASE}/functional_annotation

INPUT=${FUNC}/cluster_bin_protein_proportion_matrix.tsv
OUTDIR=${BASE}/quantification

module load Python/3.13.5-GCCcore-14.3.0 SciPy-bundle/2025.07-gfbf-2025b matplotlib/3.10.5-gfbf-2025b Seaborn/0.13.2-gfbf-2025b networkx/3.5-gfbf-2025b

echo "=========================================="
echo "Job started : $(date)"
echo "Node        : $(hostname)"
echo "Input       : ${INPUT}"
echo "Output dir  : ${OUTDIR}"
echo "=========================================="

python 19_relational_table_abundance_functionality_multigen.py all
