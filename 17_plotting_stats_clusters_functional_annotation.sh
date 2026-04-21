#!/bin/bash
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 50GB
#SBATCH -t 1:00:00
#SBATCH -J mcl_stats
##SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_stats.err
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_improve.err
##SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_stats.out
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_improve.out

BASE=/proj/naiss2024-23-57/C57_female_lineage_microbiota
FUNC=${BASE}/functional_annotation

INPUT=${FUNC}/relational_table_functional_bins_multigenerational_microbiota.txt
OUTDIR=${FUNC}/mcl

module load Python/3.13.5-GCCcore-14.3.0 SciPy-bundle/2025.07-gfbf-2025b matplotlib/3.10.5-gfbf-2025b Seaborn/0.13.2-gfbf-2025b networkx/3.5-gfbf-2025b

echo "=========================================="
echo "Job started : $(date)"
echo "Node        : $(hostname)"
echo "Input       : ${INPUT}"
echo "Output dir  : ${OUTDIR}"
echo "=========================================="

#python 17_mcl_stats_plots.py \
#    --input  ${INPUT} \
#    --outdir ${OUTDIR}

#python 17_change.py --input  ${INPUT} --outdir ${OUTDIR}
python 17_bin_02.py --input  ${INPUT} --outdir ${OUTDIR}

echo "Finished: $(date)"
