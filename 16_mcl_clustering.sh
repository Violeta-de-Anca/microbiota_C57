#!/bin/bash
#SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 250GB
#SBATCH --ntasks=4 --cpus-per-task=8
#SBATCH -t 2-00:00:00
#SBATCH -J mcl_cluster
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_clustering.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mcl_clustering.out

# =============================================================================
# 16_mcl_clustering.sh
#
# Cluster proteins into functional groups using MCL (Markov Cluster Algorithm)
# on the full merged HSSP distance matrix.
#
# Based on: Zhu et al. (2015) PLOS Comput Biol (doi:10.1371/journal.pcbi.1004472)
#   - Inflation parameter -I 1.4 (as used in the paper)
#   - Input: all proteins from F0/F1/F2 x cecum/feces as one network
#   - Output: protein clusters = functional groups
#
# Pipeline:
#   1. mcxload  — converts the tab-separated edge list into MCL binary format
#                 and builds a dictionary mapping protein IDs to integers
#   2. mcl      — runs the Markov Cluster Algorithm on the binary matrix
#   3. mcxdump  — converts binary cluster output back to readable protein IDs
# =============================================================================

set -euo pipefail

BASE=/proj/naiss2024-23-57/C57_female_lineage_microbiota
FUNC=${BASE}/functional_annotation
OUT=${FUNC}/mcl

mkdir -p ${OUT}

# Input: merged HSSP matrix (output of 15_merge_hssp.sh)
HSSP_MATRIX=${FUNC}/hssp_matrix.txt

# Intermediate files
EDGE_LIST=${OUT}/hssp_edges.txt         # 3-column edge list for mcxload

# Final output (human-readable, one cluster per line, proteins tab-separated)
MCL_CLUSTERS=${OUT}/mcl_clusters.txt

module load MCL/22.282-GCCcore-13.3.0

echo "=========================================="
echo "Job started : $(date)"
echo "Node        : $(hostname)"
echo "SLURM job   : ${SLURM_JOB_ID}"
echo "Input       : ${HSSP_MATRIX}"
echo "Output dir  : ${OUT}"
echo "=========================================="

# --- Step 0: extract the 3-column edge list (qseqid, sseqid, hssp_distance) ---
# Skip the header line, keep only columns 1, 2, 5
echo "[STEP 0] Extracting edge list from HSSP matrix ..."
tail -n +2 ${HSSP_MATRIX} | awk 'BEGIN{OFS="\t"} {print $1, $2, $5}' > ${EDGE_LIST}
echo "  Edges: $(wc -l < ${EDGE_LIST})"

# --- Step 1: mcl ---
# --abc     : read input in abc format (node1  node2  weight)
# -I 1.4  : inflation parameter from the paper; controls cluster granularity
#            lower = larger clusters, higher = smaller tighter clusters
# -te     : number of threads (use available cores)
# -o      : output binary cluster file
echo "[STEP 1] Running MCL (I=1.4) ..."
mcl \
	${EDGE_LIST} \
	--abc \
	-I  1.4 \
	-te 8 \
	-o  ${MCL_CLUSTERS}

echo "  MCL done: $(date)"

# --- Summary ---
echo ""
echo "[INFO] Number of clusters:"
wc -l < ${MCL_CLUSTERS}

echo "[INFO] First 5 clusters:"
head -5 ${MCL_CLUSTERS}

echo ""
echo "Job finished: $(date)"
echo "=========================================="
