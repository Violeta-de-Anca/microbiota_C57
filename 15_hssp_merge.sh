#!/bin/bash
##SBATCH -A uppmax2025-2-302
#SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 16GB
#SBATCH -t 1:00:00
#SBATCH -J hssp_merge
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/hssp_merge.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/hssp_merge.out

# =============================================================================
# 15_merge_hssp.sh
#
# Merge all 600 per-chunk HSSP result files into one final matrix.
# Run this AFTER all 15_submit_hssp_array.sh jobs have completed.
#
# No deduplication is needed at this stage because each query protein
# appeared in exactly one original chunk → each (qseqid, sseqid) pair
# is present in at most one chunk output file.
#
# Output columns:
#   qseqid  sseqid  qbin  sbin  hssp_distance
# =============================================================================

set -euo pipefail

BASE=/proj/naiss2024-23-57/C57_female_lineage_microbiota
HSSP_OUT=${BASE}/functional_annotation/query_chunks
FINAL=${BASE}/functional_annotation/hssp_matrix.txt

echo "=========================================="
echo "Merge started : $(date)"
echo "Input dir     : ${HSSP_OUT}"
echo "Output file   : ${FINAL}"
echo "=========================================="

# --- Sanity check: verify all 600 chunk outputs exist ---
echo "[INFO] Checking all 600 chunk outputs ..."
missing=0
for i in $(seq 0 599); do
	F=$(printf "%04d" $i)
	f="${HSSP_OUT}/hssp_chunk_${F}.txt"
	if [ ! -f "$f" ]; then
		echo "  MISSING: $f"
		((missing++))
	fi
done

if [ "${missing}" -gt 0 ]; then
	echo "[ERROR] ${missing} chunk file(s) missing. Re-run failed array jobs before merging."
	exit 1
fi
echo "[INFO] All 600 chunk files present."

# --- Merge: write header once, then append all data rows ---
echo "[INFO] Writing header ..."
echo -e "qseqid\tsseqid\tqbin\tsbin\thssp_distance" > "${FINAL}"

echo "[INFO] Concatenating chunk files (skipping per-file headers) ..."
for i in $(seq 0 599); do
	F=$(printf "%04d" $i)
	echo $F
	# tail -n +2 skips the header line of each chunk file
	tail -n +2 "${HSSP_OUT}/hssp_chunk_${F}.txt"
done >> "${FINAL}"

# --- Summary ---
echo ""
echo "[INFO] Merge complete."
echo "  Total pairs in HSSP matrix:"
tail -n +2 "${FINAL}" | wc -l

echo ""
echo "[INFO] First 5 data lines:"
head -6 "${FINAL}"

echo ""
echo "Merge finished: $(date)"
echo "=========================================="
