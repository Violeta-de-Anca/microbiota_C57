#!/bin/bash
#SBATCH -A uppmax2025-2-536
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 1:00:00
#SBATCH -J fun_table
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/functional_tables_clean_merge.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/functional_tables_clean_merge.out

#export the particular temporary directory for checkM
temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir
export TMPDIR="$temp_dir"
export TEMP="$temp_dir"
export TMP="$temp_dir"

#for i in ${1}/bin_funct_annotations/*.strict.gff; do
#	a=${i##*/}
#	b=${a/.strict.gff/}
#	c=${i%/*}
#	grep 'COG:' $i > ${c}/${b}_strict_only_COG_ID.gff
#	grep -o "COG:[A-Z0-9]*" ${c}/${b}_strict_only_COG_ID.gff | cut -f 2 -d ":" | sort | uniq > ${c}/${b}_COG_ID.txt
#done

for i in ${1}/bin_translated_genes/*.strict.faa; do
	a=${i##*/}
	b=${a/.strict.faa/}
	c=${i%/*}
	sample=${i%/*/*}
	sample=${sample##*/}
	awk 'BEGIN{p=1} /^>/{p=($0 !~ /hypothetical|putative/)} p' $i > ${c}/${b}_no_hypothetical.faa
	rm ${c}/${b}_origin_track_only_annotated.faa
	sed "s/^>/>${sample}_${b}_/" ${c}/${b}_no_hypothetical.faa >> ${c}/${b}_origin_track_only_annotated.faa
	rm ${c}/${b}_no_hypothetical.faa
done
