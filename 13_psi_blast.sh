#!/bin/bash
##SBATCH -A uppmax2025-2-536
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 100GB
#SBATCH -t 1-00:00:00
#SBATCH -J arrblast
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/blast_functional_matrix.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/blast_functional_matrix.out

#export the particular temporary directory for checkM
temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir
export TMPDIR="$temp_dir"
export TEMP="$temp_dir"
export TMP="$temp_dir"

module load  BLAST+/2.17.0-GCC-13.3.0

input_folder_query=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks
input_db=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

if [ ! -f ${input_folder}/all_proteins_multigen.faa ] || [ ! -f ${input_folder}/all_proteins_multigen_db.pdb ]; then
	for i in F1 F0 F2; do
		for f in cecum_samples last_feces; do
			bin=${input_folder}/${i}_${f}/bin_translated_genes
			for a in ${bin}/*_origin_track_only_annotated.faa; do
				cat ${a} >> ${input_folder}/all_proteins_multigen.faa
			done
		done
	done
	makeblastdb -in ${input_folder}/all_proteins_multigen.faa -dbtype prot -out ${input_folder}/all_proteins_multigen_db
else
	echo "database exist, skipping until the psi-blast"
fi

query=${input_folder_query}/chunk_${1}.faa


psiblast -db ${input_db}/all_proteins_multigen_db \
-query ${query} \
-out ${input_folder_query}/psi_results_chunck_${1}_multigen.txt \
-evalue 1e-3 \
-inclusion_ethresh 1e-10 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident" \
-num_iterations 3
# same as in: https://doi.org/10.1371/journal.pcbi.1004472 -num_iterations
