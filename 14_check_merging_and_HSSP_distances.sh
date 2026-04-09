#!/bin/bash
#SBATCH -A uppmax2025-2-536
##SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 10GB
##SBATCH --mem 100GB
#SBATCH -t 2:00:00
##SBATCH -t 1-00:00:00
#SBATCH -J QC_blast
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/merge_psi_blast_HSSP_matrix.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/merge_psi_blast_HSSP_matrix.out

module load  BLAST+/2.17.0-GCC-13.3.0
in_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation/query_chunks
out_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation

#give access to all the results
chmod -R 777 ${in_folder}

#first let's double check that all the chunks finished correctly
for i in $(seq 0 599); do
	F=$(printf "%04d" $i)
	chunk=${in_folder}/chunk_${F}.faa
	result=${in_folder}/psi_results_chunck_${F}_multigen.txt
	echo $chunk
	echo $result
	# skip if result file doesn't exist yet
	if [ ! -f $result ]; then
		echo "MISSING: chunk $F"
		continue
	fi

	expected=$(grep -c ">" $chunk)
	actual=$(cut -f1 $result | sort -u | wc -l)

	if [ "$expected" -eq "$actual" ]; then
		echo "OK: chunk $F ($actual/$expected queries)"
	else
		echo "INCOMPLETE: chunk $F ($actual/$expected queries)"
	fi
done

#before this there is this print: Search has CONVERGED!
#which we need to delete from all chunks
for i in $(seq 0 599); do
	F=$(printf "%04d" $i)
	result=${in_folder}/psi_results_chunck_${F}_multigen.txt
	if [ -f $result ]; then
	awk 'NF==13' $result
	else
		echo "WARNING, something is wrong with $result" >&2
	fi
done > ${out_folder}/psi_results_multigen.txt

#then do the HSSP distance matrix
