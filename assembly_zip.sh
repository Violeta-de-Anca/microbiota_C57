#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00
#SBATCH -J tidy_up
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/tidy_up_zip.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/tidy_up_zip.out

#this script zips the assembly so it does not take lots of space

input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes

for i in F0 F2 F1; do
	for e in cecum_samples last_feces; do
		for F in $input_fasta_path/${i}_${e}; do
			echo $F
			b=${F##*/}
			c=${F%/*}
			echo $c
			tar -czvf $F/megahit_intermediate_tar.gz $F/megahit
			rm -r $F/megahit
			gzip $F/final_assembly.fasta
		done
	done
done
