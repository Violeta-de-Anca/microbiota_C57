#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00
#SBATCH -J tidy_up
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/tidy_up_zip.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/tidy_up_zip.out


input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics

for F in $(cat $input_fasta_path/samples_multigenerational.txt); do
	echo $F
	b=${F##*/}
	c=${F%/*}
	echo $c
	#rm -r $F/work_files
	#tar -czvf $F/concoct/work_files.tar.gz $F/concoct/work_files
	#tar -czvf $F/metaBAT2/work_files/work_files.tar.gz $F/metaBAT2/work_files
        tar -czvf $F/maxbin2/work_files.tar.gz $F/work_files
	rm -r $F/work_files
	rm -r $F/concoct/work_files
	rm $F/metaBAT2/work_files/assembly.fa
	rm $F/metaBAT2/work_files/assembly.fa.amb
	rm $F/metaBAT2/work_files/assembly.fa.ann
	rm $F/metaBAT2/work_files/assembly.fa.sa
	rm $F/metaBAT2/work_files/assembly.fa.bwt
	rm $F/metaBAT2/work_files/assembly.fa.pac
	rm $F/metaBAT2/work_files/final_pure_reads.bam
done

