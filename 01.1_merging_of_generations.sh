#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 10-00:00:00
#SBATCH -J mergeLF
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_lastfeces_samples.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_lastfeces_samples.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

# I need to do last_feces when plate 4 finish the host removal
#for suffix in cecum_samples; do
#       for F in $(cat $input_main_path/$suffix/trimmed_host_removed/trimmed.files); do
#                if [[ "$F" == *"F0"* ]]; then
#		    zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
#		    zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq
#		elif [[ "$F" == *"F1"* ]]; then
#		    zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_1.fastq
#                    zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_2.fastq
#		elif [[ "$F" == *"F2"* ]]; then
#		    zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_1.fastq
#                    zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_2.fastq
#		else
#		    echo "No matching pattern found in folder name."
#		fi
#        done
#done

#gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_2.fastq

# for last_feces
for suffix in last_feces; do
       for F in $input_main_path/$suffix/trimmed_host_removed/LF_F*; do
                if [[ "$F" == *"F0"* ]]; then
                    zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_1.fastq
                    zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_2.fastq
                elif [[ "$F" == *"F1"* ]]; then
			rm $F/fastq_1.fasta_val_1.fq
                	rm $F/fastq_1.fasta_val_2.fq
			if ls $F/final_pure_reads_1.fastq.gz > /dev/null 2>&1; then
				zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_1.fastq
                		zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_2.fastq
                	fi
			if ls $F/final_pure_reads_1.fastq > /dev/null 2>&1; then
                                cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_1.fastq
                                cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_2.fastq
			fi
		elif [[ "$F" == *"F2"* ]]; then
			rm $F/fastq_1.fasta_val_1.fq
                        rm $F/fastq_1.fasta_val_2.fq
			cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_1.fastq
			cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_2.fastq
                else
                    echo "No matching pattern found in folder name."
                fi
        done
done

gzip $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_1.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_2.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_1.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_2.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_1.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_2.fastq


