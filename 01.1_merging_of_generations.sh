#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J merge
#SBATCH --mail-type=BEGIN,FAIL
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_cecum_samples.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_cecum_samples.out
##SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_lastfeces_samples.err
##SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/merging_lastfeces_samples.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
suffix=cecum_samples

# I need to do last_feces when plate 4 finish the host removal
#for suffix in cecum_samples; do
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq.gz
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq.gz
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq.gz
#	rm $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq.gz
#	touch $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq
#	touch $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq
#	touch $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq
#	touch $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq
#       for F in $(cat $input_main_path/$suffix/trimmed_host_removed/trimmed.files); do
#                if [[ "$F" == *"F0"* && "$F" == *"A"* ]]; then
#			echo $F
#			r1="$F/final_pure_reads_1.fastq.gz"
#			r2="$F/final_pure_reads_2.fastq.gz"
#			n1=$(zcat $r1 | wc -l)
#			n2=$(zcat $r2 | wc -l)
#			if [[ $n1 -eq $n2 ]]; then
#				zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq
#				zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq
#			else
#				echo "⚠️ Skipping $F because of mismatch"
#			fi
#		elif [[ "$F" == *"F0"* && "$F" == *"B"* ]]; then
#			echo $F
#			r1="$F/final_pure_reads_1.fastq.gz"
#                        r2="$F/final_pure_reads_2.fastq.gz"
#                        n1=$(zcat $r1 | wc -l)
#                        n2=$(zcat $r2 | wc -l)
#                        if [[ $n1 -eq $n2 ]]; then
#				zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq
#				zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq
#			else
#                                echo "⚠️ Skipping $F because of mismatch"
#                        fi
#		elif [[ "$F" == *"F1"* && "$F" == *"A"* ]]; then
#			echo $F
#			r1="$F/final_pure_reads_1.fastq.gz"
#			r2="$F/final_pure_reads_2.fastq.gz"
#			n1=$(zcat $r1 | wc -l)
#			n2=$(zcat $r2 | wc -l)
#			if [[ $n1 -eq $n2 ]]; then
#				zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_cecum_samples_1.fastq
#				zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_cecum_samples_2.fastq
#			else
#				echo "⚠️ Skipping $F because of mismatch"
#			fi
#		elif [[ "$F" == *"F1"* && "$F" == *"B"* ]]; then
#			echo $F
#			r1="$F/final_pure_reads_1.fastq.gz"
#                        r2="$F/final_pure_reads_2.fastq.gz"
#                        n1=$(zcat $r1 | wc -l)
#                        n2=$(zcat $r2 | wc -l)
#                        if [[ $n1 -eq $n2 ]]; then
#                                zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_control_cecum_samples_1.fastq
#                                zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_control_cecum_samples_2.fastq
#                        else
#                                echo "⚠️ Skipping $F because of mismatch"
#                        fi
#		elif [[ "$F" == *"F2"* && "$F" == *"A"* ]]; then
#		    	echo $F
			#zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_obese_cecum_samples_1.fastq
                    	#zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_obese_cecum_samples_2.fastq
#		elif [[ "$F" == *"F2"* && "$F" == *"B"* ]]; then
#                    	echo $F
			#zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_control_cecum_samples_1.fastq
                    	#zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F2_control_cecum_samples_2.fastq
#		else
#		    	echo "No matching pattern found in folder name."
#		fi
#        done
#done

#gzip $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_obese_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_obese_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_control_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_control_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_obese_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_obese_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_control_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_control_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_cecum_samples_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_cecum_samples_2.fastq

#do the merging of the F1 cecum samples into one single file to do the assembling, as in here all individuals have equal number of read in 1 and 2
zcat $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_1.fastq.gz > $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
zcat $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
zcat $input_main_path/$suffix/trimmed_host_removed/F0_obese_cecum_samples_2.fastq.gz > $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq
zcat $input_main_path/$suffix/trimmed_host_removed/F0_control_cecum_samples_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_1.fastq
gzip $input_main_path/$suffix/trimmed_host_removed/F0_cecum_samples_2.fastq

# for last_feces
#for suffix in last_feces; do
#       for F in $input_main_path/$suffix/trimmed_host_removed/LF_F*; do
#                if [[ "$F" == *"F0"* && "$F" == *"A"* ]]; then
#                    zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_obese_last_feces_1.fastq
#                    zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_obese_last_feces_2.fastq
#		elif [[ "$F" == *"F0"* && "$F" == *"B"* ]]; then
#			zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_control_last_feces_1.fastq
#			zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F0_control_last_feces_2.fastq
#                elif [[ "$F" == *"F1"* && "$F" == *"A"* ]]; then
#			rm $F/fastq_1.fasta_val_1.fq
#                	rm $F/fastq_2.fasta_val_2.fq
#			if ls $F/final_pure_reads_1.fastq.gz > /dev/null 2>&1; then
#				zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_1.fastq
#                		zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_2.fastq
#                	fi
#			if ls $F/final_pure_reads_1.fastq > /dev/null 2>&1; then
#                                cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_1.fastq
#                                cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_2.fastq
#			fi
#		elif [[ "$F" == *"F1"* && "$F" == *"B"* ]]; then
#                        rm $F/fastq_1.fasta_val_1.fq
#                        rm $F/fastq_2.fasta_val_2.fq
#                        if ls $F/final_pure_reads_1.fastq.gz > /dev/null 2>&1; then
#                                zcat $F/final_pure_reads_1.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_1.fastq
#                                zcat $F/final_pure_reads_2.fastq.gz >> $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_2.fastq
#                        fi
#                        if ls $F/final_pure_reads_1.fastq > /dev/null 2>&1; then
#                                cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_1.fastq
#                                cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_2.fastq
#                        fi
#		elif [[ "$F" == *"F2"* && "$F" == *"A"* ]]; then
#			rm $F/fastq_1.fasta_val_1.fq
#                        rm $F/fastq_2.fasta_val_2.fq
#			cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_obese_last_feces_1.fastq
#			cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_obese_last_feces_2.fastq
#		elif [[ "$F" == *"F2"* && "$F" == *"B"* ]]; then
#                        rm $F/fastq_1.fasta_val_1.fq
#                        rm $F/fastq_2.fasta_val_2.fq
#                        cat $F/final_pure_reads_1.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_control_last_feces_1.fastq
#                        cat $F/final_pure_reads_2.fastq >> $input_main_path/$suffix/trimmed_host_removed/F2_control_last_feces_2.fastq
#                else
#                    echo "No matching pattern found in folder name."
#                fi
#        done
#done

#gzip $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_control_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_obese_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_obese_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_control_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F2_control_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_control_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_control_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_obese_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F0_obese_last_feces_2.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_1.fastq
#gzip $input_main_path/$suffix/trimmed_host_removed/F1_obese_last_feces_2.fastq
