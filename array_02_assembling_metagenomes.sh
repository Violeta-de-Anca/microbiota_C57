#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -N 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_02_assembling_metagenomes.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_02_assembling_metagenomes.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes

mkdir -p output_path

# do this by generation and group separatedly
#for suffix in cecum_samples last_feces; do
#	for F in F0 F1 F2; do
#		for i in control obese; do
#			sample=$input_main_path/$suffix/trimmed_host_removed/${F}_${i}
#			echo $sample
#			sbatch --export=ALL,sample=$sample assembling_metagenomes_array.sh $sample
#		done
#	done
#done

#do it with the merge data by generation
# last_feces cecum_samples
for suffix in cecum_samples; do
	# F0 F1 F2
	for F in F0; do
		generation=$input_main_path/$suffix/trimmed_host_removed/$F
                echo $generation
                sbatch --export=ALL,sample=$generation /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/assembling_metagenomes_array.sh $generation
        done
done

#for suffix in last_feces; do
#F0 F2
#        for F in F1; do
#                generation=$input_main_path/$suffix/trimmed_host_removed/$F
#                echo $generation
#                sbatch --export=ALL,sample=$generation /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/assembling_metagenomes_array.sh $generation
#        done
#done

#do it individually

#for suffix in last_feces cecum_samples; do
#       for F in $(cat $input_main_path/$suffix/trimmed_host_removed/trimmed.files); do
#		echo $F
#                sbatch --export=ALL,sample=$F /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/assembling_metagenomes_array.sh $F
#	done
#done

#for suffix in cecum_samples last_feces; do
#       for F in $(cat $input_main_path/$suffix/trimmed_host_removed/trimmed.files); do
#		b=${F##*/}
#                c=${F%/*}
#                echo $b
#                echo $c
#                mkdir -p $output_path/${b}
#                metawrap assembly -1 $c/${b}/final_pure_reads_1.fastq  -2 $c/${b}/final_pure_reads_2.fastq \
#                -o $output_path/${b} -m 24
#        done
#done


