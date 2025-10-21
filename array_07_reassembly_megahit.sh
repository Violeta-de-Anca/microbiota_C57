#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_07_megahit.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_07_megahit.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

#for gen reassembly
input_main_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples
#last_feces cecum_samples
for suffix in last_feces cecum_samples; do
        for F in F0 F1 F2 F3 F4 F5; do
                generation=$input_main_path/$suffix/trimmed_host_removed/$F
                echo $generation
                sbatch --export=ALL,sample=$generation reassembly_per_gen.sh $generation
        done
done

#for individual reassembly
#input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

#for suffix in last_feces cecum_samples; do
#       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
#                echo $F
#                b=${F##*/}
#                c=${F%/*}
#		echo $b
#                echo $c
#                sbatch --export=ALL,a=$F /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/reassembly_megahits.sh $F
#        done
#done

