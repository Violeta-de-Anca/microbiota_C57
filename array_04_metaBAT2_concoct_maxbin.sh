#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_04_metaBAT2.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_04_metaBAT2.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

# do it by generation and type of sample

for suffix in last_feces; do
#F0 F2
        for F in F0 F2 F1; do
                generation=$input_fasta_path/$suffix/trimmed_host_removed/$F
                echo $generation
                sbatch --export=ALL,sample=$generation bin_concoct_MaxBin_metaBAT_per_gen.sh $generation
        done
done

# do it individually

#for suffix in last_feces cecum_samples; do
#       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
#                echo $F
#                echo $b
#                b=${F##*/}
#                c=${F%/*}
#                echo $c
#                sbatch --export=ALL,a=$F \
#                /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/bin_metaBAT.sh $F
#        done
#done

#for suffix in last_feces cecum_samples; do
#       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
#                echo $F
#                echo $b
#                b=${F##*/}
#                c=${F%/*}
#                echo $c
#                sbatch --export=ALL,a=$F /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/bin_MaxBin.sh $F
                #/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/bin_concoct.sh $F
#        done
#done
