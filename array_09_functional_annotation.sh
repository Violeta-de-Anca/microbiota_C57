#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_09_functional_annotation.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_09_functional_annotation.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

#for doing per generations/per sample type
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes
for suffix in cecum_samples last_feces; do
# F0 F2 F1 F3 F4 F5
        for F in F0 F1 F2 F3 F4 F5; do
                generation=$input_path/${F}_${suffix}
                echo $generation
                sbatch --export=ALL,sample=$generation functional_annotation_per_gen.sh $generation
        done
done

