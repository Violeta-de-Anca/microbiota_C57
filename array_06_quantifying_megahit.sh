#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_06_quantifying_megahit.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_06_quantifying_megahit.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools metaWRAP/1.3.2

#for doing per generations/per sample type, but to get individual information
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics
#F0 F2 F1 F3 F4 F5
#last_feces cecum_samples
for suffix in cecum_samples; do
        for F in F3 F4 F5; do
                generation=$input_path/${F}_${suffix}
                echo $generation
                sbatch --export=ALL,sample=$generation quantifying_per_gen.sh $generation
        done
done

#for doing per generations/per sample type, but to get per group information
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics
#last_feces cecum_samples
for suffix in cecum_samples; do
#F0 F2 F1 F3 F4 F5
	for F in F3 F4 F5; do
#control obese
		for e in control; do
			generation=$input_path/${F}_${e}_${suffix}
			echo $generation
			sbatch --export=ALL,sample=$generation quantifying_per_group.sh $generation
		done
	done
done

# for doing it per individual
#input_fasta_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

#for suffix in last_feces cecum_samples; do
#       for F in $(cat $input_fasta_path/$suffix/trimmed_host_removed/trimmed.files); do
#                echo $F
#                echo $b
#                b=${F##*/}
#                c=${F%/*}
#                echo $c
#                sbatch --export=ALL,a=$F /proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/quantifying_megahit.sh $F
#        done
#done
