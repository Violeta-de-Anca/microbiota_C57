#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/sanity_F0_cecum_feces.out
module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

#I want to see which of the files have the uneven number of samples
temp_fastq1=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_cecum_samples_1.fastq.gz | grep -c ^@)
temp_fastq2=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_cecum_samples_2.fastq.gz | grep -c ^@)
temp_fastqc1=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_control_cecum_samples_1.fastq.gz | grep -c ^@)
temp_fastqc2=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_control_cecum_samples_2.fastq.gz | grep -c ^@)
temp_fastqo1=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_obese_cecum_samples_1.fastq.gz | grep -c ^@)
temp_fastqo2=$(zcat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed/F0_obese_cecum_samples_2.fastq.gz | grep -c ^@)

echo "Reads in total R1: $temp_fastq1"
echo "Reads in total R2: $temp_fastq2"
echo "Reads in control R1: $temp_fastqc1"
echo "Reads in control R2: $temp_fastqc2"
echo "Reads in obese R1: $temp_fastqo1"
echo "Reads in obese R2: $temp_fastqo2"
