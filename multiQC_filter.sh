#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J multiQC
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/multiQC_filtering_step.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/multiQC_filtering_step.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

module load bioinfo-tools FastQC MultiQC

for F in $(cat /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed/fastQC_pure_reads.txt); do
	fastqc $F
done

multiqc -o /proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter -l /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed/trimmed.files
multiqc -o /proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter -l /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_feces/trimmed_host_removed/trimmed.files


