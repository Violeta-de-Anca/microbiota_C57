#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00
#SBATCH -J checkm
#SBATCH --mail-type=BEGIN
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/checkm_megahit_refined_bins.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/checkm_megahit_refined_bins.out


module load bioinfo-tools CheckM

input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics
output=/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/QC_filter/checkm_refined_bins_megahit

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

for F in $(cat $input_path/samples_multigenerational.txt); do
	a=${F##*/}
	mkdir -p $output/$a
	checkm lineage_wf -x fa --tmpdir $temp_dir --tab_table $F/refined_libraries_megahit/binsABC $output/$a
done
