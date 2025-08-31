#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J quant
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_group.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_group.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics)
gen=${a%%_*}
group_and_type=${a#*_}
group=${group_and_type%%_*}
sample=${group_and_type#*_}

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly=$temp_dir/assembly.fasta

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/${a}/per_group
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/metawrap_70_10_bins
folder_fastq=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/$gen

echo $input_megahit
echo $a
echo $group
echo $sample
echo $output_folder
echo $bins_folder

# Decompress the assembly into the temporary file
if file $input_megahit | grep -q 'gzip compressed'; then
        echo "assemly compressed in gz: $$input_megahit"
        zcat $input_megahit > $temp_assembly
else
        echo "assembly not compressed: $$input_megahit"
        cat $input_megahit > $temp_assembly
fi

mkdir -p $output_folder

#uncompress the fastq files
while IFS= read -r -d '' dir; do
	echo $dir
	sample_id=$(basename ${dir})
	for i in 1 2; do
		zcat ${dir}_${i}.fastq.gz > $temp_dir/${sample_id}_${i}.fastq
	done
done < <(find $folder_fastq/trimmed_host_removed -type d -path "*/${gen}_${group_and_type}*" -print0

#metawrap quant_bins -b $bins_folder -o $output_folder -t 8 -a $temp_assembly $temp_dir/*fastq*

