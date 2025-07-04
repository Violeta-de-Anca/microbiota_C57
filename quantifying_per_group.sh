#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J quant
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_group.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_group.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics)
c=${a%_last_feces}

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly="$temp_dir/assembly.fasta"

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}/final_assembly.fasta
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/${a}/per_group
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/metawrap_70_10_bins
folder_fastq=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed

echo $input_megahit
echo $a
echo $input_megahit
echo $output_folder
echo $bins_folder

# Decompress the assembly into the temporary file
cat "$input_megahit" > "$temp_assembly"
#zcat $folder_fastq/${c}_last_feces_1.fastq.gz > "$temp_fastq1"
#zcat $folder_fastq/${c}_last_feces_2.fastq.gz > "$temp_fastq2"

mkdir -p $output_folder

#uncompress the fastq files
while read -r gzfile; do
        a=${gzfile##*/}
        base=${a%.gz}
        dest=$temp_dir/$base
        echo "Decompressing $gzfile -> $dest"
        zcat $gzfile > $dest
done < $folder_fastq/${c}_fastq

metawrap quant_bins -b $bins_folder -o $output_folder -t 8 -a $temp_assembly $temp_dir/*fastq*

