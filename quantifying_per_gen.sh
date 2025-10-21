#!/bin/bash
#SBATCH -A uppmax2025-2-150
##SBATCH -p node
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J quant
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/quantifying/quantifying_megahit_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

a=${1##*/} #without path, generation and type of sample
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics)
gen=${a%%_*}
sample=${a#*_}

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p "$temp_dir"

# Create a temporary file on the scratch storage for the decompressed assembly
temp_assembly=$temp_dir/assembly.fasta

input_megahit=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/${a}
bins_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/metawrap_70_10_bins
folder_fastq=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/$sample

echo $input_megahit
echo $a
echo $output_folder
echo $bins_folder

# Decompress the assembly into the temporary file if it compressed
if file $input_megahit/final_assembly.fasta.gz | grep -q 'gzip compressed'; then
	echo "assemly compressed in gz: $input_megahit"
	zcat $input_megahit/final_assembly.fasta.gz > $temp_assembly
else
	echo "assembly not compressed: $input_megahit"
	cat $input_megahit/final_assembly.fasta > $temp_assembly
fi

mkdir -p $output_folder

#uncompress the fastq files - multigen codes
#while IFS= read -r -d '' file; do
#	echo $file
#	base=$(basename $file)
#	sample_id=$(basename $(dirname $file))
#	if [[ $base =~ ^final_pure_reads_([12])\.fastq\.gz$ ]]; then
#		i=${BASH_REMATCH[1]}
#		dest=$temp_dir/${sample_id}_${i}.fastq
#		echo "Decompressing: $file -> $dest"
#		gzip -cd -- $file > $dest
#	else
#		echo "Skipping unexpected file name: $file" >&2
#	fi
#done < <(
#    find $folder_fastq/trimmed_host_removed -type f -path "*/$gen*" -name 'final_pure_reads_[12].fastq.gz' -print0
#)

#uncompress the fastq files - transgen codes
awk -v gen="$gen" -F'\t' '  $0 !~ /^[[:space:]]*($|#)/ && $2 == gen { printf "%s\0", $1 }' $folder_fastq/trimmed_host_removed/transgenerational_relational_table.txt | \
while IFS= read -r -d '' sample_dir; do
	while IFS= read -r -d '' file; do
		base=$(basename "$file")
		sample_id=$(basename "$sample_dir")
		if [[ $base =~ ^final_pure_reads_([12])\.fastq\.gz$ ]]; then
			i=${BASH_REMATCH[1]}
			print $i
			dest=$temp_dir/${sample_id}_${i}.fastq
			print $dest
			gzip -cd -- $file > $dest
		else
			echo "Skipping unexpected file name: $file" >&2
		fi
	done < <(find "$sample_dir" -type f -name 'final_pure_reads_[12].fastq.gz' -print0)
done


#metawrap quant_bins -b $bins_folder -o $output_folder -t 8 -a $temp_assembly $temp_dir/*
