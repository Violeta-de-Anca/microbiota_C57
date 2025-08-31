#!/bin/bash
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
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
if file $input_megahit | grep -q 'gzip compressed'; then
	echo "assemly compressed in gz: $input_megahit"
	zcat $input_megahit > $temp_assembly
else
	echo "assembly not compressed: $input_megahit"
	cat $input_megahit > $temp_assembly
fi

mkdir -p $output_folder

#uncompress the fastq files
while IFS= read -r -d '' dir; do
	echo $dir
	sample_id=$(basename ${dir})
	for i in 1 2; do
		fasqgz=$dir/final_pure_reads_${i}.fastq.gz
		fastq=$dir/final_pure_reads_${i}.fastq
		if [[ -f $fasqgz ]]; then
			file=$fasqgz
		elif [[ -f $fastq ]]; then
			file=$fastq
		else
			echo "something went wrong in directory $dir" >&2
			continue
		fi
		dest=$temp_dir/${sample_id}_${i}.fastq
		echo "Decompressing/Copying: $file -> $dest"
		if file --brief --mime-type $file | grep -q '^application/x-gzip$'; then
			echo "compressed fastq file: $file"
			zcat -- $file > $dest
		else
			echo "not compressed fastq file: $file"
			cat -- $file > $dest
		fi
	done
done < <(
    find $folder_fastq/trimmed_host_removed -type d -path "*/$gen*" -print0
)

#metawrap quant_bins -b $bins_folder -o $output_folder -t 8 -a $temp_assembly $temp_dir/*
