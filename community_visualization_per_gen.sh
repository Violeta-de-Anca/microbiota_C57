#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p node
##SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10-00:00:00
##SBATCH --ntasks-per-core 8
#SBATCH -J community
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/community_visualization_per_gen.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/community_visualization_per_gen.out

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

#for doing it per generations/sample type
a=${1##*/} #without path (F0/1/2_last_feces)
b=${1%/*} #the path (/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/)
c=${a%_last_feces}

out_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/community_visualization
assembly_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${a}/final_assembly.fasta
refined_bins=/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/${a}/refined_libraries_megahit/metawrap_70_10_bins
reads_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed

#make the directory
mkdir -p $out_folder
mkdir -p $out_folder/${a}

temp_dir=${TMPDIR:-${SNIC_TMP:-/scratch/$SLURM_JOB_ID}}
mkdir -p $temp_dir

#uncompress the fastq files
while read -r gzfile; do
        a=${gzfile##*/}
        base=${a%.gz}
        dest=$temp_dir/$base
        echo "Decompressing $gzfile -> $dest"
        zcat $gzfile > $dest
done < $reads_folder/${c}_fastq


metawrap blobology -a $assembly_folder -o $out_folder/${a} --bins $refined_bins $temp_dir/*fastq*
