#!/bin/bash
#########!/bin/bash -x
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J kraken
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/kraken2_gens.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/kraken2_gens.out


sample=$1
echo $sample
b=${sample##*/} # file name but with obese/control
f=${b%_obese_last_feces_1.fastq.gz} #this is only the generation
temp=${sample%/*/*} #with this I will get the type of sample with the path
d=${temp##*/} #type of sample: cecum or last feces
e=${temp##*/}

echo $f
echo $b
echo $temp
echo $d
echo $e

module load bioinfo-tools metaWRAP/1.3.2 CheckM

module load biopython/1.76-py3

checkm data setRoot /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16
export CHECKM_DATA_PATH=/proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/CheckM_data/2015_01_16

assembly_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes/${b}_${d}
kraken_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy

# do it by generations
temp_dir=/scratch/${b}_${d}

mkdir -p $temp_dir
echo $temp_dir

temp_fastq1=${temp_dir}/reads_1.fastq
temp_fastq2=${temp_dir}/reads_2.fastq

echo $temp_fastq1
zcat ${sample}_${e}_1.fastq.gz > $temp_fastq1 || { echo "failed decompression";exit 1; }
zcat ${sample}_${e}_2.fastq.gz > $temp_fastq2 || { echo "failed decompression";exit 1; }

echo "Reads in 1: $(grep -c ^@ $temp_fastq1)"
echo "Reads in 2: $(grep -c ^@ $temp_fastq2)"

mkdir -p $kraken_path/${b}_${d}

echo $kraken_path/${b}_${d}

metawrap kraken2 -o $kraken_path/${b}_${d} $assembly_path $temp_fastq1 $temp_fastq2
