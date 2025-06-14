#!/bin/bash -x
#SBATCH -A uppmax2025-2-222
#SBATCH -p node
##SBATCH -p core
#SBATCH -n 2
##SBATCH --cpus-per-task 16
#SBATCH -t 10-00:00:00
#SBATCH -J metaspa_%j
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/assembly_gens_metaspades_%j.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/assembly_gens_metaspades_%j.out


sample=$1
echo $sample
b=${sample##*/} # file name
temp=${sample%/*/*} #with this I will get the type of sample with the path
d=${temp##*/} #type of sample: cecum or last feces
e=${temp##*/}

module load bioinfo-tools metaWRAP/1.3.2

output_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/assembled_metagenomes

# do it by generations
temp_dir=/scratch/${b}_${d}

mkdir -p $temp_dir
echo $temp_dir

temp_fastq1=${temp_dir}/reads_1.fastq
temp_fastq2=${temp_dir}/reads_2.fastq

echo $temp_fastq1
zcat ${1}_${e}_1.fastq.gz > $temp_fastq1 || { echo "failed decompression";exit 1; }
zcat ${1}_${e}_2.fastq.gz > $temp_fastq2 || { echo "failed decompression";exit 1; }

echo "Reads in 1: $(grep -c ^@ $temp_fastq1)"
echo "Reads in 2: $(grep -c ^@ $temp_fastq2)"

mkdir -p $output_path/${b}_${d}

echo $output_path/${b}_${d}
metawrap assembly --metaspades -1 ${temp_dir}/reads_1.fastq  -2 ${temp_dir}/reads_2.fastq  -o ${output_path}/${b}_${d} -m 128
#metawrap assembly --megahit -1 ${temp_dir}/reads_1.fastq  -2 ${temp_dir}/reads_2.fastq  -o ${output_path}/${b}_${d} -m 256 -t 16
rm -r $temp_dir

############################################################################################################

#do it individually

#metawrap assembly --metaspades -1 $c/${b}/final_pure_reads_1.fastq  -2 $c/${b}/final_pure_reads_2.fastq \
#                -o $output_path/${b} -m 64
