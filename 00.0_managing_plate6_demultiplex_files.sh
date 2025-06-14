#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J plate6
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/gzip_plate6.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/gzip_plate6.out




# for plate 6
cd /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/raw_files_F0_till_f2/transgenerational_samples

#for file in fastq_2_[0-9]*_M.fasta; do
#    num1=${file#fastq_2_}
#    num=${num1%_M.fasta}
#    if [ -d "M_${num}" ]; then
#        echo "Moviendo $file a M_${num}/"
#        mv "$file" "M_${num}/"
#    else
#        mkdir -p M_$num
#	mv "$file" "M_${num}/"
#    fi
#done

#for file in fastq_1_[0-9]*_M.fasta; do
   # Extraer la parte numérica quitando el prefijo y el sufijo:
#    num1=${file#fastq_1_}
#    num=${num1%_M.fasta}
#    if [ -d "M_${num}" ]; then
#        echo "Moviendo $file a M_${num}/"
#        mv "$file" "M_${num}/"
#    else
#        echo "No existe la carpeta M_${num} para $file"
#    fi
#done

folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/raw_files_F0_till_f2/transgenerational_samples

for i in M_*; do
#	a=${i#M_}
#	echo $i
#	echo $a
#	echo $folder
#	if ls $folder/$i/${i}_*_1.fq.gz 1> /dev/null 2>&1; then
#		zcat $folder/$i/${i}_*_1.fq.gz > $folder/$i/${a}_1.merged.fastq
#		echo $folder/$i/${a}_1.merged.fastq
#	fi
#
#	if ls $folder/$i/${i}_*_2.fq.gz 1> /dev/null 2>&1; then
#		zcat $folder/$i/${i}_*_2.fq.gz >> $folder/$i/${a}_2.merged.fastq
#		echo $folder/$i/${i}_*-1A_22TWKHLT3_L8_2.fq.gz
#	fi
#
#	if [ -f $folder/$i/fastq_1_${a}.fasta ]; then
#		cat $folder/$i/fastq_1_${a}.fasta >> $folder/$i/${a}_1.merged.fastq
#		echo $folder/$i/fastq_1_${a}.fasta
#	fi
#
#	if [ -f $folder/$i/fastq_2_${a}.fasta ]; then
#		cat $folder/$i/fastq_2_${a}.fasta >> $folder/$i/${a}_2.merged.fastq
#	fi
	gzip $folder/$i/*fasta
	gzip $folder/$i/*fastq
done
