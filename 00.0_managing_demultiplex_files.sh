#!/bin/bash
#SBATCH -A uppmax2025-2-151
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J plate1.1
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/gzip_plate1.1..err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/gzip_plate1.1.out


# for plate 1.1
cd /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/whole_families_plate1.1

#for file in fastq_2_[0-9]*.fasta; do
#    num1=${file#fastq_2_}
#    num=${num1%.fasta}
#    if [ -d "I${num}" ]; then
#        echo "Moviendo $file a I${num}/"
#        mv "$file" "I${num}/"
#    else
#        mkdir -p I$num
#	mv "$file" "I${num}/"
#    fi
#done

#for file in fastq_1_[0-9]*.fasta; do
    # Extraer la parte numérica quitando el prefijo y el sufijo:
#    num1=${file#fastq_1_}
#    num=${num1%.fasta}
    # Si existe la carpeta con el nombre "I<num>" la movemos:
#    if [ -d "I${num}" ]; then
#        echo "Moviendo $file a I${num}/"
#        mv "$file" "I${num}/"
#    else
#        echo "No existe la carpeta I${num} para $file"
#    fi
#done

folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/whole_families_plate1.1

for i in I*; do
#	a=${i#I}
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
#	rm $folder/$i/*_merged.fastq
	gzip $folder/$i/*.fasta
	gzip $folder/$i/*.fastq
done
