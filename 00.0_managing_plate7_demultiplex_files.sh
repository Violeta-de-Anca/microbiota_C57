#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 10GB
#SBATCH -t 1:00:00
#SBATCH -J raw_LF
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/managing_plate7.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/managing_plate7.out


# for plate 7, last feces from the transgenerational part
cd /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/transgenerational_samples

#in here we are going to move to a definitive folder each R1 and R2 of each individual to a folder with the definitive identifier
#for file in Undetermined_Undetermined_22TWL7LT3_L3_2_[0-9]*.fq.gz; do
#    num1=${file#Undetermined_Undetermined_22TWL7LT3_L3_2_}
#    num=${num1%.fq.gz}
#	if [ -d "L_${num}" ]; then
#		echo "Moviendo $file a L_${num}/"
#		mv "$file" "L_${num}/"
#	else
#		mkdir -p L_$num
#		mv "$file" "L_${num}/"
#	fi
#	if [ -d "LP_${num}" ]; then
#		echo "Moviendo $file a LP_${num}/"
#		mv "LP_${num}/" "L_${num}/"
#		mv "$file" "L_${num}/"
#	fi
#done

#for file in Undetermined_Undetermined_22TWL7LT3_L3_1_*.fq.gz; do
#	echo $file
#	num1=${file#Undetermined_Undetermined_22TWL7LT3_L3_1_}
#	num=${num1%.fq.gz}
#	if [ -d "LP_${num}" ]; then
#        	echo "Moviendo $file a LP_${num}/"
#        	mv "LP_${num}/" "L_${num}/"
#        	mv "$file" "L_${num}/"
#	fi
   # Extraer la parte numérica quitando el prefijo y el sufijo:
#	echo $file
#	num1=${file#Undetermined_Undetermined_22TWL7LT3_L3_1_}
#	num=${num1%.fq.gz}
#	if [ -d L_${num} ]; then
#		echo "Moviendo $file a L_${num}/"
#		mv $file L_${num}/
#	else
#		echo "No existe la carpeta L_${num} para $file"
#	fi
#done

folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/transgenerational_samples

#for i in L_*; do
for i in L_12231_L; do
	a=${i#L_}
	echo $i
	echo $a
	echo $folder
#New names for all the R1 and R2 files with the identification of each individual
	if [ -f $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz ]; then
		zcat $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz > $folder/$i/L_${a}_1.fastq
		echo $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz
	fi
	if [ -f $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_2_${a}.fq.gz ]; then
		zcat $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_2_${a}.fq.gz > $folder/$i/L_${a}_2.fastq
	fi

	gzip $folder/$i/*fasta
	gzip $folder/$i/*fastq
done

