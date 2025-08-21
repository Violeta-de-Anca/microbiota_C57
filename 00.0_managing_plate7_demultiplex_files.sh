#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J plate7
#SBATCH --error=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/managing_plate7.err
#SBATCH --output=/proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/managing_plate7.out


# for plate 7
cd /proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/transgenerational_samples

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

for i in L_*; do
	a=${i#L_}
	echo $i
	echo $a
	echo $folder
#first do forward file (1)
	if ls $folder/$i/LP_*_1.*> /dev/null 2>&1; then
		echo $folder/$i/LP_${a}_*-1A_22TWL7LT3_L3_1.fq.gz
		zcat $folder/$i/LP_${a}_*-1A_22TWL7LT3_L3_1.fq.gz > $folder/$i/L_${a}_1.merged.fastq
		echo $folder/$i/${a}_1.merged.fastq
	fi
#also for the ones that are called different (LP_)
	if ls $folder/$i/L_*_1.*> /dev/null 2>&1; then
		echo $folder/$i/L_${a}_*-1A_22TWL7LT3_L3_1.fq.gz
                zcat $folder/$i/L_${a}_*-1A_22TWL7LT3_L3_1.fq.gz > $folder/$i/L_${a}_1.merged.fastq
                echo $folder/$i/L_${a}_1.merged.fastq
        fi
#now do it for the reverse file (2)
	if ls $folder/$i/LP_${a}_*_2.* 1> /dev/null 2>&1; then
		zcat $folder/$i/LP_${a}_*-1A_22TWL7LT3_L3_2.fq.gz > $folder/$i/L_${a}_2.merged.fastq
		echo $folder/$i/LP_${a}_*-1A_22TWL7LT3_L3_2.fq.gz
	fi
#also for the ones that are called different (LP_)
	if ls $folder/$i/L_${a}_*_2.* 1> /dev/null 2>&1; then
        	zcat $folder/$i/L_${a}_*-1A_22TWL7LT3_L3_2.fq.gz > $folder/$i/L_${a}_2.merged.fastq
        	echo $folder/$i/L_${a}_*-1A_22TWL7LT3_L3_2.fq.gz
        fi

	if [ -f $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz ]; then
		zcat $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz >> $folder/$i/L_${a}_1.merged.fastq
		echo $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_1_${a}.fq.gz
	fi
	if [ -f $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_2_${a}.fq.gz ]; then
		zcat $folder/$i/Undetermined_Undetermined_22TWL7LT3_L3_2_${a}.fq.gz >> $folder/$i/L_${a}_2.merged.fastq
	fi

	gzip $folder/$i/*fasta
	gzip $folder/$i/*fastq
done

