#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 1GB
#SBATCH -t 10:00
#SBATCH -J jobarray
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_01_plate7and1.1_removal.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/array_01_plate7and1.1_removal.out

# SLURM_ARRAY_TASK_ID tells the script which iteration to run
echo $SLURM_ARRAY_TASK_ID

#pelle
module load metaWRAP/1.4-20230728-foss-2024a-Python-2.7.18

#For plate 7, which is the transgenerational part of the last feces samples
input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/transgenerational_samples

for F in $input_path/L_*; do
	echo $F
	sbatch --export=ALL,a=$F plate_7_01_removal_of_host.sh $F
done

#plate 4, which is the last feces F1 and F2
#input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/X204SC23116322_Z01_F010/X204SC23116322-Z01-F010_01/01.RawData

#for F in $(cat $input_path/plate_4_sample.list.txt); do
#        echo $F
#        sbatch --export=ALL,a=$F plate_4_01_removal_of_host.sh $F
#done

#plate 6, which is the transgenerational cecum feces
#input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/raw_files_F0_till_f2/transgenerational_samples

#for F in $input_path/*; do
#        echo $F
#        sbatch --export=ALL,a=$F plate_6_01_removal_of_host.sh $F
#done

#plate 1.1, which is the whole families of cecum feces
#input_path=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/whole_families_plate1.1

#for F in $input_path/*; do
#        echo $F
#        sbatch --export=ALL,a=$F plate_1.1_01_removal_of_host.sh $F
#done

