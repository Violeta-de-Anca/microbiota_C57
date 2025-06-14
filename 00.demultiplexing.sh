#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J demulti_7
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/demultiplex_plate7.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/demultiplex_plate7.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

# Sample folder
#input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/X204SC23116322-Z01-F005/01.RawData
#barcode_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

#demultiplex the unmultiplex files of the multigeneration part, plate 3

#./demultiplex demux -p $input_folder $barcode_folder/plate_3_barcodes.txt \
#$input_folder/Undetermined/Undetermined_Undetermined_22MG2HLT3_L8_1.fq \
#$input_folder/Undetermined/Undetermined_Undetermined_22MG2HLT3_L8_2.fq

# Sample folder, F010
input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/X204SC23116322_Z01_F010/X204SC23116322-Z01-F010_02/01.RawData
barcode_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

#demultiplex the unmultiplex files of the multigeneration part, plate 7

./demultiplex demux -p $input_folder $barcode_folder/plate_7_barcodes.txt \
$input_folder/Undetermined_Undetermined_22TWL7LT3_L3_1.fq.gz \
$input_folder/Undetermined_Undetermined_22TWL7LT3_L3_2.fq.gz

# Sample folder, F09
#input_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/X204SC23116322_Z01_F010/X204SC23116322-Z01-F010_02/01.RawData/Undetermined
#barcode_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples

#demultiplex the unmultiplex files of the multigeneration part, plate 7

#./demultiplex demux -p $input_folder $barcode_folder/plate_7_barcodes.txt \
#$input_folder/Undetermined_Undetermined_22TWL7LT3_L3_1.fq.gz \
#$input_folder/Undetermined_Undetermined_22TWL7LT3_L3_2.fq.gz

