#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10-00:00:00
#SBATCH -J rezip
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_rezip_fastq.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/01_rezip_fastq.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2

#input files:
folder_fastq=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/raw_files_F0_till_f2/X204SC23116322-Z01-F001/01.RawData
output_folder=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/cecum_samples/trimmed_host_removed

#let's start by doing the QC of the plate 1, 2 (cecum feces, generations F0 till F2)
#metawrap read_qc \
#-1 $folder_fastq/F0_103_A_female_M/F0_103_A_female_M_EKDL230054078-1A_22H37MLT3_L4_1.fq \
#-2 $folder_fastq/F0_103_A_female_M/F0_103_A_female_M_EKDL230054078-1A_22H37MLT3_L4_2.fq \
#-o $output_folder -x mm39 \
#this works, so now let's do all the individuals

#So I Need to unzip them all, otherwise it does not work
#for F in $(cat $folder_fastq/sample.list.txt); do
#	b=${F##*/}
#	c=${F%/*}
#	echo $b
#	echo $c
#	gzip -d $c/${b}/${b}*_1*.gz
#	gzip -d $c/${b}/${b}*_2*.gz
#done

#for F in $(cat $folder_fastq/sample.list.txt); do
#        b=${F##*/}
#	c=${F%/*}
#	mkdir -p $output_folder/${b}
#	metawrap read_qc -1 $c/${b}/${b}*_1.fq \
#	-2 $c/${b}/${b}*_2.fq \
#	-o $output_folder/${b} -x mm39
#done

# Rezip all the fastq files, otherwise it takes so much space
#So I Need to unzip them all, otherwise it does not work
#for F in $(cat $folder_fastq/sample.list.txt); do
#       b=${F##*/}
#       c=${F%/*}
#       echo $b
#       echo $c
#       gzip $c/${b}/${b}*_1*
#       gzip $c/${b}/${b}*_2*
#done

#now for the last feces, from F0 to F1
#input files:
folder_fastq_LF=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/raw_files_F0_till_F1/X204SC23116322-Z01-F005/01.RawData
output_folder_LF=/proj/naiss2024-23-57/C57_female_lineage_microbiota/samples/last_feces/trimmed_host_removed

#for F in $(cat $folder_fastq_LF/samples.list.txt); do
#        b=${F##*/}
#        c=${F%/*}
#	echo $b
#	gzip -d $c/${b}/${b}*_1*.gz
#	gzip -d $c/${b}/${b}*_2*.gz
#        mkdir -p $output_folder_LF/${b}
#        metawrap read_qc -1 $c/${b}/${b}*_1.fq \
#        -2 $c/${b}/${b}*_2.fq \
#        -o $output_folder_LF/${b} -x mm39
#done

#Also zip these files
for F in $(cat $folder_fastq_LF/samples.list.txt); do
       b=${F##*/}
       c=${F%/*}
       echo $b
       echo $c
       gzip $c/${b}/${b}*_1*
       gzip $c/${b}/${b}*_2*
done
