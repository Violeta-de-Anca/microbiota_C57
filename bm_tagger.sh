#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 5-00:00:00
#SBATCH -J bm_tagger
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bm_tagger.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/bm_tagger.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools metaWRAP/1.3.2 SRPRISM

#conda activate microbiota


cd /proj/naiss2024-23-57/C57_female_lineage_microbiota/databases/BMTAGGER_DB
#wget ftp://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/*fa.gz
#gunzip *fa.gz
#cat *fa > mm39.fa
#rm chr*.fa

#bmtool -d mm39.fa -o mm39.bitmask
srprism mkindex -i mm39.fa -o mm39.srprism

#conda deactivate
