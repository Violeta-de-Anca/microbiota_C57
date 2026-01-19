#!/bin/bash
#SBATCH -A uppmax2025-2-302
#SBATCH -p pelle
#SBATCH --mem 100GB
#SBATCH -t 10:00:00
#SBATCH -J srprism
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/srprism.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/srprism.out

#for metawrap removal of host
#srprism
#/sw/arch/eb/software/SRPRISM/3.3.2-GCCcore-13.3.0-Java-17/bin/srprism mkindex ../../reference_genomes/mus_musculus/uscs_ref/mm39.fa -o ../../reference_genomes/mus_musculus/uscs_ref/mm39.srprism

#bmtool
/sw/arch/eb/software/bmtagger/3.101-gompi-2024a-Java-17/bin/bmtool -d ../../reference_genomes/mus_musculus/uscs_ref/mm39.fa -o ../../reference_genomes/mus_musculus/uscs_ref/mm39.bitmask
