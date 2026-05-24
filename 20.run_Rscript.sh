#!/bin/bash -l
##SBATCH -A uppmax2025-2-302
##SBATCH -A uppmax2025-2-536
#SBATCH -A uppmax2026-1-34
#SBATCH -p pelle
#SBATCH --mem 250GB
#SBATCH -t 24:00:00
#SBATCH -J mef_sex
#SBATCH --error /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mefisto_microbiotafunctional_and_adipocites.err
#SBATCH --output /proj/naiss2024-23-57/C57_female_lineage_microbiota/log_files/mefisto_microbiotafunctional_and_adipocites.out

module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load R-bundle-CRAN/2024.11-foss-2024a

R --no-save --quiet < 20_integration_multigen_microbio_adipocytes.R
