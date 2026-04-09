library(BiocManager)
library(basilisk)
# BiocManager::install("MOFA2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
Sys.setenv(BASILISK_HOME="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/basilisk")
library(MOFA2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(tidyverse)
library(plyr)
library(pheatmap)
library(purrr)
library(tibble)
library(data.table)
# BiocManager::install("karyoploteR", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(GenomicRanges)
library(karyoploteR, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
require(devtools)
# devtools::install_github("MiguelCastresana/anubix", lib="/crex/proj/naiss2024-23-57/reference_genomes")
library(ANUBIX , lib.loc = "/crex/proj/naiss2024-23-57/reference_genomes")
# install.packages("neat", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(neat, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
# install.packages("gprofiler2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(gprofiler2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(KEGGgraph)
library(KEGGREST)
# devtools::install_github("noriakis/ggkegg", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggkegg, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
# install.packages("ggfx", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggfx, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(ggraph)
library(igraph)
library(clusterProfiler)
library(tidygraph)
library(msigdbr)
library(org.Mm.eg.db)
# install.packages("pathfindR", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(pathfindR, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(biomaRt)
library(UpSetR)
setwd("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration")
set.seed(123)

#load the relational tables with the methylation data ####
relational_table_C57=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)

#load the microbiota tables ####
F0_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F0_cecum_samples/bin_abundance_table_F0_cecum_per_individual.tab",header = T)
F1_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F1_cecum_samples/bin_abundance_table_F1_cecum_per_individual.tab",header = T)
F2_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F2_cecum_samples/bin_abundance_table_F2_cecum_per_individual.tab",header = T)
cecum_relational=data.frame(file_name_M=c(names(F0_cecum[,-1]),names(F1_cecum[,-1]),names(F2_cecum[,-1])))
cecum_relational=cecum_relational%>%separate(file_name_M,into = c("generation","sample_name","group","sex","view"),sep = "_",remove = F)

F0_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F0_last_feces/bin_abundance_table_F0_last_feces_individuals.tab", header = T)
F1_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F1_last_feces/bin_abundance_table_F1_last_feces_individuals.tab",header = T)
F2_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F2_last_feces/bin_abundance_table_F2_last_feces_individuals.tab",header = T)
lf_relational=data.frame(file_name_LF=c(names(F0_lf[,-1]),names(F1_lf[,-1]),names(F2_lf[,-1])))
relational_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/last_feces_relational_table.txt",header = T)
relational_lf=inner_join(relational_lf,lf_relational,by = c("file_name"="file_name_LF"))

relational_microbiota_multigen=inner_join(cecum_relational,relational_lf,by = c("sample_name","generation","group","sex"))
relational_microbiota_multigen=fread(file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/relational_microbiota_multigen.txt",header = T)

total_relational_table=inner_join(relational_microbiota_multigen,relational_table_C57,by = c("sample_name","generation","sex"))

#functional annotation of bins from each gen, each type of sample ####
bin_dir="/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation"
generations=c("F0","F1","F2")
sample_type=c("cecum_samples","last_feces")
gff_col=c("COG")
all_COG=rbindlist(lapply(generations,function(gen){
  rbindlist(lapply(sample_type,function(st){
    gff_path=file.path(bin_dir,paste0(gen,"_",st),"bin_funct_annotations")
    gff_files=list.files(gff_path,pattern = "_COG_ID\\.txt$",full.names = T)
    if(length(gff_files)==0) return(NULL)
    rbindlist(lapply(gff_files,function(f){
      dt=fread(f,header = F,col.names = gff_col)
      dt[, bin := sub("_COG_ID\\.txt$","",basename(f))]
      dt[, generation := gen]
      dt[, sample_type := st]
      dt
    }))
  }))
}))
















