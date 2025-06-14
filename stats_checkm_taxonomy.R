setwd("/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/megahit/")
library(data.table)

#import files #####
files=list.files(path = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/megahit/", pattern = "\\_megahit_checkm_taxonomy_abundance.txt$",full.names = T)

data_list=lapply(files, fread)
names(data_list)=tools::file_path_sans_ext(basename(files))
data_list=Filter(function(df) nrow(df)>0, data_list)

#we need to merge all the taxas and do the stats on them




















