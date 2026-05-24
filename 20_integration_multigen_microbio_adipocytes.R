library(ggplot2)
library(BiocManager)
library(basilisk)
# BiocManager::install("MOFA2", lib="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
Sys.setenv(BASILISK_HOME="/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/basilisk")
library(MOFA2, lib.loc = "/crex/proj/naiss2024-23-57/C57_female_lineage_adipocytes/bin/")
library(tidyverse)
library(stringr)
library(dplyr)
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

# #load the microbiota tables ####
# F0_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F0_cecum_samples/bin_abundance_table_F0_cecum_per_individual.tab",header = T)
# F1_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F1_cecum_samples/bin_abundance_table_F1_cecum_per_individual.tab",header = T)
# F2_cecum=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F2_cecum_samples/bin_abundance_table_F2_cecum_per_individual.tab",header = T)
# cecum_relational=data.frame(file_name_M=c(names(F0_cecum[,-1]),names(F1_cecum[,-1]),names(F2_cecum[,-1])))
# cecum_relational=cecum_relational%>%separate(file_name_M,into = c("generation","sample_name","group","sex","view"),sep = "_",remove = F)
# write.table(cecum_relational,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/cecum_relational_table.txt",quote = F,row.names = F,sep = "\t")
# 
# F0_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F0_last_feces/bin_abundance_table_F0_last_feces_individuals.tab", header = T)
# F1_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F1_last_feces/bin_abundance_table_F1_last_feces_individuals.tab",header = T)
# F2_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F2_last_feces/bin_abundance_table_F2_last_feces_individuals.tab",header = T)
# lf_relational=data.frame(file_name_LF=c(names(F0_lf[,-1]),names(F1_lf[,-1]),names(F2_lf[,-1])))
# relational_lf=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/last_feces_relational_table.txt",header = T)
# relational_lf=inner_join(relational_lf,lf_relational,by = c("file_name"="file_name_LF"))
# write.table(relational_lf,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/last_feces_C57_relational_table.txt",quote = F,row.names = F,sep = "\t")

# relational_microbiota_multigen=inner_join(cecum_relational,relational_lf,by = c("sample_name","generation","group","sex"))
# relational_microbiota_multigen=fread(file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/relational_microbiota_multigen.txt",header = T)
# 
# total_relational_table=inner_join(relational_microbiota_multigen,relational_table_C57,by = c("sample_name","generation","sex"))
relational_LF_multigen=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/last_feces_C57_relational_table.txt",header = T)
relational_CS_multigen=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/cecum_relational_table.txt",header = T)
relational_multiomics_C57=full_join(relational_CS_multigen,relational_LF_multigen, by = c("sample_name","generation","group","sex"))
relational_multiomics_C57$group[relational_multiomics_C57$group ==  "B"] = "C"
relational_multiomics_C57$group[relational_multiomics_C57$group ==  "A"] = "SL"
relational_multiomics_C57=full_join(relational_multiomics_C57,relational_table_C57,by = c("sample_name","generation","sex","group"))

# #creating the long dataframe for mefisto c(sample, feature, view, value, group)####
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/methylation_long_dataframe.rda")
# microbiota_functional_C57=fread("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/cluster_individual_per_bin_weighted_abundance_matrix.tsv",header = T)
# setDT(microbiota_functional_C57)
# microbiota_functional_C57_long=melt(
#   microbiota_functional_C57,
#   id.vars = "cluster",
#   variable.name = "col_name",
#   value.name = "value",
#   variable.factor = F
# )
# microbiota_functional_C57_long[, `:=`(
#   generation = sub("_.*","",col_name),
#   bin = stringr::str_extract(col_name,"bin\\.\\d+"),
#   view = fcase(
#     grepl("LF",col_name), "LF",
#     grepl("_M$",col_name), "CC",
#     default = NA_character_
#   ),
#   group = fcase(
#     grepl("_B_",col_name), "C",
#     grepl("_A_",col_name), "SL",
#     default = NA_character_
#   ),
#   sex = fcase(
#     grepl("female",col_name), "female",
#     grepl("male",col_name), "male",
#     default = NA_character_
#   )
# )]
# microbiota_functional_C57_long[,feature := paste(cluster,bin, sep = "_")]
# microbiota_functional_C57_long=microbiota_functional_C57_long[!is.na(view)& !is.na(sex) & !is.na(bin)]
# save(microbiota_functional_C57_long,
#      file="/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/microbiota_functional_C57_long.rda")
# # load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/microbiota_functional_C57_long.rda")
# # load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/methylation_long_dataframe.rda")
# 
# # #get the sample name
# microbiota_functional_C57_long$sample= sub(".*__","",microbiota_functional_C57_long$col_name)
# lookup=relational_multiomics_C57 |> as_tibble() |>
#   dplyr::select(sample_name,file_name_M,file_name.x) |>
#   pivot_longer(c(file_name_M,file_name.x), values_to = "sample") |>
#   filter(!is.na(sample)) |>
#   distinct(sample, sample_name)
# microbiota_functional_C57_long=microbiota_functional_C57_long |>
#   left_join(lookup,by = "sample")
# microbiota_functional_C57_long=microbiota_functional_C57_long |>
#   filter(!is.na(sample_name))
# lookup2=relational_table_C57%>%as.tibble() %>%dplyr::select(sample_name,sex)
# long_rmp_counts_C571=long_rmp_counts_C571%>%left_join(lookup2,by = c("sample"="sample_name"))
# final_microbiota_C57=tibble(sample = microbiota_functional_C57_long$sample_name,
#                                 feature = microbiota_functional_C57_long$feature,
#                                 value = microbiota_functional_C57_long$value,
#                                 group = microbiota_functional_C57_long$group,
#                                 view = microbiota_functional_C57_long$view,
#                                 generation = microbiota_functional_C57_long$generation,
#                                 sex = microbiota_functional_C57_long$sex)
# save(final_microbiota_C57,
#      file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/final_microbiota_C57.rda")
# save(long_rmp_counts_C571,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/final_adipocytes_C57.rda")
# 
# rm(list = ls())

# #merge methylation and abundance data
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/final_microbiota_C57.rda")
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/final_adipocytes_C57.rda")
# data_integration_C57_maternal_line=data.table::rbindlist(list(long_rmp_counts_C571,final_microbiota_C57),use.names = T, fill = T)
# rm(long_rmp_counts_C571)
# rm(final_microbiota_C57)
# #we also need to filter by MAD (median absolute deviation, which is more resistant to outliers), as there is too many rows and mofa cannot work with them
# top_features= function(dt,N){
#   v=dt[!is.na(value), 
#        .(spread = mad(value, na.rm = T)), 
#        by = feature][order(-spread)]
#   v$feature[seq_len(min(N, nrow(v)))]
# }
# #and we apply the same but by sexes also to the methylation data
# keep_meth=top_features(data_integration_C57_maternal_line[view == "methylation"], 50000)
# keep_cc=top_features(data_integration_C57_maternal_line[view == "CC"], 50000)
# keep_lf=top_features(data_integration_C57_maternal_line[view == "LF"], 50000)
# save(keep_meth,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/meth_keep.rda")
# save(keep_cc,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/keep_cc.rda")
# save(keep_lf,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/keep_lf.rda")
# 
# data_integration_C57_maternal_line=data_integration_C57_maternal_line[
#   (view == "methylation" & feature %in% keep_meth) |
#     (view == "CC" & feature %in% keep_cc) |
#     (view == "LF" & feature %in% keep_lf)
# ]
# data_integration_C57_maternal_line[, treatment := group]
# data_integration_C57_maternal_line[, group := sex]
# data_integration_C57_maternal_line[, generation_num := as.integer(sub("^F", "",generation))]
# 
# save(data_integration_C57_maternal_line,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/meth_microbiota_c57_maternal_line.rda")
# rm(list = ls())

#create de mofa object ####
load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/meth_microbiota_c57_maternal_line.rda")

#quality checks
# setDT(data_integration_C57_maternal_line)
# print(data_integration_C57_maternal_line[, .(
#   n_rows= .N,
#   n_keys= uniqueN(.SD),
#   n_sex = uniqueN(sex)
# ), by = view, .SDcols = c("sample", "feature", "sex")])
# 
# print(data_integration_C57_maternal_line[view == "LF", .N, by = "sex"])
# 
# ind_problem=data_integration_C57_maternal_line[, uniqueN(group), by=sample][V1>1, sample]
# 
# data_integration_C57_maternal_line[sample == "4_2_4", .(view,group,sex,generation)] |> unique()
# data_integration_C57_maternal_line[sample == "4_2_4" & view == "CC" & sex == "female"]

combination_C57=create_mofa(data=data_integration_C57_maternal_line,groups = "group")
combination_C57=set_covariates(combination_C57,covariates = "generation_num")
save(combination_C57,file="/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
rm(list = ls())

# now let's put the options ####
load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
# #now we need to set the options
options_data=get_default_data_options(combination_C57)
# gaussian if you are using normalized counts, poisson for raw counts
opts_model=get_default_model_options(combination_C57)
# opts_model$likelihoods[]="poisson"
# #number of components, or factors as they call them here to calculate
# opts_model$num_factors=10
# opts_training=get_default_training_options(combination_C57)
# opts_training$maxiter=1000
# opts_training$convergence_mode="slow"
# opts_training$seed=123
# opts_mefisto=get_default_mefisto_options(combination_C57)
# opts_mefisto$model_groups=T
# #
# combination_C57=prepare_mofa(combination_C57,model_options = opts_model,
#                           mefisto_options = opts_mefisto,
#                           training_options = opts_training,
#                           data_options = options_data)
# #
# combination_C57=run_mofa(combination_C57,
#                       outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto_run_ouput.hdf5",
#                       use_basilisk = T)
# save(combination_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto_ran.rda")

# #number of components, or factors as they call them here to calculate, we needed to rerun with 5 factors as there were warnings ####
# load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
# 
# options_data=get_default_data_options(combination_C57)
# # gaussian if you are using normalized counts, poisson for raw counts
# opts_model=get_default_model_options(combination_C57)
# # opts_model$likelihoods[]="poisson"
# opts_model$num_factors=5
# opts_training=get_default_training_options(combination_C57)
# opts_training$maxiter=1000
# opts_training$convergence_mode="slow"
# opts_training$seed=123
# opts_mefisto=get_default_mefisto_options(combination_C57)
# opts_mefisto$model_groups=T
# #
# combination_C57=prepare_mofa(combination_C57,model_options = opts_model,
#                              mefisto_options = opts_mefisto,
#                              training_options = opts_training,
#                              data_options = options_data)
# #
# combination_C57=run_mofa(combination_C57,
#                          outfile = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_run_ouput.hdf5",
#                          use_basilisk = T)
# save(combination_C57,file = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_ran.rda")
load("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_ran.rda")
#check if we have any type of bias with sequencing depth in factor 1, we need to import the amount of reads per individual
factors=get_factors(combination_C57,as.data.frame = T)
factor_1=subset(factors,factor=="Factor1")
#get the data overlapp per individual
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/data_overlapp_all_omics_C57.tiff",
     height = 2000, width = 3000, res = 150)
plot_data_overview(combination_C57)
dev.off()
#get the treatment differences ####
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/factors_visualization_by_sex.tiff",
     height = 2000, width = 2500, res = 150)
plot_variance_explained(combination_C57)
dev.off()

#get the variance that each factor explains ####
var_factors_sex_grouping=get_variance_explained(combination_C57)

#see if we have a correlation with any factor and the time/treatment
plot_factor_cor(combination_C57)

#check if the components variation also varies during the different generations
get_scales(combination_C57)
plot_factors_vs_cov(combination_C57,color_by = "generation_num",factors = c(4),warped = F)

plot_factors_vs_cov(female_57,color_by = "treatment",factors = "all")
plot_factors_vs_cov(female_57,color_by = "treatment",factors = c(4))
get_scales(female_57)

plot_factors_vs_cov(male_57,color_by = "treatment",factors = "all")
plot_factors_vs_cov(male_57,color_by = "treatment",factors = c(5))
get_scales(male_57)

#distribution of the factors per sex ####
female_57=subset_groups(combination_C57,groups = "female")
male_57=subset_groups(combination_C57,groups = "male")
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/factors_evolution_sex_grouping_color_treatment_onlyfemales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(female_57,color_by = "treatment",factors = 1:5)
dev.off()
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/factors_evolution_sex_grouping_color_treatmentallfactors_onlymales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(male_57,color_by = "treatment",factors = 1:5)
dev.off()

#get the scores of the factors so we can do our own density plots and test ####
md=samples_metadata(combination_C57)
factors_sex_grouping=get_factors(combination_C57)
female_factors_sex_grouping=as.data.frame(factors_sex_grouping["female"])
female_factors_sex_grouping=female_factors_sex_grouping%>%mutate(sample_name=row.names(female_factors_sex_grouping))
female_factors_sex_grouping=left_join(female_factors_sex_grouping,md,by = c("sample_name"="sample"))
male_factors_sex_grouping=as.data.frame(factors_sex_grouping["male"])
male_factors_sex_grouping=male_factors_sex_grouping%>%mutate(sample_name=row.names(male_factors_sex_grouping))
male_factors_sex_grouping=left_join(male_factors_sex_grouping,md,by = c("sample_name"="sample"))

tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/factor_4_male_density.tiff", height = 1500, width = 2500, res = 150)
male_factor_density=ggplot2::ggplot(male_factors_sex_grouping,mapping =aes(x=male.Factor4, color=treatment,fill = treatment))+
  geom_density(alpha=0.4,linewidth= 3)+
  theme_classic()+theme(
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28)
  )+
  scale_colour_manual(values = c("#85bc37","#e74269"),
                      breaks = c("C","SL"),
                      labels=c("Control","Small Litter"))+
  labs(x="Male factor 4",
       y="Density of the score",
       color="Treatment")+
  scale_fill_manual(values = c("#85bc37","#e74269"),
                    breaks = c("C","SL"),
                    labels=c("Control","Small Litter"),
                    guide="none")
male_factor_density
dev.off()

#now also plot female factors that seems different, 2 and 4
tiff(filename = "/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/mofa_integration/factor_4_female_density.tiff", height = 1500, width = 2500, res = 150)
female_factor_density=ggplot(female_factors_sex_grouping,aes(x=female.Factor4, color=treatment,fill = treatment))+
  geom_density(alpha=0.4,linewidth= 3)+
  theme_classic()+theme(
    axis.title = element_text(size = 30),
    axis.text  = element_text(size = 28),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 28)
  )+
  scale_colour_manual(values = c("#85bc37","#e74269"),
                      breaks = c("C","SL"),
                      labels=c("Control","Small Litter"))+
  labs(x="Female factor 4",
       y="Density of the score",
       color="Treatment")+
  scale_fill_manual(values = c("#85bc37","#e74269"),
                    breaks = c("C","SL"),
                    labels=c("Control","Small Litter"),
                    guide="none")
female_factor_density
dev.off()

#now we are going to see which of the factor, for each of the sex, have different distributions, with the kolmogorov-smirnov test ####
ks_dist=function(df,group_col = "treatment", g1= "C", g2= "SL",
                 factor_pattern= "\\.Factor\\d+$"){
  factor_cols=grep(factor_pattern,names(df),value = T)
  map_dfr(factor_cols,function(col){
    x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
    y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
    kt=suppressWarnings(stats::ks.test(x,y,exact=F))
    tibble(
      factor=col,
      n_C=length(x),
      n_SL=length(y),
      D= unname(kt$statistic),
      p_value=kt$p.value
    )
  })
}
female_ks_dist=ks_dist(female_factors_sex_grouping)
male_ks_dist=ks_dist(male_factors_sex_grouping)

#now we are going to see which of the factor, for each of the sex, have the same distributions, but different medians, with the mann whitney test ####
mw_dist=function(df,group_col = "treatment", g1= "C", g2= "SL",
                 factor_pattern= "\\.Factor\\d+$"){
  factor_cols=grep(factor_pattern,names(df),value = T)
  map_dfr(factor_cols,function(col){
    x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
    y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
    mw=suppressWarnings(stats::wilcox.test(x,y))
    tibble(
      factor=col,
      n_C=length(x),
      n_SL=length(y),
      D= unname(mw$statistic),
      p_value=mw$p.value
    )
  })
}
female_mw_dist=mw_dist(female_factors_sex_grouping)
male_mw_dist=mw_dist(male_factors_sex_grouping)









# #functional annotation of bins from each gen, each type of sample ####
# bin_dir="/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/functional_annotation"
# generations=c("F0","F1","F2")
# sample_type=c("cecum_samples","last_feces")
# gff_col=c("COG")
# all_COG=rbindlist(lapply(generations,function(gen){
#   rbindlist(lapply(sample_type,function(st){
#     gff_path=file.path(bin_dir,paste0(gen,"_",st),"bin_funct_annotations")
#     gff_files=list.files(gff_path,pattern = "_COG_ID\\.txt$",full.names = T)
#     if(length(gff_files)==0) return(NULL)
#     rbindlist(lapply(gff_files,function(f){
#       dt=fread(f,header = F,col.names = gff_col)
#       dt[, bin := sub("_COG_ID\\.txt$","",basename(f))]
#       dt[, generation := gen]
#       dt[, sample_type := st]
#       dt
#     }))
#   }))
# }))
