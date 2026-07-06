library(ggplot2)
library(BiocManager)
library(basilisk)
# BiocManager::install("MOFA2", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
Sys.setenv(BASILISK_HOME="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/basilisk")
library(MOFA2, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(tidyverse)
library(stringr)
library(dplyr)
library(pheatmap)
library(purrr)
library(tibble)
library(data.table)
# BiocManager::install("karyoploteR", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(GenomicRanges)
library(karyoploteR, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
require(devtools)
# devtools::install_github("MiguelCastresana/anubix", lib="/gorilla/proj/microbiota_project/reference_genomes")
library(ANUBIX , lib.loc = "/gorilla/proj/microbiota_project/reference_genomes")
# install.packages("neat", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(neat, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
# install.packages("gprofiler2", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(gprofiler2, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(KEGGgraph)
library(KEGGREST)
# devtools::install_github("noriakis/ggkegg", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(ggkegg, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
# install.packages("ggfx", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(ggfx, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(ggraph)
library(igraph)
library(clusterProfiler)
library(tidygraph)
library(msigdbr)
library(org.Mm.eg.db)
# install.packages("pathfindR", lib="/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(pathfindR, lib.loc = "/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/bin/")
library(biomaRt)
library(UpSetR)
setwd("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration")
set.seed(123)

#load the relational tables with the methylation data ####
relational_table_C57=fread("/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
relational_table_C57[,sex := tolower(sex)]
relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)

#load the microbiota tables ####
F0_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F0_cecum_samples/bin_abundance_table_F0_cecum_per_individual.tab",header = T)
F1_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F1_cecum_samples/bin_abundance_table_F1_cecum_per_individual.tab",header = T)
F2_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F2_cecum_samples/bin_abundance_table_F2_cecum_per_individual.tab",header = T)
cecum_relational=data.frame(file_name_M=c(names(F0_cecum[,-1]),names(F1_cecum[,-1]),names(F2_cecum[,-1])))
cecum_relational=cecum_relational%>%separate(file_name_M,into = c("generation","sample_name","group","sex","view"),sep = "_",remove = F)
write.table(cecum_relational,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/cecum_relational_table.txt",quote = F,row.names = F,sep = "\t")

F0_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F0_last_feces/bin_abundance_table_F0_last_feces_individuals.tab", header = T)
F1_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F1_last_feces/bin_abundance_table_F1_last_feces_individuals.tab",header = T)
F2_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F2_last_feces/bin_abundance_table_F2_last_feces_individuals.tab",header = T)
lf_relational=data.frame(file_name_LF=c(names(F0_lf[,-1]),names(F1_lf[,-1]),names(F2_lf[,-1])))
relational_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/last_feces_relational_table.txt",header = T)
relational_lf=inner_join(relational_lf,lf_relational,by = c("file_name"="file_name_LF"))
write.table(relational_lf,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/last_feces_C57_relational_table.txt",quote = F,row.names = F,sep = "\t")

# relational_microbiota_multigen=inner_join(cecum_relational,relational_lf,by = c("sample_name","generation","group","sex"))
# relational_microbiota_multigen=fread(file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/relational_microbiota_multigen.txt",header = T)
# 
# total_relational_table=inner_join(relational_microbiota_multigen,relational_table_C57,by = c("sample_name","generation","sex"))
relational_LF_multigen=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/last_feces_C57_relational_table.txt",header = T)
relational_CS_multigen=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/cecum_relational_table.txt",header = T)
relational_multiomics_C57=full_join(relational_CS_multigen,relational_LF_multigen, by = c("sample_name","generation","group","sex"))
relational_multiomics_C57$group[relational_multiomics_C57$group ==  "B"] = "C"
relational_multiomics_C57$group[relational_multiomics_C57$group ==  "A"] = "SL"
relational_multiomics_C57=full_join(relational_multiomics_C57,relational_table_C57,by = c("sample_name","generation","sex","group"))

# #creating the long dataframe for mefisto c(sample, feature, view, value, group)####
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/methylation_long_dataframe.rda")
# microbiota_functional_C57=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/cluster_individual_per_bin_weighted_abundance_matrix.tsv",header = T)
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
#      file="/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/microbiota_functional_C57_long.rda")
# # load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/microbiota_functional_C57_long.rda")
# # load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/methylation_long_dataframe.rda")
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
#      file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/final_microbiota_C57.rda")
# save(long_rmp_counts_C571,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/final_adipocytes_C57.rda")
# 
# rm(list = ls())

# #merge methylation and abundance data #####
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/final_microbiota_C57.rda")
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/final_adipocytes_C57.rda")
# data_integration_C57_maternal_line=data.table::rbindlist(list(long_rmp_counts_C571,final_microbiota_C57),use.names = T, fill = T)
# rm(long_rmp_counts_C571)
# rm(final_microbiota_C57)
# #we also need to filter by MAD (median absolute deviation, which is more resistant to outliers), as there is too many rows and mofa cannot work with them #####
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
# save(keep_meth,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/meth_keep.rda")
# save(keep_cc,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/keep_cc.rda")
# save(keep_lf,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/keep_lf.rda")
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
# save(data_integration_C57_maternal_line,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/meth_microbiota_c57_maternal_line.rda")
# rm(list = ls())

#create de mofa object ####
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/meth_microbiota_c57_maternal_line.rda")

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
# 
# combination_C57=create_mofa(data=data_integration_C57_maternal_line,groups = "group")
# combination_C57=set_covariates(combination_C57,covariates = "generation_num")
# save(combination_C57,file="/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
# rm(list = ls())
# 
# # now let's put the options ####
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
# # #now we need to set the options
# options_data=get_default_data_options(combination_C57)
# # gaussian if you are using normalized counts, poisson for raw counts
# opts_model=get_default_model_options(combination_C57)
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
#                       outfile = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto_run_ouput.hdf5",
#                       use_basilisk = T)
# save(combination_C57,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto_ran.rda")

# #number of components, or factors as they call them here to calculate, we needed to rerun with 5 factors as there were warnings ####
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M.rda")
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
#                          outfile = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_run_ouput.hdf5",
#                          use_basilisk = T)
# save(combination_C57,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_ran.rda")
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_ran.rda")

# #check if we have any type of bias with sequencing depth in factor 1, we need to import the amount of reads per individual
# factors=get_factors(combination_C57,as.data.frame = T)
# factor_1=subset(factors,factor=="Factor1")

# #get the data overlapp per individual
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/data_overlapp_all_omics_C57.tiff",
#      height = 2000, width = 3000, res = 150)
# plot_data_overview(combination_C57)
# dev.off()

# #get the treatment differences ####
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_visualization_by_sex.tiff",
#      height = 2000, width = 2500, res = 150)
# plot_variance_explained(combination_C57)
# dev.off()
# 
# #get the variance that each factor explains ####
# var_factors_sex_grouping=get_variance_explained(combination_C57)
# 
# #see if we have a correlation with any factor and the time/treatment
# plot_factor_cor(combination_C57)
# 
# #check if the components variation also varies during the different generations
# get_scales(combination_C57)
# plot_factors_vs_cov(combination_C57,color_by = "generation_num",factors = c(4),warped = F)
# 
# plot_factors_vs_cov(female_57,color_by = "treatment",factors = "all")
# plot_factors_vs_cov(female_57,color_by = "treatment",factors = c(4))
# get_scales(female_57)
# 
# plot_factors_vs_cov(male_57,color_by = "treatment",factors = "all")
# plot_factors_vs_cov(male_57,color_by = "treatment",factors = c(5))
# get_scales(male_57)
# 
# #distribution of the factors per sex ####
# female_57=subset_groups(combination_C57,groups = "female")
# male_57=subset_groups(combination_C57,groups = "male")
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_evolution_sex_grouping_color_treatment_onlyfemales.tiff", height = 2000, width = 2500, res = 150)
# plot_factors(female_57,color_by = "treatment",factors = 1:5)
# dev.off()
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_evolution_sex_grouping_color_treatmentallfactors_onlymales.tiff", height = 2000, width = 2500, res = 150)
# plot_factors(male_57,color_by = "treatment",factors = 1:5)
# dev.off()
# 
# #get the scores of the factors so we can do our own density plots and test ####
# md=samples_metadata(combination_C57)
# factors_sex_grouping=get_factors(combination_C57)
# female_factors_sex_grouping=as.data.frame(factors_sex_grouping["female"])
# female_factors_sex_grouping=female_factors_sex_grouping%>%mutate(sample_name=row.names(female_factors_sex_grouping))
# female_factors_sex_grouping=left_join(female_factors_sex_grouping,md,by = c("sample_name"="sample"))
# male_factors_sex_grouping=as.data.frame(factors_sex_grouping["male"])
# male_factors_sex_grouping=male_factors_sex_grouping%>%mutate(sample_name=row.names(male_factors_sex_grouping))
# male_factors_sex_grouping=left_join(male_factors_sex_grouping,md,by = c("sample_name"="sample"))
# 
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factor_4_male_density.tiff", height = 1500, width = 2500, res = 150)
# male_factor_density=ggplot2::ggplot(male_factors_sex_grouping,mapping =aes(x=male.Factor4, color=treatment,fill = treatment))+
#   geom_density(alpha=0.4,linewidth= 3)+
#   theme_classic()+theme(
#     axis.title = element_text(size = 30),
#     axis.text  = element_text(size = 28),
#     legend.title = element_text(size = 30),
#     legend.text = element_text(size = 28)
#   )+
#   scale_colour_manual(values = c("#85bc37","#e74269"),
#                       breaks = c("C","SL"),
#                       labels=c("Control","Small Litter"))+
#   labs(x="Male factor 4",
#        y="Density of the score",
#        color="Treatment")+
#   scale_fill_manual(values = c("#85bc37","#e74269"),
#                     breaks = c("C","SL"),
#                     labels=c("Control","Small Litter"),
#                     guide="none")
# male_factor_density
# dev.off()
# 
# #now also plot female factors that seems different, 2 and 4
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factor_4_female_density.tiff", height = 1500, width = 2500, res = 150)
# female_factor_density=ggplot(female_factors_sex_grouping,aes(x=female.Factor4, color=treatment,fill = treatment))+
#   geom_density(alpha=0.4,linewidth= 3)+
#   theme_classic()+theme(
#     axis.title = element_text(size = 30),
#     axis.text  = element_text(size = 28),
#     legend.title = element_text(size = 30),
#     legend.text = element_text(size = 28)
#   )+
#   scale_colour_manual(values = c("#85bc37","#e74269"),
#                       breaks = c("C","SL"),
#                       labels=c("Control","Small Litter"))+
#   labs(x="Female factor 4",
#        y="Density of the score",
#        color="Treatment")+
#   scale_fill_manual(values = c("#85bc37","#e74269"),
#                     breaks = c("C","SL"),
#                     labels=c("Control","Small Litter"),
#                     guide="none")
# female_factor_density
# dev.off()
# 
# #now we are going to see which of the factor, for each of the sex, have different distributions, with the kolmogorov-smirnov test ####
# ks_dist=function(df,group_col = "treatment", g1= "C", g2= "SL",
#                  factor_pattern= "\\.Factor\\d+$"){
#   factor_cols=grep(factor_pattern,names(df),value = T)
#   map_dfr(factor_cols,function(col){
#     x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
#     y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
#     kt=suppressWarnings(stats::ks.test(x,y,exact=F))
#     tibble(
#       factor=col,
#       n_C=length(x),
#       n_SL=length(y),
#       D= unname(kt$statistic),
#       p_value=kt$p.value
#     )
#   })
# }
# female_ks_dist=ks_dist(female_factors_sex_grouping)
# male_ks_dist=ks_dist(male_factors_sex_grouping)
# 
# #now we are going to see which of the factor, for each of the sex, have the same distributions, but different medians, with the mann whitney test ####
# mw_dist=function(df,group_col = "treatment", g1= "C", g2= "SL",
#                  factor_pattern= "\\.Factor\\d+$"){
#   factor_cols=grep(factor_pattern,names(df),value = T)
#   map_dfr(factor_cols,function(col){
#     x=df%>% filter(.data[[group_col]]==g1)%>%pull(.data[[col]])%>%na.omit()
#     y=df%>% filter(.data[[group_col]]==g2)%>%pull(.data[[col]])%>%na.omit()
#     mw=suppressWarnings(stats::wilcox.test(x,y))
#     tibble(
#       factor=col,
#       n_C=length(x),
#       n_SL=length(y),
#       D= unname(mw$statistic),
#       p_value=mw$p.value
#     )
#   })
# }
# female_mw_dist=mw_dist(female_factors_sex_grouping)
# male_mw_dist=mw_dist(male_factors_sex_grouping)
# 
# #get the top 100 windows which contributes the most to the relevant factors ####
# M_M_features=get_weights(male_57,as.data.frame = T,factors = c(4,5))
# M_M_features=M_M_features%>%pivot_wider(names_from = "factor",values_from = "value")
# factor4=M_M_features[,-4]
# top_factor4=factor4%>%arrange(desc(abs(Factor4)))%>%slice_head(n=100)
# factor5=M_M_features[,-3]
# top_factor5=factor5%>%arrange(desc(abs(Factor5)))%>%slice_head(n=100)
# save(top_factor5,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_5_M_M.rda")
# save(top_factor4,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_4_M_M.rda")
# 
# 
# #let's plot the top weights of methylation to see if we have more incidence in a chromosome #####
# methylation_features_M_M=M_M_features%>%filter(view=="methylation")
# factor4=as.data.frame(methylation_features_M_M$Factor4)
# 
# factor4=factor4%>%mutate(feature=as.character(feature))%>%
#   separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
#   separate(pos, into = c("start","end"),sep = "-",convert = T)
# factor4$factor="factor_8"
# gr4=GRanges(seqnames = factor4$chr,ranges = IRanges(start = factor4$start,end = factor4$end),
#             factor=factor4$factor,value=factor4$value)
# #guardalo como tabla para anotar
# bed_4=data.frame(
#   chrom=factor4$chr,
#   chromStart=factor4$start,
#   chromEnd=factor4$end,
#   score=0,
#   strand=".",
#   value=factor4$value,
#   factor=factor4$factor
# )
# write.table(bed_4,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_meth_factor_4.bed",
#             sep = "\t",
#             quote = F,
#             row.names = F,
#             col.names = F)
# 
# ##################################################################################################################
# # let's do the mefisto but only with the abundance from the bins, not with the weights of the clusterization #####
# ##################################################################################################################
# load the microbiota tables ####
# F0_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F0_cecum_samples/bin_abundance_table_F0_cecum_per_individual.tab",header = T)
# F1_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F1_cecum_samples/bin_abundance_table_F1_cecum_per_individual.tab",header = T)
# F2_cecum=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F2_cecum_samples/bin_abundance_table_F2_cecum_per_individual.tab",header = T)
# F0_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F0_last_feces/bin_abundance_table_F0_last_feces_individuals.tab", header = T)
# F1_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F1_last_feces/bin_abundance_table_F1_last_feces_individuals.tab",header = T)
# F2_lf=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/quantification/F2_last_feces/bin_abundance_table_F2_last_feces_individuals.tab",header = T)
# 
# # load the metadata #####
# relational_table_C57=fread("/gorilla/proj/microbiota_project/C57_female_lineage_adipocytes/qsea_results/qsea_with_file_path.txt")
# relational_table_C57[,sex := tolower(sex)]
# relational_table_C57=as.data.frame(relational_table_C57,stringsAsFactors=F)
# relational_table_C57$family=sub("_.*","", relational_table_C57$sample_name)
# relational_LF_multigen=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/last_feces_C57_relational_table.txt",header = T)
# relational_CS_multigen=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/cecum_relational_table.txt",header = T)
# relational_multiomics_C57=bind_rows(relational_CS_multigen%>%rename(sample=file_name_M),
#                                     relational_LF_multigen%>%rename(sample=file_name))
# relational_multiomics_C57$group[relational_multiomics_C57$group ==  "B"] = "C"
# relational_multiomics_C57$group[relational_multiomics_C57$group ==  "A"] = "SL"
# 
# # creating the long dataframe for mefisto c(sample, feature, view, value, group)####
# all_microbiota=list(F0_cecum=F0_cecum, F1_cecum=F1_cecum, F2_cecum=F2_cecum, F0_lf=F0_lf, F1_lf=F1_lf, F2_lf=F2_lf)
# long_all_microbiota=lapply(all_microbiota, function(x){
#   pivot_longer(x, cols = -"Genomic bins", names_to = "sample", values_to = "value")
# })
# long_all_microbiota=bind_rows(long_all_microbiota,.id = "source")
# long_all_microbiota=inner_join(long_all_microbiota,relational_multiomics_C57,by = c("sample"))
# long_all_microbiota=long_all_microbiota%>%filter(sample!="F2_mistery_female_m")
# 
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/final_adipocytes_C57.rda")
# #filter by at least 1 count in at least 80% of the total individuals
# long_rmp_counts_C57=long_rmp_counts_C571%>%group_by(feature)%>%filter(sum(value > 0) / n() >= 0.8)%>% ungroup()
# long_all_microbiota=long_all_microbiota%>% mutate(feature= paste(source, `Genomic bins`,sep = "_"))
# long_all_microbiota=long_all_microbiota[,c(6,10, 4,7,9,5,8)]
# top_features= function(dt,N){
#     v=dt[!is.na(value),
#          .(spread = mad(value, na.rm = T)),
#          by = feature][order(-spread)]
#     v$feature[seq_len(min(N, nrow(v)))]
#   }
# # and we apply the same but by sexes also to the methylation data
# long_rmp_counts_C571=as.data.table(long_rmp_counts_C57)
# keep_meth=top_features(long_rmp_counts_C571[view == "methylation"], 50000)
# long_rmp_counts_C57=long_rmp_counts_C571[
#   (view == "methylation" & feature %in% keep_meth)
# ]
# 
# long_all_microbiota=as.data.table(long_all_microbiota)
# setnames(long_all_microbiota, "sample_name","sample")
# full_M_M=rbind(long_all_microbiota,long_rmp_counts_C57)
# full_M_M[, treatment := group]
# full_M_M[, group := sex]
# full_M_M[, generation_num := as.integer(sub("^F", "",generation))]
# full_M_M$value=log2(full_M_M$value + 1)
# save(full_M_M,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/M_M_taxo_abun_final.rda")
# rm(list = ls())
# 
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/M_M_taxo_abun_final.rda")
# combination_C57=create_mofa(data=full_M_M,groups = "group")
# combination_C57=set_covariates(combination_C57,covariates = "generation_num")
# save(combination_C57,file="/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M_taxo_abun.rda")
# rm(list = ls())
# 
# load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/first_mefisto_C57_M_M_taxo_abun.rda")
# options_data=get_default_data_options(combination_C57)
# 
# # gaussian if you are using normalized counts, poisson for raw counts
# opts_model=get_default_model_options(combination_C57)
# 
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
#                          outfile = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_taxo_run_ouput.hdf5",
#                          use_basilisk = T)
# save(combination_C57,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_taxo_ran.rda")

load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/c57_M_M_mefisto5_taxo_ran.rda")

#is factor 1 due to technical bias? or is it biology? ####
factors=get_factors(combination_C57)
data=get_data(combination_C57)
nonzero=lapply(names(data),function(view){
  combined=do.call(cbind,data[[view]])
  n_nonzero=colSums(combined != 0 & !is.na(combined))
  data.frame(
    sample= names(n_nonzero),
    n_nonzero = as.numeric(n_nonzero),
    view = view
  )
})
non_zero=do.call(rbind, nonzero)
factor1=get_factors(combination_C57,factors = 1)
factor_df=do.call(rbind, lapply(names(factor1),function(g){
  data.frame(
    sample = names(factor1[[g]][,1]),
    factor1=as.numeric(factor1[[g]][,1]),
    group = g
  )
}))
merged=merge(non_zero,factor_df,by = "sample")

for(v in unique(merged$view)){
  sub=merged[merged$view==v,]
  r=cor(sub$factor1, sub$n_nonzero,use = "complete.obs")
  cat(v, ": r =",round(r,3), "\n")
}

#let's get the variance explained ####
# get the treatment differences ####
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_visualization_taxo__by_sex.tiff",
     height = 2000, width = 2500, res = 150)
plot_variance_explained(combination_C57)
dev.off()
# 
# #get the variance that each factor explains ####
var_factors_sex_grouping=get_variance_explained(combination_C57)

# #see if we have a correlation with any factor and the time/treatment
plot_factor_cor(combination_C57)
# 
# #check if the components variation also varies during the different generations
get_scales(combination_C57)
MOFA2::plot_factors_vs_cov(combination_C57,color_by = "generation_num",factors = c(4),warped = F)
# 
female_57=subset_groups(combination_C57,groups = "female")
male_57=subset_groups(combination_C57,groups = "male")
MOFA2::plot_factors_vs_cov(female_57,color_by = "treatment",factors = "all")
MOFA2::plot_factors_vs_cov(female_57,color_by = "treatment",factors = c(4))
get_scales(female_57)
# 
plot_factors_vs_cov(male_57,color_by = "treatment",factors = "all")
plot_factors_vs_cov(male_57,color_by = "treatment",factors = c(1))
# get_scales(male_57)
# 
#distribution of the factors per sex ####
female_57=subset_groups(combination_C57,groups = "female")
male_57=subset_groups(combination_C57,groups = "male")
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_taxo_evolution_sex_grouping_color_treatment_onlyfemales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(female_57,color_by = "treatment",factors = 1:5)
dev.off()
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factors_taxo_evolution_sex_grouping_color_treatmentallfactors_onlymales.tiff", height = 2000, width = 2500, res = 150)
plot_factors(male_57,color_by = "treatment",factors = 1)
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

tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factor_1_taxo_male_density.tiff", height = 1500, width = 2500, res = 150)
male_factor_density=ggplot2::ggplot(male_factors_sex_grouping,mapping =aes(x=male.Factor5, color=treatment,fill = treatment))+
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
  labs(x="Male factor 1",
       y="Density of the score",
       color="Treatment")+
  scale_fill_manual(values = c("#85bc37","#e74269"),
                    breaks = c("C","SL"),
                    labels=c("Control","Small Litter"),
                    guide="none")
male_factor_density
dev.off()
# 
# for females
# tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/factor_4_female_density.tiff", height = 1500, width = 2500, res = 150)
# female_factor_density=ggplot(female_factors_sex_grouping,aes(x=female.Factor5, color=treatment,fill = treatment))+
#   geom_density(alpha=0.4,linewidth= 3)+
#   theme_classic()+theme(
#     axis.title = element_text(size = 30),
#     axis.text  = element_text(size = 28),
#     legend.title = element_text(size = 30),
#     legend.text = element_text(size = 28)
#   )+
#   scale_colour_manual(values = c("#85bc37","#e74269"),
#                       breaks = c("C","SL"),
#                       labels=c("Control","Small Litter"))+
#   labs(x="Female factor 4",
#        y="Density of the score",
#        color="Treatment")+
#   scale_fill_manual(values = c("#85bc37","#e74269"),
#                     breaks = c("C","SL"),
#                     labels=c("Control","Small Litter"),
#                     guide="none")
# female_factor_density
# dev.off()
# 
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

#get the top 100 windows which contributes the most to the relevant factors, and extract by omic!! ####
M_M_features=get_weights(male_57,as.data.frame = T,factors = c(3))
M_M_features=M_M_features%>%pivot_wider(names_from = "factor",values_from = "value")
top_factor3=M_M_features%>%group_by(view)%>%slice_max(abs(Factor3),n = 100)%>%ungroup()
save(top_factor3,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_3_taxo_M_M.rda")

factor_1_M_M=get_weights(male_57,as.data.frame = T,factors = c(1))
factor_1_M_M=factor_1_M_M%>%pivot_wider(names_from = "factor",values_from = "value")
top_factor1=factor_1_M_M%>%group_by(view)%>%slice_max(abs(Factor1),n = 100)%>%ungroup()
save(top_factor1,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_1_taxo_M_M.rda")


load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_3_taxo_M_M.rda")
load("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_factor_1_taxo_M_M.rda")

#let's plot the top weights of methylation to see if we have more incidence in a chromosome #####
methylation_features_M_M=top_factor3%>%filter(view=="methylation")
meth_1_M_M=top_factor1%>%filter(view=="methylation")

factor3=methylation_features_M_M%>%mutate(feature=as.character(feature))%>%
  separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
  separate(pos, into = c("start","end"),sep = "-",convert = T)
factor3$factor="factor_3"
gr3=GRanges(seqnames = factor3$chr,ranges = IRanges(start = factor3$start,end = factor3$end),
            factor=factor3$factor,value=factor3$Factor3)
factor1=meth_1_M_M%>%mutate(feature=as.character(feature))%>%
  separate(feature,into = c("chr","pos"),sep = ":",remove = F)%>%
  separate(pos, into = c("start","end"),sep = "-",convert = T)
factor1$factor="factor_1"
gr1=GRanges(seqnames = factor1$chr,ranges = IRanges(start = factor1$start,end = factor1$end),
            factor=factor1$factor,value=factor1$Factor1)

#guardalo como tabla para anotar
bed_3=data.frame(
  chrom=factor3$chr,
  chromStart=factor3$start,
  chromEnd=factor3$end,
  score=0,
  strand=".",
  value=factor3$Factor3,
  factor=factor3$factor
)
write.table(bed_3,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_meth_factor_3_taxo.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

bed_1=data.frame(
  chrom=factor1$chr,
  chromStart=factor1$start,
  chromEnd=factor1$end,
  score=0,
  strand=".",
  value=factor1$Factor1,
  factor=factor1$factor
)
write.table(bed_1,file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_meth_factor_1_taxo.bed",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = F)

#plot them in a graph the incidence in chromosomes ####
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_weight_1_3_meth_loc_male_taxo.tiff",
     height = 3000, width = 1500, res = 150)
cyto=getCytobands(genome = "mm39")
cyto$gieStain[cyto$name=="p"]="pArm"
cyto$gieStain[cyto$name=="q"]="qArm"
pp=getDefaultPlotParams(plot.type = 2)
pp$data2height=30
bw=c(
  gneg= "white",
  qArm= "grey80",
  gpos50 = "grey60",
  gpos75 = "grey40",
  gpos100 = "black",
  gvar = "grey70",
  stalk = "grey50",
  acen = "black"
)
mm10=getCytobands(genome = "mm10")
kp=plotKaryotype(genome = "mm39",plot.type = 1,plot.params = pp,
                 ideogram.plotter = NULL,cex=2)
kpAddCytobands(kp,color.table = bw,color.schema = "biovizbase")
kpAddBaseNumbers(kp)
kpPlotRegions(kp,data = gr1,col = "#D55E00",border = "#D55E00",avoid.overlapping = T)
kpPlotRegions(kp,data = gr3,col = "#009E73",border = "#009E73",avoid.overlapping = T)
legend("right",
       legend = c("Factor 1","Factor 3"),
       fill = c("#D55E00","#009E73"),
       cex = 2,
       bty = "n")

dev.off()

#let's get the top weights of microbiota #####
m_features_M_M=top_factor3%>%filter(view=="M")
lf_features_M_M=top_factor3%>%filter(view=="last_feces")
microbiota_top_features=rbind(m_features_M_M,lf_features_M_M)
write.table(microbiota_top_features,
            file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_microbiota.txt",
            quote = F,row.names = F,sep = "\t")

m1_features_M_M=top_factor1%>%filter(view=="M")
lf1_features_M_M=top_factor1%>%filter(view=="last_feces")
microbiota_top_features1=rbind(m1_features_M_M,lf1_features_M_M)
write.table(microbiota_top_features1,
            file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_factor_1_microbiota.txt",
            quote = F,row.names = F,sep = "\t")

# get the factor correlation against time ####
corr_scores_time=plot_factors_vs_cov(combination_C57,color_by = "treatment",factors = c(1,3),shape_by = "group",scale = T,return_data = T) #

tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/correlation_sig_factors_time.tiff",
     height = 1500, width = 2000, res = 150)
ggplot(corr_scores_time,aes(x=value.covariate,y=value.factor,colour = group,shape = factor))+geom_jitter(size=6)+
  scale_color_manual(values = c("#CC79A7","#009E73"),labels=c("Female","Male"))+
  scale_shape_manual(values = c(8,19),labels=c("Factor 1","Factor 3"))+
  labs(x="Generations",y="Scores from the significant factors",colour="Sex",shape="Factor")+
  geom_smooth(method = "auto",se = T)+theme_minimal(base_size = 30)
dev.off()

# heatmaps of the relevant bins #####
male_dat=get_data(male_57,views = "methylation")[["methylation"]][[1]]
male_dat=as.matrix(male_dat)
w=get_weights(male_57,views = "methylation",factors = c(1))[["methylation"]]
w3=get_weights(male_57,views = "methylation",factors = c(3))[["methylation"]]
top1=names(sort(abs(w[,1]),decreasing = T))[1:100]
top3=names(sort(abs(w3[,1]),decreasing = T))[1:100]
top=c(top1,top3)
male_mat=male_dat[top,]
meta_male=male_57@samples_metadata
meta_male=meta_male[match(colnames(male_mat),meta_male$sample),]
ord=order(meta_male$treatment)
male_mat=male_mat[,ord]
meta_male=meta_male[ord,]



# start to work with the relevant bins of the microbiota #####
M_f1=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_factor_1_microbiota.txt",header = T)
M_f3=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/top_windows_microbiota.txt",header = T)
clusters_M=fread("/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/functional_annotation/cluster_bin_protein_proportion_matrix.tsv",header = T)
cc_f3=M_f3[view=="M",]
lf_f3=M_f3[view=="last_feces",]

# but only work with the cc as that is the one that explains more variability #####
cc_f3[,feature_col := gsub("_cecum_bin\\.","_cecum_samples_bin.",feature)]
matched=cc_f3$feature_col[cc_f3$feature_col %in%names(clusters_M)]
f3_clusters=clusters_M[, c("cluster",matched),with = F]
#all clusters
f3_clusters=f3_clusters[rowSums(!is.na(f3_clusters[,-"cluster"]))>0,]
# clusters that are in at least 2 groups
f3_clusters_min2=f3_clusters[rowSums(!is.na(f3_clusters[,-"cluster"]))>1,]
# clusters that are always affected in all groups
f3_clusters_all=f3_clusters[rowSums(is.na(f3_clusters[,-"cluster"]))==0,]

#now do this per generation
gen_cols=list(
  F0=grep("^F0_", matched, value = TRUE),
  F1 = grep("^F1_", matched, value = TRUE),
  F2 = grep("^F2_", matched, value = TRUE)
)

# ── Per-cluster stats for each generation ────────────────────────────────────
gen_stats <- rbindlist(lapply(names(gen_cols), function(gen) {
  cols <- gen_cols[[gen]]
  if (length(cols) == 0) return(NULL)
  mat  <- as.matrix(f3_clusters[, cols, with = FALSE])
  data.table(
    generation   = gen,
    cluster      = f3_clusters$cluster,
    n_bins       = length(cols),
    n_nonNA      = rowSums(!is.na(mat)),
    has_any      = rowSums(!is.na(mat)) > 0,
    all_complete = rowSums(is.na(mat)) == 0
  )
}))

# ── Summary table ─────────────────────────────────────────────────────────────
summary_gen <- gen_stats[, .(
  n_bins              = unique(n_bins),
  n_clusters_total    = .N,
  n_clusters_any_nonNA = sum(has_any),
  n_clusters_complete = sum(all_complete),
  mean_nonNA_bins     = round(mean(n_nonNA), 2)
), by = generation]

print(summary_gen)

# ── Classify each cluster×generation into a category ─────────────────────────
gen_stats[, status := fcase(
  all_complete, "Complete (no NAs)",
  has_any,      "Partial (some NAs)",
  default =     "All NA"
)]

plot_data <- gen_stats[, .N, by = .(generation, status)]

# ── Plot ──────────────────────────────────────────────────────────────────────
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/clusters_in_relevant_features_factor_3.tiff",
     height = 1500, width = 2500, res = 200)
ggplot(plot_data[status !="All NA"], aes(x = generation, y = N, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(
    "Complete (no NAs)"  = "#2ecc71",
    "Partial (some NAs)" = "#f39c12"
  ),
  labels=c(
    "Complete (no NAs)"  = "Clusters with incidence in all\nCC features in factor 3",
    "Partial (some NAs)" = "Clusters with partial incidence in\nCC features in factor 3"
  )) +
  labs(
    x     = "Generation",
    y     = "Number of clusters",
    fill  = "Status"
  ) +
  theme_minimal(base_size = 13)
dev.off()

#first do an upset to also see which is the overlapp between clusters by generation
library(UpSetR)
write.table(gen_stats,file = "summary_clusters_in_factor3.txt",quote = F,sep = "\t",row.names = F)
gens <- c("F0", "F1", "F2")

# --- long -> wide binary membership matrix -------------------------------
# swap value.var to "has_any" for detection, "all_complete" for completeness
membership_var <- "has_any"
gen_stats[, .(present=sum(has_any),n=.N),by=generation]

wide <- dcast(gen_stats, cluster ~ generation, value.var = membership_var)

# logical -> 0/1; a cluster absent from a generation becomes NA -> 0
wide[, (gens) := lapply(.SD, \(x) as.integer(fifelse(is.na(x), FALSE, x))),
     .SDcols = gens]

df <- as.data.frame(wide)          # UpSetR uses only the columns named in `sets`

# --- plot ----------------------------------------------------------------
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/upset_overlapp_clusters_factor_3.tiff",
     height = 1500, width = 2500, res = 200)
upset(df,
      sets          = rev(gens),   # rev() puts F0 at the top of the matrix
      keep.order    = TRUE,        # don't reorder sets by size
      order.by      = "freq",      # intersections largest-first
      main.bar.color = "#37474F",
      sets.bar.color = "#90A4AE",
      text.scale    = 1.3,
      mainbar.y.label = "Clusters in intersection",
      sets.x.label    = "Clusters per generation")
dev.off()
#let's also do a heatmap to see how many bins are per cluster that is affected
library(data.table)
library(pheatmap)

gens <- c("F0", "F1", "F2")

gen_stats[, pct := 100 * n_nonNA / n_bins]   # percentage of bins with data
keep <- gen_stats[n_nonNA > 0]               # drop zero-coverage observations
keep = keep[pct > 50]
w <- dcast(keep, cluster ~ generation, value.var = "pct")
m <- as.matrix(w[, ..gens])
rownames(m) <- w$cluster                     # kept for export, not shown

m <- m[order(-rowMeans(m, na.rm = TRUE)), ]  # NA-safe row ordering
tiff(filename = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/mofa_integration/heatmap_percentage50atleast_bins_per_cluster_factor_3.tiff",
     height = 2500, width = 1500, res = 200)
pheatmap(m,
         cluster_rows  = FALSE, cluster_cols = FALSE,   # keep F0→F1→F2 order
         show_rownames = FALSE,
         na_col        = "grey90",                      # cells deleted / absent
         color         = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
         breaks = seq(0,100,length.out=101),
         main          = "",
         legend_breaks = c(0,25,50,75,100),
         legend_labels = c("0%","25%","50%","75%","100%"))

dev.off()

#now load the functional annotations of the clusters ###
functanno=fread(file = "/gorilla/proj/microbiota_project/C57_female_lineage_microbiota/")






