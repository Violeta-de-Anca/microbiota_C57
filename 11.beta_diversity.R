library(vegan)
library(phyloseq)
library(tidyverse)
library(patchwork)
library(agricolae)
library(data.table)
library(readr)
# library(FSA)
# library(rcompanion)
setwd("/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin")
set.seed(123)
#first load all the abundance tables for all the datasets
path_quan="/crex/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification"
list_quantification=list.files(path_quan,pattern = "\\.tab$",recursive = T,full.names = T)
df_quantification=lapply(list_quantification,fread)
names(df_quantification)=basename(list_quantification)
df_quantification=lapply(df_quantification, function(x) {subset(x, select= -`Genomic bins`)
  })
#also we cannot have 0s as then the bray curtis method goes bonkers
df_quantification=lapply(df_quantification, function(x_num) {
  x_num=as.data.frame(x_num)
  x_num= x_num[rowSums(x_num,na.rm = T)>0,,drop=F]
  x_num= x_num[,colSums(x_num,na.rm = T)>0,drop=F]
  x_num
})

# now calculate the distance matrix using the bray curtis method
distance_quantification=lapply(df_quantification,vegdist)
distance_quantification=lapply(df_quantification,as.matrix)
#we need to have the matrices in a phyloseq object
otu_df_quantification=lapply(df_quantification,otu_table,taxa_are_rows = T)
phylo_quant=lapply(otu_df_quantification,phyloseq)
names_samples_quanti=lapply(phylo_quant,sample_names)
#subset the multigenerational samples
phylo_quant_multigen=phylo_quant[grepl("F[0-2]",names(phylo_quant))]

phylo_quant_multigen=lapply(phylo_quant_multigen,function(x){
  sam=sample_names(x)
  group=ifelse(grepl("A",sam),"Overnutrition",
               ifelse(grepl("B",sam),"Control",NA))
  meta=data.frame(group=factor(group),
                  row.names = sam,
                  check.names = F)
  x=phyloseq(x,sample_data(meta))
})

# now let's do the PCoA, principal coordinate analysis
pcoa_quant_multigen=lapply(phylo_quant_multigen,ordinate,method = "PCoA",distance="bray")
example=plot_ordination(phylo_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab,
                pcoa_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab,type = "samples",
                label = sample_variables(phylo_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab),
                axes=2:3)
  geom_point(size=3)

#plotting F0 cecum
F0_cecum=plot_ordination(phylo_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab,
                pcoa_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab,type = "samples",
                justDF = T,axes = 3:4)

F0_cf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F0_cecum_per_individual.tab))
F0_cecum$group=F0_cf_names$group[match(rownames(F0_cecum),rownames(F0_cf_names))]
ggplot(F0_cecum,aes(x=Axis.3,y=Axis.4,colour = group))+geom_point(size=3)

#plotting F1 cecum, podria valer pero hay q representar familias para ver si se hacen centroides x familia y tratamiento
F1_cecum=plot_ordination(phylo_quant_multigen$bin_abundance_table_F1_cecum_per_individual.tab,
                         pcoa_quant_multigen$bin_abundance_table_F1_cecum_per_individual.tab,type = "samples",
                         justDF = T,axes = 3:4)

F1_cf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F1_cecum_per_individual.tab))
F1_cecum$group=F1_cf_names$group[match(rownames(F1_cecum),rownames(F1_cf_names))]
ggplot(F1_cecum,aes(x=Axis.3,y=Axis.4,colour = group))+geom_point(size=3)
F1_cecum=F1_cecum[rownames(F1_cecum)!="fastq",]
#centroids
centroids_F1_cf=F1_cecum%>%group_by(group)%>%summarise_at(vars(matches("Axis")),mean)
ggplot()+
  geom_point(data = F1_cecum,aes(x=Axis.3,y=Axis.4,colour = group,shape = "Samples"),size=3)+
  geom_point(data = centroids_F1_cf,aes(x=Axis.3,y=Axis.4,colour=group, shape="Centroids"),size=5)+
  scale_color_discrete(name="Cecum samples - F1")+
  scale_shape_manual(name="Point type",
                     values = c(Samples=16,Centroids=8))+
  theme_bw()


#plotting F2 cecum
F2_cecum=plot_ordination(phylo_quant_multigen$bin_abundance_table_F2_cecum_per_individual.tab,
                         pcoa_quant_multigen$bin_abundance_table_F2_cecum_per_individual.tab,type = "samples",
                         justDF = T,axes = 1:2)

F2_cf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F2_cecum_per_individual.tab))
F2_cecum$group=F2_cf_names$group[match(rownames(F2_cecum),rownames(F2_cf_names))]
ggplot(F2_cecum,aes(x=Axis.1,y=Axis.2,colour = group))+geom_point(size=2)

#plotting F0 last feces, podria valer el 2 y 3, necesitaria hacer centroides
F0_lf=plot_ordination(phylo_quant_multigen$bin_abundance_table_F0_last_feces_individuals.tab,
                         pcoa_quant_multigen$bin_abundance_table_F0_last_feces_individuals.tab,type = "samples",justDF = T,axes = 2:3)

F0_lf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F0_last_feces_individuals.tab))
F0_lf$group=F0_lf_names$group[match(rownames(F0_lf),rownames(F0_lf_names))]
F0_lf=F0_lf[rownames(F0_lf)!="fastq",]
#centroids
centroids_F0_lf=F0_lf%>%group_by(group)%>%summarise_at(vars(matches("Axis")),mean)
ggplot()+
  geom_point(data = F0_lf,aes(x=Axis.2,y=Axis.3,colour = group,shape = "Samples"),size=3)+
  geom_point(data = centroids_F0_lf,aes(x=Axis.2,y=Axis.3,colour=group, shape="Centroids"),size=5)+
  scale_color_discrete(name="Last feces - F0")+
  scale_shape_manual(name="Point type",
                     values = c(Samples=16,Centroids=8))+
  theme_bw()


#plotting F1 last feces, podria valer pero hay q representar familias para ver si se hacen centroides x familia y tratamiento
F1_lf=plot_ordination(phylo_quant_multigen$bin_abundance_table_F1_last_feces_individuals.tab,
                      pcoa_quant_multigen$bin_abundance_table_F1_last_feces_individuals.tab,type = "samples",
                      justDF = T,axes = 3:4)

F1_lf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F1_last_feces_individuals.tab))
F1_lf$group=F1_lf_names$group[match(rownames(F1_lf),rownames(F1_lf_names))]
F1_lf=F1_lf[rownames(F1_lf)!="fastq",]
#centroids
centroids_F1_lf=F1_lf%>%group_by(group)%>%summarise_at(vars(matches("Axis")),mean)
ggplot()+
  geom_point(data = F1_lf,aes(x=Axis.3,y=Axis.4,colour = group,shape = "Samples"),size=3)+
  geom_point(data = centroids_F1_lf,aes(x=Axis.3,y=Axis.4,colour=group, shape="Centroids"),size=5)+
  scale_color_discrete(name="Last feces - F1")+
  scale_shape_manual(name="Point type",
                     values = c(Samples=16,Centroids=8))+
  theme_bw()

#plotting F2 last feces, podria valer el 2 y 3, necesitaria hacer centroides
F2_lf=plot_ordination(phylo_quant_multigen$bin_abundance_table_F2_last_feces_individuals.tab,
                      pcoa_quant_multigen$bin_abundance_table_F2_last_feces_individuals.tab,type = "samples",justDF = T,axes = 2:3)

F2_lf_names=data.frame(group=sample_data(phylo_quant_multigen$bin_abundance_table_F2_last_feces_individuals.tab))
F2_lf$group=F2_lf_names$group[match(rownames(F2_lf),rownames(F2_lf_names))]
F2_lf=F2_lf[rownames(F2_lf)!="fastq",]
#centroids
centroids_F2_lf=F2_lf%>%group_by(group)%>%summarise_at(vars(matches("Axis")),mean)
ggplot()+
  geom_point(data = F2_lf,aes(x=Axis.2,y=Axis.3,colour = group,shape = "Samples"),size=3)+
  geom_point(data = centroids_F2_lf,aes(x=Axis.2,y=Axis.3,colour=group, shape="Centroids"),size=5)+
  scale_color_discrete(name="Last feces - F2")+
  scale_shape_manual(name="Point type",
                     values = c(Samples=16,Centroids=8))+
  theme_bw()

#permanova to test significance of different beta diversities, done with the out table
#first we need to transpose the dataframe, rows as individuals and comunity as columns, but it can be the distance matrix
test_otus=lapply(df_quantification,t)
relation_list=lapply(phylo_quant_multigen,function(ps){
  sd=as.data.frame(ps@sam_data)
  data.frame(individual=row.names(sd),
             group=sd$group,
             row.names = NULL)
})

test_otus_multigen=test_otus[grepl("F[0-2]",names(test_otus))]
test_otus_multigen=test_otus_multigen[grepl("indi",names(test_otus_multigen))]
relation_list=relation_list[grepl("indi",names(relation_list))]

#test is the order of the indiiduals is the same to do the test
order_ok=sapply(names(test_otus_multigen),function(nm){
  out=test_otus_multigen[[nm]]
  rel=relation_list[[nm]]
  ind_otu=row.names(out)
  ind_rel=rel$individual
  identical(ind_rel,ind_otu)
})

#check for NAs in the relation table, we have individuals that do not have assigned a group
sapply(relation_list,function(rel) any(is.na(rel$group))  )
#delete the individuals that do not have a group assigned while also doing the 



permanova_all_out_tables=Map(
  f=function(out,rel){
    keep=!is.na(rel$group)
    out2=out[keep, ,drop=F]
    rel2=rel[keep, , drop=F]
    rel2$group=droplevels(rel2$group)
    adonis2(out2~group,data = rel2,
            method = "bray",
            permutations = 9999) 
  },
  out=test_otus_multigen,
  rel=relation_list
)






