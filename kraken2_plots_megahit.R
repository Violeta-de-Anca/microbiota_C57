setwd("/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/")
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(ggrepel)
library(stringr)
library(forcats)
# library(miaViz)
# library(scater)
# library(SingleCellExperiment)

#for the color ####
complementary_col=function(n){
  if(n%%2==0){
    n_par=n/2
    hues=seq(0,360-360/n_par,length.out=n_par)
    color1=hcl(h=hues,c=100,l=65)
    color2=hcl(h=(hues+100)%%360,c=100,l=65)
    palette1=c(rbind(color1,color2))
  } else{
    n_par=floor(n/2)
    hues=seq(0,360-360/n_par,length.out=n_par)
    color1=hcl(h=hues,c=100,l=65)
    color2=hcl(h=(hues+100)%%360,c=100,l=65)
    palette1=c(rbind(color1,color2))
    palette1=c(palette1,hcl(h=0,c=0,l=50))
  }
  return(palette1)
}
#read the datasets for the cecum feces ####
SL_cecum=fread("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/small_megahit_litter_report_cecum_feces.txt")
C_cecum=fread("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/control_megahit_report_cecum_feces.txt")

names(SL_cecum)=gsub("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/","",names(SL_cecum))
names(C_cecum)=gsub("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/","",names(C_cecum))

SL_cecum$key <- gsub("Candidatus", "candidatus", SL_cecum$key)
C_cecum$key <- gsub("Candidatus", "candidatus", C_cecum$key)

#separate the key into the different classification of taxa ####
SL_cecum=SL_cecum%>%mutate(split=strsplit(key, " (?=[A-Z])",perl = T))%>%unnest_wider(split,names_sep = "_")
names(SL_cecum)[names(SL_cecum) %in% c("split_2","split_3","split_4","split_5","split_6","split_7","split_8")] = c("domain","kingdom","phylum","order","family","genus","specie")
C_cecum=C_cecum%>%mutate(split=strsplit(key, " (?=[A-Z])",perl = T))%>%unnest_wider(split,names_sep = "_")
names(C_cecum)[names(C_cecum) %in% c("split_2","split_3","split_4","split_5","split_6","split_7","split_8")] = c("domain","kingdom","phylum","order","family","genus","specie")

#filter the rows that have always a NA ####
SL_cecum=SL_cecum[rowSums(!is.na(SL_cecum[,-c(1,34:45)]))>0,]
C_cecum=C_cecum[rowSums(!is.na(C_cecum[,-c(1,78:89)]))>0,]

#separate the generations ####
F0_SL=names(SL_cecum)[grepl("F0",names(SL_cecum))]
F0_SL=c(names(SL_cecum)[c(1,34:45)],F0_SL)
F0_sl_cecum=SL_cecum[, F0_SL]


F0_c=names(C_cecum)[grepl("F0",names(C_cecum))]
F0_c=c(names(C_cecum)[c(1,78:89)],F0_c)
F0_c_cecum=C_cecum[, F0_c]

#separate the generations - F1 ####
F1_SL=names(SL_cecum)[grepl("F1",names(SL_cecum))]
F1_SL=c(names(SL_cecum)[c(1,34:45)],F1_SL)
F1_sl_cecum=SL_cecum[, F1_SL]

F1_c=names(C_cecum)[grepl("F1",names(C_cecum))]
F1_c=c(names(C_cecum)[c(1,78:89)],F1_c)
F1_c_cecum=C_cecum[, F1_c]

#separate the generations - F2 ####
F2_SL=names(SL_cecum)[grepl("F2",names(SL_cecum))]
F2_SL=c(names(SL_cecum)[c(1,34:45)],F2_SL)
F2_sl_cecum=SL_cecum[, F2_SL]

F2_c=names(C_cecum)[grepl("F2",names(C_cecum))]
F2_c=c(names(C_cecum)[c(1,78:89)],F2_c)
F2_c_cecum=C_cecum[, F2_c]

# Plot the abundances of the phylum in a stacked plot #####
#first convert data
F0_SL_cecum_long=pivot_longer(F0_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F0_SL_cecum_summary=F0_SL_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F0_SL_cecum_summary_bacteria=F0_SL_cecum_summary[!grepl("irae$",F0_SL_cecum_summary$phylum),]
F0_SL_cecum_summary_virus=F0_SL_cecum_summary[grepl("irae$",F0_SL_cecum_summary$phylum),]
factor(F0_SL_cecum_summary$phylum)
lvl_bac=length(levels(factor(F0_SL_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F0_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F0_SL_cecum_summary_bacteria[complete.cases(F0_SL_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   scale_y_continuous(labels = scales::percent_format())+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   theme_minimal()+
#   labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F0 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F0_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F0_SL_cecum_summary_virus[complete.cases(F0_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F0 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

F0_c_cecum_long=pivot_longer(F0_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F0_c_cecum_summary=F0_c_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F0_c_cecum_summary_bacteria=F0_c_cecum_summary[!grepl("irae$",F0_c_cecum_summary$phylum),]
F0_c_cecum_summary_virus=F0_c_cecum_summary[grepl("irae$",F0_c_cecum_summary$phylum),]
factor(F0_c_cecum_summary$phylum)
lvl_bac=length(levels(factor(F0_c_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F0_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F0_c_cecum_summary_bacteria[complete.cases(F0_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F0 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F0_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F0_c_cecum_summary_virus[complete.cases(F0_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F0 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

#merge the bacteria datasets ####
F0_c_cecum_summary_bacteria$group="Control"
F0_SL_cecum_summary_bacteria$group="Small_litter"
F0_cecum_summary_bacteria=rbind(F0_SL_cecum_summary_bacteria,F0_c_cecum_summary_bacteria)

#merge the virus datasets ####
F0_c_cecum_summary_virus$group="Control"
F0_SL_cecum_summary_virus$group="Small_litter"
F0_cecum_summary_virus=rbind(F0_SL_cecum_summary_virus,F0_c_cecum_summary_virus)


#F1 generation #####
# Plot the abundances of the phylum in a stacked plot #####
#first convert data
F1_SL_cecum_long=pivot_longer(F1_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F1_SL_cecum_summary=F1_SL_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F1_SL_cecum_summary_bacteria=F1_SL_cecum_summary[!grepl("irae$",F1_SL_cecum_summary$phylum),]
F1_SL_cecum_summary_virus=F1_SL_cecum_summary[grepl("irae$",F1_SL_cecum_summary$phylum),]
factor(F1_SL_cecum_summary$phylum)
lvl_bac=length(levels(factor(F1_SL_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F1_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F1_SL_cecum_summary_bacteria[complete.cases(F1_SL_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   scale_y_continuous(labels = scales::percent_format())+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   theme_minimal()+
#   labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F1 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F1_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F1_SL_cecum_summary_virus[complete.cases(F1_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F1 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

F1_c_cecum_long=pivot_longer(F1_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F1_c_cecum_summary=F1_c_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F1_c_cecum_summary_bacteria=F1_c_cecum_summary[!grepl("irae$",F1_c_cecum_summary$phylum),]
F1_c_cecum_summary_virus=F1_c_cecum_summary[grepl("irae$",F1_c_cecum_summary$phylum),]
factor(F1_c_cecum_summary$phylum)
lvl_bac=length(levels(factor(F1_c_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F1_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F1_c_cecum_summary_bacteria[complete.cases(F1_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F1 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F1_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F1_c_cecum_summary_virus[complete.cases(F1_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F1 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

#merge the bacteria datasets ####
F1_c_cecum_summary_bacteria$group="Control"
F1_SL_cecum_summary_bacteria$group="Small_litter"
F1_cecum_summary_bacteria=rbind(F1_SL_cecum_summary_bacteria,F1_c_cecum_summary_bacteria)

#merge the virus datasets ####
F1_c_cecum_summary_virus$group="Control"
F1_SL_cecum_summary_virus$group="Small_litter"
F1_cecum_summary_virus=rbind(F1_SL_cecum_summary_virus,F1_c_cecum_summary_virus)


#F2 generation #####
# Plot the abundances of the phylum in a stacked plot #####
#first convert data
F2_SL_cecum_long=pivot_longer(F2_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F2_SL_cecum_summary=F2_SL_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F2_SL_cecum_summary_bacteria=F2_SL_cecum_summary[!grepl("irae$",F2_SL_cecum_summary$phylum),]
F2_SL_cecum_summary_virus=F2_SL_cecum_summary[grepl("irae$",F2_SL_cecum_summary$phylum),]
factor(F2_SL_cecum_summary$phylum)
lvl_bac=length(levels(factor(F2_SL_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F2_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F2_SL_cecum_summary_bacteria[complete.cases(F2_SL_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   scale_y_continuous(labels = scales::percent_format())+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   theme_minimal()+
#   labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F2 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/SL_F2_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F2_SL_cecum_summary_virus[complete.cases(F2_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F2 small litter - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

F2_c_cecum_long=pivot_longer(F2_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F2_c_cecum_summary=F2_c_cecum_long%>%group_by(phylum, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F2_c_cecum_summary_bacteria=F2_c_cecum_summary[!grepl("irae$",F2_c_cecum_summary$phylum),]
F2_c_cecum_summary_virus=F2_c_cecum_summary[grepl("irae$",F2_c_cecum_summary$phylum),]
factor(F2_c_cecum_summary$phylum)
lvl_bac=length(levels(factor(F2_c_cecum_summary$phylum)))
color_bac=complementary_col(lvl_bac)
#now plot the data
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F2_cecum_bacteria.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F2_c_cecum_summary_bacteria[complete.cases(F2_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "F2 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+theme(legend.position = "none")
# dev.off()
# 
# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/control/C_F2_cecum_virus.tiff",width = 2000,height = 1000,units = "px",res = 150)
# ggplot(F2_c_cecum_summary_virus[complete.cases(F2_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=phylum),
#             position = position_fill(vjust = 0.5),
#             color="white", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+labs(
#     title = "Relative abundance of virus phylum",
#     x= "F2 control - cecum samples",
#     y= "Relative abundance",
#     fill= "Virus phylum"
#   )+theme(legend.position = "none")
# dev.off()

#merge the bacteria datasets ####
F2_c_cecum_summary_bacteria$group="Control"
F2_SL_cecum_summary_bacteria$group="Small_litter"
F2_cecum_summary_bacteria=rbind(F2_SL_cecum_summary_bacteria,F2_c_cecum_summary_bacteria)

#merge the virus datasets ####
F2_c_cecum_summary_virus$group="Control"
F2_SL_cecum_summary_virus$group="Small_litter"
F2_cecum_summary_virus=rbind(F2_SL_cecum_summary_virus,F2_c_cecum_summary_virus)


#Now merge all the cecum bacteria ####
F0_cecum_summary_bacteria$generation="F0"
F1_cecum_summary_bacteria$generation="F1"
F2_cecum_summary_bacteria$generation="F2"
#labels_x=interaction(cecum_summary_bacteria$group,cecum_summary_bacteria$generation,sep = " ")
cecum_summary_bacteria=rbind(F0_cecum_summary_bacteria,F1_cecum_summary_bacteria,F2_cecum_summary_bacteria)
cecum_summary_bacteria=cecum_summary_bacteria%>%
  mutate(generation=factor(generation,levels=c("F0","F1","F2")),
         group=factor(group,levels=c("Control","Small_litter")))%>%
  arrange(generation,group)%>%
  mutate(individual=factor(individual,levels=unique(individual)))
cecum_summary_bacteria=cecum_summary_bacteria%>% mutate(phylum=fct_reorder(phylum,total_value, .fun=mean,.desc = T))
lvl_bac=length(levels(factor(cecum_summary_bacteria$phylum)))
color_bac=complementary_col(lvl_bac)
colors <- c(
  "#FF0000", "#00FFFF", "#DFFF00", "#00DFFF", "#FF4000", "#00BFFF", "#FF6000", "#009FFF",
  "#FF8000", "#0080FF", "#FF9F00", "#0060FF", "#FFBF00", "#0040FF", "#FFDF00", "#0020FF",
  "#FFFF00", "#0000FF", "#FF2000", "#DF00FF", "#BFFF00", "#BF00FF", "#9FFF00", "#9F00FF",
  "#7FFF00", "#7F00FF", "#60FF00", "#6000FF", "#40FF00", "#4000FF", "#20FF00", "#2000FF",
  "#00FF00", "#FF00FF", "#00FFDF", "#FF00DF", "#00FFBF", "#FF00BF", "#00FF9F", "#FF009F",
  "#00FF80", "#FF007F", "#00FF60", "#FF005F", "#00FF40", "#FF003F", "#00FF20", "#FF001F"
)

# tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/cecum_bacteria_megahit_assembly.tiff",
#      width = 4000,height = 1500,units = "px",res = 150)
# ggplot(cecum_summary_bacteria,aes(x=individual,y=total_value,fill=phylum))+
#   geom_bar(stat = "identity",position = "fill")+
#   geom_text(aes(label=str_wrap(as.character(phylum),width = 5)),
#             position = position_fill(vjust = 0.5),
#             color="black", size=3, check_overlap = T)+
#   scale_y_continuous(labels = scales::percent_format())+
#   theme_minimal()+scale_fill_manual(values = colors)+labs(
#     title = "Relative abundance of bacteria phylum",
#     x= "Cecum samples",
#     y= "Relative abundance",
#     fill= "Bacteria phylum"
#   )+
#   theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))
# dev.off()

#merge and plot all the virus ####
#Now merge all the cecum bacteria ####
F0_cecum_summary_virus$generation="F0"
F1_cecum_summary_virus$generation="F1"
F2_cecum_summary_virus$generation="F2"
#labels_x=interaction(cecum_summary_virus$group,cecum_summary_virus$generation,sep = " ")
cecum_summary_virus=rbind(F0_cecum_summary_virus,F1_cecum_summary_virus,F2_cecum_summary_virus)
cecum_summary_virus=cecum_summary_virus%>%
  mutate(generation=factor(generation,levels=c("F0","F1","F2")),
         group=factor(group,levels=c("Control","Small_litter")))%>%
  arrange(generation,group)%>%
  mutate(individual=factor(individual,levels=unique(individual)))
cecum_summary_virus=cecum_summary_virus%>% mutate(phylum=fct_reorder(phylum,total_value, .fun=mean,.desc = T))
lvl_vir=length(levels(factor(cecum_summary_virus$phylum)))
color_vir=complementary_col(lvl_vir)

tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/cecum_virus_megahit_assembly.tiff",
     width = 4000,height = 1500,units = "px",res = 150)
ggplot(cecum_summary_virus,aes(x=individual,y=total_value,fill=phylum))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=str_wrap(as.character(phylum),width = 5)),
            position = position_fill(vjust = 0.5),
            color="black", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+scale_fill_manual(values = color_vir)+labs(
    title = "Relative abundance of virus phylum",
    x= "Cecum samples",
    y= "Relative abundance",
    fill= "Virus phylum"
  )+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

################################
#Now plot at the order level####
################################
# Plot the abundances of the order in a stacked plot #####
#first convert data
F0_SL_cecum_long=pivot_longer(F0_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F0_SL_cecum_summary=F0_SL_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F0_SL_cecum_summary_bacteria=F0_SL_cecum_summary[!grepl("ricota$",F0_SL_cecum_summary$order),]
F0_SL_cecum_summary_virus=F0_SL_cecum_summary[grepl("ricota$",F0_SL_cecum_summary$order),]
factor(F0_SL_cecum_summary$order)
lvl_bac=length(levels(factor(F0_SL_cecum_summary_bacteria$order)))
color_bac=complementary_col(lvl_bac)
#now plot the data
tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/F0_SL_order_lvl_cecum_bacteria_megahit_assembly.tiff",width = 2000,height = 1000,units = "px",res = 150)
ggplot(F0_SL_cecum_summary_bacteria,
       aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  scale_y_continuous(labels = scales::percent_format())+
  geom_text(aes(label= str_wrap(as.character(order),width = 10)),
            color="black", size=3, position = position_fill(vjust = 0.5),
            check_overlap = T)+
  theme_minimal()+
  scale_fill_manual(values = color_bac)+
  labs(
    title = "Relative abundance of bacteria order",
    x= "F0 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")
dev.off()

tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/small_litter/F0_SL_order_lvl_cecum_bacteria_megahit_assembly.tiff",width = 2000,height = 1000,units = "px",res = 150)

ggplot(F0_SL_cecum_summary_virus[complete.cases(F0_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F0 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")


F0_c_cecum_long=pivot_longer(F0_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F0_c_cecum_summary=F0_c_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F0_c_cecum_summary_bacteria=F0_c_cecum_summary[!grepl("ricota$",F0_c_cecum_summary$order),]
F0_c_cecum_summary_virus=F0_c_cecum_summary[grepl("ricota$",F0_c_cecum_summary$order),]
factor(F0_c_cecum_summary$order)
#now plot the data
ggplot(F0_c_cecum_summary_bacteria[complete.cases(F0_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of bacteria order",
    x= "F0 control - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")

ggplot(F0_c_cecum_summary_virus[complete.cases(F0_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F0 control - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")

#merge the bacteria datasets ####
F0_c_cecum_summary_bacteria$group="Control"
F0_SL_cecum_summary_bacteria$group="Small_litter"
F0_cecum_summary_bacteria=rbind(F0_SL_cecum_summary_bacteria,F0_c_cecum_summary_bacteria)

#merge the virus datasets ####
F0_c_cecum_summary_virus$group="Control"
F0_SL_cecum_summary_virus$group="Small_litter"
F0_cecum_summary_virus=rbind(F0_SL_cecum_summary_virus,F0_c_cecum_summary_virus)

#F1 generation #####
# Plot the abundances of the order in a stacked plot #####
#first convert data
F1_SL_cecum_long=pivot_longer(F1_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F1_SL_cecum_summary=F1_SL_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F1_SL_cecum_summary_bacteria=F1_SL_cecum_summary[!grepl("ricota$",F1_SL_cecum_summary$order),]
F1_SL_cecum_summary_virus=F1_SL_cecum_summary[grepl("ricota$",F1_SL_cecum_summary$order),]
factor(F1_SL_cecum_summary$order)
#now plot the data
ggplot(F1_SL_cecum_summary_bacteria[complete.cases(F1_SL_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  scale_y_continuous(labels = scales::percent_format())+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  theme_minimal()+
  labs(
    title = "Relative abundance of bacteria order",
    x= "F1 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")

ggplot(F1_SL_cecum_summary_virus[complete.cases(F1_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F1 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")


F1_c_cecum_long=pivot_longer(F1_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F1_c_cecum_summary=F1_c_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F1_c_cecum_summary_bacteria=F1_c_cecum_summary[!grepl("ricota$",F1_c_cecum_summary$order),]
F1_c_cecum_summary_virus=F1_c_cecum_summary[grepl("ricota$",F1_c_cecum_summary$order),]
factor(F1_c_cecum_summary$order)
#now plot the data
ggplot(F1_c_cecum_summary_bacteria[complete.cases(F1_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of bacteria order",
    x= "F1 control - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")

ggplot(F1_c_cecum_summary_virus[complete.cases(F1_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F1 control - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")

#merge the bacteria datasets ####
F1_c_cecum_summary_bacteria$group="Control"
F1_SL_cecum_summary_bacteria$group="Small_litter"
F1_cecum_summary_bacteria=rbind(F1_SL_cecum_summary_bacteria,F1_c_cecum_summary_bacteria)

#merge the virus datasets ####
F1_c_cecum_summary_virus$group="Control"
F1_SL_cecum_summary_virus$group="Small_litter"
F1_cecum_summary_virus=rbind(F1_SL_cecum_summary_virus,F1_c_cecum_summary_virus)


#F2 generation #####
# Plot the abundances of the order in a stacked plot #####
#first convert data
F2_SL_cecum_long=pivot_longer(F2_sl_cecum,cols = 14:22,names_to = "individual",values_to = "abundance",values_drop_na = T)
F2_SL_cecum_summary=F2_SL_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F2_SL_cecum_summary_bacteria=F2_SL_cecum_summary[!grepl("ricota$",F2_SL_cecum_summary$order),]
F2_SL_cecum_summary_virus=F2_SL_cecum_summary[grepl("ricota$",F2_SL_cecum_summary$order),]
factor(F2_SL_cecum_summary$order)
#now plot the data
ggplot(F2_SL_cecum_summary_bacteria[complete.cases(F2_SL_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  scale_y_continuous(labels = scales::percent_format())+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  theme_minimal()+
  labs(
    title = "Relative abundance of bacteria order",
    x= "F2 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")

ggplot(F2_SL_cecum_summary_virus[complete.cases(F2_SL_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F2 small litter - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")


F2_c_cecum_long=pivot_longer(F2_c_cecum,cols = 14:31,names_to = "individual",values_to = "abundance",values_drop_na = T)
F2_c_cecum_summary=F2_c_cecum_long%>%group_by(order, individual)%>%summarise(total_value=sum(abundance,na.rm = T), .groups = "drop")
F2_c_cecum_summary_bacteria=F2_c_cecum_summary[!grepl("ricota$",F2_c_cecum_summary$order),]
F2_c_cecum_summary_virus=F2_c_cecum_summary[grepl("ricota$",F2_c_cecum_summary$order),]
factor(F2_c_cecum_summary$order)
#now plot the data
ggplot(F2_c_cecum_summary_bacteria[complete.cases(F2_c_cecum_summary_bacteria),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of bacteria order",
    x= "F2 control - cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+theme(legend.position = "none")

ggplot(F2_c_cecum_summary_virus[complete.cases(F2_c_cecum_summary_virus),],aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=order),
            position = position_fill(vjust = 0.5),
            color="white", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+labs(
    title = "Relative abundance of virus order",
    x= "F2 control - cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+theme(legend.position = "none")

#merge the bacteria datasets ####
F2_c_cecum_summary_bacteria$group="Control"
F2_SL_cecum_summary_bacteria$group="Small_litter"
F2_cecum_summary_bacteria=rbind(F2_SL_cecum_summary_bacteria,F2_c_cecum_summary_bacteria)

#merge the virus datasets ####
F2_c_cecum_summary_virus$group="Control"
F2_SL_cecum_summary_virus$group="Small_litter"
F2_cecum_summary_virus=rbind(F2_SL_cecum_summary_virus,F2_c_cecum_summary_virus)


#Now merge all the cecum bacteria ####
F0_cecum_summary_bacteria$generation="F0"
F1_cecum_summary_bacteria$generation="F1"
F2_cecum_summary_bacteria$generation="F2"
cecum_summary_bacteria=rbind(F0_cecum_summary_bacteria,F1_cecum_summary_bacteria,F2_cecum_summary_bacteria)
cecum_summary_bacteria=cecum_summary_bacteria%>%
  mutate(generation=factor(generation,levels=c("F0","F1","F2")),
         group=factor(group,levels=c("Control","Small_litter")))%>%
  arrange(generation,group)%>%
  mutate(individual=factor(individual,levels=unique(individual)))
cecum_summary_bacteria=cecum_summary_bacteria%>% mutate(order=fct_reorder(order,total_value, .fun=mean,.desc = T))

lvl_bac=length(levels(factor(cecum_summary_bacteria$order)))
color_bac=complementary_col(lvl_bac)

tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/cecum_order_lvl_bacteria_megahit_assembly.tiff",width = 4000,height = 1500,units = "px",res = 150)
ggplot(cecum_summary_bacteria,aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=str_wrap(as.character(order),width = 5)),
            position = position_fill(vjust = 0.5),
            color="black", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+scale_fill_manual(values = color_bac)+labs(
    title = "Relative abundance of bacteria orders",
    x= "Cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

#let's do this but by whole groups and generations #####
ggplot(cecum_summary_bacteria,aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=str_wrap(as.character(order),width = 5)),
            position = position_fill(vjust = 0.5),
            color="black", size=3, check_overlap = T)+
  theme_minimal()+scale_fill_manual(values = color_bac)+labs(
    title = "Relative abundance of bacteria orders",
    x= "Cecum samples",
    y= "Relative abundance",
    fill= "Bacteria order"
  )+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1))

#merge and plot all the virus ####
#Now merge all the cecum bacteria ####
F0_cecum_summary_virus$generation="F0"
F1_cecum_summary_virus$generation="F1"
F2_cecum_summary_virus$generation="F2"
#labels_x=interaction(cecum_summary_virus$group,cecum_summary_virus$generation,sep = " ")
cecum_summary_virus=rbind(F0_cecum_summary_virus,F1_cecum_summary_virus,F2_cecum_summary_virus)
cecum_summary_virus=cecum_summary_virus%>%
  mutate(generation=factor(generation,levels=c("F0","F1","F2")),
         group=factor(group,levels=c("Control","Small_litter")))%>%
  arrange(generation,group)%>%
  mutate(individual=factor(individual,levels=unique(individual)))
cecum_summary_virus=cecum_summary_virus%>% mutate(phylum=fct_reorder(order,total_value, .fun=mean,.desc = T))
lvl_vir=length(levels(factor(cecum_summary_virus$order)))
color_vir=complementary_col(lvl_vir)

tiff("/proj/naiss2024-23-57/C57_female_lineage_microbiota/kraken2_taxonomy/cecum_feces/cecum_order_lvl_virus_megahit_assembly.tiff",
     width = 4000,height = 1500,units = "px",res = 150)
ggplot(cecum_summary_virus,aes(x=individual,y=total_value,fill=order))+
  geom_bar(stat = "identity",position = "fill")+
  geom_text(aes(label=str_wrap(as.character(order),width = 5)),
            position = position_fill(vjust = 0.5),
            color="black", size=3, check_overlap = T)+
  scale_y_continuous(labels = scales::percent_format())+
  theme_minimal()+scale_fill_manual(values = color_vir)+labs(
    title = "Abundance of virus order",
    x= "Cecum samples",
    y= "Relative abundance",
    fill= "Virus order"
  )+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

###########################################
#read the datasets for the last feces######
###########################################







