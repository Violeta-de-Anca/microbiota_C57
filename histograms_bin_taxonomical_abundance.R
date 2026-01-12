setwd("/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin/")
library(tidyverse)
library(data.table)
library(dplyr)
#load quantification of bins ####
direc="/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification"
subfolder=c("F0_last_feces","F1_last_feces","F2_last_feces")
paths=file.path(direc,subfolder,paste0("bin_abundance_table_",subfolder,"_individuals.tab"))
quantification=paths%>%set_names(subfolder)%>%map_dfr(read_tsv,.id = "generation")
colnames(quantification)[colnames(quantification)=="Genomic bins"]="bin"
quantification=quantification%>%pivot_longer(cols = starts_with("LF_"),
                                             names_to = "individuals",
                                             values_to = "relative_abundance")%>%
  mutate(treatment=case_when(str_detect(individuals,"_A_")~"small_litter",
                             str_detect(individuals,"_B_")~"control_litter", TRUE~NA_character_))

#load taxonomy of bins ####
direc="/proj/naiss2024-23-57/C57_female_lineage_microbiota/bin_metagenomics/"
subfolder=c("F0_last_feces","F1_last_feces","F2_last_feces")
paths=file.path(direc,subfolder,paste0("refined_libraries_megahit/metawrap_70_10_bins.stats"))
taxonomy=paths%>%set_names(subfolder)%>%map_dfr(read_tsv,.id = "generation")

#merge the two dataframes by generation and dataframe #####
last_feces=inner_join(taxonomy,quantification,by=c("generation","bin"))

# now do summary #####
summary_last_feces=last_feces%>%group_by(generation, lineage,treatment)%>%
  summarise(
    avg_relative_abundance=mean(relative_abundance,na.rm = T),
    sd_relative_abundance=sd(relative_abundance,na.rm=T),
    n=n()
  )
summary_last_feces=summary_last_feces%>%filter(!is.na(treatment))
summary_last_feces=summary_last_feces%>%filter(lineage!="Bacteria")
#now do the plotting ####
summary_last_feces_F0=summary_last_feces%>%filter(generation=="F0_last_feces")
tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F0_last_feces/relative_abundance_F0.tiff",
     width = 1000,height = 1000,units = "px",res = 150)
ggplot(summary_last_feces_F0,aes(x=lineage,y=avg_relative_abundance,fill=treatment))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_errorbar(aes(ymin=avg_relative_abundance-sd_relative_abundance,
                    ymax=avg_relative_abundance+sd_relative_abundance),
                width=0.2,
                position = position_dodge(0.7))+
  labs(
    title = "Average abundance of bins in F0 generation",
    x="Lineage",
    y="Average relative abundance (+-SD)",
    fill="Litter size"
  )+
  theme_minimal()+
  scale_fill_manual(values = c("control_litter"="#66C2A5","small_litter"="#FCB862"),
                    labels=c("control_litter"="Control litter","small_litter"="Small litter"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
        plot.margin = margin(t=10,r=10,b=10,l=20))
dev.off()

tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F1_last_feces/relative_abundance_F1.tiff",
     width = 1000,height = 1000,units = "px",res = 150)
summary_last_feces_F1=summary_last_feces%>%filter(generation=="F1_last_feces")
ggplot(summary_last_feces_F1,aes(x=lineage,y=avg_relative_abundance,fill=treatment))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_errorbar(aes(ymin=avg_relative_abundance-sd_relative_abundance,
                    ymax=avg_relative_abundance+sd_relative_abundance),
                width=0.2,
                position = position_dodge(0.7))+
  labs(
    title = "Average abundance of bins in F1 generation",
    x="Lineage",
    y="Average relative abundance (+-SD)",
    fill="Litter size"
  )+
  theme_minimal()+
  scale_fill_manual(values = c("control_litter"="#66C2A5","small_litter"="#FCB862"),
                    labels=c("control_litter"="Control litter","small_litter"="Small litter"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
        plot.margin = margin(t=10,r=10,b=10,l=20))
dev.off()

tiff(filename = "/proj/naiss2024-23-57/C57_female_lineage_microbiota/quantification/F2_last_feces/relative_abundance_F2.tiff",
     width = 1000,height = 1000,units = "px",res = 150)
summary_last_feces_F2=summary_last_feces%>%filter(generation=="F2_last_feces")
ggplot(summary_last_feces_F2,aes(x=lineage,y=avg_relative_abundance,fill=treatment))+
  geom_col(position = position_dodge(),width = 0.7)+
  geom_errorbar(aes(ymin=avg_relative_abundance-sd_relative_abundance,
                    ymax=avg_relative_abundance+sd_relative_abundance),
                width=0.2,
                position = position_dodge(0.7))+
  labs(
    title = "Average abundance of bins in F2 generation",
    x="Lineage",
    y="Average relative abundance (+-SD)",
    fill="Litter size"
  )+
  theme_minimal()+
  scale_fill_manual(values = c("control_litter"="#66C2A5","small_litter"="#FCB862"),
                    labels=c("control_litter"="Control litter","small_litter"="Small litter"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 14),
        plot.margin = margin(t=10,r=10,b=10,l=20))
dev.off()




