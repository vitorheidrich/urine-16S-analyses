library(ggplot2)
library(decontam)
library(qiime2R)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(ggpubr)

#load data for each amplicon as phyloseq objects
metadata<-read.csv("Results/Dadaed-chimremoved-nonbactremoved/decontam-metadata.tsv", sep = "\t") %>% column_to_rownames(var = "sampleid")

load_physeq_region<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  
  physeq_region<-phyloseq(sample_data(metadata),
                        otu_table(table,taxa_are_rows = T),
                        tax_table(as.matrix(taxa)))
  return(physeq_region)
}

physeq_V1V2<-load_physeq_region("V1V2")
physeq_V2V3<-load_physeq_region("V2V3")
physeq_V3V4<-load_physeq_region("V3V4")
physeq_V4V5<-load_physeq_region("V4V5")
physeq_V5V7<-load_physeq_region("V5V7")
physeq_V7V9<-load_physeq_region("V7V9")

######plot library size vs sample type (control vs. true sample)
plot_lib_size<-function(ps, plottitle){
  df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  plot<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control))+geom_text(aes(label=rownames(df)),nudge_y = 2000) + 
    geom_point() + labs (title = plottitle) + guides(color=F)
return(plot)
}

ggarrange(
plot_lib_size(physeq_V1V2,"V1V2"),
plot_lib_size(physeq_V2V3,"V2V3"),
plot_lib_size(physeq_V3V4,"V3V4"),
plot_lib_size(physeq_V4V5,"V4V5"),
plot_lib_size(physeq_V5V7,"V5V7"),
plot_lib_size(physeq_V7V9,"V7V9"))

######identify.contaminants by frequency and prevalence methods (method = "combined") using decontam
sample_data(physeq_V1V2)$is.neg <- sample_data(physeq_V1V2)$Sample_or_Control == "Control Sample"
contamdf.V1V2 <- merge(isContaminant(physeq_V1V2, method="combined",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V1V2)),by=0)

sample_data(physeq_V2V3)$is.neg <- sample_data(physeq_V2V3)$Sample_or_Control == "Control Sample"
contamdf.V2V3 <- merge(isContaminant(physeq_V2V3, method="combined",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V2V3)),by=0)

sample_data(physeq_V3V4)$is.neg <- sample_data(physeq_V3V4)$Sample_or_Control == "Control Sample"
contamdf.V3V4 <- merge(isContaminant(physeq_V3V4, method="combined",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V3V4)),by=0)

sample_data(physeq_V4V5)$is.neg <- sample_data(physeq_V4V5)$Sample_or_Control == "Control Sample"
contamdf.V4V5 <- merge(isContaminant(physeq_V4V5, method="both",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V4V5)),by=0)

sample_data(physeq_V5V7)$is.neg <- sample_data(physeq_V5V7)$Sample_or_Control == "Control Sample"
contamdf.V5V7 <- merge(isContaminant(physeq_V5V7, method="combined",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V5V7)),by=0)

sample_data(physeq_V7V9)$is.neg <- sample_data(physeq_V7V9)$Sample_or_Control == "Control Sample"
contamdf.V7V9 <- merge(isContaminant(physeq_V7V9, method="combined",conc="quant",neg="is.neg",threshold=0.25),
                       as.data.frame(tax_table(physeq_V7V9)),by=0)

#merge results
contamdf.all<-do.call("rbind",list(
subset(contamdf.V1V2,contaminant==T) %>% mutate(region="V1V2"),
subset(contamdf.V2V3,contaminant==T) %>% mutate(region="V2V3"),
subset(contamdf.V3V4,contaminant==T) %>% mutate(region="V3V4"),
subset(contamdf.V4V5,contaminant==T) %>% mutate(region="V4V5"),
subset(contamdf.V5V7,contaminant==T) %>% mutate(region="V5V7"),
subset(contamdf.V7V9,contaminant==T) %>% mutate(region="V7V9")))

#visualize unique results
unique(contamdf.all$Genus)
