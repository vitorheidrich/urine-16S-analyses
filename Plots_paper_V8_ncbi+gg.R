library(ggplot2)
library(RColorBrewer)
library(plyr)
library(readr)
library(btools)
library(qiime2R)
library(QsRutils)
library(gridExtra)
library(gplots)
library(pals)
library(metagMisc)
library(beanplot)
library(picante)
library(ggpubr)
library(ComplexHeatmap)
library(gridGraphics)
library(cowplot)
library(UpSetR)
library(facetscales)
library(forcats)
library(ggforce)
library(vegan)
library(maditr)
library(SRS)
library(phyloseq)
library(ggbeeswarm)
library(tidyverse)
library(rstatix)
library(metacoder)
source("plotDistances.R")
library(DECIPHER)
library(Bios2cor)
library(ggpattern)
theme_set(theme_bw())

#function to replace generic proxies in taxonomy to NA
dropna_taxa_names<-function(physeq){
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))][tax_table(physeq)[, colnames(tax_table(physeq))] == ""] <- NA 
  return(physeq)
}

#loading amplicon data as phyloseq objects
metadata<-read.csv("Results/Dadaed-chimremoved-nonbactremoved/decontam-metadata.tsv", sep = "\t") %>% column_to_rownames(var = "sampleid")

###GREENGENES###
load_physeq_region_gg<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/gg_",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  tree<-read_qza(paste0("Results/Tree/",region,"_rooted-tree.qza"))[["data"]]
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)),
                          tree=tree)
  return(physeq_region)
}

physeq_V1V2_gg<-load_physeq_region_gg("V1V2")
physeq_V2V3_gg<-load_physeq_region_gg("V2V3")
physeq_V3V4_gg<-load_physeq_region_gg("V3V4")
physeq_V4V5_gg<-load_physeq_region_gg("V4V5")
physeq_V5V7_gg<-load_physeq_region_gg("V5V7")
physeq_V7V9_gg<-load_physeq_region_gg("V7V9")

#removing libraries 3, 4, 5 (low-quality). ps: different names in the publication
physeq_V1V2_gg<-prune_samples(!(sample_names(physeq_V1V2_gg) %in% c("3","4","5")), physeq_V1V2_gg)
physeq_V2V3_gg<-prune_samples(!(sample_names(physeq_V2V3_gg) %in% c("3","4","5")), physeq_V2V3_gg)
physeq_V3V4_gg<-prune_samples(!(sample_names(physeq_V3V4_gg) %in% c("3","4","5")), physeq_V3V4_gg)
physeq_V4V5_gg<-prune_samples(!(sample_names(physeq_V4V5_gg) %in% c("3","4","5")), physeq_V4V5_gg)
physeq_V5V7_gg<-prune_samples(!(sample_names(physeq_V5V7_gg) %in% c("3","4","5")), physeq_V5V7_gg)
physeq_V7V9_gg<-prune_samples(!(sample_names(physeq_V7V9_gg) %in% c("3","4","5")), physeq_V7V9_gg)

#removing ASVs which had non-zero counts only in the removed libraries
physeq_V1V2_gg<-prune_taxa(taxa_sums(physeq_V1V2_gg) > 0, physeq_V1V2_gg) 
physeq_V2V3_gg<-prune_taxa(taxa_sums(physeq_V2V3_gg) > 0, physeq_V2V3_gg) 
physeq_V3V4_gg<-prune_taxa(taxa_sums(physeq_V3V4_gg) > 0, physeq_V3V4_gg) 
physeq_V4V5_gg<-prune_taxa(taxa_sums(physeq_V4V5_gg) > 0, physeq_V4V5_gg) 
physeq_V5V7_gg<-prune_taxa(taxa_sums(physeq_V5V7_gg) > 0, physeq_V5V7_gg) 
physeq_V7V9_gg<-prune_taxa(taxa_sums(physeq_V7V9_gg) > 0, physeq_V7V9_gg)

####Figure S3A####
#20/02/22 - % of assigned sequences per level per dataset
assigned_seqs<-function(physeq,rank){
  physeq_treated<-merge_samples(tax_glom(dropna_taxa_names(physeq), taxrank = rank, NArm = F), group = "Sample_or_Control")
  data<-merge(as.data.frame(tax_table(physeq_treated)),t(as.data.frame(otu_table(physeq_treated))), by = 0)

  result<-100*(sum(subset(data, is.na(data[[rank]])==F)[["True Sample"]])/sum(data[["True Sample"]]))
  return(result)
}


assigned_seqs_perregion_perlevel_gg<-data.frame(region=c(rep("V1V2",6),rep("V2V3",6),rep("V3V4",6),rep("V4V5",6),rep("V5V7",6),rep("V7V9",6)),
                                                rank=rep(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),6),
                                                assigned=c(assigned_seqs(physeq_V1V2_gg, rank = "Phylum"), assigned_seqs(physeq_V1V2_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V1V2_gg, rank = "Order"),assigned_seqs(physeq_V1V2_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V1V2_gg, rank = "Genus"),assigned_seqs(physeq_V1V2_gg, rank = "Species"),
                                                           assigned_seqs(physeq_V2V3_gg, rank = "Phylum"), assigned_seqs(physeq_V2V3_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V2V3_gg, rank = "Order"),assigned_seqs(physeq_V2V3_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V2V3_gg, rank = "Genus"),assigned_seqs(physeq_V2V3_gg, rank = "Species"),
                                                           assigned_seqs(physeq_V3V4_gg, rank = "Phylum"), assigned_seqs(physeq_V3V4_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V3V4_gg, rank = "Order"),assigned_seqs(physeq_V3V4_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V3V4_gg, rank = "Genus"),assigned_seqs(physeq_V3V4_gg, rank = "Species"),
                                                           assigned_seqs(physeq_V4V5_gg, rank = "Phylum"), assigned_seqs(physeq_V4V5_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V4V5_gg, rank = "Order"),assigned_seqs(physeq_V4V5_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V4V5_gg, rank = "Genus"),assigned_seqs(physeq_V4V5_gg, rank = "Species"),
                                                           assigned_seqs(physeq_V5V7_gg, rank = "Phylum"), assigned_seqs(physeq_V5V7_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V5V7_gg, rank = "Order"),assigned_seqs(physeq_V5V7_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V5V7_gg, rank = "Genus"),assigned_seqs(physeq_V5V7_gg, rank = "Species"),
                                                           assigned_seqs(physeq_V7V9_gg, rank = "Phylum"), assigned_seqs(physeq_V7V9_gg, rank = "Class"),
                                                           assigned_seqs(physeq_V7V9_gg, rank = "Order"),assigned_seqs(physeq_V7V9_gg, rank = "Family"),
                                                           assigned_seqs(physeq_V7V9_gg, rank = "Genus"),assigned_seqs(physeq_V7V9_gg, rank = "Species")))

assigned_seqs_perregion_perlevel_gg$rank = factor(assigned_seqs_perregion_perlevel_gg$rank, 
                                                  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

FS3A<-ggplot(assigned_seqs_perregion_perlevel_gg,aes(x=region,y=assigned, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+scale_y_continuous(breaks=c(0,25,50,75,100),limits=c(0,100))+
  labs(y="% of assigned sequences")+facet_wrap(.~rank, nrow = 1)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###NCBI###
load_physeq_region_ncbi<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/ncbi_",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  tree<-read_qza(paste0("Results/Tree/",region,"_rooted-tree.qza"))[["data"]]
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)),
                          tree=tree)
  return(physeq_region)
}

physeq_V1V2_ncbi<-load_physeq_region_ncbi("V1V2")
physeq_V2V3_ncbi<-load_physeq_region_ncbi("V2V3")
physeq_V3V4_ncbi<-load_physeq_region_ncbi("V3V4")
physeq_V4V5_ncbi<-load_physeq_region_ncbi("V4V5")
physeq_V5V7_ncbi<-load_physeq_region_ncbi("V5V7")
physeq_V7V9_ncbi<-load_physeq_region_ncbi("V7V9")

#removing libraries 3, 4, 5 (low-quality). ps: different names in the publication
physeq_V1V2_ncbi<-prune_samples(!(sample_names(physeq_V1V2_ncbi) %in% c("3","4","5")), physeq_V1V2_ncbi)
physeq_V2V3_ncbi<-prune_samples(!(sample_names(physeq_V2V3_ncbi) %in% c("3","4","5")), physeq_V2V3_ncbi)
physeq_V3V4_ncbi<-prune_samples(!(sample_names(physeq_V3V4_ncbi) %in% c("3","4","5")), physeq_V3V4_ncbi)
physeq_V4V5_ncbi<-prune_samples(!(sample_names(physeq_V4V5_ncbi) %in% c("3","4","5")), physeq_V4V5_ncbi)
physeq_V5V7_ncbi<-prune_samples(!(sample_names(physeq_V5V7_ncbi) %in% c("3","4","5")), physeq_V5V7_ncbi)
physeq_V7V9_ncbi<-prune_samples(!(sample_names(physeq_V7V9_ncbi) %in% c("3","4","5")), physeq_V7V9_ncbi)

#removing ASVs which had non-zero counts only in the removed libraries
physeq_V1V2_ncbi<-prune_taxa(taxa_sums(physeq_V1V2_ncbi) > 0, physeq_V1V2_ncbi) 
physeq_V2V3_ncbi<-prune_taxa(taxa_sums(physeq_V2V3_ncbi) > 0, physeq_V2V3_ncbi) 
physeq_V3V4_ncbi<-prune_taxa(taxa_sums(physeq_V3V4_ncbi) > 0, physeq_V3V4_ncbi) 
physeq_V4V5_ncbi<-prune_taxa(taxa_sums(physeq_V4V5_ncbi) > 0, physeq_V4V5_ncbi) 
physeq_V5V7_ncbi<-prune_taxa(taxa_sums(physeq_V5V7_ncbi) > 0, physeq_V5V7_ncbi) 
physeq_V7V9_ncbi<-prune_taxa(taxa_sums(physeq_V7V9_ncbi) > 0, physeq_V7V9_ncbi)

assigned_seqs_perregion_perlevel_ncbi<-data.frame(region=c(rep("V1V2",6),rep("V2V3",6),rep("V3V4",6),rep("V4V5",6),rep("V5V7",6),rep("V7V9",6)),
                                                rank=rep(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),6),
                                                assigned=c(assigned_seqs(physeq_V1V2_ncbi, rank = "Phylum"), assigned_seqs(physeq_V1V2_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V1V2_ncbi, rank = "Order"),assigned_seqs(physeq_V1V2_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V1V2_ncbi, rank = "Genus"),assigned_seqs(physeq_V1V2_ncbi, rank = "Species"),
                                                           assigned_seqs(physeq_V2V3_ncbi, rank = "Phylum"), assigned_seqs(physeq_V2V3_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V2V3_ncbi, rank = "Order"),assigned_seqs(physeq_V2V3_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V2V3_ncbi, rank = "Genus"),assigned_seqs(physeq_V2V3_ncbi, rank = "Species"),
                                                           assigned_seqs(physeq_V3V4_ncbi, rank = "Phylum"), assigned_seqs(physeq_V3V4_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V3V4_ncbi, rank = "Order"),assigned_seqs(physeq_V3V4_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V3V4_ncbi, rank = "Genus"),assigned_seqs(physeq_V3V4_ncbi, rank = "Species"),
                                                           assigned_seqs(physeq_V4V5_ncbi, rank = "Phylum"), assigned_seqs(physeq_V4V5_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V4V5_ncbi, rank = "Order"),assigned_seqs(physeq_V4V5_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V4V5_ncbi, rank = "Genus"),assigned_seqs(physeq_V4V5_ncbi, rank = "Species"),
                                                           assigned_seqs(physeq_V5V7_ncbi, rank = "Phylum"), assigned_seqs(physeq_V5V7_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V5V7_ncbi, rank = "Order"),assigned_seqs(physeq_V5V7_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V5V7_ncbi, rank = "Genus"),assigned_seqs(physeq_V5V7_ncbi, rank = "Species"),
                                                           assigned_seqs(physeq_V7V9_ncbi, rank = "Phylum"), assigned_seqs(physeq_V7V9_ncbi, rank = "Class"),
                                                           assigned_seqs(physeq_V7V9_ncbi, rank = "Order"),assigned_seqs(physeq_V7V9_ncbi, rank = "Family"),
                                                           assigned_seqs(physeq_V7V9_ncbi, rank = "Genus"),assigned_seqs(physeq_V7V9_ncbi, rank = "Species")))

assigned_seqs_perregion_perlevel_ncbi$rank = factor(assigned_seqs_perregion_perlevel_ncbi$rank, 
                                                  levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

FS3B<-ggplot(assigned_seqs_perregion_perlevel_ncbi,aes(x=region,y=assigned, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+scale_y_continuous(breaks=c(0,25,50,75,100),limits=c(0,100))+
  labs(y="% of assigned sequences")+facet_wrap(.~rank, nrow = 1)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


################################Arranging figures############################
####Figure 1####
ggarrange(annotate_figure(FS3A+labs(title="Greengenes"),fig.lab = "A",fig.lab.size = 16),
          annotate_figure(FS3B+labs(title="NCBI 16S RefSeq"),fig.lab = "B",fig.lab.size = 16),
          nrow=2, heights = c(1,1))
