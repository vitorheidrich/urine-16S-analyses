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

sample_list<-c('#3','#4','#5','#6','#7','#9','#10','#12','#13','#14','#16','#19','#21','#23','#25','#37','#39','#42')
relabeling_samples<-c('#3'='#1','#4'='#2','#5'='#3','#6'='#4','#7'='#5','#9'='#6','#10'='#7.1',
                      '#12'='#7.2','#13'='#7.3','#14'='#8.1','#16'='#8.2','#19'='#9.1','#21'='#9.2',
                      '#23'='#10.1','#25'='#10.2','#37'='#13.1','#39'='#13.2','#42'='#14') #relabeling samples to match the convention of the paper
sample_labeller <- function(variable,value){
  return(tail(relabeling_samples,n=15)[value])
}

#function to replace generic proxies in taxonomy to NA
dropna_taxa_names<-function(physeq){
  if ("2ddfe48ea4c5e20e648355a2fe987839" %in% rownames(tax_table(physeq))){tax_table(physeq)["2ddfe48ea4c5e20e648355a2fe987839", "Genus"]<-"Allisonella"}
  if ("baedbf0feaa251e4838b99bc66cf86f8" %in% rownames(tax_table(physeq))){tax_table(physeq)["baedbf0feaa251e4838b99bc66cf86f8", "Genus"]<-"Globicatella"}
  if ("d92d64056b1e18ee026378b57a977bfb" %in% rownames(tax_table(physeq))){tax_table(physeq)["d92d64056b1e18ee026378b57a977bfb", "Genus"]<-"Globicatella"}
   
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "d__", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "uncultured", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "metagenome", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "human_gut", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "_sp.", replacement = NA)
  
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  return(physeq)
}

#function to make taxonomic nomenclature homogeneous, preventing generic proxies and empty entries
coherent_taxa_names<-function(physeq){
  if ("2ddfe48ea4c5e20e648355a2fe987839" %in% rownames(tax_table(physeq))){tax_table(physeq)["2ddfe48ea4c5e20e648355a2fe987839", "Genus"]<-"Allisonella"}
  if ("baedbf0feaa251e4838b99bc66cf86f8" %in% rownames(tax_table(physeq))){tax_table(physeq)["baedbf0feaa251e4838b99bc66cf86f8", "Genus"]<-"Globicatella"}
  if ("d92d64056b1e18ee026378b57a977bfb" %in% rownames(tax_table(physeq))){tax_table(physeq)["d92d64056b1e18ee026378b57a977bfb", "Genus"]<-"Globicatella"}
  
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "d__", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "uncultured", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "metagenome", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "human_gut", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "_sp.", replacement = NA)
  
  tax.clean <- data.frame(tax_table(physeq))
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
      tax.clean[i, 2:7] <- kingdom
    } else if (tax.clean[i,3] == ""){
      phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
      tax.clean[i, 3:7] <- phylum
    } else if (tax.clean[i,4] == ""){
      class <- paste0("Class_", tax.clean[i,3], sep = "")
      tax.clean[i, 4:7] <- class
    } else if (tax.clean[i,5] == ""){
      order <- paste("Order_", tax.clean[i,4], sep = "")
      tax.clean[i, 5:7] <- order
    } else if (tax.clean[i,6] == ""){
      family <- paste("Family_", tax.clean[i,5], sep = "")
      tax.clean[i, 6:7] <- family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
    }
  }
  tax_table(physeq) <- as.matrix(tax.clean)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  return(physeq)
}

#function to perform SRS normalization
SRS_normalize<-function(physeq){ #Cmin = minimum sequencing depth
  physeq_copy<-physeq
  SRS_output <- SRS(data = as.data.frame(otu_table(physeq_copy)), Cmin = min(sample_sums(physeq))) #running SRS
  rownames(SRS_output)<-rownames(as.data.frame(otu_table(physeq_copy))) #reassigning features names as rownames
  otu_table(physeq_copy)<-otu_table(SRS_output,taxa_are_rows = T)
  return(physeq_copy)
}

#loading amplicon data as phyloseq objects
metadata<-read.csv("Results/Dadaed-chimremoved-nonbactremoved/decontam-metadata.tsv", sep = "\t") %>% column_to_rownames(var = "sampleid")

load_physeq_region<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  tree<-read_qza(paste0("Results/Tree/",region,"_rooted-tree.qza"))[["data"]]
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)),
                          tree=tree)
  return(physeq_region)
}

physeq_V1V2<-load_physeq_region("V1V2")
physeq_V2V3<-load_physeq_region("V2V3")
physeq_V3V4<-load_physeq_region("V3V4")
physeq_V4V5<-load_physeq_region("V4V5")
physeq_V5V7<-load_physeq_region("V5V7")
physeq_V7V9<-load_physeq_region("V7V9")

#removing libraries 3, 4, 5 (low-quality). ps: different names in the publication
physeq_V1V2<-prune_samples(!(sample_names(physeq_V1V2) %in% c("3","4","5")), physeq_V1V2)
physeq_V2V3<-prune_samples(!(sample_names(physeq_V2V3) %in% c("3","4","5")), physeq_V2V3)
physeq_V3V4<-prune_samples(!(sample_names(physeq_V3V4) %in% c("3","4","5")), physeq_V3V4)
physeq_V4V5<-prune_samples(!(sample_names(physeq_V4V5) %in% c("3","4","5")), physeq_V4V5)
physeq_V5V7<-prune_samples(!(sample_names(physeq_V5V7) %in% c("3","4","5")), physeq_V5V7)
physeq_V7V9<-prune_samples(!(sample_names(physeq_V7V9) %in% c("3","4","5")), physeq_V7V9)

#removing ASVs which had non-zero counts only in the removed libraries
physeq_V1V2<-prune_taxa(taxa_sums(physeq_V1V2) > 0, physeq_V1V2) 
physeq_V2V3<-prune_taxa(taxa_sums(physeq_V2V3) > 0, physeq_V2V3) 
physeq_V3V4<-prune_taxa(taxa_sums(physeq_V3V4) > 0, physeq_V3V4) 
physeq_V4V5<-prune_taxa(taxa_sums(physeq_V4V5) > 0, physeq_V4V5) 
physeq_V5V7<-prune_taxa(taxa_sums(physeq_V5V7) > 0, physeq_V5V7) 
physeq_V7V9<-prune_taxa(taxa_sums(physeq_V7V9) > 0, physeq_V7V9)

####Figure S1####
preprocess_length<-function(region_R){
  df<-read.csv(paste0("Results/Read_length/",region_R,"_read_length.txt"), sep=" ", row.names = NULL)[,-1]
  df$freq<-apply(df, MARGIN = 1, function(x) {x[!is.na(x)][1]})
  df$size<-apply(df, MARGIN = 1, function(x) {x[!is.na(x)][2]})
  df<-df[c("freq","size")]
  df<-aggregate(freq~size, data = df, FUN=sum)
  df$region_R<-region_R
  return(df)
}

lengths<-do.call("rbind",list(preprocess_length("V1V2_R1"), preprocess_length("V1V2_R2"),
               preprocess_length("V2V3_R1"), preprocess_length("V2V3_R2"),
               preprocess_length("V3V4_R1"), preprocess_length("V3V4_R2"),
               preprocess_length("V4V5_R1"), preprocess_length("V4V5_R2"),
               preprocess_length("V5V7_R1"), preprocess_length("V5V7_R2"),
               preprocess_length("V7V9_R1"), preprocess_length("V7V9_R2"),
               #preprocess_length("ITS_R1"), preprocess_length("ITS_R2"),
               preprocess_length("unknown_R1"), preprocess_length("unknown_R2")))

lengths$region_R<- str_replace(lengths$region_R,"_"," - ")

FS1<-ggplot(lengths, aes(x=size,y=freq))+geom_col(color = "black", fill = "black")+labs(x="Read length in bp",y="Frequency") +
  facet_wrap(.~region_R,nrow = 8)+scale_x_continuous(breaks = c(1,50,100,150,200,250,276), limits = c(0,280))+
  scale_y_continuous(breaks = c(0, 500000, 1000000), limits=c(0,1000000))+
  theme(panel.grid = element_blank())#, axis.text.y = element_blank(), axis.ticks.y = element_blank())
FS1

####Figure 1A####
#02/06/21 number of reads generated for each amplicon
preprocess_nreads<-function(region){
  df<-read.csv(paste0("Results/Amplicons_demultiplexed/",region,"_per-sample-fastq-counts.tsv"), sep ="\t")[,1:2]
  colnames(df)<-c("sampleid","count")
  df$region<-region
  df$sampleid<-as.factor(df$sampleid)
  return(df)
}

nreads<-do.call("rbind", list(preprocess_nreads("V1V2"),preprocess_nreads("V2V3"),preprocess_nreads("V3V4"),preprocess_nreads("V4V5"),
              preprocess_nreads("V5V7"),preprocess_nreads("V7V9"),preprocess_nreads("ITS"),preprocess_nreads("unknown")))

#median read number per library
median(as.data.frame(subset(nreads, region!="ITS" & sampleid!='56' & sampleid!='57') %>% group_by(sampleid) %>% summarise(sum=sum(count)) %>% arrange(sum))$sum)

#not considering its from now on  
nreads_totais_uteis_summary<-data.frame(region = c(#'ITS',
                                                          'V1V2','V2V3','V3V4','V4V5','V5V7','V7V9'),
                                   total=c(#14275,
                                                  750924,801072,668509,1674525,881612,749534),
                                   after_filtering=c(#NA,
                                                        614973,624939,488560,477043,766525,634974))
reshape2::melt(nreads_totais_uteis_summary)

F1A<-ggplot(reshape2::melt(nreads_totais_uteis_summary), aes(x=region,y=value,fill=region,color=variable))+
  geom_bar_pattern(stat = "identity", 
                   pattern = c(#"stripe", 
                                        "stripe", "stripe", "stripe","stripe","stripe","stripe", # 1st col
                               #"none", 
                                        "none", "none", "none", "none", "none","none" # 2nd col
                   ),
                   pattern_angle = c(rep(45, 6), rep(0, 6)), #6 instead of 7
                   pattern_density = .1,
                   pattern_spacing = .05,
                   pattern_fill = 'white',
                   aes(fill = region),
                   position = "dodge")+ 
  labs(y="Number of reads")+
  scale_y_continuous(breaks=c(0,500000,1000000,1500000), limits = c(0,1674525),
                     labels = function(x) format(x, scientific = TRUE),expand = c(0.01,0))+
  scale_fill_manual(values = c(#"gray55", 
                                          brewer.pal(6, "Dark2")))+
  scale_color_manual(values=c("white","white"), labels = c("Total", "After filtering"))+
  guides(color = guide_legend(override.aes =  list(fill = "white",color = "black",pattern = c("stripe", "none"), pattern_spacing = .05, pattern_angle = c(45, 0))), fill = "none")+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(), axis.title.x = element_blank(), legend.title = element_blank())
F1A

####Figure S3####
#10/06/21 histogram with ASV lengths in bp per amplicon
preprocess_length_asv<-function(region){
  df<-read.csv(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_ASV_lengths.txt"), header = F)
  df$region<-region
  colnames(df)<-c("size","region")
  return(df)
}

lengths_asv<-do.call("rbind",list(preprocess_length_asv("V1V2"),
                              preprocess_length_asv("V2V3"),
                              preprocess_length_asv("V3V4"),
                              preprocess_length_asv("V4V5"),
                              preprocess_length_asv("V5V7"),
                              preprocess_length_asv("V7V9")))

lengths_asv_melt<-reshape2::melt(table(lengths_asv),c("size"))
colnames(lengths_asv_melt)<-c("size","region","freq")

lengths_asv_median <- lengths_asv %>% group_by(region) %>%
  summarise(median_size = median(size))

FS3<-ggplot(lengths_asv_melt, aes(x=size,y=freq))+
  geom_vline(data = lengths_asv_median, aes(xintercept = median_size), size = 0.6, color = "red", linetype="dotted")+
  geom_col(color = "black", fill = "black")+labs(x="ASV length in bp",y="Frequency") +
  facet_wrap(.~region, nrow = 6)+scale_x_continuous(breaks = seq(0,450,50),limits = c(0,450))+
  scale_y_continuous(limits = c(0,125), breaks = seq(0,120,60))+
  theme(panel.grid = element_blank())
FS3
####Figure 1B,1E####
#10/06/21 number of useful (good quality) reads per amplicon dataset
preprocess_nreads_uteis<-function(region){
  df<-read.csv(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_sample-frequency-detail.csv"))
  colnames(df)<-c("sampleid","count")
  df$region<-region
  df$sampleid<-as.factor(df$sampleid)
  return(df)
}

nreads_uteis<-do.call("rbind", list(preprocess_nreads_uteis("V1V2"),preprocess_nreads_uteis("V2V3"),preprocess_nreads_uteis("V3V4"),preprocess_nreads_uteis("V4V5"),
                              preprocess_nreads_uteis("V5V7"),preprocess_nreads_uteis("V7V9")))

#barplot com total de reads gerados por região (todas as amostras)
# F1B<-ggplot(nreads_uteis%>%
#          group_by(region)%>%
#                    summarise(total=sum(count)), aes(x=region,y=total,fill=region))+geom_bar(stat="identity")+
#             scale_fill_manual(values = brewer.pal(6, "Dark2"))+guides(fill="none")+labs(y = "Reads after filtering")+
#   scale_y_continuous(breaks=c(0,500000,1000000,1500000), limits = c(0,1750000),labels = function(x) format(x, scientific = TRUE),expand = c(0.01,0))+
#   theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank())
# 
# F1B

#nreads_uteis$sampleid<-paste0("#",nreads_uteis$sampleid)

# F1E<-ggplot(nreads_uteis, aes(x=reorder(sampleid,-count),y=count,fill=region,group=region))+geom_bar(stat="identity")+
#   labs(x = "Sample", fill = "Amplicon", y = "Reads after filtering" )+scale_y_continuous(expand=c(0.01,0))+
#   scale_fill_manual(values = brewer.pal(6, "Dark2"))+theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank())
# 
# F1E

####Figure 1B,1C####
#calculating relative frequencies
nreads_uteis_unmelted<-dcast(data = nreads_uteis, formula = sampleid~region,fun.aggregate = sum,value.var = "count")
nreads_uteis_relfreq<-as.data.frame(cbind(as.character(nreads_uteis_unmelted$sampleid),t(as.data.frame(apply(nreads_uteis_unmelted[,-1],MARGIN = 1, function(x) {x/sum(x)})))))
colnames(nreads_uteis_relfreq)<-colnames(nreads_uteis_unmelted)
nreads_uteis_relfreq_melted<-reshape2::melt(nreads_uteis_relfreq, "sampleid")
colnames(nreads_uteis_relfreq_melted)<-c("sampleid","region","relfreq")
nreads_uteis_relfreq_melted$relfreq<-as.numeric(nreads_uteis_relfreq_melted$relfreq)
nreads_uteis_relfreq_melted$sampleid<-paste0("#",nreads_uteis_relfreq_melted$sampleid)

#lvls <- subset(nreads_uteis_relfreq_melted,region=="V5V7")[order(-subset(nreads_uteis_relfreq_melted,region=="V5V7")$relfreq),]$sampleid

F1C<-ggplot(nreads_uteis_relfreq_melted, aes(x=factor(sampleid, levels = sample_list),y=relfreq,fill=region))+geom_bar(stat="identity")+
  labs(x = "Sample", y = "Relative frequency", fill = "Amplicon")+scale_y_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = brewer.pal(8, "Dark2"))+scale_x_discrete(labels = relabeling_samples)
F1C

#recalculating relative frequences merging all samples (whole dataset evaluation) 
nreads_uteis_unmelted<-dcast(data = nreads_uteis, formula = sampleid~region,fun.aggregate = sum,value.var = "count")
nreads_uteis_total<-as.data.frame(apply(nreads_uteis_unmelted[,-1], MARGIN=2, FUN=sum))
colnames(nreads_uteis_total)<-"count"
nreads_uteis_total$relfreq<-nreads_uteis_total$count/sum(nreads_uteis_total)

F1B<-ggplot(nreads_uteis_total, aes(x=0,y=relfreq,fill=rownames(nreads_uteis_total)))+
  geom_bar(stat="identity")+labs(y ="Relative frequency",fill="Amplicon", x = "All")+
  scale_fill_manual(values = brewer.pal(8, "Dark2"))+scale_y_continuous(expand = c(0.01,0))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position = "right")

F1B


####Figure S2####
#rarefaction curves
par(mfrow = c(3, 2))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V1V2)),
         step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V1V2", label = F, xlim = c(1,100000),ylim=c(0,85))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V2V3)),
         step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V2V3", label = F, xlim = c(1,100000),ylim=c(0,85))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V3V4)),
          step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V3V4", label = F, xlim = c(1,100000),ylim=c(0,85))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V4V5)),
          step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V4V5", label = F, xlim = c(1,100000),ylim=c(0,85))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V5V7)),
          step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V5V7", label = F, xlim = c(1,100000))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))
par(mar = c(5, 5, 5, 5))
rarecurve(t(otu_table(physeq_V7V9)),
          step=50, ylab = "ASVs", xlab = "Sequencing depth", main = "V7V9", label = F, xlim = c(1,100000))
axis(1,at = seq(0,100000,20000))
axis(2,at = seq(0,80,20))

####Figure 2A####
#15/06/21 - plot taxa richness for whole dataset per region

taxa_richness<-function(physeq, region){
  df<-data.frame(region = rep(region,7), rank = c(rank_names(physeq)[-1],"ASV"),
                 richness = c(phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[2], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[3], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[4], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[5], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[6], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[7], NArm = T)),
                              phyloseq::ntaxa(physeq)))
  return (df)
}

taxa_richness_perregion<-do.call("rbind",list(taxa_richness(dropna_taxa_names(physeq_V1V2),"V1V2"),taxa_richness(dropna_taxa_names(physeq_V2V3),"V2V3"),
                                            taxa_richness(dropna_taxa_names(physeq_V3V4),"V3V4"),taxa_richness(dropna_taxa_names(physeq_V4V5),"V4V5"),
                                            taxa_richness(dropna_taxa_names(physeq_V5V7),"V5V7"),taxa_richness(dropna_taxa_names(physeq_V7V9),"V7V9")))

taxa_richness_perregion$rank = factor(taxa_richness_perregion$rank, levels = c(rank_names(physeq_V1V2)[-1],"ASV"))

F2A<-ggplot(subset(taxa_richness_perregion,rank=="ASV"),aes(x=region,y=richness, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+scale_y_continuous(breaks=c(0,100,200),limits=c(0,200),expand=c(0.01,0))+
  facet_wrap(.~rank, nrow = 1)+theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                     legend.position = "none",axis.title.x = element_blank())+labs(y="Richness")
F2B<-ggplot(subset(taxa_richness_perregion,rank!="ASV"),aes(x=region,y=richness, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+
  scale_y_continuous(breaks=c(0,35,70),limits=c(0,75),expand=c(0.01,0))+labs(y="Richness")+
  facet_wrap(.~rank, nrow = 1)+theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                        legend.position = "none",axis.title = element_blank())

F2A
F2B

####Figure 2C####
#08/09/21 - sequence variability (entropy) per region
calc_entropy<-function(region, median_asv_size){
  fa<-readDNAStringSet(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_dna-sequences.fasta"))
  fa.al<-AlignSeqs(fa)
  entropy<-entropy(fa.al, gap_ratio = 1)
  entropy.rollavg<-zoo::rollmean(entropy, k = 20)
  entropy.df<-as.data.frame(entropy.rollavg)
  colnames(entropy.df)<-"entropy"
  entropy.df$pos<-as.numeric(rownames(entropy.df))
  entropy.df$region<-region
  entropy.df<-subset(entropy.df,pos<=median_asv_size)
  return(entropy.df)
}

entropies<-do.call("rbind",list(calc_entropy("V1V2", 322),calc_entropy("V2V3", 380),calc_entropy("V3V4", 409),
                                calc_entropy("V4V5", 373),calc_entropy("V5V7", 369),calc_entropy("V7V9", 377)))
F2D<-ggplot(entropies, aes(x=pos, y=entropy))+geom_line()+labs(y= "Entropy", x = "Nucleotide position")+
  facet_wrap(.~region, nrow = 1)+scale_x_continuous(breaks=c(10,100,200,300,400),limits = c(10,400))+
  scale_y_continuous(expand = c(0.01,0),limits = c(0,0.4))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
F2D


####Figure 2B####
#12/08/21 - faith's pd of the whole dataset
trace(estimate_richness,edit = T) #tracing with...
# if( "FaithPD" %in% measures){
#   outlist <- c(outlist, list(FaithPD = t(picante::pd(samp = OTU, tree = phy_tree(physeq), include.root = F))[1,] ))
# }
asv_pd<-data.frame(region = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"),
                   fpd=c(estimate_richness(merge_samples(physeq_V1V2, group = "Sample_or_Control"),measures="FaithPD")[,1],
                         estimate_richness(merge_samples(physeq_V2V3, group = "Sample_or_Control"),measures="FaithPD")[,1],
                         estimate_richness(merge_samples(physeq_V3V4, group = "Sample_or_Control"),measures="FaithPD")[,1],
                         estimate_richness(merge_samples(physeq_V4V5, group = "Sample_or_Control"),measures="FaithPD")[,1],
                         estimate_richness(merge_samples(physeq_V5V7, group = "Sample_or_Control"),measures="FaithPD")[,1],
                         estimate_richness(merge_samples(physeq_V7V9, group = "Sample_or_Control"),measures="FaithPD")[,1]))

F2C<-ggplot(asv_pd,aes(x=region,y=fpd, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+scale_y_continuous(breaks=c(0,5,10,15),limits=c(0,15),expand = c(0.01,0))+
  labs(y="Faith's PD")+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",axis.title.x = element_blank())

F2C


#correlation between median asv size and fpd
cor.test(asv_pd$fpd,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.9194; rho=-0.08571429 

#correlation between median asv and taxonomic richness
cor.test(subset(taxa_richness_perregion,rank=="Phylum")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.1583; rho=0.6546537
cor.test(subset(taxa_richness_perregion,rank=="Class")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.1583; rho=0.6546537 
cor.test(subset(taxa_richness_perregion,rank=="Order")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.4247; rho=-0.4058397
cor.test(subset(taxa_richness_perregion,rank=="Family")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.6175; rho=-0.260897
cor.test(subset(taxa_richness_perregion,rank=="Genus")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.6583; rho=-0.2571429
cor.test(subset(taxa_richness_perregion,rank=="Species")$richness,c(322,380,409,373,369,377), method = "spearman") #p-value = 0.2307; rho=-0.5768179

####Figure 1F####
#16/08/21 - % of assigned sequences per level per dataset
assigned_seqs<-function(physeq,rank){
  physeq_treated<-merge_samples(tax_glom(coherent_taxa_names(physeq), taxrank = rank), group = "Sample_or_Control")
  data<-merge(as.data.frame(tax_table(physeq_treated)),t(as.data.frame(otu_table(physeq_treated))), by = 0)
  
  data <- data %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Phylum_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Class_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Order_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Family_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Genus_.*", "none"))
  
  result<-100*(sum(subset(data, data[[rank]]!="none")[["True Sample"]])/sum(data[["True Sample"]]))
  return(result)
}

assigned_seqs_perregion_perlevel<-data.frame(region=c(rep("V1V2",6),rep("V2V3",6),rep("V3V4",6),rep("V4V5",6),rep("V5V7",6),rep("V7V9",6)),
                                             rank=rep(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),6),
                                             assigned=c(assigned_seqs(physeq_V1V2, rank = "Phylum"), assigned_seqs(physeq_V1V2, rank = "Class"),
                                                        assigned_seqs(physeq_V1V2, rank = "Order"),assigned_seqs(physeq_V1V2, rank = "Family"),
                                                        assigned_seqs(physeq_V1V2, rank = "Genus"),assigned_seqs(physeq_V1V2, rank = "Species"),
                                                        assigned_seqs(physeq_V2V3, rank = "Phylum"), assigned_seqs(physeq_V2V3, rank = "Class"),
                                                        assigned_seqs(physeq_V2V3, rank = "Order"),assigned_seqs(physeq_V2V3, rank = "Family"),
                                                        assigned_seqs(physeq_V2V3, rank = "Genus"),assigned_seqs(physeq_V2V3, rank = "Species"),
                                                        assigned_seqs(physeq_V3V4, rank = "Phylum"), assigned_seqs(physeq_V3V4, rank = "Class"),
                                                        assigned_seqs(physeq_V3V4, rank = "Order"),assigned_seqs(physeq_V3V4, rank = "Family"),
                                                        assigned_seqs(physeq_V3V4, rank = "Genus"),assigned_seqs(physeq_V3V4, rank = "Species"),
                                                        assigned_seqs(physeq_V4V5, rank = "Phylum"), assigned_seqs(physeq_V4V5, rank = "Class"),
                                                        assigned_seqs(physeq_V4V5, rank = "Order"),assigned_seqs(physeq_V4V5, rank = "Family"),
                                                        assigned_seqs(physeq_V4V5, rank = "Genus"),assigned_seqs(physeq_V4V5, rank = "Species"),
                                                        assigned_seqs(physeq_V5V7, rank = "Phylum"), assigned_seqs(physeq_V5V7, rank = "Class"),
                                                        assigned_seqs(physeq_V5V7, rank = "Order"),assigned_seqs(physeq_V5V7, rank = "Family"),
                                                        assigned_seqs(physeq_V5V7, rank = "Genus"),assigned_seqs(physeq_V5V7, rank = "Species"),
                                                        assigned_seqs(physeq_V7V9, rank = "Phylum"), assigned_seqs(physeq_V7V9, rank = "Class"),
                                                        assigned_seqs(physeq_V7V9, rank = "Order"),assigned_seqs(physeq_V7V9, rank = "Family"),
                                                        assigned_seqs(physeq_V7V9, rank = "Genus"),assigned_seqs(physeq_V7V9, rank = "Species")))

assigned_seqs_perregion_perlevel$rank = factor(assigned_seqs_perregion_perlevel$rank, 
                                               levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))

F1D<-ggplot(assigned_seqs_perregion_perlevel,aes(x=region,y=assigned, fill = region))+
  geom_col()+scale_fill_manual(values = brewer.pal(6, "Dark2"))+scale_y_continuous(breaks=c(0,25,50,75,100),limits=c(0,100))+
  labs(y="% of assigned sequences")+facet_wrap(.~rank, nrow = 1)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
F1D

####Figure S6####
#16/01/21 - relative abundances
  
plot_physeq_bar_all<-function(rank,min_ra=0,in_n_samples=1){
physeq_all_ra<-transform_sample_counts(coherent_taxa_names(tax_glom(merge_phyloseq(
  physeq_V1V2_tomerge,physeq_V2V3_tomerge,physeq_V3V4_tomerge,physeq_V4V5_tomerge,physeq_V5V7_tomerge,physeq_V7V9_tomerge),rank,NArm = F)),
  function(x) x/sum(x))
metadata<-data.frame(sample_data(physeq_all_ra))
metadata$sampleid<-factor(paste0("#",data.frame(sample_data(physeq_all_ra))$sampleid),levels = sample_list)
sample_data(physeq_all_ra)<-metadata
plot_bar(prune_taxa(genefilter_sample(physeq_all_ra, filterfun_sample(function(x) x >= min_ra), 
                                      A = in_n_samples),physeq_all_ra), "region", fill=rank)+facet_grid(.~sampleid, labeller = sample_labeller)+
  labs(y="Relative abundance")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values = glasbey())+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1),legend.text = element_text(size=8.5), legend.key.size = unit(0.4, 'cm'))
}

FS6<-ggarrange(
plot_physeq_bar_all("Phylum")+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=2,byrow = T)),
plot_physeq_bar_all("Class")+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=2,byrow = T)),
plot_physeq_bar_all("Order",0.01,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
plot_physeq_bar_all("Family",0.0724,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
plot_physeq_bar_all("Genus",0.139,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=5,byrow = T)),
plot_physeq_bar_all("Species",0.18,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=5,byrow = T)), 
nrow=6,heights = c(0.8,0.8,0.95,0.95,1.1,1.1)) #export 1440x2160
FS6

####Figure S5####
#18/06/21 - plot trees  (taxonomic)
plot_taxonomic_tree<- function(physeq,region){
set.seed(1)
heat_tree(parse_phyloseq(tax_glom(dropna_taxa_names(physeq),"Genus"), class_regex = "(.*)", class_key = "taxon_name"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs,
          layout = "reingold-tilford",
          title = region,
          title_size = 0.05,
          node_label_size_range = c(0.015, 0.05),
          node_size_range = c(0.01,0.05),
          edge_size_range = c(0.005,0.025),
          node_color_interval = c(1,75),
          edge_color_interval = c(1,75),
          node_color_axis_label = "observed genera")+theme_nothing()
}

ggarrange(#FS5
plot_taxonomic_tree(physeq_V1V2,"V1V2"),
plot_taxonomic_tree(physeq_V2V3,"V2V3"),
plot_taxonomic_tree(physeq_V3V4,"V3V4"),
plot_taxonomic_tree(physeq_V4V5,"V4V5"),
plot_taxonomic_tree(physeq_V5V7,"V5V7"),
plot_taxonomic_tree(physeq_V7V9,"V7V9"),ncol = 2, nrow = 3) #export 1600x2000

####Figure 3A####
#19/01/21 - intersections between identified genera
plot_intersect<-function(min_ra = 0, in_n_samples = 1){

observed_genera_filtered<-list(V1V2 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                               V2V3 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                               V3V4 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                               V4V5 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                               V5V7 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                               V7V9 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                  filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Genus",NArm = T), function(x) x/sum(x)))))$Genus))

p<-upset(fromList(observed_genera_filtered), sets = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"),
      order.by = "freq", nintersects = NA, sets.x.label = "Genus richness", mainbar.y.label = "Intersection size", 
      point.size = 2.5, line.size = 1.5, text.scale = c(1.3, 1.3, 1, 1, 1.3, 1), set_size.scale_max = 80, mainbar.y.max = 30, color.pal = "black")

#print(paste("n of genera:",nrow(fromList(observed_genera_filtered))))
return(p)

}
upset_plot_1<-plot_intersect()
F3A <- cowplot::plot_grid(NULL, NULL, upset_plot_1$Main_bar,
                          NULL, NULL, NULL,
                          upset_plot_1$Sizes, NULL, upset_plot_1$Matrix,
                          nrow=3, align='hv', rel_heights = c(2,-0.25,1),
                          rel_widths = c(1,-0.05,3))
F3A

####Figure S4B####
#21/01/21 - hierarchical clustering
fromList <- function(input){
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){x <- as.vector(match(elements, x))}))
  data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  rownames(data)<-elements
  return(data)
}

min_ra=0
in_n_samples=1

observed_genera<-list(V1V2 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                                 V2V3 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                                 V3V4 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                                 V4V5 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                                 V5V7 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Genus",NArm = T), function(x) x/sum(x)))))$Genus),
                                 V7V9 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Genus",NArm = T), function(x) x/sum(x)), 
                                                                                                    filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                                  transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Genus",NArm = T), function(x) x/sum(x)))))$Genus))
  
observed_genera<-as.data.frame(fromList(observed_genera))
observed_genera<-observed_genera[order(rowSums(observed_genera),decreasing=T),] #reordering dataframe based on sum of rows
observed_genera<-observed_genera[with(observed_genera, order(-V1V2, -V2V3, -V3V4, -V4V5, -V5V7, -V7V9)), ]
ordered_genera<-rownames(observed_genera)
observed_genera$taxa<-ordered_genera
observed_genera_melt<-reshape2::melt(observed_genera,"taxa")
observed_genera_melt$taxa<-factor(observed_genera_melt$taxa, levels=ordered_genera)

FS4B<-ggplot(observed_genera_melt, aes(y=variable,x=taxa))+geom_tile(aes(fill=as.factor(value)),color="black",size=1)+
  scale_fill_manual(values=c("lightyellow","firebrick3"),labels=c("0"="Undetected","1"="Detected"))+
  scale_y_discrete(expand = c(0,0))+labs(x = "Genus",fill = "")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1, size = 7.5),legend.position = "top")
FS4B

####Figure S4A####
observed_phyla<-list(V1V2 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V1V2),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum),
                     V2V3 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V2V3),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum),
                     V3V4 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V3V4),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum),
                     V4V5 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V4V5),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum),
                     V5V7 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V5V7),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum),
                     V7V9 = unique(as.data.frame(tax_table(prune_taxa(genefilter_sample(transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Phylum",NArm = T), function(x) x/sum(x)), 
                                                                                        filterfun_sample(function(x) x >= min_ra), A = in_n_samples),
                                                                      transform_sample_counts(tax_glom(dropna_taxa_names(physeq_V7V9),"Phylum",NArm = T), function(x) x/sum(x)))))$Phylum))
observed_phyla<-as.data.frame(fromList(observed_phyla))
observed_phyla<-observed_phyla[order(rowSums(observed_phyla),decreasing=T),] #reordering dataframe based on sum of rows
observed_phyla<-observed_phyla[with(observed_phyla, order(-V1V2, -V2V3, -V3V4, -V4V5, -V5V7, -V7V9)), ]
ordered_phyla<-rownames(observed_phyla)
observed_phyla$taxa<-ordered_phyla
observed_phyla_melt<-reshape2::melt(observed_phyla,"taxa")
observed_phyla_melt$taxa<-factor(observed_phyla_melt$taxa, levels=ordered_phyla)

FS4A<-ggplot(observed_phyla_melt, aes(y=variable,x=taxa))+geom_tile(aes(fill=as.factor(value)),color="black",size=1)+
  scale_fill_manual(values=c("lightyellow","firebrick3"),labels=c("0"="Undetected","1"="Detected"))+
  scale_y_discrete(expand = c(0,0))+labs(x = "Phylum",fill = "")+
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),legend.position = "top")

FS4A

####Figure 3C####
#03/08/21 - average relative abundance per amplicon
preprocess_physeq_tomerge<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)))
  physeq_region<-prune_samples(!(sample_names(physeq_region) %in% c("3","4","5")), physeq_region)
  physeq_region<-prune_taxa(taxa_sums(physeq_region) > 0, physeq_region) 
  
  physeq_pcoa<-physeq_region
  sample_names(physeq_pcoa)<-paste0(sample_names(physeq_region),"_",region)
  sample_data(physeq_pcoa)<-data.frame(sample_data(physeq_pcoa)) %>% mutate(region=region) %>% mutate(sampleid=sample_names(physeq_region))
  return(physeq_pcoa)
}

physeq_V1V2_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V1V2"))
physeq_V2V3_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V2V3"))
physeq_V3V4_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V3V4"))
physeq_V4V5_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V4V5"))
physeq_V5V7_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V5V7"))
physeq_V7V9_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V7V9"))

plot_physeq_bar_avg<-function(rank,min_ra=0,in_n_samples=1){
  physeq_all_ra<-transform_sample_counts(coherent_taxa_names(tax_glom(merge_phyloseq(
    physeq_V1V2_tomerge,physeq_V2V3_tomerge,physeq_V3V4_tomerge,physeq_V4V5_tomerge,physeq_V5V7_tomerge,physeq_V7V9_tomerge),rank,NArm = F)),
    function(x) x/sum(x))
  metadata<-data.frame(sample_data(physeq_all_ra))
  metadata$sampleid<-factor(data.frame(sample_data(physeq_all_ra))$sampleid,levels = c('6','7','9','10','12','13','14','16','19','21','23','25','37','39','42'))
  sample_data(physeq_all_ra)<-metadata
  physeq_avg<-phyloseq_average(physeq_all_ra,group = "region",avg_type = "arithmetic")
  physeq_avg_filtered<-prune_taxa(genefilter_sample(physeq_avg, filterfun_sample(function(x) x >= min_ra), 
                                                       A = in_n_samples),physeq_avg)
  
  plot_bar(physeq_avg_filtered, fill=rank)+
    labs(y="Relative abundance")+scale_y_continuous(expand=c(0,0),limits = c(0,1))+scale_fill_manual(values = glasbey())+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1))
}


F3C<-plot_physeq_bar_avg("Genus",0.01,1)+theme(axis.title.x = element_blank(),panel.grid = element_blank())+
  guides(fill=guide_legend(ncol=3,byrow = T))
F3C

####Figure 3B####
#10/08/21 - beta diversity boxplots (distance within and between samples)
plot_physeq_dist_all_wb<-function(rank,metric="bray", dont_consider_abundance = F){
  physeq_all_SRS<-SRS_normalize(coherent_taxa_names(tax_glom(merge_phyloseq(
    physeq_V1V2_tomerge,physeq_V2V3_tomerge,physeq_V3V4_tomerge,physeq_V4V5_tomerge,
    physeq_V5V7_tomerge,physeq_V7V9_tomerge), rank, NArm = F)))
  
  metadata<-as(sample_data(physeq_all_SRS),"data.frame")
  metadata$sampleid_r<-rownames(metadata)
  sample_data(physeq_all_SRS)<-metadata
  
  distances<-plotDistances(physeq_all_SRS,m=metric,plot=F,s="sampleid_r",d="sampleid", dont_consider_abundance = dont_consider_abundance)
  
  distances$group<-apply(distances, 1, function(x) ifelse(x[4]==x[5], "Within", "Between"))
  
  ggplot(distances, aes(x = group, y = value, fill = group))+geom_boxplot(color="black", outlier.alpha = 0.2)+
    labs(y = 'Dissimilarity (J)')+facet_wrap(eval(rank)~.)+
    scale_fill_manual(values = c("lightyellow","firebrick3"))+scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2))+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none",axis.title.x = element_blank())+
    stat_compare_means(label.y.npc = 0.95,label.x.npc = 0.5, method = "wilcox.test", label = "p.signif")
}

F3Ba<-ggarrange(
  plot_physeq_dist_all_wb("Phylum"),
  plot_physeq_dist_all_wb("Class")+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Order")+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Family")+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Genus")+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Species")+theme(axis.title.y = element_blank()),
  nrow=1,widths = c(1.1,1,1,1,1,1))

#MODIFY YLAB IN THE FUNCTION CALL BEFORE PLOTTING
F3Bb<-ggarrange(
  plot_physeq_dist_all_wb("Phylum",metric = "jaccard",dont_consider_abundance = T),
  plot_physeq_dist_all_wb("Class",metric = "jaccard",dont_consider_abundance = T)+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Order",metric = "jaccard",dont_consider_abundance = T)+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Family",metric = "jaccard",dont_consider_abundance = T)+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Genus",metric = "jaccard",dont_consider_abundance = T)+theme(axis.title.y = element_blank()),
  plot_physeq_dist_all_wb("Species",metric = "jaccard",dont_consider_abundance = T)+theme(axis.title.y = element_blank()),
  nrow=1,widths = c(1.1,1,1,1,1,1))

F3B<-ggarrange(F3Ba,F3Bb,nrow=2)
F3B


##############################SIDLE#########################################
#loading sidle datasets as phyloseq objects
load_physeq_region_sidle<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Sidle/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Sidle/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)))
  return(physeq_region)
}

physeq_V1V2_sidle<-load_physeq_region_sidle("V1V2")
physeq_V2V3_sidle<-load_physeq_region_sidle("V2V3")
physeq_V3V4_sidle<-load_physeq_region_sidle("V3V4")
physeq_V4V5_sidle<-load_physeq_region_sidle("V4V5")
physeq_V5V7_sidle<-load_physeq_region_sidle("V5V7")
physeq_V7V9_sidle<-load_physeq_region_sidle("V7V9")

physeq_V1V2V4V5_sidle<-load_physeq_region_sidle("V1V2-V4V5")
physeq_V1V2V5V7_sidle<-load_physeq_region_sidle("V1V2-V5V7")
physeq_V1V2V7V9_sidle<-load_physeq_region_sidle("V1V2-V7V9")
physeq_V2V3V5V7_sidle<-load_physeq_region_sidle("V2V3-V5V7")
physeq_V2V3V7V9_sidle<-load_physeq_region_sidle("V2V3-V7V9")
physeq_V3V4V7V9_sidle<-load_physeq_region_sidle("V3V4-V7V9")

physeq_full_sidle<-load_physeq_region_sidle("full")

#removing libraries 3, 4, 5. ps: other names in the publication
physeq_V1V2_sidle<-prune_samples(!(sample_names(physeq_V1V2_sidle) %in% c("3","4","5")), physeq_V1V2_sidle)
physeq_V2V3_sidle<-prune_samples(!(sample_names(physeq_V2V3_sidle) %in% c("3","4","5")), physeq_V2V3_sidle)
physeq_V3V4_sidle<-prune_samples(!(sample_names(physeq_V3V4_sidle) %in% c("3","4","5")), physeq_V3V4_sidle)
physeq_V4V5_sidle<-prune_samples(!(sample_names(physeq_V4V5_sidle) %in% c("3","4","5")), physeq_V4V5_sidle)
physeq_V5V7_sidle<-prune_samples(!(sample_names(physeq_V5V7_sidle) %in% c("3","4","5")), physeq_V5V7_sidle)
physeq_V7V9_sidle<-prune_samples(!(sample_names(physeq_V7V9_sidle) %in% c("3","4","5")), physeq_V7V9_sidle)

#removing ASVs which had non-zero counts only in the removed libraries
physeq_V1V2_sidle<-prune_taxa(taxa_sums(physeq_V1V2_sidle) > 0, physeq_V1V2_sidle) 
physeq_V2V3_sidle<-prune_taxa(taxa_sums(physeq_V2V3_sidle) > 0, physeq_V2V3_sidle) 
physeq_V3V4_sidle<-prune_taxa(taxa_sums(physeq_V3V4_sidle) > 0, physeq_V3V4_sidle) 
physeq_V4V5_sidle<-prune_taxa(taxa_sums(physeq_V4V5_sidle) > 0, physeq_V4V5_sidle) 
physeq_V5V7_sidle<-prune_taxa(taxa_sums(physeq_V5V7_sidle) > 0, physeq_V5V7_sidle) 
physeq_V7V9_sidle<-prune_taxa(taxa_sums(physeq_V7V9_sidle) > 0, physeq_V7V9_sidle)

#same
physeq_V1V2V4V5_sidle<-prune_samples(!(sample_names(physeq_V1V2V4V5_sidle) %in% c("3","4","5")), physeq_V1V2V4V5_sidle)
physeq_V1V2V5V7_sidle<-prune_samples(!(sample_names(physeq_V1V2V5V7_sidle) %in% c("3","4","5")), physeq_V1V2V5V7_sidle)
physeq_V1V2V7V9_sidle<-prune_samples(!(sample_names(physeq_V1V2V7V9_sidle) %in% c("3","4","5")), physeq_V1V2V7V9_sidle)
physeq_V2V3V5V7_sidle<-prune_samples(!(sample_names(physeq_V2V3V5V7_sidle) %in% c("3","4","5")), physeq_V2V3V5V7_sidle)
physeq_V2V3V7V9_sidle<-prune_samples(!(sample_names(physeq_V2V3V7V9_sidle) %in% c("3","4","5")), physeq_V2V3V7V9_sidle)
physeq_V3V4V7V9_sidle<-prune_samples(!(sample_names(physeq_V3V4V7V9_sidle) %in% c("3","4","5")), physeq_V3V4V7V9_sidle)

#same
physeq_V1V2V4V5_sidle<-prune_taxa(taxa_sums(physeq_V1V2V4V5_sidle) > 0, physeq_V1V2V4V5_sidle) 
physeq_V1V2V5V7_sidle<-prune_taxa(taxa_sums(physeq_V1V2V5V7_sidle) > 0, physeq_V1V2V5V7_sidle) 
physeq_V1V2V7V9_sidle<-prune_taxa(taxa_sums(physeq_V1V2V7V9_sidle) > 0, physeq_V1V2V7V9_sidle) 
physeq_V2V3V5V7_sidle<-prune_taxa(taxa_sums(physeq_V2V3V5V7_sidle) > 0, physeq_V2V3V5V7_sidle) 
physeq_V2V3V7V9_sidle<-prune_taxa(taxa_sums(physeq_V2V3V7V9_sidle) > 0, physeq_V2V3V7V9_sidle) 
physeq_V3V4V7V9_sidle<-prune_taxa(taxa_sums(physeq_V3V4V7V9_sidle) > 0, physeq_V3V4V7V9_sidle)

physeq_full_sidle<-prune_samples(!(sample_names(physeq_full_sidle) %in% c("3","4","5")), physeq_full_sidle)
physeq_full_sidle<-prune_taxa(taxa_sums(physeq_full_sidle) > 0, physeq_full_sidle) 

#function to replace generic proxies in taxonomy to NA (for sidle)
dropna_taxa_names_sidle<-function(physeq){
  if ("MH364377.1.1545" %in% rownames(tax_table(physeq))){tax_table(physeq)["MH364377.1.1545", "Genus"]<-"Allisonella"}
  if ("JN713272.1.1538" %in% rownames(tax_table(physeq))){tax_table(physeq)["JN713272.1.1538", "Genus"]<-"Globicatella"}
  
  #tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "uncultured", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "metagenome", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "human_gut", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "_sp.", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "\\|", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "d__Bacteria", replacement = "Bacteria")
  
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  return(physeq)
}

#function to make taxonomic nomenclature homogeneous, preventing generic proxies and empty entries (for sidle)
coherent_taxa_names_sidle<-function(physeq){
  if ("MH364377.1.1545" %in% rownames(tax_table(physeq))){tax_table(physeq)["MH364377.1.1545", "Genus"]<-"Allisonella"}
  if ("JN713272.1.1538" %in% rownames(tax_table(physeq))){tax_table(physeq)["JN713272.1.1538", "Genus"]<-"Globicatella"}
  
  #tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "uncultured", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "metagenome", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "human_gut", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "_sp.", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "\\|", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "d__Bacteria", replacement = "Bacteria")
  
  tax.clean <- data.frame(tax_table(physeq))
  
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
      tax.clean[i, 2:7] <- kingdom
    } else if (tax.clean[i,3] == ""){
      phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
      tax.clean[i, 3:7] <- phylum
    } else if (tax.clean[i,4] == ""){
      class <- paste0("Class_", tax.clean[i,3], sep = "")
      tax.clean[i, 4:7] <- class
    } else if (tax.clean[i,5] == ""){
      order <- paste("Order_", tax.clean[i,4], sep = "")
      tax.clean[i, 5:7] <- order
    } else if (tax.clean[i,6] == ""){
      family <- paste("Family_", tax.clean[i,5], sep = "")
      tax.clean[i, 6:7] <- family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
    }
  }
  tax_table(physeq) <- as.matrix(tax.clean)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  return(physeq)
}

#function to remove non-bacterial and with unindetified phylum from phyloseq objects
remove_unassigned<-function(physeq){
  return(subset_taxa(physeq, Species != "Kingdom_Bacteria" & Family != "Mitochondria" & Family!= "Chloroplast"))
}
#function to remove non-bacterial from phyloseq objects
dropna_unassigned<-function(physeq){
  return(subset_taxa(physeq, Family != "Mitochondria" & Family!= "Chloroplast"))
}

####Figure 4A####

#05/08/21 - sidle richness
taxa_richness_sidle<-function(physeq, region){
  df<-data.frame(region = rep(region,6), rank = c(rank_names(physeq)[-1]),
                 richness = c(phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[2], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[3], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[4], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[5], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[6], NArm = T)),
                              phyloseq::ntaxa(tax_glom(physeq,taxrank = rank_names(physeq)[7], NArm = T))))
  return (df)
}


taxa_richness_perregion_sidle<-do.call("rbind",list(taxa_richness_sidle(dropna_taxa_names(physeq_V1V2),"V1V2"),taxa_richness_sidle(dropna_taxa_names(physeq_V2V3),"V2V3"),
                                              taxa_richness_sidle(dropna_taxa_names(physeq_V3V4),"V3V4"),taxa_richness_sidle(dropna_taxa_names(physeq_V4V5),"V4V5"),
                                              taxa_richness_sidle(dropna_taxa_names(physeq_V5V7),"V5V7"),taxa_richness_sidle(dropna_taxa_names(physeq_V7V9),"V7V9"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V1V2V4V5_sidle)),"V1V2-V4V5"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V1V2V5V7_sidle)),"V1V2-V5V7"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V1V2V7V9_sidle)),"V1V2-V7V9"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V2V3V5V7_sidle)),"V2V3-V5V7"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V2V3V7V9_sidle)),"V2V3-V7V9"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_V3V4V7V9_sidle)),"V3V4-V7V9"),
                                              taxa_richness_sidle(dropna_unassigned(dropna_taxa_names_sidle(physeq_full_sidle)),"full")))

taxa_richness_perregion_sidle$rank = factor(taxa_richness_perregion_sidle$rank, levels = c(rank_names(physeq_V1V2)[-1]))
taxa_richness_perregion_sidle$region= factor(taxa_richness_perregion_sidle$region, levels = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9",
                                                                                              "V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                                                                              "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"))

F4A<-ggplot(subset(taxa_richness_perregion_sidle),aes(x=region,y=richness, fill = region))+
    geom_col()+scale_fill_manual(values = c(brewer.pal(6, "Dark2"),"gray70","gray65","gray60","gray55","gray50","gray45","black"))+
  scale_y_continuous(breaks=c(0,35,70,105,140),limits=c(0,140),expand = c(0.01,0))+labs(y="Richness")+
    facet_wrap(.~rank, nrow = 1)+theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                       legend.position = "none",axis.title.x = element_blank())
F4A

####Figure S7A####
#05/08/21 - sidle richness fold change
taxa_richness_perregion_sidle_fc<-data.frame(amplicons = rep(c("two","six (all)"),6), 
                                             rank = c(rep('Phylum',2),rep('Class',2),rep('Order',2),
                                                      rep('Family',2),rep('Genus',2),rep('Species',2)))
taxa_richness_perregion_sidle_fc$avg_richness_fc<-c(mean(subset(taxa_richness_perregion_sidle,rank == "Phylum"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9","V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Phylum"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
  mean(subset(taxa_richness_perregion_sidle,rank == "Phylum"&region=="full")$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Phylum"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Class"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9","V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Class"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Class"&region=="full")$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Class"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
  mean(subset(taxa_richness_perregion_sidle,rank == "Order"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9","V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Order"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
  mean(subset(taxa_richness_perregion_sidle,rank == "Order"&region=="full")$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Order"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Family"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9", "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Family"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Family"&region=="full")$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Family"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
  mean(subset(taxa_richness_perregion_sidle,rank == "Genus"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9","V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Genus"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
  mean(subset(taxa_richness_perregion_sidle,rank == "Genus"&region=="full")$richness)/
  mean(subset(taxa_richness_perregion_sidle,rank == "Genus"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Species"&region%in%c("V1V2-V4V5","V1V2-V5V7","V1V2-V7V9","V2V3-V5V7","V2V3-V7V9","V3V4-V7V9"))$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Species"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness),
mean(subset(taxa_richness_perregion_sidle,rank == "Species"&region=="full")$richness)/
mean(subset(taxa_richness_perregion_sidle,rank == "Species"&region%in%c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9"))$richness))

taxa_richness_perregion_sidle_fc$amplicons = factor(taxa_richness_perregion_sidle_fc$amplicons, levels = c("two","six (all)"))
taxa_richness_perregion_sidle_fc$rank = factor(taxa_richness_perregion_sidle_fc$rank, levels = c(rank_names(physeq_V1V2)[-1]))

FS7A<-ggplot(subset(taxa_richness_perregion_sidle_fc),aes(x=amplicons,y=avg_richness_fc, fill = amplicons))+
  geom_col()+scale_fill_manual(values = c("gray55","black"))+
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(0,4),expand=c(0.01,0))+
  facet_wrap(.~rank, nrow = 1)+theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                     legend.position = "none")+
  geom_hline(yintercept = 1,linetype="dashed")+
  labs(x="N of amplicons",y="Richness fold-change")
FS7A

####Figure S8####
#06/08 - sidle - taxa relative abundance
preprocess_physeq_tomerge<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Dadaed-chimremoved-nonbactremoved-contamremoved/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Taxonomy/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)))
  physeq_region<-prune_samples(!(sample_names(physeq_region) %in% c("3","4","5")), physeq_region)
  physeq_region<-prune_taxa(taxa_sums(physeq_region) > 0, physeq_region) 
  
  physeq_pcoa<-physeq_region
  sample_names(physeq_pcoa)<-paste0(sample_names(physeq_region),"_",region)
  sample_data(physeq_pcoa)<-data.frame(sample_data(physeq_pcoa)) %>% mutate(region=region) %>% mutate(sampleid=sample_names(physeq_region))
  return(physeq_pcoa)
}

physeq_V1V2_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V1V2"))
physeq_V2V3_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V2V3"))
physeq_V3V4_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V3V4"))
physeq_V4V5_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V4V5"))
physeq_V5V7_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V5V7"))
physeq_V7V9_tomerge<-coherent_taxa_names(preprocess_physeq_tomerge("V7V9"))

preprocess_physeq_tomerge_sidle<-function(region){
  table<-as.data.frame(read_qza(paste0("Results/Sidle/",region,"_table.qza"))[["data"]])
  taxa<-parse_taxonomy(read_qza(paste0("Results/Sidle/",region,"_taxonomy.qza"))[["data"]][,1:2], tax_sep=";")
  taxa<-subset(taxa,rownames(taxa)%in%rownames(table))
  physeq_region<-phyloseq(sample_data(metadata),
                          otu_table(table,taxa_are_rows = T),
                          tax_table(as.matrix(taxa)))
  physeq_region<-prune_samples(!(sample_names(physeq_region) %in% c("3","4","5")), physeq_region)
  physeq_region<-prune_taxa(taxa_sums(physeq_region) > 0, physeq_region) 
  
  physeq_pcoa<-physeq_region
  sample_names(physeq_pcoa)<-paste0(sample_names(physeq_region),"_",region)
  sample_data(physeq_pcoa)<-data.frame(sample_data(physeq_pcoa)) %>% mutate(region=region) %>% mutate(sampleid=sample_names(physeq_region))
  return(physeq_pcoa)
}

physeq_V1V2V4V5_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V1V2-V4V5")))
physeq_V1V2V5V7_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V1V2-V5V7")))
physeq_V1V2V7V9_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V1V2-V7V9")))
physeq_V2V3V5V7_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V2V3-V5V7")))
physeq_V2V3V7V9_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V2V3-V7V9")))
physeq_V3V4V7V9_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("V3V4-V7V9")))
physeq_full_tomerge_sidle<-remove_unassigned(coherent_taxa_names_sidle(preprocess_physeq_tomerge_sidle("full")))

plot_physeq_bar_all_sidle<-function(rank,min_ra=0,in_n_samples=1){
  physeq_all_ra<-transform_sample_counts(tax_glom(merge_phyloseq(
    physeq_V1V2_tomerge,physeq_V2V3_tomerge,physeq_V3V4_tomerge,physeq_V4V5_tomerge,physeq_V5V7_tomerge,physeq_V7V9_tomerge,
    physeq_V1V2V4V5_tomerge_sidle,physeq_V1V2V5V7_tomerge_sidle,physeq_V1V2V7V9_tomerge_sidle,
    physeq_V2V3V5V7_tomerge_sidle,physeq_V2V3V7V9_tomerge_sidle,physeq_V3V4V7V9_tomerge_sidle,
    physeq_full_tomerge_sidle),rank,NArm = F),
    function(x) x/sum(x))
  metadata<-data.frame(sample_data(physeq_all_ra))
  metadata$sampleid<-factor(paste0("#",data.frame(sample_data(physeq_all_ra))$sampleid),levels = sample_list)
 
  metadata$region<-factor(data.frame(sample_data(physeq_all_ra))$region,levels = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9",
                                                                                     "V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                                                                     "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"))
  sample_data(physeq_all_ra)<-metadata
  plot_bar(prune_taxa(genefilter_sample(physeq_all_ra, filterfun_sample(function(x) x >= min_ra), 
                                        A = in_n_samples),physeq_all_ra), "region", fill=rank)+facet_wrap(.~sampleid, nrow = 2, labeller = sample_labeller)+
    labs(y="Relative abundance")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values = glasbey())+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0),legend.text = element_text(size=8.5), legend.key.size = unit(0.4, 'cm'))
}

FS8<-ggarrange(plot_physeq_bar_all_sidle("Phylum")+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=2,byrow = T)),
plot_physeq_bar_all_sidle("Class")+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=2,byrow = T)),
plot_physeq_bar_all_sidle("Order",0.0195,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
plot_physeq_bar_all_sidle("Family",0.1,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
plot_physeq_bar_all_sidle("Genus",0.2,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
plot_physeq_bar_all_sidle("Species",0.245,1)+theme(legend.position = "bottom", axis.title.x = element_blank())+guides(fill=guide_legend(nrow=4,byrow = T)),
heights=c(0.9,0.9,1.1,1.1,1.1,1.1),nrow = 6) #export as 1650x3000
FS8

####Figure 4C, S7C####
#06/08 - sidle - average taxa relative abundance
plot_physeq_bar_avg_sidle<-function(rank,min_ra=0,in_n_samples=1){
  physeq_all_ra<-transform_sample_counts(tax_glom(merge_phyloseq(
    physeq_V1V2_tomerge,physeq_V2V3_tomerge,physeq_V3V4_tomerge,physeq_V4V5_tomerge,physeq_V5V7_tomerge,physeq_V7V9_tomerge,
    physeq_V1V2V4V5_tomerge_sidle,physeq_V1V2V5V7_tomerge_sidle,physeq_V1V2V7V9_tomerge_sidle,
    physeq_V2V3V5V7_tomerge_sidle,physeq_V2V3V7V9_tomerge_sidle,physeq_V3V4V7V9_tomerge_sidle,
    physeq_full_tomerge_sidle),rank,NArm = F),
    function(x) x/sum(x))
  metadata<-data.frame(sample_data(physeq_all_ra))
  metadata$sampleid<-factor(data.frame(sample_data(physeq_all_ra))$sampleid,levels = c('6','7','9','10','12','13','14','16','19','21','23','25','37','39','42'))
  
  sample_data(physeq_all_ra)<-metadata
  physeq_avg<-phyloseq_average(physeq_all_ra,group = "region",avg_type = "arithmetic")
  physeq_avg_filtered<-prune_taxa(genefilter_sample(physeq_avg, filterfun_sample(function(x) x >= min_ra), 
                                                    A = in_n_samples),physeq_avg)

  p<-plot_bar(physeq_avg_filtered, fill=rank)+
    labs(y="Relative abundance")+scale_y_continuous(expand=c(0,0), limits = c(0,1))+scale_fill_manual(values = glasbey())+
    theme(axis.text.x = element_text(vjust = 0.5, hjust=0),legend.text = element_text(size=10), legend.key.size = unit(0.5, 'cm'))
  
  p$data$Sample <- factor(p$data$Sample, levels = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9",
                                                    "V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                                    "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"))
  print(p)
}

F4C<-plot_physeq_bar_avg_sidle("Genus",0.015,1)+theme(axis.title.x = element_blank())+guides(fill=guide_legend(ncol=2,byrow = T))

FS7C<-plot_physeq_bar_avg_sidle("Species",0.02,1)+theme(axis.title.x = element_blank())+guides(fill=guide_legend(ncol=2,byrow = T))

####Figure 4B####
#10/08 - sidle - ambiguity
calc_ambiguity<-function(physeq,region){
  taxa<-as.data.frame(tax_table(remove_unassigned(coherent_taxa_names_sidle(physeq))))
  result<-as.data.frame(matrix(nrow = length(unique(taxa$Species)), ncol = 4))
  colnames(result)<-c("region","taxon","16s_mapped", "total_ambiguity")
  result$region<-region
  result$taxon<-unique(taxa$Species)
  pos=1
  for (taxon in unique(taxa$Species)){
    taxa_subset<-subset(taxa, Species == taxon)
    result[pos,3]<-sum(str_count(rownames(taxa_subset),pattern="\\|"))+nrow(taxa_subset)
    pos=pos+1
  }
  result$total_ambiguity<-sum(log(result[["16s_mapped"]]))/nrow(result)#log(sum(result[["16s_mapped"]])/nrow(result))
  #result$total_ambiguity<-log(sum(result[["16s_mapped"]])/nrow(result))
  return(result)
}

ambiguity<-data.frame(region = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9",
                                 "V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                 "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"),
ambiguity=c(calc_ambiguity(physeq_V1V2_sidle,region="V1V2")[1,4],
            calc_ambiguity(physeq_V2V3_sidle,region="V2V3")[1,4],
            calc_ambiguity(physeq_V3V4_sidle,region="V3V4")[1,4],
            calc_ambiguity(physeq_V4V5_sidle,region="V4V5")[1,4],
            calc_ambiguity(physeq_V5V7_sidle,region="V5V7")[1,4],
            calc_ambiguity(physeq_V7V9_sidle,region="V7V9")[1,4],
            calc_ambiguity(physeq_V1V2V4V5_sidle, region = "V1V2-V4V5")[1,4],
            calc_ambiguity(physeq_V1V2V5V7_sidle, region = "V1V2-V5V7")[1,4],
            calc_ambiguity(physeq_V1V2V7V9_sidle, region = "V1V2-V7V9")[1,4],
            calc_ambiguity(physeq_V2V3V5V7_sidle, region = "V2V3-V5V7")[1,4],
            calc_ambiguity(physeq_V2V3V7V9_sidle, region = "V2V3-V7V9")[1,4],
            calc_ambiguity(physeq_V3V4V7V9_sidle, region = "V3V4-V7V9")[1,4],
            calc_ambiguity(physeq_full_sidle,region="full")[1,4]))

ambiguity$region= factor(ambiguity$region, levels = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9","V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                                       "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"))

F4B<-ggplot(ambiguity,aes(x=region,y=ambiguity, fill = region))+labs(y="Ambiguity")+
  geom_col()+scale_fill_manual(values = c(brewer.pal(6, "Dark2"),"gray70","gray65","gray60","gray55","gray50","gray45","black"))+
  scale_y_continuous(breaks=c(0,1,2,3),limits=c(0,3),expand = c(0.01,0))+theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
                                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                     legend.position = "none",axis.title.x = element_blank())

F4B
####Figure S7B####
#16/08/21 - sidle - % of assigned sequences per level per dataset

assigned_seqs_sidle<-function(physeq,rank){
  physeq_treated<-merge_samples(tax_glom(coherent_taxa_names_sidle(physeq), taxrank = rank), group = "Sample_or_Control")
  data<-merge(as.data.frame(tax_table(physeq_treated)),t(as.data.frame(otu_table(physeq_treated))), by = 0)
  
  data <- data %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Phylum_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Class_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Order_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Family_.*", "none")) %>% 
    mutate(!!rank := str_replace(!!sym(rank), "Genus_.*", "none"))
  
  result<-100*(sum(subset(data, data[[rank]]!="none")[["True Sample"]])/sum(data[["True Sample"]]))
  return(result)
}

assigned_seqs_perregion_perlevel_sidle<-data.frame(region=c(rep("V1V2",6),rep("V2V3",6),rep("V3V4",6),rep("V4V5",6),rep("V5V7",6),rep("V7V9",6),
                                                      rep("V1V2-V4V5",6),rep("V1V2-V5V7",6),rep("V1V2-V7V9",6),rep("V2V3-V5V7",6),rep("V2V3-V7V9",6),rep("V3V4-V7V9",6),
                                                      rep("full",6)),
                                             rank=rep(c("Phylum", "Class", "Order", "Family", "Genus", "Species"),13),
                                             assigned=c(assigned_seqs(physeq_V1V2, rank = "Phylum"), assigned_seqs(physeq_V1V2, rank = "Class"),
                                                        assigned_seqs(physeq_V1V2, rank = "Order"),assigned_seqs(physeq_V1V2, rank = "Family"),
                                                        assigned_seqs(physeq_V1V2, rank = "Genus"),assigned_seqs(physeq_V1V2, rank = "Species"),
                                                        assigned_seqs(physeq_V2V3, rank = "Phylum"), assigned_seqs(physeq_V2V3, rank = "Class"),
                                                        assigned_seqs(physeq_V2V3, rank = "Order"),assigned_seqs(physeq_V2V3, rank = "Family"),
                                                        assigned_seqs(physeq_V2V3, rank = "Genus"),assigned_seqs(physeq_V2V3, rank = "Species"),
                                                        assigned_seqs(physeq_V3V4, rank = "Phylum"), assigned_seqs(physeq_V3V4, rank = "Class"),
                                                        assigned_seqs(physeq_V3V4, rank = "Order"),assigned_seqs(physeq_V3V4, rank = "Family"),
                                                        assigned_seqs(physeq_V3V4, rank = "Genus"),assigned_seqs(physeq_V3V4, rank = "Species"),
                                                        assigned_seqs(physeq_V4V5, rank = "Phylum"), assigned_seqs(physeq_V4V5, rank = "Class"),
                                                        assigned_seqs(physeq_V4V5, rank = "Order"),assigned_seqs(physeq_V4V5, rank = "Family"),
                                                        assigned_seqs(physeq_V4V5, rank = "Genus"),assigned_seqs(physeq_V4V5, rank = "Species"),
                                                        assigned_seqs(physeq_V5V7, rank = "Phylum"), assigned_seqs(physeq_V5V7, rank = "Class"),
                                                        assigned_seqs(physeq_V5V7, rank = "Order"),assigned_seqs(physeq_V5V7, rank = "Family"),
                                                        assigned_seqs(physeq_V5V7, rank = "Genus"),assigned_seqs(physeq_V5V7, rank = "Species"),
                                                        assigned_seqs(physeq_V7V9, rank = "Phylum"), assigned_seqs(physeq_V7V9, rank = "Class"),
                                                        assigned_seqs(physeq_V7V9, rank = "Order"),assigned_seqs(physeq_V7V9, rank = "Family"),
                                                        assigned_seqs(physeq_V7V9, rank = "Genus"),assigned_seqs(physeq_V7V9, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V1V2V4V5_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V1V2V5V7_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V1V2V7V9_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V2V3V5V7_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V2V3V7V9_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Order"),assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_V3V4V7V9_sidle, rank = "Species"),
                                                        assigned_seqs_sidle(physeq_full_sidle, rank = "Phylum"), assigned_seqs_sidle(physeq_full_sidle, rank = "Class"),
                                                        assigned_seqs_sidle(physeq_full_sidle, rank = "Order"),assigned_seqs_sidle(physeq_full_sidle, rank = "Family"),
                                                        assigned_seqs_sidle(physeq_full_sidle, rank = "Genus"),assigned_seqs_sidle(physeq_full_sidle, rank = "Species")))

assigned_seqs_perregion_perlevel_sidle$rank = factor(assigned_seqs_perregion_perlevel_sidle$rank, 
                                               levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
assigned_seqs_perregion_perlevel_sidle$region = factor(assigned_seqs_perregion_perlevel_sidle$region, 
                                                       levels = c("V1V2","V2V3","V3V4","V4V5","V5V7","V7V9","V1V2-V4V5","V1V2-V5V7","V1V2-V7V9",
                                                                  "V2V3-V5V7","V2V3-V7V9","V3V4-V7V9","full"))

FS7B<-ggplot(assigned_seqs_perregion_perlevel_sidle,aes(x=region,y=assigned, fill = region))+
  geom_col()+scale_fill_manual(values = c(brewer.pal(6, "Dark2"),"gray70","gray65","gray60","gray55","gray50","gray45","black"))+
  scale_y_continuous(breaks=c(0,25,50,75,100),limits=c(0,100),expand = c(0.01,0))+
  labs(y="% of assigned sequences")+facet_wrap(.~rank, nrow = 1)+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

FS7B

################################Arranging figures############################
####Figure 1####
ggarrange(annotate_figure(F1A,fig.lab = "A",fig.lab.size = 16),
          ggarrange(annotate_figure(F1B,fig.lab ="B",fig.lab.size = 16),
                    annotate_figure(F1C+labs(x="Library")+theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
                                    fig.lab ="C",fig.lab.size = 16),nrow=1,widths = c(0.35,1)),
          annotate_figure(F1D,fig.lab = "D",fig.lab.size = 16),
          nrow=3, heights = c(1.2,1,1))

####Figure S4####
ggarrange(annotate_figure(FS4A+coord_equal(), fig.lab = "A",fig.lab.size = 16),
          annotate_figure(FS4B+theme(legend.position = "none")+coord_equal(), fig.lab = "B", fig.lab.size = 16), nrow=2,heights = c(0.6,1))

####Figure 3####
ggarrange(annotate_figure(F3A, fig.lab = "A", fig.lab.size = 16),
          annotate_figure(F3B, fig.lab = "B", fig.lab.size = 16),
          annotate_figure(F3C+theme(legend.text = element_text(size=8.5), legend.key.size = unit(0.4, 'cm')), fig.lab = "C", fig.lab.size = 16),
          nrow = 3, widths = c(0.8,0.9,1)) #900x1200

####Figure 2####
ggarrange(ggarrange(ggarrange(annotate_figure(F2A, fig.lab = "A", fig.lab.size = 16),
                              annotate_figure(F2B, fig.lab = "B", fig.lab.size = 16), nrow = 1, widths = c(1/5,1)),
          annotate_figure(F2C, fig.lab = "C", fig.lab.size = 16),nrow=1,widths=c(1,0.25)),
          annotate_figure(F2D, fig.lab = "D", fig.lab.size = 16),heights = c(1,1),nrow = 2)

####Figure 4####
ggarrange(annotate_figure(F4A, fig.lab = "A", fig.lab.size = 16),
          annotate_figure(F4B, fig.lab = "B", fig.lab.size = 16),
          annotate_figure(F4C, fig.lab = "C", fig.lab.size = 16),
          nrow = 3, heights = c(0.55,0.45,1)) #1050x1050

####Figure S7####
ggarrange(annotate_figure(FS7A, fig.lab = "A", fig.lab.size = 16),
          annotate_figure(FS7B, fig.lab = "B", fig.lab.size = 16),
          annotate_figure(FS7C, fig.lab = "C", fig.lab.size = 16),
          nrow = 3, heights = c(0.45,0.55,1)) #1000x1000

#########################TESTS - 21/10/21#########################
#generating some results to write V1V2-specific section

plot_bar(phyloseq_average(tax_glom(coherent_taxa_names(physeq_V1V2),taxrank = "Phylum"), avg_type = "arithmetic"), fill = "Phylum")

otu_table(phyloseq_average(tax_glom(coherent_taxa_names(physeq_V1V2),taxrank = "Phylum"), avg_type = "arithmetic"))
tax_table(phyloseq_average(tax_glom(coherent_taxa_names(physeq_V1V2),taxrank = "Phylum"), avg_type = "arithmetic"))

eval_exclusive_genera<-function(physeq, list_of_genera){
  
physeq_avg<-phyloseq_average(tax_glom(coherent_taxa_names(physeq),taxrank = "Genus"), avg_type = "arithmetic")
otu_exc<-as.data.frame(otu_table(physeq_avg))
tax_exc<-as.data.frame(tax_table(physeq_avg))
tax_exc<-subset(tax_exc, Genus %in% list_of_genera)
otu_exc<-subset(otu_exc,rownames(otu_exc)%in%rownames(tax_exc))
rownames(otu_exc)<-tax_exc$Genus

print(paste("sum:",sum(otu_exc)))
print(paste("min:",min(otu_exc)))
print(paste("max:",max(otu_exc)))

return(otu_exc)
}

v1v2_exc_g<-eval_exclusive_genera(physeq_V1V2, c('Alkalibacterium','Arcanobacterium','Chromohalobacter','Chryseobacterium',
                                    'Clostridium_sensu_stricto_1','Comamonas','Ezakiella','Facklamia','Herbaspirillum',
                                    'Jeotgalibaca','Lactococcus','Mycoplasma','Oribacterium',
                                    'Salipaludibacillus','Tepidimonas'))

# [1] "sum: 0.0398540145796059"
# [1] "min: 5.47660122128207e-06"
# [1] "max: 0.0285581826611034"

eval_exclusive_genera(physeq_V2V3, c('Azoarcus','Flavonifractor','Fournierella','Luteimonas','oc32','Polaromonas',
                                     'Pseudonocardia','[Ruminococcus]_torques_group','0319-6G20'))

# [1] "sum: 0.0148098151268408"
# [1] "min: 1.62881924178464e-05"
# [1] "max: 0.00929539699091315"

eval_exclusive_genera(physeq_V3V4, c('Brachybacterium','Chthoniobacter','Microbacterium','Nannocystis',
                                     'Porphyromonas','S5-A14a'))

# [1] "sum: 0.0644935802052112"
# [1] "min: 1.28802352789644e-05"
# [1] "max: 0.0278227984873042"

eval_exclusive_genera(physeq_V4V5, c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium','Azovibrio',
                                     'Bifidobacterium','Delftia','Schlegella'))

# [1] "sum: 0.00406277606997556"
# [1] "min: 1.80836497972048e-05"
# [1] "max: 0.00395909851188937"

eval_exclusive_genera(physeq_V5V7, c('Butyricicoccus','Desulfovermiculus','Eikenella','Veillonella'))

# [1] "sum: 0.00118704594787457"
# [1] "min: 2.46263717658648e-05"
# [1] "max: 0.000713315217391304"

eval_exclusive_genera(physeq_V7V9, c('Acidovorax','Alcanivorax','Cruoricaptor','Desulfovibrio',
                                     'Moraxella','Sphingobacterium','Sphingomonas'))
# [1] "sum: 0.0202571735250034"
# [1] "min: 1.73456578797216e-05"
# [1] "max: 0.012931422907092"
