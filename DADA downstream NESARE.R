#NESARE microbiome data analysis 
#Got ASVs from DADA2 pipeline using default parameters.
#All sequences from the 3 objectives were analyzed together

#Last updated: MLR 9/30/2022

#Set working directory
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/NESARE/ASV")

#Attach libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(zCompositions)
library(compositions)
library(viridis)
library(svglite)
library(pairwiseAdonis)
library(decontam)

#### Import data ####
asvs_all<-read.csv('ASV_all.csv', header = TRUE, row.names = 1)
taxon_all<-as.data.frame(read.csv('Taxon_all.csv', header = TRUE, row.names = 1))
metadata_all<-read.csv('metadata_all.csv', header = TRUE, row.names = 1)

#Add '_unclassified' marker to NAs in taxon table
taxon_all$Phylum<-ifelse(is.na(taxon_all$Phylum), paste(taxon_all$Kingdom, "unclassified", sep = '_'), taxon_all$Phylum)
taxon_all$Class<-ifelse(is.na(taxon_all$Class), paste(taxon_all$Phylum, "unclassified", sep = '_'), taxon_all$Class)
taxon_all$Order<-ifelse(is.na(taxon_all$Order), paste(taxon_all$Class, "unclassified", sep = '_'), taxon_all$Order)
taxon_all$Family<-ifelse(is.na(taxon_all$Family), paste(taxon_all$Order, "unclassified", sep = '_'), taxon_all$Family)
taxon_all$Genus<-ifelse(is.na(taxon_all$Genus), paste(taxon_all$Family, "unclassified", sep = '_'), taxon_all$Genus)
taxon_all$Species<-ifelse(is.na(taxon_all$Species), paste(taxon_all$Genus, "unclassified", sep = '_'), taxon_all$Species)

#Remove extra _unclassified
taxon_all$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Class)
taxon_all$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Species)

#Convert asv and taxon tables to matrix
asvs_all<-as.matrix(asvs_all)
taxon_all<-as.matrix(taxon_all)

#Make phyloseq object
phyloseq_16s<-phyloseq(otu_table(asvs_all, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_all))

#Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- phyloseq_16s %>% subset_taxa( Order!="Chloroplast" )
physeq_16s <- physeq_16s %>% subset_taxa( Family!="Mitochondria" )


#Get ASV table from phyloseq object
asv.16s<-as.data.frame(t(otu_table(physeq_16s)))
tail(rowSums(asv.16s))

#Remove ASVs with zero counts in all samples
asv.16s<-asv.16s[ which(rowSums(asv.16s)>0),]
asv.16s<-t(asv.16s)

taxon_16s_clean<-as.matrix(tax_table(physeq_16s))


#### Objective 1: Characterization of the microbiiota of ice cream processing facilities ####
#Split phyloseq by Objective
ps_Obj1<-subset_samples(physeq_16s, Obj == 1)

#Get ASV ad metadata table from phyloseq object
asv_obj1<-as.data.frame(t(otu_table(ps_Obj1)))
meta_obj1<-sample_data(ps_Obj1)

#Remove ASVs with zero counts in all samples
asv_obj1<-asv_obj1[ which(rowSums(asv_obj1)>0),]
asv_obj1<-as.data.frame(asv_obj1)

## STACKED BARPLOT 
#Convert ASV table to appropriate format. Following step requires samples on rows and ASVs in columns
head(t(asv_obj1)) 

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_obj1<-t(cmultRepl(t(asv_obj1), label=0, method="CZM", output="p-counts")) #331,770 corrected values

head(asv.n0_obj1) #output table needs to have samples in columns and ASVs in rows

#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions
asv.n0.acomp_obj1<-as.data.frame(acomp(t(asv.n0_obj1)), total=1)
rowSums(asv.n0.acomp_obj1) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_Obj1 = phyloseq(otu_table(asv.n0.acomp_obj1, taxa_are_rows = FALSE), tax_table(taxon_all),sample_data(meta_obj1))

#Make tables with RA at ASV, family and genus level
asv_CoDa_obj1 <- phyloseq_Obj1 %>%  
  transform_sample_counts(function(x) {x *100} ) %>% #Changes RA to % 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

family_CoDa_obj1 <- phyloseqObj1 %>%  
  tax_glom(taxrank = "Family") %>% #Agglomerates all ASVs that have the same Family
  transform_sample_counts(function(x) {x*100} ) %>% #Changes RA to %  
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

genus_CoDa_obj1 <- phyloseqObj1 %>%  
  tax_glom(taxrank = "Genus") %>% #Agglomerates all ASVs that have the same Family
  transform_sample_counts(function(x) {x*100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

write.csv(ASV_CoDa_obj1, "ASV Obj1.csv")
write.csv(family_CoDa_obj1, "Family Obj1.csv")
write.csv(genus_CoDa_obj1, "Genus Obj1.csv")

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_over2abund_obj1 <- filter(asv_CoDa_obj1, Abundance >2)

#Calculate total abundance of each ASV by Facility
asv_total_obj1<-as.data.frame(xtabs(Abundance ~ Facility, asv_over2abund_obj1))

#Calculate "Others" total
asv_total_obj1$Diff<-100-asv_total_obj1$Freq

#Add Others Category to dataframe
FA_asv_others_Obj1<-c("Others","FA",asv_total_obj1[1,3],"FA","A","Bacteria","","","","","","","","","","Others (less than 2% relative abundance)","Others (less than 2% relative abundance)")
FB_asv_others_Obj1<-c("Others","FB",asv_total_obj1[2,3],"FB","B","Bacteria","","","","","","","","","","Others (less than 2% relative abundance)","Others (less than 2% relative abundance)")
FC_asv_others_Obj1<-c("Others","FC",asv_total_obj1[3,3],"FC","C","Bacteria","","","","","","","","","","Others (less than 2% relative abundance)","Others (less than 2% relative abundance)")


#Bind rows
asv_over2abund_final_obj1<-rbind(asv_over2abund_obj1, FA_asv_others_Obj1, FB_asv_others_Obj1, FC_asv_others_Obj1)
asv_over2abund_final_obj1$Abundance<-as.numeric(asv_over2abund_final_obj1$Abundance)

#Calculate % RA by genera per facility
genera_perc<-asv_over2abund_final_obj1%>%
  group_by(Sample, Genus)%>%
  tally(Abundance)
  
#Use ggplot2 to make the stacked barplot 
#Genus with Family fill color
barplot_obj1<- ggplot(asv_over2abund_final_obj1, aes(x = Facility, y = Abundance , fill = Genus))  + 
  geom_bar(stat = "identity", color='black') + 
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="right")+
  ggtitle("Microbiota of ice cream processing facilities", subtitle = "Genus taxonomic level")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  scale_fill_discrete(name="Genus")
ggsave("Barplot Obj1.png", plot=barplot_obj1, device="png", width=8, height=13, units="in", dpi=1000)
ggsave("Barplot Obj1.svg", plot=barplot_obj1, device="svg", width=8, height=13, units="in", dpi=1000)




#### Obj2: Attachment of lactic acid bacteria in the presence of microbiomes ####
#Split phyloseq by Objective
ps_Obj2<-subset_samples(physeq_16s, Obj == 2)

#Get ASV ad metadata table from phyloseq object
asv_obj2<-as.data.frame(t(otu_table(ps_Obj2)))
meta_obj2<-sample_data(ps_Obj2)

#Remove ASVs with zero counts in all samples
asv_obj2<-asv_obj2[ which(rowSums(asv_obj2)>0),]

#CHECH NEGATIVE CONTROLS & DECONTAMINATE

#Make Phyloseq
phyloseqObj2<-phyloseq(otu_table(asv_obj2, taxa_are_rows = TRUE),tax_table(taxon_all), sample_data(meta_obj2))

#Subset Controls
phyloseq_Obj2NC<-subset_samples(phyloseqObj2, SampleID =="NegContr")
phyloseq_Obj2PC<-subset_samples(phyloseqObj2, SampleID =="PosContr")

#Make long table
phyloseq_Obj2_NC<-psmelt(phyloseq_Obj2NC)
phyloseq_Obj2_PC<-psmelt(phyloseq_Obj2PC)

#Remove rows with less than 100 reads
phyloseq_Obj2NC<-subset(phyloseq_Obj2_NC, Abundance>100)
phyloseq_Obj2PC<-subset(phyloseq_Obj2_PC, Abundance>100)

#Plot reads of control by  NC
barplot_Obj2_NC<-ggplot(phyloseq_Obj2NC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black', fill='#00C08D')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Obj2 - NC")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_Obj2_NC.png", plot =barplot_Obj2_NC, device="png", width=20, height=10, units="in",dpi=600)

#Plot reads of control by  PC
barplot_Obj2_PC<-ggplot(phyloseq_Obj2PC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black', fill='#00C08D')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Obj2 - PC")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_Obj2_PC.png", plot =barplot_Obj2_PC, device="png", width=20, height=10, units="in",dpi=600)

#Calculate relative abundance
Obj2PC_RA<-phyloseq_Obj2PC%>%
  filter_taxa(function(x) x >50, TRUE)%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()%>%
  arrange(desc(Abundance))



#Decontaminate sequencing data

#Detect contaminants with prevalence method stringent method (threshold=0.5)
sample_data(phyloseqObj2)$is.neg <- sample_data(phyloseqObj2)$SampleID == "NegContr"
contamdf.prev <- isContaminant(phyloseqObj2, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) 

#Remove identified contaminants (prevalence method) from ASV table
phyloseq_obj2.noncontam <- prune_taxa(!contamdf.prev$contaminant, phyloseqObj2)
phyloseq_obj2.noncontam

#Remove negative and positive controls from ASV table
phyloseq_obj2_clean<-subset_samples(phyloseq_obj2.noncontam, Facility!="Control" )

#Get ASV table from phyloseq object
asv_obj2_clean<-as.data.frame(otu_table(phyloseq_obj2_clean))
tail(rowSums(asv_obj2_clean))

#Remove ASVs with zero counts in all samples
asv_obj2_clean<-asv_obj2_clean[ which(rowSums(asv_obj2_clean)>0),]

#Get metadata 
metadata_obj2_clean<-subset(metadata_all, Obj==2 & Facility != "Control")


#Compositional analysis of microbiome data -based on Microbiome Analysis in R. Chap 10.
#Step 1: Convert ASV table to appropriate format. Following step requires samples on rows and ASVs in columns
head(t(asv_obj2_clean)) 

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_obj2<-t(cmultRepl(t(asv_obj2_clean), label=0, method="CZM", output="p-counts")) 

head(asv.n0_obj2) #output table needs to have samples in columns and ASVs in rows

#Step 3: Convert data to proportions
asv.n0_obj2_prop<-apply(asv.n0_obj2, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_obj2_prop_f<-asv.n0_obj2[apply(asv.n0_obj2_prop, 1, min) > 0.0000001, ]
head(asv.n0_obj2_prop_f) #Check that samples are on columns and asvs in rows


#Step 5: perform CLR transformation
asv.n0.clr_obj2<-t(apply(asv.n0_obj2_prop_f, 2, function(x){log(x)-mean(log(x))}))
head(asv.n0.clr_obj2) #Check output table. Samples should be in rows and asvs in columns

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_obj2<-prcomp(asv.n0.clr_obj2)

png("Screeplot - PCA - Obj 2.png")
par(mar=c(2,2,2,2))
screeplot(pc.clr_obj2, type='lines', main="Obj 2")
dev.off()

#Calculate total variance of the data
mvar.clr_obj2<-mvar(asv.n0.clr_obj2)

#Display results
row_obj2<-rownames(asv.n0.clr_obj2) #Make vector with sample names
pc_out_obj2<-as.data.frame(pc.clr_obj2$x[,1:2]) #Get PC1 and PC2
pc_out_meta_obj2<-as.data.frame(bind_cols(pc_out_obj2,metadata_obj2_clean)) #Add metadata information
row.names(pc_out_meta_obj2)<-row_obj2 #Add rownames to dataframe

# Make PCA plot
PCA_Obj2 <- ggplot(pc_out_meta_obj2, aes(x=PC1,y=PC2, color=Facility, shape=Treatment))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_obj2$sdev[1]^2/mvar.clr_obj2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_obj2$sdev[2]^2/mvar.clr_obj2*100, digits=1), "%", sep="")) +
  ggtitle("PCA - Obj 2 by facility")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_Obj2
ggsave("PCA_Obj2.png", plot=PCA_Obj2, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_Obj2.svg", plot=PCA_Obj2, device="svg", width=8, height=8, units="in", dpi=600)

# PERMANOVA
#Calculate Aitchinson distance
dist_obj2<-dist(asv.n0.clr_obj2, method='euclidean')

#By Facility
permanova_obj2_fac<-pairwise.adonis2(dist_obj2~Facility, data=metadata_obj2_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_obj2_fac


#By Treatment
permanova_obj2_trt<-pairwise.adonis2(dist_obj2~Treatment, data=metadata_obj2_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_obj2_trt


# STACKED BARPLOT 
#Note: used compositional approach to transform the sample counts to compositions. 


#Transform sample counts into compositions
asv.n0.acomp_obj2<-as.data.frame(acomp(t(asv.n0_obj2)), total=1)
rowSums(asv.n0.acomp_obj2) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_obj2_RA <- phyloseq(otu_table(asv.n0.acomp_obj2, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_obj2_clean))

#Make long format table from Phyloseq object - only tretment samples
asv_Obj2_long <- subset_samples(phyloseq_obj2_RA, Facility!="Control") %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment for each ASV
asv_Obj2_mean<-asv_Obj2_long %>%
  group_by(OTU, Treatment, Facility, Family, Genus, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_Obj2_over2abund <- filter(asv_Obj2_mean, Mean>2)
write.csv(asv_Obj2_over2abund, file = 'mean RA Obj2.csv')

#Stacked barplot by treatment at the ASV level
barplot_obj2<-ggplot(asv_Obj2_over2abund, aes(x=reorder(Treatment,SampleOrder), y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+facet_grid(.~Facility, scales = "free", space = 'free')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms - Obj 2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_obj2.png", plot=barplot_obj2, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_obj2.svg", plot=barplot_obj2, device="svg", width=13, height=10, units="in", dpi=600)

#Make long format table from Phyloseq object - only postive control
asv_Obj2_PC <- subset_samples(phyloseq_obj2_RA, Facility=="Control") %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment for each ASV
asv_Obj2_mean<-asv_Obj2_PC %>%
  group_by(OTU, Treatment, Facility, Family, Genus, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_Obj2_over2abund_PC <- filter(asv_Obj2_mean, Mean>0.5)


#Stacked barplot by treatment at the ASV level
barplot_obj2_PC<-ggplot(asv_Obj2_over2abund_PC, aes(x=reorder(Treatment,SampleOrder), y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms - Obj 2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_obj2.png", plot=barplot_obj2, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_obj2.svg", plot=barplot_obj2, device="svg", width=13, height=10, units="in", dpi=600)



#### Obj3: Inhibition of L. monocytogenes by lactic acid bacteria in the presence of microbiomes ####
#Split phyloseq by Objective
ps_Obj3<-subset_samples(physeq_16s, Obj == 3)

#Get ASV ad metadata table from phyloseq object
asv_obj3<-as.data.frame(t(otu_table(ps_Obj3)))
meta_obj3<-sample_data(ps_Obj3)

#Remove ASVs with zero counts in all samples
asv_obj3<-asv_obj3[ which(rowSums(asv_obj3)>0),]


# CHECH CONTROLS NAGATIVE CONTROL AND DECONTAMINATE READS
#Make Phyloseq
phyloseqObj3<-phyloseq(otu_table(asv_obj3, taxa_are_rows = TRUE),tax_table(taxon_all), sample_data(meta_obj3))

#Subset Controls
phyloseq_Obj3NC<-subset_samples(phyloseqObj3, SampleID =="NC")
phyloseq_Obj3PC<-subset_samples(phyloseqObj3, SampleID =="PC")
phyloseq_Obj3PCDNA<-subset_samples(phyloseqObj3, SampleID=="PCDNA")

#Make long table
phyloseq_Obj3_NC<-psmelt(phyloseq_Obj3NC)
phyloseq_Obj3_PC<-psmelt(phyloseq_Obj3PC)
phyloseq_Obj3_PCDNA<-psmelt(phyloseq_Obj3PCDNA)

#Remove rows with less than 100 reads
phyloseq_Obj3NC<-subset(phyloseq_Obj3_NC, Abundance>100)
phyloseq_Obj3PC<-subset(phyloseq_Obj3_PC, Abundance>100)
phyloseq_Obj3PCDNA<-subset(phyloseq_Obj3_PCDNA, Abundance>100)

#Plot reads of control by  NC
barplot_Obj3_NC<-ggplot(phyloseq_Obj3NC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black', fill='#00C08D')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Obj3 - NC")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_Obj3_NC.png", plot =barplot_Obj3_NC, device="png", width=20, height=10, units="in",dpi=600)

#Plot reads of control by  PC
barplot_Obj3_PC<-ggplot(phyloseq_Obj3PC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black', fill='#00C08D')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Obj3 - PC")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_Obj3_PC.png", plot =barplot_Obj3_PC, device="png", width=20, height=10, units="in",dpi=600)

#Plot reads of control by  PCDNA
barplot_Obj3_PCDNA<-ggplot(phyloseq_Obj3PCDNA, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=Facility))+
  geom_bar(stat='identity', color='black', fill='#00C08D')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Obj3 - PCDNA")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_Obj3_PCDNA.png", plot =barplot_Obj3_PCDNA, device="png", width=20, height=10, units="in",dpi=600)


#Calculate relative abundance for positive controls
Obj3PC_RA<-phyloseq_Obj3PC%>%
  psmelt()%>%
  arrange(desc(Abundance))

ASV_Obj3_PC<-subset(Obj3PC_RA, OTU=="ASV27" | OTU=="ASV1" | OTU=="ASV36" | OTU=="ASV41" | OTU=="ASV8" | OTU=="ASV74" | OTU=="ASV38" | OTU=="ASV29")
ASV_Obj3_PC$RA<-ASV_Obj3_PC$Abundance/sum(ASV_Obj3_PC$Abundance)*100


Obj3PCDNA_RA<-phyloseq_Obj3PCDNA%>%
  psmelt()%>%
  arrange(desc(Abundance))

ASV_Obj3_PCDNA<-subset(Obj3PCDNA_RA, OTU=="ASV27" | OTU=="ASV1" | OTU=="ASV36" | OTU=="ASV41" | OTU=="ASV8" | OTU=="ASV74" | OTU=="ASV38" | OTU=="ASV29")
ASV_Obj3_PCDNA$RA<-ASV_Obj3_PCDNA$Abundance/sum(ASV_Obj3_PCDNA$Abundance)*100


#Decontaminate sequencing data

#Detect contaminants with prevalence method stringent method (threshold=0.5)
sample_data(phyloseqObj3)$is.neg <- sample_data(phyloseqObj3)$Treatment == "NegControl"
contamdf.prev <- isContaminant(phyloseqObj3, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) 

#Remove identified contaminants (prevalence method) from OTU table
phyloseq_obj3.noncontam <- prune_taxa(!contamdf.prev$contaminant, phyloseqObj3)
phyloseq_obj3.noncontam

#Remove controls from my ASV table
phyloseq_obj3_clean<-subset_samples(phyloseq_obj3.noncontam, Facility != "Control")

#Get ASV table from phyloseq object
asv_obj3_clean<-as.data.frame(otu_table(phyloseq_obj3_clean))
tail(rowSums(asv_obj3_clean))

#Remove ASVs with zero counts in all samples
asv_obj3_clean<-asv_obj3_clean[ which(rowSums(asv_obj3_clean)>0),]

#Get metadata 
metadata_obj3_clean<-subset(metadata_all, Obj==3 & Facility != "Control")


#Compositional analysis of microbiome data -based on Microbiome Analysis in R. Chap 10.
#Step 1: Convert ASV table to appropriate format. Following step requires samples on rows and ASVs in columns
head(t(asv_obj3_clean)) 

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_obj3<-t(cmultRepl(t(asv_obj3_clean), label=0, method="CZM", output="p-counts")) 

head(asv.n0_obj3) #output table needs to have samples in columns and ASVs in rows

#Step 3: Convert data to proportions
asv.n0_obj3_prop<-apply(asv.n0_obj3, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_obj3_prop_f<-asv.n0_obj3[apply(asv.n0_obj3_prop, 1, min) > 0.0000001, ]
head(asv.n0_obj3_prop_f) #Check that samples are on columns and asvs in rows


#Step 5: perform CLR transformation
asv.n0.clr_obj3<-t(apply(asv.n0_obj3_prop_f, 2, function(x){log(x)-mean(log(x))}))
head(asv.n0.clr_obj3) #Check output table. Samples should be in rows and asvs in columns

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_obj3<-prcomp(asv.n0.clr_obj3)

png("Screeplot - PCA - Obj 3.png")
par(mar=c(2,2,2,2))
screeplot(pc.clr_obj3, type='lines', main="Obj 3")
dev.off()

#Calculate total variance of the data
mvar.clr_obj3<-mvar(asv.n0.clr_obj3)

#Display results
row_obj3<-rownames(asv.n0.clr_obj3) #Make vector with sample names
pc_out_obj3<-as.data.frame(pc.clr_obj3$x[,1:2]) #Get PC1 and PC2
pc_out_meta_obj3<-as.data.frame(bind_cols(pc_out_obj3,metadata_obj3_clean)) #Add metadata information
row.names(pc_out_meta_obj3)<-row_obj3 #Add rownames to dataframe

# Make PCA plot
PCA_Obj3 <- ggplot(pc_out_meta_obj3, aes(x=PC1,y=PC2, color=Facility, shape=Treatment))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_obj3$sdev[1]^2/mvar.clr_obj3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_obj3$sdev[2]^2/mvar.clr_obj3*100, digits=1), "%", sep="")) +
  ggtitle("PCA - Obj 3 by facility")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_Obj3
ggsave("PCA_Obj3.png", plot=PCA_Obj3, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_Obj3.svg", plot=PCA_Obj3, device="svg", width=8, height=8, units="in", dpi=600)

# PERMANOVA
#Calculate Aitchinson distance
dist_obj3<-dist(asv.n0.clr_obj3, method='euclidean')

#By Facility
permanova_obj3_fac<-pairwise.adonis2(dist_obj3~Facility, data=metadata_obj3_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_obj3_fac


#By Treatment
permanova_obj3_trt<-pairwise.adonis2(dist_obj3~Treatment, data=metadata_obj3_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_obj3_trt


# STACKED BARPLOT 
#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions
asv.n0.acomp_obj3<-as.data.frame(acomp(t(asv.n0_obj3)), total=1)
rowSums(asv.n0.acomp_obj3) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_obj3_RA <- phyloseq(otu_table(asv.n0.acomp_obj3, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_obj3_clean))

#Make long format table from Phyloseq object - only treatment samples
asv_Obj3_long <- subset_samples(phyloseq_obj3_RA, Facility!="Control") %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment for each ASV
asv_Obj3_mean<-asv_Obj3_long %>%
  group_by(OTU, Treatment, Facility, Family, Genus, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_Obj3_over2abund <- filter(asv_Obj3_mean, Mean>2)
write.csv(asv_Obj3_over2abund, file = 'mean RA Obj3.csv')

#Stacked barplot by treatment at the ASV level
barplot_obj3<-ggplot(asv_Obj3_over2abund, aes(x=reorder(Treatment,SampleOrder), y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+facet_grid(.~Facility, scales = "free", space = 'free')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms - Obj 3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_obj3.png", plot=barplot_obj3, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_obj3.svg", plot=barplot_obj3, device="svg", width=13, height=10, units="in", dpi=600)

#Make long format table from Phyloseq object - only positive controls
asv_Obj3_long_PC <- subset_samples(phyloseq_obj3_RA, Facility=="Control") %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_Obj3_over2abund <- filter(asv_Obj3_long_PC, Abundance>2)


#Stacked barplot by treatment at the ASV level
barplot_obj3_PC<-ggplot(asv_Obj3_over2abund, aes(x=Treatment, y=Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms - Obj 3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_obj3.png", plot=barplot_obj3, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_obj3.svg", plot=barplot_obj3, device="svg", width=13, height=10, units="in", dpi=600)
