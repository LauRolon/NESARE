#DADA2 pipeline applied on NESARE data
#From https://benjjneb.github.io/dada2/tutorial.html

#Last updated: MLR 12/2/21

### Obj1 ####
setwd('/storage/home/mlr355/work/NESARE/3.DADA2/Obj1')

#Install DADA2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

#Attach libraries
library(dada2)

#Make path for fastq files
path <- '/storage/home/mlr355/work/NESARE/3.DADA2/Obj1' 
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality for forward and reverse reads of first thee files
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Trim and filter for Forward at 200bp, for Reverse at 180 - Check based on graph
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]


#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove non-target sequences from the sequence table
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline - check that there is no great drop in reads throughout any of the previous steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/storage/home/mlr355/work/NESARE/3.DADA2/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#Add species column
taxa <- addSpecies(taxa, "/storage/home/mlr355/work/NESARE/3.DADA2/tax/silva_species_assignment_v132.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa.print

#Make phyloseq
library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

#Change name from sequence to ASV#
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#Save ASV table as csv
asvs<-as.data.frame(otu_table(ps))
Taxon<-as.data.frame(tax_table(ps))

write.csv(asvs, file = "ASV.csv")
write.csv(Taxon, file = "Taxon.csv")


### All sequences from project ####
setwd('/storage/home/mlr355/work/NESARE/3.DADA2/All/')

#Attach libraries
library(dada2)

#Make path for fastq files
path_all <- '/storage/home/mlr355/work/NESARE/3.DADA2/All/' 
list.files(path_all)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs_all <- sort(list.files(path_all, pattern="_1.fq", full.names = TRUE))
fnRs_all <- sort(list.files(path_all, pattern="_2.fq", full.names = TRUE))

# Extract sample names
sample.names_all <- sapply(strsplit(basename(fnFs_all), "_"), `[`, 1)

#Inspect read quality for forward and reverse reads of first thee files
plotQualityProfile(fnFs_all[1:40]) #Looks like the samples sequenced for Obj1 had lower sequence quality
plotQualityProfile(fnFs_all[41:93])
plotQualityProfile(fnRs_all[1:40]) #Looks like the samples sequenced for Obj1 had lower sequence quality
plotQualityProfile(fnRs_all[41:93])

#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs_all <- file.path(path_all, "filtered", paste0(sample.names_all, "_F_filt.fastq.gz"))
filtRs_all <- file.path(path_all, "filtered", paste0(sample.names_all, "_R_filt.fastq.gz"))
names(filtFs_all) <- sample.names_all
names(filtRs_all) <- sample.names_all

#Trim and filter for Forward at 200bp, for Reverse at 180 - Check based on graph
out_all <- filterAndTrim(fnFs_all, filtFs_all, fnRs_all, filtRs_all, truncLen=c(210,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out_all

#Learn error rates
errF_all <- learnErrors(filtFs_all, multithread=TRUE)
errR_all <- learnErrors(filtRs_all, multithread=TRUE)
plotErrors(errF_all, nominalQ=TRUE)


#Sample inference
dadaFs_all <- dada(filtFs_all, err=errF_all, multithread=TRUE)
dadaRs_all <- dada(filtRs_all, err=errR_all, multithread=TRUE)

dadaFs_all[[1]]


#Merge paired reads
mergers_all <- mergePairs(dadaFs_all, filtFs_all, dadaRs_all, filtRs_all, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_all[[1]])

#construct sequence table
seqtab_all <- makeSequenceTable(mergers_all)
dim(seqtab_all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_all)))

#Remove non-target sequences from the sequence table
seqtab2_all <- seqtab_all[,nchar(colnames(seqtab_all)) %in% 250:256]

#Remove chimeras
seqtab.nochim_all <- removeBimeraDenovo(seqtab2_all, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_all)
sum(seqtab.nochim_all)/sum(seqtab_all)

#Track reads through the pipeline - check that there is no great drop in reads throughout any of the previous steps
getN <- function(x) sum(getUniques(x))
track_all <- cbind(out_all, sapply(dadaFs_all, getN), sapply(dadaRs_all, getN), sapply(mergers_all, getN), rowSums(seqtab.nochim_all))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_all) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_all) <- sample.names_all
track_all

#Assign taxonomy
taxa_all <- assignTaxonomy(seqtab.nochim_all, "/storage/home/mlr355/work/NESARE/3.DADA2/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#Add species column
taxa_all <- addSpecies(taxa_all, "/storage/home/mlr355/work/NESARE/3.DADA2/tax/silva_species_assignment_v132.fa.gz")
taxa.print_all <- taxa # Removing sequence rownames for display only
rownames(taxa.print_all) <- NULL
taxa.print_all

#Make phyloseq
library(phyloseq)
ps_all <- phyloseq(otu_table(seqtab.nochim_all, taxa_are_rows=FALSE), 
               tax_table(taxa_all))

#Change name from sequence to ASV#
dna_all <- Biostrings::DNAStringSet(taxa_names(ps_all))
names(dna_all) <- taxa_names(ps_all)
ps_all <- merge_phyloseq(ps_all, dna_all)
taxa_names(ps_all) <- paste0("ASV", seq(ntaxa(ps_all)))
ps_all


#Save ASV table as csv
asvs_all<-as.data.frame(otu_table(ps_all))
Taxon_all<-as.data.frame(tax_table(ps_all))

write.csv(asvs_all, file = "ASV_all.csv")
write.csv(Taxon_all, file = "Taxon_all.csv")