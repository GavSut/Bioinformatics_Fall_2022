#Bioinformtics Research Project R Script
#Gavin Sutter and Anothony Kieffaber

#installation using biocmanager
#double check the most current version of dada2 and update version = ""
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.11")

#loading packages
library(dada2); packageVersion("dada2") #loads the dada2 package into R session
library(ggplot2)
library(phyloseq)
library(ggfortify)


#data import - after forward and reverse reads have been demultiplexed with cutadpt
path <- "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/trimmed" 
path_filtered = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data"
metadata =  read.csv("/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/arachnid_meta_data.csv", header = TRUE)
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



#### Inspect Read Qulaity Profiles ####
plotQualityProfile(fnFs[1:2]) #helps identify where to trim reads due to lack of quality at end of reads
#plotQualityProfile(fnRs[1:2]) #reverse read quality analysis
ggsave("quality_scores_forward_reads.jpeg", dpi = 300, path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs")

#### Filter and Trim ####
# Place filtered files in filtered/ subdirectory - creates directory for filtered reads
filtFs <- file.path(path_filtered, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(path_filtered, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trim
out <- filterAndTrim(fnFs, filtFs, truncLen=15,
                     maxN=1, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 15) # On Windows set multithread=FALSE
head(out)

#### Learn the Error Rates ####

#function for learning the error rates stored in errF for forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
#function for learning the error rates stored in errR for reverse reads
#errR <- learnErrors(filtRs, multithread=TRUE)

#plotting errors for the forward reads - as long as trend is following the predicted at all we assume our data is good. 
#severe changes from the red line are to be examined further - something is wrong with the reads
plotErrors(errF, nominalQ=TRUE)
#plotting errors for the reverse reads
#plotErrors(errR, nominalQ=TRUE)

#### Dereplication ####
#Dereplication combines all identical sequencing reads into into ???unique sequences??? with a corresponding ???abundance??? equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE) #forward reads dereplicating step
#derepRs <- derepFastq(filtRs, verbose=TRUE) #reverse reads dereplicating step
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
#names(derepRs) <- sample.names

#### Sample Inference ####
#APPLYING THE CORE DADA2 ALGORITHM
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32) #forward reads
#dadaRs <- dada(derepRs, err=errR, multithread=TRUE) #reverse reads
#inspecting results
dadaFs[[1]]
#dadaRs[[1]]

#### Merge Paired Reads ####
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]]) #shows head of output - check for appearance

#### Construct Sequence Table ####
#constructing an ASV table from the merged reads
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab) #shows the length of shortest and longest sequences in the file
table(nchar(getSequences(seqtab)))

#### Remove Chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #removes the chimeras from our sequence table
dim(seqtab.nochim) #shows the new lengths of the shortest and longest sequences

sum(seqtab.nochim)/sum(seqtab) #from the tutorial "The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on...

#### Tracks Reads through the Pipeline ####
#Shows a comparison of the input reads, filtered reads, denoisedF and denoisedR reads, merged reads, and nonchimera reads. 
#allows you to see where you are losing reads to filters/and or if there is a specific sample that isn't meeting parameters as well as others. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
track
#shows the range of the values that passed the nochim step - allows for removing allow read samples
summary(track[,4])

#removing samples with less than 1000 reads
st <- seqtab[rowSums(seqtab.nochim) >= 1000,]

#### Assigning Taxonomy ####
#you can use multiple ways to assign taxonomy to a ASV table, in this case we use the Silva database that can be downloaded locally. 
#make sure to change file pathway to where the silva taxonoy file is. 
taxa <- assignTaxonomy(st, "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#looking at taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL #removes NULLS to be able to better visually look at the data without long lists of NULLS in the results. 
head(taxa.print)

#saving output
saveRDS(taxa.print, "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/R Objects/taxa.print.RDS")
saveRDS(seqtab.nochim, "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/R Objects/seqtab.nochim.RDS")

#### Phyloseq ####
#manual installation of phyloseq - needed because of version update error in current version of R
#manually by downloading the tar file from the Phyloseq github, manually updating the location of the tar file download and unpacking with:
#install.packages("/Users/gav/Downloads/joey711-phyloseq-9b211d9.tar.gz", repos = NULL, type="source")

library(phyloseq) #load phyloseq
library(ggplot2) #load ggplot2
library(vegan)
library(dplyr)#load vegan

#setting up a data table - this would look different if not running this tutorial data set. 
samples.out <- rownames(st)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
organism = metadata$taxa_id
spider = metadata$spider
subject <- substr(subject,1,8)
samdf <- data.frame(Subject=subject, Organism=organism, Spider=spider)
rownames(samdf) <- samples.out
samdf
organism
spider
st
#creating a phyloseq object out of dada2 data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps
#ps <- prune_samples(sample_sums(ps) >= 1000, ps)

#Example of plotting alpha diversity using the Shannon and Simpson Alpha Diversity metrics
plot_richness(ps, x="Spider", measures=c("Shannon", "Simpson")) #curious why there isn't any singletons in this data set, phyloseq shows a warning that this is highly suspicious given the data. 
plot_richness(ps, x="Organism", measures=c("Shannon", "Simpson"), color = "Spider")+ 
              theme(panel.background = element_rect(
                fill = "white",
                colour = "black",
                linewidth = 1, 
                linetype = 1), 
                panel.grid.major.y = element_line(color = "grey90",
                                                  size = 0.5,
                                                  linetype = 3), 
                panel.grid.major.x = element_line(color = "grey90",
                                                  size = 0.5,
                                                  linetype = 3)) #curious why there isn't any singletons in this data set, phyloseq shows a warning that this is highly suspicious given the data. 
ggsave("shannon_simpson_div_by_sample.jpeg", path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs", dpi = 300)

#diversity by spider
plot_richness(ps, x="Spider", measures=c("Shannon", "Simpson"), color = "Spider")+ 
  theme(panel.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 1, 
    linetype = 1), 
    panel.grid.major.y = element_line(color = "grey90",
                                      size = 0.5,
                                      linetype = 3), 
    panel.grid.major.x = element_line(color = "grey90",
                                      size = 0.5,
                                      linetype = 3))


#shannon_div_samples = estimate_richness(ps, split = TRUE, measures = "Shannon")
#shannon_div_samples
#metadata
#shannon_div = left_join(samdf,shannon_div_samples, by = c("Subject" == "experiment_accession"))

shannon_div = read.csv("/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/data/shannon_spider_data.csv", header = TRUE)
shannon_div

ggplot(data = shannon_div, aes(spider, shannon))+
  geom_boxplot()+
  xlab("")+
  ylab("Shannon Index")+
  theme(axis.text.x = element_text(size = 20),
        text = element_text(size = 20),
        panel.background = element_rect(fill = 'white', color = 'black', linetype = "solid", linewidth = 1))
ggsave("shannon_diversity.jpeg", path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs", dpi = 300 )
        

#Ordinate
#ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu)) - normalizing is causing all ordination values to be incredibly similar. 
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

#ordination with vegan\
# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
ps_vegan = psotu2veg(ps)
ord = metaMDS(ps_vegan, distance = "bray")
plot_ordination(ps, ord, color="Spider", title="Bray NMDS")+
  theme_classic()
ord
adonis2(ord~samples, data = ord)


#Plotting Ordination
plot_ordination(ps, ord.nmds.bray, color="Spider", title="Bray NMDS")+
  theme(axis.text.x = element_text(size = 20),
        text = element_text(size = 20),
        panel.background = element_rect(fill = 'white', color = 'white'), 
        axis.line.x = element_line(color = 'black'), 
        axis.line.y = element_line(color = 'black'))
ggsave("ord_by_spider.jpeg", path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs", dpi = 900)







#plotting abundances in a bar chart by sample (filtered by the top 20 taxa in each sample)
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Organism", fill="Phylum")+
  theme(panel.background = element_rect(fill = 'white', color = 'black', linetype = "solid", linewidth = 1))
ggsave("taxa_sample_family.jpeg", path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs", dpi = 300)
ggsave("taxa_sample_family_legend.jpeg", path = "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Research_Project/working_dir/scripts/Outputs", dpi = 1200)

ps.top20

plot_bar(ps.top20, x="Spider", fill="Phylum") 

ggplot(ps.top20, aes( x= "Spider", Fill = "Family"))+
  geom_bar()

#calculating shannon diversity for each sample in datatable - needs work
#samples_shannon =  diversity("Genus", index = "shannon", exp = FALSE) 

#more ordinations 
ps.ord <- ordinate(ps, "NMDS", "bray")
p1 = plot_ordination(ps, ps.ord, type="taxa", color="Class", title="taxa")
print(p1)

