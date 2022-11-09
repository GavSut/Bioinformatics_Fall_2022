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
path <- "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Bioinformatics_DADA2_Pipeline/working_dir2/data" 
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)



#### Inspect Read Qulaity Profiles ####
plotQualityProfile(fnFs[1:5]) #helps identify where to trim reads due to lack of quality at end of reads
plotQualityProfile(fnRs[1:2]) #reverse read quality analysis

#### Filter and Trim ####
# Place filtered files in filtered/ subdirectory - creates directory for filtered reads
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#filter and trim
out <- filterAndTrim(fnFs, filtFs, truncLen=120,
                     maxN=1, maxEE=c(8), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#### Learn the Error Rates ####

#function for learning the error rates stored in errF for forward reads
errF <- learnErrors(filtFs, multithread=TRUE)
#function for learning the error rates stored in errR for reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

#plotting errors for the forward reads - as long as trend is following the predicted at all we assume our data is good. 
#severe changes from the red line are to be examined further - something is wrong with the reads
plotErrors(errF, nominalQ=TRUE)
#plotting errors for the reverse reads
plotErrors(errR, nominalQ=TRUE)

#### Dereplication ####
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE) #forward reads dereplicating step
derepRs <- derepFastq(filtRs, verbose=TRUE) #reverse reads dereplicating step
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### Sample Inference ####
#APPLYING THE CORE DADA2 ALGORITHM
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #forward reads
dadaRs <- dada(derepRs, err=errR, multithread=TRUE) #reverse reads
#inspecting results
dadaFs[[1]]
dadaRs[[1]]

#### Merge Paired Reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]]) #shows head of output - check for appearance

#### Construct Sequence Table ####
#constructing an ASV table from the merged reads
seqtab <- makeSequenceTable(dadaRs)
dim(seqtab) #shows the length of shortest and longest sequences in the file
table(nchar(getSequences(seqtab)))

#### Remove Chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #removes the chimeras from our sequence table
dim(seqtab.nochim) #shows the new lengths of the shortest and longest sequences

sum(seqtab.nochim)/sum(seqtab) #from the tutorial "The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on...
#...factors including experimental procedures and sample complexity. Here chimeras make up about 21% of the merged sequence variants, but when we account for the abundances of those variants we see they account for only about 4% of the merged sequence reads.

#### Tracks Reads through the Pipeline ####
#Shows a comparison of the input reads, filtered reads, denoisedF and denoisedR reads, merged reads, and nonchimera reads. 
#allows you to see where you are losing reads to filters/and or if there is a specific sample that isn't meeting parameters as well as others. 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#### Assigning Taxonomy ####
#you can use multiple ways to assign taxonomy to a ASV table, in this case we use the Silva database that can be downloaded locally. 
#make sure to change file pathway to where the silva taxonoy file is. 
taxa <- assignTaxonomy(seqtab.nochim, "/Users/gav/Library/CloudStorage/GoogleDrive-gavin.sutter0@gmail.com/My Drive/2022/HSU/Fall 2022/Bioinformatics/Bioinformatics_DADA2_Pipeline/Tutorial/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#looking at taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL #removes NULLS to be able to better visually look at the data without long lists of NULLS in the results. 
head(taxa.print)

#### Evaluating Accuracy ####
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#### Phyloseq ####
#manual installation of phyloseq - needed because of version update error in current version of R
#manually by downloading the tar file from the Phyloseq github, manually updating the location of the tar file download and unpacking with:
#install.packages("/Users/gavinsutter/Downloads/joey711-phyloseq-9b211d9.tar.gz", repos = NULL, type="source")

library(phyloseq) #load phyloseq
library(ggplot2) #load ggplot2

#setting up a data table - this would look different if not running this tutorial data set. 
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

#creating a phyloseq object out of dada2 data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

#Example of plotting alpha diversity using the Shannon and Simpson Alpha Diversity metrics
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") #curious why there isn't any singletons in this data set, phyloseq shows a warning that this is highly suspicious given the data. 

#Ordinate
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

#Plotting Ordination
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
#clear separation between the early and late samples

#plotting abundances in a bar chart by sample (filtered by the top 20 taxa in each sample)
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")




