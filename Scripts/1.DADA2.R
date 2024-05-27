####Introduction####

# This code was used for the analysis of ITS2 gene sequencing of clinical ascites samples
# from the INTeGRATE clinical study. The analysis starts from raw Illumina reads where used
# primers get trimmed using cutadapt and further analysis is carried out using DADA2 to
# end up with amplicon sequence variant (ASV) table. Taxonomy analysis combines DADA2 and
# Blast approaches. Phyloseq is used for downstream analysis.

#### End ####

#### Package setup ####

library(dada2); packageVersion("dada2")
library(ShortRead)
library(Biostrings)
library(vegan)
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(openxlsx)
library(microViz)
theme_set(theme_bw())

#### End ####

#### Set Directory ####
setwd("path/to/working/directory") 
path = "Fastq" # Change to the directory containing the fastq files
list.files(path)

# Create directory to save results of the analysis and subdirectory for quality control

dir.create("Results/1.QC", recursive = T)

#### End ####

#### Quality Control and "Pre-Filtering" ####
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz

fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

plotQualityProfile(fnFs)

plotQualityProfile(fnRs)

# Identify primers

FWD <- "GTGARTCATCGARTCTTTG"   ## Change to your forward primer sequence
REV <- "TTCCTSCGCTTATTGATATGC" ## Change to your reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# “Pre-filter” the sequences to remove those with Ns, but perform no other filtering.

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, truncQ=0, multithread = F)

# Count the number of times the primers appear in the forward and reverse read
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fn=fnFs.filtN[[1]]),
      FWD.ReverseReads=sapply(FWD.orients,primerHits,fn=fnRs.filtN[[1]]),
      REV.ForwardReads=sapply(REV.orients,primerHits,fn=fnFs.filtN[[1]]),
      REV.ReverseReads=sapply(REV.orients,primerHits,fn=fnRs.filtN[[1]]))

#### End ####

#### Remove primers with cutadapt ####
# Performed in Ubuntu environment. See script 1.1.Primers_trim

#### Check that primers have been removed ####
# Forward and reverse fastq filenames have the format:
path.cut <- file.path(path, "filt_cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Count the number of times the primers appear in the forward and reverse read
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#### End ####

#### Quality Filter ####
fnFs.filtQ <- file.path(path, "filtQ", basename(fnFs.cut))
fnRs.filtQ <- file.path(path, "filtQ", basename(fnRs.cut))

out2 <- filterAndTrim(fnFs.cut,fnFs.filtQ,fnRs.cut,fnRs.filtQ,minLen=20,maxEE=8,truncQ=8,compress=TRUE,multithread=F)  

#### End ####

####Learn the Error Rates####

errF <- learnErrors(fnFs.filtQ,multithread=F)
errR <- learnErrors(fnRs.filtQ,multithread=F)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#### End ####

#### Dada ####
# Apply the core sample inference algorithm to the the filtered and trimmed sequence data.

dadaFs <- dada(fnFs.filtQ, err=errF, multithread=FALSE)
dadaRs <- dada(fnRs.filtQ, err=errR, multithread=FALSE)

# Inspecting the returned dada-class object

dadaFs[[1]]
dadaRs[[1]]

#Merge paired reads
# Some samples lost more than 50% of the reads after merging with default min overlap,
# Adjusting it to 8 did not improve the merging in the most affected samples

mergers <- mergePairs(dadaFs, fnFs.filtQ, dadaRs, fnRs.filtQ, verbose=TRUE)

#Construct sequence table for mergers and only for forward reads

seqtab <- makeSequenceTable(mergers)  
seqtabF <- makeSequenceTable(dadaFs)  

dim(seqtab)
dim(seqtabF)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtabF)))

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
seqtab.nochimF <- removeBimeraDenovo(seqtabF, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)
dim(seqtab.nochimF)

# Percentage of not chimeric reads

sum(seqtab.nochim)/sum(seqtab)
sum(seqtab.nochimF)/sum(seqtabF)

#### End ####

#### Track Pipeline ####
#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, out2[,2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
trackF <- cbind(out, out2[,2], sapply(dadaFs, getN), rowSums(seqtab.nochimF))

colnames(track) <- c("input", "filteredN", "filteredQ", "denoisedF", "denoisedR", "merged", "nonchim")
colnames(trackF) <- c("input", "filteredN", "filteredQ", "denoisedF", "nonchimF")

head(track)
head(trackF)

# Save the track

write.csv(track, file = "Results/1.QC/Pipeline_Track_final_merged_reads.csv")
write.csv(trackF, file = "Results/1.QC/Pipeline_Track_final_forward_reads.csv")

#### End ####

#### Proceed only with forward reads ####
# Further steps were carried out only on the forward reads to preserve the read number

#### Make Count and Fasta Tables ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochimF)

asv_headers <- vector(dim(seqtab.nochimF)[2], mode="character")

for (i in 1:dim(seqtab.nochimF)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab <- c(rbind(asv_headers, asv_seqs))
write(fasta_tab, "Results/fasta_tab.fa")

# count table:

count_tab <- t(seqtab.nochimF)
row.names(count_tab) <- sub(">", "", asv_headers)
colnames(count_tab) <- gsub("_1.fastq.gz", "", colnames(count_tab))
write.table(count_tab, "Results/count_tab.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

#### Assign taxonomy with blast ####
# Performed in Ubuntu environment. See script 1.2.Taxonomy.assignment and 1.3.Filter_blast

# Read taxa table created with blast:
tax_tab <- as.matrix(read.xlsx("Results/2.blast/blast_final_95_nt_unite_customized.xlsx", rowNames = T))
tax_tab <- tax_tab[,c("Kingdom_unite", "Phylum_unite", "Class_unite", "Order_unite", "Family_unite", "Genus_custom", "Species_custom")]
colnames(tax_tab) <- gsub("_unite|_custom", "", colnames(tax_tab))
  
# Add ASVs that are not present, for phyloseq object
nASV <- rownames(count_tab)[!rownames(count_tab) %in% rownames(tax_tab)]
tax_tab <- rbind(tax_tab, matrix(, nrow = length(nASV), ncol = dim(tax_tab)[2], dimnames = list(nASV,colnames(tax_tab))))

#### End ####

#### Add samples metadata and match to count table####

sample_info_tab <- read.xlsx("00-INTeGRATE_Fungome_paper_alldata_30102023.xlsx", sheet = "0-metadata", rowNames = T)
sample_info_tab$sample = rownames(sample_info_tab)

# Examine all tables
head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

# Examine consistancy in order between count_tab colnames and coldata rownames 

all(colnames(count_tab) %in% rownames(sample_info_tab))

all(rownames(tax_tab) %in% rownames(count_tab))
all(rownames(count_tab) %in% rownames(tax_tab))

gplots::venn(list(taxonomy=rownames(tax_tab), featuretable=rownames(count_tab)))

#### End ####

#### Making our phyloseq object ####
# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

physeq <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                   sample_data(sample_info_tab), 
                   tax_table(tax_tab))

#### End ####

#### Count total reads ####

sample_sums(physeq) #Nr of reads per Sample
sample_data(physeq)$ITS2_Read_Nr <- sample_sums(physeq) 

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

dir.create("Results/RData")

save.image("Results/RData/1.physeq.original.RData")

#### End ####

## Filter low Abundance Taxa and count table normalization ##

#### Remove NA Phyla ####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum))

table(tax_table(physeq_oNA)[, "Phylum"], exclude = NULL)

physeq 
physeq_oNA
physeq = physeq_oNA

#### End ####

####Define prevalence of each taxa (in how many samples did each taxa appear at least once) ####

prev0 = apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))

#save ASV Prevalence and Abundance table before filtering

write.table(prevdf, "Results/asv_prevdf.tsv", sep="\t", quote=F, col.names=NA)

#Plot Taxa prevalence v. total counts. Each point is a different taxa. 

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Species)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Species) + theme(legend.position="none")

#### End ####

#### Get reads and ASVs corresponding to fungi ####

length(taxa_sums(physeq))

track_fungal <- cbind (trackF, sample_sums(physeq))
colnames(track_fungal) <- c(colnames(trackF),"fungal")
write.xlsx(as.data.frame(track_fungal), "Results/1.QC/Pipeline_track_fungal.xlsx", rowNames = T)

sample_data(physeq)$Fungal_ASV_counts <- sample_sums(physeq) 

#### End ####

#### Remove taxa not seen more than 10 times in at least 1 sample #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperasv = 10

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperasv) >= 1, TRUE)

colSums(otu_table(physeq))
colSums(otu_table(physeq_filtered))
 
physeq
physeq_filtered

#### End ####

#### Remove taxa in each sample with less than 0.1% ####

# # Define a threshold for filtering ASVs 
threshold <- 0.001

# Get the ASV table
asv_table <- otu_table(physeq_filtered)

# Calculate the total read count for each sample
sample_total_counts <- colSums(asv_table)

# Iterate through each sample
for (sample_name in sample_names(physeq_filtered)) {
  # Get the ASV counts for the current sample
  sample_counts <- asv_table[, sample_name]
  
  # Calculate the total count for the current sample
  total_count <- sample_total_counts[sample_name]
  
  # Identify ASVs with counts less than the threshold
  asvs_to_zero <- sample_counts / total_count < threshold
  
  # Set the counts of those ASVs to 0 for the current sample
  sample_counts[asvs_to_zero] <- 0
  
  # Update the ASV counts for the current sample in the phyloseq object
  asv_table[, sample_name] <- sample_counts
}


# Update the phyloseq object with the modified ASV counts
otu_table(physeq_filtered) <- asv_table

physeq = physeq_filtered

#### End ####

#### Assign upper taxonomic ranks to unassigned taxa ####

physeq <- tax_fix(physeq)

#### End ####

### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq, standf)

# Transform to relative abundance. Save as new object.

physeq_re = transform_sample_counts(physeq, function(x){x / sum(x)})

#### End ####

#### Exploratory plots after filtering and normalization ####
# Check individual phylum Abundance
# Abundance value transformation function

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Species",
                          Color = "Genus",
                          n = NULL){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "category_fungi", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
    scale_y_log10()
}


#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####


#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image("Results/RData/2.phyloseq.filtered.RData")

#### End ####
