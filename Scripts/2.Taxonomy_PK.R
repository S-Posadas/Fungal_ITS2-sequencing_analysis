
# 2.Taxonomy analysis from positive control #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(openxlsx)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

#load("Results/RData/2.phyloseq.filtered.RData")
load("Results/RData/2.phyloseq.filtered.RData")

# Create directories for results
dir.create("Results/3.Taxonomy/PK", recursive = T)

#### End ####

#### Keep only controls and remove tax IDs that are not present ####
PC = c("PC-A-Mock_community", "PC-B-Mock_community", "PC-Candida_albicans", "PC-Saccharomyces_cerevisiae")

physeq_K = subset_samples(physeq, sample_names(physeq) %in% c(PC))
physeq_K = subset_taxa(physeq_K, rowSums(otu_table(physeq_K)) > 0)

physeq_K_re = subset_samples(physeq_re, sample_names(physeq_re) %in% c(PC))
physeq_K_re = subset_taxa(physeq_K_re, rowSums(otu_table(physeq_K_re)) > 0)

#Select only sample data for controls

sample_data(physeq_K)$sample <- gsub("PC-", "", sample_data(physeq_K)$sample)
sample_data(physeq_K)$Control_name <- gsub("A-|B-", "", sample_data(physeq_K)$sample)

sample_data(physeq_K_re)$sample <- gsub("PC-", "", sample_data(physeq_K_re)$sample)
sample_data(physeq_K_re)$Control_name <- gsub("A-|B-", "", sample_data(physeq_K_re)$sample)

#### End ####

#### Bar plots ####

# Absolute abundance

plot_bar(physeq_K, y =  "Abundance", title = "Abundance based on Genus", fill="Genus")+ facet_wrap(~Control_name, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot_bar(physeq_K, y =  "Abundance", title = "Abundance based on Species", fill="Species")+ facet_wrap(~Control_name, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Relative abundance

# Agglomerate at Species level
physeq_K_re <- tax_glom(physeq_K_re, taxrank = "Species")

# Remove taxa in each sample with less than 3% 

# Get the ASV table
asv_table <- otu_table(physeq_K_re)

# Calculate the total read count for each sample
sample_total_counts <- colSums(asv_table)

# Iterate through each sample
for (sample_name in sample_names(physeq_K_re)) {
  # Get the counts for the current sample
  sample_counts <- asv_table[, sample_name]

  # Identify Species with counts less than the threshold
  asvs_to_zero <- sample_counts < 0.03
  
  # Set the counts of those Species to 0 for the current sample
  sample_counts[asvs_to_zero] <- 0
  
  # Update the counts for the current sample in the phyloseq object
  asv_table[, sample_name] <- sample_counts
}


# Update the phyloseq object with the modified ASV counts
otu_table(physeq_K_re) <- asv_table
physeq_K_re = subset_taxa(physeq_K_re, rowSums(otu_table(physeq_K_re)) > 0)

plot_bar(physeq_K_re, y =  "Abundance", title = "Relative abundance based on Species", fill="Species")+ facet_wrap(~Control_name, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL)


#### End ####
