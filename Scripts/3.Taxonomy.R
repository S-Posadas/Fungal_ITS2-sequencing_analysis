
# 3.Taxonomy analysis #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(cowplot)
library(grid)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

load("Results/RData/2.phyloseq.filtered.RData")

# Create directories for results
dir.create("Results/3.Taxonomy/Bar_plots", recursive = T)
dir.create("Results/4.Stats", recursive = T)

#### End ####

#### Keep only patient samples and remove tax IDs that are not present ####
PC = c("PC-A-Mock_community", "PC-B-Mock_community", "PC-Candida_albicans", "PC-Saccharomyces_cerevisiae")

physeq = subset_samples(physeq, !sample_names(physeq) %in% PC)
physeq = subset_taxa(physeq, rowSums(otu_table(physeq)) > 0)

physeq_re = subset_samples(physeq_re, !sample_names(physeq_re) %in% c(PC,  "INT-43")) # Remove sample 43 because it has no reads
physeq_re = subset_taxa(physeq_re, rowSums(otu_table(physeq_re)) > 0)

#### End ####

#### Plot relative abundance ####

# Genus

plot_bar(physeq_re, y =  "Abundance", title = "Relative abundance based on Genus", fill="Genus")+ facet_wrap(~sample, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15)) + ylab("Relative abundance") +
  guides(fill = guide_legend("Genus"), color = guide_legend("Genus"))

# Reformat taxa names for plotting until species level

tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "s__Pyronemataceae_sp"] = "f__Pyronemataceae"
tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "o__Hypocreales Order"] = "o__Hypocreales"
tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "o__Pleosporales Order"] = "o__Pleosporales"
tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "s__Inocybe_sp"] = "g__Inocybe"
tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "s__Otidea_sp"] = "g__Otidea"
tax_table(physeq_re)[,"Species"][tax_table(physeq_re)[,"Species"] == "s__Nakaseomyces_sp"] = "g__Nakaseomyces"

# Species

plot_bar(physeq_re, y =  "Abundance", title = "Relative abundance", fill="Species")+ facet_wrap(~sample, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15)) + ylab("Relative abundance") +
  guides(fill = guide_legend("Taxa"), color = guide_legend("Taxa"))

#ggsave("Results/3.Taxonomy/Bar_plots/Figure_2.tiff", dpi=600, compression = 'lzw')

#### Plot absolute abundance ####

# Genus

plot_bar(physeq, y =  "Abundance", title = "Absolute abundance based on Genus", fill="Genus")+ facet_wrap(~sample, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))  +
  guides(fill = guide_legend("Genus"), color = guide_legend("Genus"))

# Reformat taxa names for plotting until species level

tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "s__Pyronemataceae_sp"] = "f__Pyronemataceae"
tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "o__Hypocreales Order"] = "o__Hypocreales"
tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "o__Pleosporales Order"] = "o__Pleosporales"
tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "s__Inocybe_sp"] = "g__Inocybe"
tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "s__Otidea_sp"] = "g__Otidea"
tax_table(physeq)[,"Species"][tax_table(physeq)[,"Species"] == "s__Nakaseomyces_sp"] = "g__Nakaseomyces"

# Species

plot_bar(physeq, y =  "Abundance", title = "Absolute abundance based on Species", fill="Species")+ facet_wrap(~sample, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))  +
  guides(fill = guide_legend("Taxa"), color = guide_legend("Taxa")) 

# Plot absolute from high and low count samples

High = c("INT-04", "INT-05", "INT-11", "INT-13", "INT-36", "INT-47")
Low= c("INT-09", "INT-31", "INT-40", "INT-41", "INT-43", "INT-46")

physeq_high = subset_samples(physeq, sample_names(physeq) %in% High)
physeq_high = subset_taxa(physeq_high, rowSums(otu_table(physeq_high)) > 0)

physeq_low = subset_samples(physeq, sample_names(physeq) %in% Low)
physeq_low = subset_taxa(physeq_low, rowSums(otu_table(physeq_low)) > 0)


H = plot_bar(physeq_high, y =  "Abundance", fill="Species")+ facet_wrap(~sample, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + xlab (NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))  +
  guides(fill = guide_legend("Taxa"), color = guide_legend("Taxa"))

L = plot_bar(physeq_low, y =  "Abundance", fill="Species")+ facet_wrap(~sample, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + xlab (NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))  +
  guides(fill = guide_legend("Taxa"), color = guide_legend("Taxa"))


plot_grid(top = textGrob("Absolute abundance based on Species",gp=gpar(fontsize=20)),
          H, L, 
          nrow=3, rel_heights = c(1,8,8))

#### Plot absolute abundance according to fungal growth ####
# Absolute
# Transform all variables to factors just in case...
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)

ps.culture <- merge_samples(physeq, "Fungal_growth")  # Merging 
sample_data(ps.culture)$Fungal_growth <- rownames(sample_data(ps.culture))

ps.plot = plot_bar(ps.culture, y =  "Abundance", fill="Species") +
  facet_wrap(~Fungal_growth, scales="free_x")

ps.plot + geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  ylab("Absolute abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.03, "npc")) +
  guides(fill = guide_legend("Taxa"), color = guide_legend("Taxa"))

#### End ####

#### Create matrix for correlation ####

# Agglomerate taxa
physeq <- tax_glom(physeq, taxrank = "Species")

counts <- otu_table(physeq)
rownames(counts) <- tax_table(physeq)[,"Species"]

counts <- merge(sample_data(physeq)[,c("ITS2_Read_Nr", "Fungal_ASV_counts")], t(counts), by =0)
counts$ITS2_Read_Nr <- as.numeric(counts$ITS2_Read_Nr)
counts$Fungal_ASV_counts <- as.numeric(counts$Fungal_ASV_counts)

corr_data <- merge(sample_info_tab, counts, all.x = T, by.x=0, by.y = "Row.names")
corr_data[!corr_data$Row.names %in% counts$Row.names, colnames(counts)[-1]] <- 0

rownames(corr_data) <- corr_data$Row.names ; corr_data$Row.names <- NULL

write.xlsx(corr_data, "Results/4.Stats/correlations_tab.xlsx", rowNames = T)


#### End ####

