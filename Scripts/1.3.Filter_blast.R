library(dplyr)
library(openxlsx)

#### Format blast table to get a single taxa per ASV UNITE ####
# Read alignment table from blast and arrange column names
blastu <- read.table("Results/2.blast/tax_tab_95_10align_unite", row.names = NULL)
colnames(blastu) = blastu[1,] ; blastu = blastu[blastu$seqid != "seqid",]

# Split the stitle column using "|" and ";"
blastu_split <- data.frame(do.call(rbind, strsplit(as.character(blastu$stitle), "[|;]")))

# Rename the columns
colnames(blastu_split) <- c("ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "X9")

# Combine the original dataframe and the split columns
blastu <- cbind(blastu, blastu_split[,2:8])

# Remove unnecessary columns
blastu$sseqid <- NULL ; blastu$stitle <- NULL 

# Make sure pident and evalue are numeric
blastu$pident = as.numeric(blastu$pident)
blastu$evalue = as.numeric(blastu$evalue)
head(blastu)
dim(blastu)

# Filter alingments with percentage of identical matches higher than n 
# We already used 95% in blast, but just in case we want to be more strict
# In our case we are not really filtering
blastu = blastu[blastu$pident >= 95,] 
dim(blastu)

# Check how many unique genera were obtained
unique(blastu$Genus)

# Check how many unique genera were obtained
unique(blastu$Genus)

# Keep alingments with the highest percentage of identical matches #and lowest e value
# For ASVs with alignment to a single species, keep only one row
filt = list()
for(asv in unique(blastu$seqid)){
  ASV_n = blastu[blastu$seqid == asv,]
  ASV_n = ASV_n[ASV_n$pident == max(ASV_n$pident),] # highest identity
  ASV_n = ASV_n[ASV_n$evalue == min(ASV_n$evalue),] # lowest e value
  if(length(unique(ASV_n$Genus)) == 1 & length(unique(ASV_n$Species)) == 1){
    ASV_n = ASV_n[1,]
  }
  filt[[asv]] = ASV_n
}
filtu <- as.data.frame(do.call(rbind, filt)) 
dim(filtu)

# Check if Genus is not unique within each ASV
non_unique_genusu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)
length(unique(non_unique_genusu$seqid)) 
# "ASV_2"   "ASV_14"  "ASV_23"  "ASV_36"  "ASV_43"  "ASV_52"  "ASV_53"  "ASV_54"  "ASV_59"  "ASV_64"  "ASV_69"  "ASV_71"  "ASV_101"

# Manually revise not unique ASVs and keep upper taxon which is unique
filtu[filtu$seqid == "ASV_2",]$Family = NA ; filtu[filtu$seqid == "ASV_2",]$Genus = NA ; filtu[filtu$seqid == "ASV_2",]$Species = NA 
filtu[filtu$seqid == "ASV_14",]$Family = NA ; filtu[filtu$seqid == "ASV_14",]$Genus = NA ; filtu[filtu$seqid == "ASV_14",]$Species = NA 
filtu[filtu$seqid == "ASV_23", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_36",]$Family = NA ; filtu[filtu$seqid == "ASV_36",]$Genus = NA ; filtu[filtu$seqid == "ASV_36",]$Species = NA 
filtu[filtu$seqid == "ASV_43", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_52",]$Family = NA ; filtu[filtu$seqid == "ASV_52",]$Genus = NA ; filtu[filtu$seqid == "ASV_52",]$Species = NA 
filtu[filtu$seqid == "ASV_53", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_54", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_59", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_64", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_69", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_71", c("Phylum", "Class", "Order", "Family", "Genus", "Species")] = NA
filtu[filtu$seqid == "ASV_101",] = filtu[filtu$seqid == "ASV_101",][1,]

# Check if Genus is now unique within each ASV
non_unique_genusu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)
length(unique(non_unique_genusu$seqid)) #now it should be 0

# Check if Species is not unique within each ASV
non_unique_speciesu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)
length(unique(non_unique_speciesu$seqid)) 
# "ASV_6"   "ASV_11"  "ASV_27"  "ASV_51"  "ASV_60"  "ASV_75"  "ASV_86"  "ASV_90"  "ASV_102" "ASV_108" "ASV_122"

# Manually revise each ASV
filtu[filtu$seqid == "ASV_6",] = filtu[filtu$seqid == "ASV_6",][1,]
filtu[filtu$seqid == "ASV_11",]$Species = NA # gattii DD neoformans
filtu[filtu$seqid == "ASV_27",] = filtu[filtu$seqid == "ASV_27",][3,]
filtu[filtu$seqid == "ASV_51",] = filtu[filtu$seqid == "ASV_51",][1,]
filtu[filtu$seqid == "ASV_60",] = filtu[filtu$seqid == "ASV_60",][1,]
filtu[filtu$seqid == "ASV_75",] = filtu[filtu$seqid == "ASV_75",][1,]
filtu[filtu$seqid == "ASV_86",] = filtu[filtu$seqid == "ASV_86",][1,]
filtu[filtu$seqid == "ASV_90",] = filtu[filtu$seqid == "ASV_90",][1,]
filtu[filtu$seqid == "ASV_102",] = filtu[filtu$seqid == "ASV_102",][1,]
filtu[filtu$seqid == "ASV_108",] = filtu[filtu$seqid == "ASV_108",][1,] # s__flavisporus
filtu[filtu$seqid == "ASV_122",] = filtu[filtu$seqid == "ASV_122",][3,] # s__hirsutum DD s__gausapatum

# Check if Species is not unique within each ASV
non_unique_speciesu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)
length(unique(non_unique_speciesu$seqid)) #now it should be 0

# Keep only one row per ASV
length(unique(filtu$seqid)) # Check how many unique ASVs are there

final = list()
for(asv in unique(filtu$seqid)){
  ASV_n = filtu[filtu$seqid == asv,]
  ASV_n = ASV_n[1,]
  final[[asv]] = ASV_n
}
finalu <- as.data.frame(do.call(rbind, final)) 
dim(finalu)

min(as.numeric(finalu$pident))

#Save data in desired format
write.table(finalu, "Results/2.blast/blast_final_95_unite", sep = "\t", row.names = FALSE)
write.xlsx(finalu, "Results/2.blast/blast_final_95_unite.xlsx")

#### End ####

#### Format blast table to get a single taxa per ASV NT ####
# Read alignment table from blast and arrange column names
blastnt <- read.table("Results/2.blast/tax_tab_95_10align_nt", row.names = NULL)
colnames(blastnt) = c(colnames(blastnt)[2:9],"Genus", "Species")
head(blastnt)
length(unique(blastnt$seqid))

# Filter alingments with percentage of identical matches higher than n 
# We already used 95% in blast, but just in case we want to be more strict
# In our case we are not really filtering
blastnt = blastnt[blastnt$pident >= 95,] 
dim(blastnt)

# Check how many unique genera were obtained
unique(blastnt$Genus)

# Unify taxonomy
blastnt[blastnt$Genus %in% c("homo","Human"),]$Genus = "Homo"
blastnt[blastnt$Species == "DNA",]$Species = "sapiens"
blastnt[blastnt$Genus == "uncultured",]$Genus = "Uncultured"
blastnt[blastnt$Genus == "[Candida]",]$Genus = "Candida"

# Check how many unique genera were obtained
unique(blastnt$Genus)

# Keep alingments with the highest percentage of identical matches and lowest e value
# For ASVs with alignment to a single species, keep only one row
filt = list()
for(asv in unique(blastnt$seqid)){
  ASV_n = blastnt[blastnt$seqid == asv,]
  ASV_n = ASV_n[ASV_n$pident == max(ASV_n$pident),] # highest identity
  ASV_n = ASV_n[ASV_n$evalue == min(ASV_n$evalue),] # lowest e value
  if(length(unique(ASV_n$Genus)) == 1 & length(unique(ASV_n$Species)) == 1){
    ASV_n = ASV_n[1,]
  }
  filt[[asv]] = ASV_n
}
filtnt <- as.data.frame(do.call(rbind, filt)) 
dim(filtnt)

# The following steps are carried out to have a better overview in the table, since we are not interested
# in any DNA that is not fungal, and some human DNA can align to gorilla DNA
# Change gorilla and monkeys dna to human -> supervise changes!!
filtnt[filtnt$seqid %in% filtnt[filtnt$Genus %in% c("Pongo","PREDICTED:"),"seqid"],]$Species = "sapiens"
filtnt[filtnt$seqid %in% filtnt[filtnt$Genus %in% c("Pongo","PREDICTED:"),"seqid"],]$Genus = "Homo"

# Change eukaryotic synthetik dna to human or NA -> supervise changes!!
filtnt[filtnt$seqid %in% filtnt[filtnt$Genus == "Eukaryotic","seqid"],]$Species = "sapiens"
filtnt[filtnt$seqid %in% filtnt[filtnt$Genus == "Eukaryotic","seqid"],]$Genus = "Homo"

# Check if Genus is not unique within each ASV
non_unique_genusnt <- filtnt %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)
length(unique(non_unique_genusnt$seqid)) 

# Manually revise each ASV
filtnt[filtnt$seqid == "ASV_2",c("Genus", "Species")] = NA #o__Pleosporales
filtnt[filtnt$seqid == "ASV_6",]$Genus = "Aureobasidium" ; filtnt[filtnt$seqid == "ASV_6",]$Species = NA # pullulans DD proteae
filtnt[filtnt$seqid == "ASV_14",c("Genus", "Species")] = NA #o__Hypocreales
filtnt[filtnt$seqid == "ASV_15", ] = filtnt[filtnt$seqid == "ASV_15",][1,]
filtnt[filtnt$seqid == "ASV_19", ]$Genus = "Nakaseomyces" ; filtnt[filtnt$seqid == "ASV_19", ]$Species = NA
filtnt[filtnt$seqid == "ASV_21", ]$Genus = "Nakaseomyces" ; filtnt[filtnt$seqid == "ASV_21", ]$Species = NA
filtnt[filtnt$seqid == "ASV_25", ]$Genus = "Nakaseomyces" ; filtnt[filtnt$seqid == "ASV_25", ]$Species = NA
filtnt[filtnt$seqid == "ASV_27", ] = filtnt[filtnt$seqid == "ASV_27",][1,]
filtnt[filtnt$seqid == "ASV_29", ] = filtnt[filtnt$seqid == "ASV_29",][2,]
filtnt[filtnt$seqid == "ASV_30", ] = filtnt[filtnt$seqid == "ASV_30",][1,]
filtnt[filtnt$seqid == "ASV_31", ] = filtnt[filtnt$seqid == "ASV_31",][1,]
filtnt[filtnt$seqid == "ASV_43",c("Genus", "Species")] = NA #o__Pezizales DDJuglans DD Dothideomycetes
filtnt[filtnt$seqid == "ASV_47",c("Genus", "Species")] = NA # DDJuglans DD Dothideomycetes
filtnt[filtnt$seqid == "ASV_51",]$Genus = "Aureobasidium" ; filtnt[filtnt$seqid == "ASV_51",]$Species = NA # pullulans DD proteae
filtnt[filtnt$seqid == "ASV_52",c("Genus", "Species")] = NA #o__Pleosporales
filtnt[filtnt$seqid == "ASV_53",c("Genus", "Species")] = NA # DD Spinacia  DD fungus
filtnt[filtnt$seqid == "ASV_54",c("Genus", "Species")] = NA # DD Rubus
filtnt[filtnt$seqid == "ASV_56",c("Genus", "Species")] = NA # DD Spinacia  DD fungus
filtnt[filtnt$seqid == "ASV_60",]$Genus = "Aureobasidium" ; filtnt[filtnt$seqid == "ASV_60",]$Species = NA # pullulans DD proteae
filtnt[filtnt$seqid == "ASV_66",c("Genus", "Species")] = NA # DD Rubus
filtnt[filtnt$seqid == "ASV_74", ] = filtnt[filtnt$seqid == "ASV_74",][1,]
filtnt[filtnt$seqid == "ASV_76", ] = filtnt[filtnt$seqid == "ASV_76",][1,]
filtnt[filtnt$seqid == "ASV_80",c("Genus", "Species")] = NA # DD Vaccinium    
filtnt[filtnt$seqid == "ASV_83",c("Genus", "Species")] = NA # DD Magnoliophyta DD Coriandrum       
filtnt[filtnt$seqid == "ASV_88", ] = filtnt[filtnt$seqid == "ASV_88",][1,]
filtnt[filtnt$seqid == "ASV_90", ] = filtnt[filtnt$seqid == "ASV_90",][1,]
filtnt[filtnt$seqid == "ASV_93", ] = filtnt[filtnt$seqid == "ASV_93",][2,]
filtnt[filtnt$seqid == "ASV_100",c("Genus", "Species")] = NA # DD Daucus        
filtnt[filtnt$seqid == "ASV_101",] = filtnt[filtnt$seqid == "ASV_101",][1,]
filtnt[filtnt$seqid == "ASV_102", ] = filtnt[filtnt$seqid == "ASV_102",][1,]
filtnt[filtnt$seqid == "ASV_103",c("Genus", "Species")] = NA # DD Petroselinum   
filtnt[filtnt$seqid == "ASV_108", ]$Genus = "Xylodon" ; filtnt[filtnt$seqid == "ASV_108",]$Species = NA # flaviporus DD ovisporus
filtnt[filtnt$seqid == "ASV_110",c("Genus", "Species")] = NA # DD Rubus
filtnt[filtnt$seqid == "ASV_112",c("Genus", "Species")] = NA # DD Triticum DD Helminthosporium
filtnt[filtnt$seqid == "ASV_113",c("Genus", "Species")] = NA # DD Epicoccum DD Glonium
filtnt[filtnt$seqid == "ASV_122", ] = filtnt[filtnt$seqid == "ASV_122",][1,]

# Check if Genus is not unique within each ASV
non_unique_genusnt <- filtnt %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)
length(unique(non_unique_genusnt$seqid)) #now it should be 0

# Check if Species is not unique within each ASV
non_unique_speciesnt <- filtnt %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)
length(unique(non_unique_speciesnt$seqid)) 
# "ASV_49"  "ASV_55"  "ASV_72"  "ASV_75"  "ASV_114"

# Manually revise each ASV
filtnt[filtnt$seqid == "ASV_49",] = filtnt[filtnt$seqid == "ASV_49",][1,]
filtnt[filtnt$seqid == "ASV_55",] = filtnt[filtnt$seqid == "ASV_55",][1,]
filtnt[filtnt$seqid == "ASV_72",] = filtnt[filtnt$seqid == "ASV_72",][1,] # albicans DD africana in unite albicans
filtnt[filtnt$seqid == "ASV_75",] = filtnt[filtnt$seqid == "ASV_75",][1,] 
filtnt[filtnt$seqid == "ASV_114",]$Species = NA # Rubus

# Check if Species is not unique within each ASV
non_unique_speciesnt <- filtnt %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)
length(unique(non_unique_speciesnt$seqid))  #now it should be 0

# Keep only one row per ASV
length(unique(filtnt$seqid)) # Check how many unique ASVs are there

final = list()
for(asv in unique(filtnt$seqid)){
  ASV_n = filtnt[filtnt$seqid == asv,]
  ASV_n = ASV_n[1,]
  final[[asv]] = ASV_n
}
finalnt <- as.data.frame(do.call(rbind, final)) 
dim(finalnt)

# Check minimum percentage of identical matches
min(finalnt$pident)

write.table(finalnt, "Results/2.blast/blast_final_95_nt", sep = "\t", row.names = FALSE)
write.xlsx(finalnt, "Results/2.blast/blast_final_95_nt.xlsx")

#### End ####

#### Merge UNITE and NT results ####

# Merge data frames
colnames(finalnt)[2:10] = paste0(colnames(finalnt)[2:10], "_nt")
colnames(finalu)[2:14] = paste0(colnames(finalu)[2:14], "_unite")

final = merge(finalnt, finalu, by = "seqid", all = T)
dim(final)
head(final)
write.table(final, "Results/2.blast/blast_final_95_nt_unite", sep = "\t", row.names = FALSE)
write.xlsx(final, "Results/2.blast/blast_final_95_nt_unite.xlsx")

# Create column for customized Genus and Species
final$Genus_custom <- final$Genus_unite
final$Species_custom <- final$Species_unite

# Manually revise ASVs with no genus identification in GTDB
final[is.na(final$Genus_unite),c("seqid","Genus_nt", "Species_nt")]

# Only one ASV was assigned to a fungus (ASV_125)
# Keep this assignment in our Genus_custom
final[final$seqid == "ASV_125",c("Genus_custom", "Species_custom")] = final[final$seqid == "ASV_125",c("Genus_nt", "Species_nt")]

# Manually revise ASVs with no species identification in GTDB
final[is.na(final$Species_unite),c("seqid","Genus_nt", "Species_nt")]
final[grep("_sp",final$Species_unite),c("seqid","Genus_unite", "Genus_nt", "Species_nt")]

# Three ASVs have a better resolution with nt (ASV_11, ASV_15, ASV_30, ASV_31, ASV_93)
final[final$seqid == "ASV_11",c("Genus_custom", "Species_custom")] = final[final$seqid == "ASV_11",c("Genus_nt", "Species_nt")]
final[final$seqid %in% c("ASV_29","ASV_15","ASV_30", "ASV_31", "ASV_74", "ASV_93"), "Species_custom"] =
  final[final$seqid %in% c("ASV_29","ASV_15","ASV_30", "ASV_31", "ASV_74", "ASV_93"), "Species_nt"]

# Save customized taxa table

write.table(final, "Results/2.blast/blast_final_95_nt_unite_customized", sep = "\t", row.names = FALSE)
write.xlsx(final, "Results/2.blast/blast_final_95_nt_unite_customized2.xlsx")

#### End ####
