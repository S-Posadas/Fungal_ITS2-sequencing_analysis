# Fungal_ITS2-sequencing_analysis
Copyright Sara Posadas Cantera, Date: 21 February 2021 University of Freiburg Germany

This code was used for the analysis of ITS2 gene sequencing of clinical ascites samples from the INTeGRATE clinical study. The analysis starts from raw fastq files where used primers get trimmed using cutadapt and further analysis is carried out using DADA2 to end up with amplicon sequence variant (ASV) table. Taxonomy analysis combines DADA2 and Blast approaches. Phyloseq is used for downstream analysis.