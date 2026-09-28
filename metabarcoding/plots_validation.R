# ==============================================================================
# PERCIVAL - Rhizosphere Metabarcoding Technical Validation (16S & ITS)
# Script to generate Read-Depth Distributions and ASV Rarefaction Curves
# ==============================================================================

# 1. Load required libraries
library("phyloseq")
library("ggplot2")
library("vegan")

# NOTE: Set your working directory to the folder containing your BIOM, 
# tree, and FASTA files before running the script.
# Example: setwd("path/to/your/data")

# ==============================================================================
# SECTION 1: FUNGAL DATASET (ITS)
# ==============================================================================

# Load fungal biom object, rooted tree, and representative sequences
ps_its <- import_biom("output_its.biom", 
                      treefilename="tree_rooted_its.tree",
                      refseqfilename="otus_its.fasta")

# Export High-Resolution PNG for Fungal Technical Validation
png("Supplementary_Figure_S1_Fungi_Validation.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2)) # Split plotting area into 1 row and 2 columns

# Panel A: Read-depth distribution (including negative control at 0)
plot(sample_sums(ps_its), 
     xlab = "Samples", 
     ylab = "Nr. Reads", 
     main = "A: Read-depth distribution (ITS)", 
     pch = 19, col = "darkred")

# Panel B: Rarefaction curves (excluding negative control to avoid execution errors)
otu_its_matrix <- as(otu_table(ps_its), "matrix")
if(taxa_are_rows(ps_its)) {
  otu_its_matrix <- t(otu_its_matrix)
}
# Assuming the negative control is the first sample (index 1), exclude it with [-1, ]
rarecurve(otu_its_matrix[-1, ], 
          step = 1000, 
          xlab = "Nr. Reads", 
          ylab = "Observed ASVs", 
          main = "B: Rarefaction curves (ITS)", 
          col = "black")

dev.off() # Close graphic device and save file

# ==============================================================================
# SECTION 2: PROKARYOTIC DATASET (16S)
# ==============================================================================

# Load prokaryotic biom object, rooted tree, and representative sequences
ps_16s <- import_biom("output_16s.biom", 
                      treefilename="tree_rooted_16s.tree",
                      refseqfilename="otus_16s.fasta")


# Export High-Resolution PNG for Bacterial Technical Validation
png("Supplementary_Figure_S2_Bacteria_Validation.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2)) # Split plotting area into 1 row and 2 columns

# Panel A: Read-depth distribution (including negative control at 0)
plot(sample_sums(ps_16s), 
     xlab = "Samples", 
     ylab = "Nr. Reads", 
     main = "A: Read-depth distribution (16S)", 
     pch = 19, col = "darkblue")

# Panel B: Rarefaction curves (excluding negative control to avoid execution errors)
otu_16s_matrix <- as(otu_table(ps_16s), "matrix")
if(taxa_are_rows(ps_16s)) {
  otu_16s_matrix <- t(otu_16s_matrix)
}
# Assuming the negative control is the first sample (index 1), exclude it with [-1, ]
rarecurve(otu_16s_matrix[-1, ], 
          step = 1000, 
          xlab = "Nr. Reads", 
          ylab = "Observed ASVs", 
          main = "B: Rarefaction curves (16S)", 
          col = "black")

dev.off() # Close graphic device and save file
