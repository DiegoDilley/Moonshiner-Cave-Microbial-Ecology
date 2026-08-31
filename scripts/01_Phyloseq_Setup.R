# =============================================================================
# 01_Phyloseq_Setup.R
# Purpose: Build, filter, and save phyloseq objects from raw QIIME2 outputs.
#          This script must be run first — all downstream scripts depend on its
#          .rds outputs.
#
# Inputs (place in Inputs/):
#   feature-table.tsv        — QIIME2 ASV count table
#   taxonomy.tsv             — QIIME2 taxonomy with 7 rank columns
#   Karst_Input_Metadata.tsv — sample metadata
#
# Outputs (written to project root and Outputs/):
#   phyloseq-DNA.rds                  — raw phyloseq object
#   phyloseq_filtered_counts_DNA.rds  — filtered count phyloseq (primary input for most scripts)
#   phyloseq_filtered_RA_DNA.rds      — filtered relative-abundance phyloseq
#   Outputs/Unfiltered_Sample_Check.csv
#   Outputs/Rarefaction_Diagnostics.pdf
#   Outputs/Rarefaction_One_Per_Page.pdf
# =============================================================================

library(phyloseq)
library(data.table)
library(dplyr)
library(ggplot2)
library(vegan)

# Create output directory
if (!dir.exists("Outputs")) dir.create("Outputs", recursive = TRUE)

# =============================================================================
# STEP 1 — READ INPUT FILES
# =============================================================================

otu <- read.table("Inputs/feature-table.tsv",
                  header = TRUE, row.names = 1,
                  sep = "\t", skip = 1, comment.char = "")
otu_table_ps <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
colnames(otu_table_ps) <- gsub("\\.", "-", colnames(otu_table_ps))

tax <- read.table("Inputs/taxonomy.tsv",
                  header = TRUE, row.names = 1, sep = "\t",
                  fill = TRUE, comment.char = "")
tax_split <- tax[, 1:7]
colnames(tax_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_table_ps <- tax_table(as.matrix(tax_split))

metadata <- read.table("Inputs/Karst_Input_Metadata.tsv",
                       header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, fill = TRUE)
metadata <- metadata[, apply(metadata, 2, function(x) !all(is.na(x)))]
sample_data_ps <- sample_data(metadata)

# =============================================================================
# STEP 2 — BUILD RAW PHYLOSEQ AND EXPORT UNFILTERED CHECK
# =============================================================================

phy.DNA <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
saveRDS(phy.DNA, "phyloseq-DNA.rds")

sample_df <- as(sample_data(phy.DNA), "data.frame")
seq.dt <- data.table(SampleID   = rownames(sample_df),
                     sample_df,
                     TotalReads = sample_sums(phy.DNA))
fwrite(seq.dt, "Outputs/Unfiltered_Sample_Check.csv")

# =============================================================================
# STEP 3 — RAREFACTION CURVES
# =============================================================================

phy.DNA <- prune_samples(sample_sums(phy.DNA) >= 100, phy.DNA)

otu_tab <- as(otu_table(phy.DNA), "matrix")
if (taxa_are_rows(phy.DNA)) otu_tab <- t(otu_tab)

pdf("Outputs/Rarefaction_Diagnostics.pdf", width = 12, height = 10)
rarecurve(otu_tab,
          step  = 100,
          col   = "darkblue",
          label = TRUE,
          cex   = 0.6,
          main  = "Rarefaction Curves: Check for Plateau & Sample Names",
          xlab  = "Sequencing Depth (Reads)",
          ylab  = "Number of ASVs")
abline(v = 100,  col = "red",     lty = 2)
abline(v = 1000, col = "darkred", lty = 3)
while (!is.null(dev.list())) dev.off()

pdf("Outputs/Rarefaction_One_Per_Page.pdf", width = 8, height = 8)
sample_names_vec <- rownames(otu_tab)
for (i in seq_len(nrow(otu_tab))) {
  current_sample <- otu_tab[i, , drop = FALSE]
  rarecurve(current_sample,
            step  = 100,
            col   = "darkblue",
            label = FALSE,
            main  = paste("Rarefaction Curve:", sample_names_vec[i]),
            xlab  = "Sequencing Depth (Reads)",
            ylab  = "Number of ASVs")
  abline(v = 100,  col = "red",     lty = 2)
  abline(v = 1000, col = "darkred", lty = 3)
}
dev.off()

# =============================================================================
# STEP 4 — FILTER AND SAVE
# =============================================================================

# Remove samples below read depth threshold
phy.DNA <- prune_samples(sample_sums(phy.DNA) >= 1000, phy.DNA)

# Remove chloroplast and mitochondria
phy.DNA <- subset_taxa(
  phy.DNA,
  !tax_table(phy.DNA)[, "Family"] %in% "Mitochondria" &
    !tax_table(phy.DNA)[, "Order"]  %in% "Chloroplast"
)
message("Removed chloroplast and mitochondria.")

# Remove low-abundance taxa
keep_taxa <- taxa_sums(phy.DNA) > 10
phy.DNA   <- prune_taxa(keep_taxa, phy.DNA)
message("Removed low-abundance taxa (<=10 total reads).")

# Save count and relative-abundance objects
phy.filtered.DNA.rel <- transform_sample_counts(phy.DNA, function(x) x / sum(x))

saveRDS(phy.DNA,              "phyloseq_filtered_counts_DNA.rds")
saveRDS(phy.filtered.DNA.rel, "phyloseq_filtered_RA_DNA.rds")

cat("\n============================\n")
cat("SETUP COMPLETE\n")
cat("Outputs: phyloseq_filtered_counts_DNA.rds\n")
cat("         phyloseq_filtered_RA_DNA.rds\n")
cat("         Outputs/Unfiltered_Sample_Check.csv\n")
cat("         Outputs/Rarefaction_Diagnostics.pdf\n")
cat("         Outputs/Rarefaction_One_Per_Page.pdf\n")
cat("============================\n")
