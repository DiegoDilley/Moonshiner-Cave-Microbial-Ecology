#This is the first R script you run to prep data for command line stuff
#Inputs are pruned raw counts and sequence
#output is asvs_for_tree.fasta

library(phyloseq)
library(Biostrings)

phy <- readRDS("Inputs/phyloseq_filtered_counts_DNA.rds")

# Read the sequences fasta
rep_seqs <- readDNAStringSet("Inputs/sequences.fasta")

# Keep only ASVs in your filtered phyloseq object
filtered_asvs <- taxa_names(phy)
rep_seqs_filtered <- rep_seqs[names(rep_seqs) %in% filtered_asvs]

cat("ASVs in phyloseq:", length(filtered_asvs), "\n")
cat("ASVs matched in sequences.fasta:", length(rep_seqs_filtered), "\n")

# Export filtered sequences for tree building
writeXStringSet(rep_seqs_filtered, filepath = "Inputs/asvs_for_tree.fasta")

#####Below section is used to generate a subset of my ASVS for a test SLURM job, make sure i dont waste (as much) OSC time
rep_seqs <- readDNAStringSet("Inputs/asvs_for_tree.fasta")
test_seqs <- rep_seqs[1:250]
writeXStringSet(test_seqs, filepath = "Inputs/asvs_for_tree_TEST.fasta")
cat("Test fasta written:", length(test_seqs), "sequences\n")
