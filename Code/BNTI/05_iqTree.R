library(phyloseq)
library(ape)

install.packages("phangorn")
library(phangorn)
tree_rooted <- midpoint(tree)

# Load filtered phyloseq object
phy <- readRDS("Inputs/phyloseq_filtered_counts_DNA.rds")

# Load the tree
tree <- read.tree("Inputs/final_tree.treefile")

# Check tree
cat("Number of tips:", length(tree$tip.label), "\n")
cat("Is rooted:", is.rooted(tree), "\n")

# Midpoint root
tree_rooted <- midpoint(tree)
cat("Is rooted after midpoint:", is.rooted(tree_rooted), "\n")

# Merge into phyloseq
phy_tree <- merge_phyloseq(phy, phy_tree(tree_rooted))

# Verify
cat("Phyloseq taxa:", ntaxa(phy_tree), "\n")

# Save
saveRDS(phy_tree, "phyloseq_filtered_counts_DNA_tree.rds")
