# =============================================================================
# 06_FEAST_Resident_Analysis.R
# Purpose: Characterise the "resident cave" fraction identified by FEAST.
#          Identifies ASVs positively correlated with the Unknown (resident)
#          fraction, examines their taxonomy, and compares them to the full
#          cave community using beta-diversity and FAPROTAX functional annotation.
#
# Requires:
#   phyloseq_filtered_counts_DNA.rds               (from 01_Phyloseq_Setup.R)
#   FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx  (from 05_FEAST_Source_Tracking.R)
#
# Outputs (written to Outputs/FEAST_Resident/):
#   Resident_ASV_Correlations.xlsx
#   Resident_Community_Analysis.xlsx
# =============================================================================

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)
library(readxl)
library(writexl)
library(microeco)
library(tidyr)

# Output folder
out_dir <- "Outputs/FEAST_Resident"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
phy.filtered.DNA <- readRDS("phyloseq_filtered_counts_DNA.rds")

feast <- read_excel("FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx") %>%
  mutate(SampleID = sub("_Cave$", "", SampleID))

cat("FEAST samples loaded:", nrow(feast), "\n")
cat("Unknown fraction range:", round(min(feast$Unknown), 3),
    "-", round(max(feast$Unknown), 3), "\n")
cat("Unknown fraction mean:", round(mean(feast$Unknown), 3), "\n")

# =============================================================================
# STEP 1 — EXTRACT CAVE OTU TABLE AND METADATA
# =============================================================================

meta         <- as.data.frame(as.matrix(sample_data(phy.filtered.DNA)))
cave_samples <- rownames(meta)[meta$Site_Type == "Cave"]

otu_raw <- as.data.frame(otu_table(phy.filtered.DNA))
if (!taxa_are_rows(phy.filtered.DNA)) otu_raw <- t(otu_raw)
otu_raw  <- as.data.frame(otu_raw)
otu_cave <- otu_raw[, cave_samples]

feast_cave <- feast %>%
  mutate(SampleID_short = sub("-S.*", "", SampleID)) %>%
  filter(SampleID_short %in% sub("-S.*", "", cave_samples))

sample_unknown <- setNames(feast_cave$Unknown,
                           cave_samples[match(feast_cave$SampleID_short,
                                              sub("-S.*", "", cave_samples))])

cat(sprintf("\nCave samples matched to FEAST Unknown: %d\n", sum(!is.na(sample_unknown))))

# =============================================================================
# STEP 2 — PER-ASV SPEARMAN CORRELATION VS UNKNOWN FRACTION
# =============================================================================

cat("\nRunning per-ASV Spearman correlations vs Unknown fraction...\n")

otu_cave_rel <- sweep(otu_cave, 2, colSums(otu_cave), "/")

min_presence <- 3
present_in   <- rowSums(otu_cave > 0)
otu_cave_rel <- otu_cave_rel[present_in >= min_presence, ]

cat(sprintf("ASVs present in >= %d cave samples: %d\n", min_presence, nrow(otu_cave_rel)))

# Align sample order
common_samples  <- intersect(colnames(otu_cave_rel), names(sample_unknown))
otu_cave_rel    <- otu_cave_rel[, common_samples]
unknown_aligned <- sample_unknown[common_samples]

asv_results <- lapply(seq_len(nrow(otu_cave_rel)), function(i) {
  asv_abund <- as.numeric(otu_cave_rel[i, ])
  test      <- suppressWarnings(
    cor.test(asv_abund, unknown_aligned, method = "spearman", exact = FALSE)
  )
  data.frame(
    ASV_ID  = rownames(otu_cave_rel)[i],
    rho     = round(test$estimate, 4),
    P_value = round(test$p.value, 4)
  )
}) %>% bind_rows()

asv_results$P_adjusted  <- round(p.adjust(asv_results$P_value, method = "fdr"), 4)
asv_results$Significant <- ifelse(asv_results$P_adjusted < 0.05, "YES ***",
                                  ifelse(asv_results$P_adjusted < 0.10, "TRENDING ~", "NO"))
asv_results <- asv_results %>% arrange(P_value)

cat(sprintf("\nASVs positively correlated (FDR < 0.05): %d\n",
            sum(asv_results$P_adjusted < 0.05 & asv_results$rho > 0)))
cat(sprintf("ASVs negatively correlated (FDR < 0.05): %d\n",
            sum(asv_results$P_adjusted < 0.05 & asv_results$rho < 0)))

# =============================================================================
# STEP 3 — ATTACH TAXONOMY TO CORRELATED ASVs
# =============================================================================

tax_full <- as.data.frame(tax_table(phy.filtered.DNA))
for (col in c("Kingdom","Phylum","Class","Order","Family","Genus")) {
  tax_full[[col]] <- gsub("^[a-z]__", "", tax_full[[col]])
  tax_full[[col]][tax_full[[col]] == ""] <- NA
}

asv_results_tax <- asv_results %>%
  left_join(tax_full %>% rownames_to_column("ASV_ID"), by = "ASV_ID")

resident_asvs <- asv_results_tax %>% filter(P_adjusted < 0.05, rho > 0)
cat("\nTop 20 resident-associated ASVs (ranked by rho):\n")
print(head(resident_asvs %>% arrange(desc(rho)) %>%
             select(ASV_ID, rho, P_adjusted, Phylum, Class, Order, Genus), 20))

# Taxonomic composition of resident ASVs
if (nrow(resident_asvs) > 0) {
  cat("\nPhylum breakdown of resident-associated ASVs:\n")
  print(resident_asvs %>%
          count(Phylum, name = "N") %>%
          mutate(Pct = round(100 * N / sum(N), 2)) %>%
          arrange(desc(N)))
}

# =============================================================================
# STEP 4 — EXPORT
# =============================================================================

write_xlsx(
  list(
    All_ASV_Correlations    = as.data.frame(asv_results_tax),
    Resident_ASVs           = as.data.frame(resident_asvs),
    Resident_Phylum_Summary = as.data.frame(
      resident_asvs %>%
        count(Phylum, name = "N") %>%
        mutate(Pct = round(100 * N / sum(N), 2)) %>%
        arrange(desc(N))
    )
  ),
  file.path(out_dir, "Resident_ASV_Correlations.xlsx")
)
cat("\nSaved: Resident_ASV_Correlations.xlsx\n")

# =============================================================================
# STEP 5 — BETA-DIVERSITY OF RESIDENT vs NON-RESIDENT ASVs (CAVE SAMPLES)
# =============================================================================

cat("\n=== Beta-diversity: Resident vs Non-resident ASVs ===\n")

resident_ids     <- resident_asvs$ASV_ID
non_resident_ids <- asv_results_tax %>% filter(P_adjusted < 0.05, rho < 0) %>% pull(ASV_ID)

if (length(resident_ids) >= 3 && length(non_resident_ids) >= 3) {
  otu_res     <- as.data.frame(t(otu_cave[resident_ids,     cave_samples]))
  otu_non_res <- as.data.frame(t(otu_cave[non_resident_ids, cave_samples]))

  bray_res     <- vegdist(otu_res,     method = "bray")
  bray_non_res <- vegdist(otu_non_res, method = "bray")

  cat("Mean Bray-Curtis dissimilarity:\n")
  cat(sprintf("  Resident ASVs:     %.3f\n", mean(bray_res)))
  cat(sprintf("  Non-resident ASVs: %.3f\n", mean(bray_non_res)))

  community_summary <- data.frame(
    Community        = c("Resident-associated", "Non-resident-associated"),
    N_ASVs           = c(length(resident_ids), length(non_resident_ids)),
    Mean_BrayCurtis  = round(c(mean(bray_res), mean(bray_non_res)), 3)
  )

  write_xlsx(
    list(Community_Summary = community_summary),
    file.path(out_dir, "Resident_Community_Analysis.xlsx")
  )
  cat("Saved: Resident_Community_Analysis.xlsx\n")
} else {
  cat("Insufficient resident or non-resident ASVs for beta-diversity comparison.\n")
}

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("============================\n")
