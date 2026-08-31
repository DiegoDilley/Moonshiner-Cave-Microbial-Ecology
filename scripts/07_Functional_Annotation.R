# =============================================================================
# 07_Functional_Annotation.R
# Purpose: Predict metabolic functional groups for the full community and for
#          cave-exclusive ASVs using FAPROTAX via the microeco package.
#          Tests how functional profiles differ across site types and change
#          with cave depth. Also tests whether the FEAST Unknown fraction
#          directly predicts predatory activity.
#
# Requires:
#   phyloseq_filtered_counts_DNA.rds               (from 01_Phyloseq_Setup.R)
#   FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx  (from 05_FEAST_Source_Tracking.R)
#   FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx  (for predatory correlation test)
#
# Outputs (written to Outputs/Functional_Annotation/):
#   FAPROTAX_Results.xlsx
#   FAPROTAX_Statistical_Summary.xlsx
#   FAPROTAX_DepthCorrelations.xlsx
#   FAPROTAX_CaveExclusive.xlsx
#   Predatory_vs_FEAST_Unknown.xlsx
#   Plots/Plot_DepthCorrelations.png
# =============================================================================

library(phyloseq)
library(writexl)
library(readxl)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(microeco)
library(stringr)

# Output folder
out_dir <- "Outputs/Functional_Annotation"
dir.create(file.path(out_dir, "Plots"), recursive = TRUE, showWarnings = FALSE)

# Load data
phy.filtered.DNA <- readRDS("phyloseq_filtered_counts_DNA.rds")

sig <- function(p) ifelse(p < 0.05, "YES ***", ifelse(p < 0.10, "TRENDING ~", "NO"))

fix_names  <- function(x) gsub("\\.", "-", x)

clean_tax_df <- function(df) {
  for (col in c("Kingdom","Phylum","Class","Order","Family","Genus")) {
    df[[col]] <- gsub("^[a-z]__", "", df[[col]])
    df[[col]][df[[col]] == ""] <- NA
  }
  df
}

# =============================================================================
# SECTION 1 — FAPROTAX: FULL COMMUNITY (ALL SAMPLES)
# =============================================================================

cat("\n=== SECTION 1: Full Community FAPROTAX ===\n")

# Export tables for microeco
counts_full <- as.data.frame(otu_table(phy.filtered.DNA))
counts_full[is.na(counts_full)] <- 0
counts_full[] <- lapply(counts_full, as.numeric)

tax_exp <- as.data.frame(tax_table(phy.filtered.DNA))
tax_exp <- clean_tax_df(tax_exp)
tax_exp <- tax_exp[, c("Kingdom","Phylum","Class","Order","Family","Genus")]

metadata <- as.data.frame(sample_data(phy.filtered.DNA))

# Align ASVs and samples
common_asvs    <- intersect(rownames(counts_full), rownames(tax_exp))
common_samples <- intersect(colnames(counts_full), rownames(metadata))
counts_full    <- counts_full[common_asvs, common_samples]
tax_exp        <- tax_exp[common_asvs, ]
metadata       <- metadata[common_samples, ]

cat(sprintf("Matched: %d ASVs, %d samples\n", length(common_asvs), length(common_samples)))

dataset <- microtable$new(sample_table = metadata, otu_table = counts_full, tax_table = tax_exp)
dataset$tidy_dataset()
dataset$cal_abund()

faprotax <- trans_func$new(dataset)
faprotax$cal_spe_func(prok_database = "FAPROTAX")
faprotax$cal_spe_func_perc(abundance_weighted = TRUE)

results <- as.data.frame(faprotax$res_spe_func_perc)
results <- cbind(
  SampleID     = rownames(results),
  Site_Type    = metadata[rownames(results), "Site_Type",    drop = TRUE],
  Site_Number  = metadata[rownames(results), "Site_Number",  drop = TRUE],
  Distance_M   = as.numeric(metadata[["Distance_Datum_M"]][match(rownames(results), rownames(metadata))]),
  WaterPercent = as.numeric(metadata[["Water_Percantage"]][match(rownames(results), rownames(metadata))]),
  results
)

cat(sprintf("FAPROTAX: %d samples, %d functional groups\n",
            nrow(results), ncol(results) - 5))

write_xlsx(
  list(
    Functional_Abundance  = results,
    Functional_Groups_Raw = as.data.frame(faprotax$res_spe_func)
  ),
  file.path(out_dir, "FAPROTAX_Results.xlsx")
)
cat("Saved: FAPROTAX_Results.xlsx\n")

# --- Kruskal-Wallis by site type ---
df_kw <- results
meta_cols <- c("SampleID","Site_Type","Site_Number","Distance_M","WaterPercent",
               "Water_Percantage","Distance_Datum_M","SiteReplicate")
func_cols <- setdiff(colnames(df_kw), meta_cols)
df_kw[func_cols] <- lapply(df_kw[func_cols], as.numeric)
nonzero_cols <- func_cols[colSums(df_kw[func_cols], na.rm = TRUE) > 0]

kw_results <- lapply(nonzero_cols, function(f) {
  test <- kruskal.test(df_kw[[f]] ~ df_kw$Site_Type)
  data.frame(Function = f, H_statistic = round(test$statistic, 3),
             P_value = round(test$p.value, 4), df = test$parameter)
}) %>% bind_rows()

kw_results$P_adjusted  <- round(p.adjust(kw_results$P_value, method = "fdr"), 4)
kw_results$Significant <- ifelse(kw_results$P_adjusted < 0.05, "YES ***",
                                 ifelse(kw_results$P_adjusted < 0.10, "TRENDING ~", "NO"))
kw_results <- kw_results %>% arrange(P_adjusted)

cat(sprintf("\nSignificant functions (FDR < 0.05): %d\n",
            sum(kw_results$P_adjusted < 0.05, na.rm = TRUE)))

sig_funcs <- kw_results %>% filter(P_adjusted < 0.05) %>% pull(Function)

pairwise_results <- lapply(sig_funcs, function(f) {
  pairs <- combn(unique(df_kw$Site_Type), 2, simplify = FALSE)
  lapply(pairs, function(p) {
    g1   <- df_kw[[f]][df_kw$Site_Type == p[1]]
    g2   <- df_kw[[f]][df_kw$Site_Type == p[2]]
    test <- wilcox.test(g1, g2)
    data.frame(
      Function = f, Group1 = p[1], Group2 = p[2],
      Mean_G1 = round(mean(g1, na.rm = TRUE), 4), Mean_G2 = round(mean(g2, na.rm = TRUE), 4),
      W_statistic = test$statistic, P_value = round(test$p.value, 4),
      Significant = sig(test$p.value)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

mean_by_site <- df_kw %>%
  group_by(Site_Type) %>%
  summarise(across(all_of(nonzero_cols), mean, na.rm = TRUE)) %>%
  pivot_longer(-Site_Type, names_to = "Function", values_to = "Mean_Abundance") %>%
  pivot_wider(names_from = Site_Type, values_from = Mean_Abundance) %>%
  left_join(kw_results %>% select(Function, P_adjusted, Significant), by = "Function") %>%
  arrange(P_adjusted)

write_xlsx(
  list(
    KruskalWallis_All    = as.data.frame(kw_results),
    Pairwise_Significant = as.data.frame(pairwise_results),
    Mean_by_SiteType     = as.data.frame(mean_by_site),
    Raw_Data             = as.data.frame(df_kw)
  ),
  file.path(out_dir, "FAPROTAX_Statistical_Summary.xlsx")
)
cat("Saved: FAPROTAX_Statistical_Summary.xlsx\n")

# --- Depth correlations (cave samples only) ---
cave_df <- df_kw %>%
  filter(Site_Type == "Cave") %>%
  mutate(Distance_M = as.numeric(Distance_M))

sig_funcs_depth <- kw_results %>% filter(P_adjusted < 0.05) %>% pull(Function)
sig_funcs_depth <- sig_funcs_depth[sig_funcs_depth %in% colnames(cave_df)]

cor_results <- lapply(sig_funcs_depth, function(f) {
  vals <- as.numeric(cave_df[[f]])
  test <- cor.test(cave_df$Distance_M, vals, method = "spearman", exact = FALSE)
  data.frame(Function = f, rho = round(test$estimate, 3),
             P_value = round(test$p.value, 4),
             Direction = ifelse(test$estimate > 0, "Increases with depth", "Decreases with depth"))
}) %>% bind_rows()

cor_results$P_adjusted  <- round(p.adjust(cor_results$P_value, method = "fdr"), 4)
cor_results$Significant <- ifelse(cor_results$P_adjusted < 0.05, "YES ***",
                                  ifelse(cor_results$P_adjusted < 0.10, "TRENDING ~", "NO"))
cor_results <- cor_results %>% arrange(P_value)

cat("\nFunctions correlated with cave depth (p < 0.05):\n")
print(cor_results %>% filter(P_value < 0.05))

write_xlsx(cor_results, file.path(out_dir, "FAPROTAX_DepthCorrelations.xlsx"))
cat("Saved: FAPROTAX_DepthCorrelations.xlsx\n")

# --- Depth correlation plot ---
plot_funcs <- cor_results %>% filter(P_value < 0.10) %>% pull(Function)

if (length(plot_funcs) > 0) {
  clean_name <- function(x) str_replace_all(x, "_", " ") %>% str_to_sentence()

  plot_data <- cave_df %>%
    select(SampleID, Distance_M, all_of(plot_funcs)) %>%
    pivot_longer(cols = all_of(plot_funcs), names_to = "Function", values_to = "Abundance") %>%
    mutate(Function = clean_name(Function)) %>%
    left_join(cor_results %>% mutate(Function = clean_name(Function)) %>%
                select(Function, rho, P_value), by = "Function") %>%
    mutate(Label = sprintf("%s\nρ = %.3f, p = %.3f", Function, rho, P_value))

  p_depth <- ggplot(plot_data, aes(x = Distance_M, y = Abundance)) +
    geom_point(color = "#4A4A8A", size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#c0392b",
                fill = "#fadbd8", linewidth = 0.8) +
    facet_wrap(~Label, scales = "free_y") +
    labs(title    = "Functional Groups Changing with Distance into Cave",
         subtitle = "Spearman correlation, cave samples only",
         x        = "Distance from Entrance (m)",
         y        = "Relative Abundance") +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#4A4A8A"),
      strip.text       = element_text(color = "white", face = "bold", size = 8),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(color = "grey40")
    )

  ggsave(file.path(out_dir, "Plots", "Plot_DepthCorrelations.png"),
         p_depth, width = 12, height = 8, dpi = 300, limitsize = FALSE)
  cat("Saved: Plots/Plot_DepthCorrelations.png\n")
}

# =============================================================================
# SECTION 2 — FAPROTAX: CAVE-EXCLUSIVE ASVs
# =============================================================================

cat("\n=== SECTION 2: Cave-Exclusive ASV FAPROTAX ===\n")

otu_raw <- as.data.frame(otu_table(phy.filtered.DNA))
if (!taxa_are_rows(phy.filtered.DNA)) otu_raw <- t(otu_raw)
otu_raw           <- as.data.frame(otu_raw)
colnames(otu_raw) <- fix_names(colnames(otu_raw))

meta_ce           <- as.data.frame(as.matrix(sample_data(phy.filtered.DNA)))
rownames(meta_ce) <- fix_names(rownames(meta_ce))

cave_samples   <- rownames(meta_ce)[meta_ce$Site_Type == "Cave"]
source_samples <- rownames(meta_ce)[meta_ce$Site_Type != "Cave"]

cave_present        <- rowSums(otu_raw[, cave_samples])   > 0
source_present      <- rowSums(otu_raw[, source_samples]) > 0
cave_exclusive_asvs <- rownames(otu_raw)[cave_present & !source_present]

cat(sprintf("Cave-exclusive ASVs: %d\n", length(cave_exclusive_asvs)))

make_numeric <- function(df) {
  df <- as.data.frame(lapply(df, as.numeric), row.names = rownames(df))
  colnames(df) <- fix_names(colnames(df))
  df
}
safe_t <- function(df) {
  out <- as.data.frame(t(df))
  rownames(out) <- fix_names(rownames(out))
  colnames(out) <- fix_names(colnames(out))
  as.data.frame(lapply(out, as.numeric), row.names = rownames(out))
}

otu_excl_raw    <- make_numeric(otu_raw[cave_exclusive_asvs, cave_samples])
colnames(otu_excl_raw) <- fix_names(colnames(otu_excl_raw))

tax_ce        <- as.data.frame(tax_table(phy.filtered.DNA))
tax_ce        <- clean_tax_df(tax_ce)
tax_exclusive <- tax_ce[cave_exclusive_asvs, c("Kingdom","Phylum","Class","Order","Family","Genus")]
meta_mc       <- meta_ce[cave_samples, ]

# Load FEAST unknown for annotation
feast2 <- read_excel("FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx") %>%
  mutate(SampleID = fix_names(sub("_Cave$", "", SampleID)),
         SampleID_short = sub("-S.*", "", SampleID))

common_mc   <- intersect(colnames(otu_excl_raw), rownames(meta_mc))
otu_excl_raw <- otu_excl_raw[, common_mc]
meta_mc      <- meta_mc[common_mc, ]

cat(sprintf("Cave-exclusive FAPROTAX — ASVs: %d, Samples: %d\n",
            nrow(otu_excl_raw), ncol(otu_excl_raw)))
cat("Names match:", all(colnames(otu_excl_raw) == rownames(meta_mc)), "\n")

dataset_excl <- microtable$new(sample_table = meta_mc,
                               otu_table    = otu_excl_raw,
                               tax_table    = tax_exclusive)
dataset_excl$tidy_dataset()
dataset_excl$cal_abund()

faprotax_excl <- trans_func$new(dataset_excl)
faprotax_excl$cal_spe_func(prok_database = "FAPROTAX")
faprotax_excl$cal_spe_func_perc(abundance_weighted = TRUE)

results_excl <- as.data.frame(faprotax_excl$res_spe_func_perc)
results_excl <- cbind(
  SampleID     = fix_names(rownames(results_excl)),
  Distance_M   = as.numeric(meta_mc[rownames(results_excl), "Distance_Datum_M"]),
  WaterPercent = as.numeric(meta_mc[rownames(results_excl), "Water_Percantage"]),
  results_excl
) %>%
  mutate(SampleID_short = sub("-S.*", "", SampleID)) %>%
  left_join(feast2 %>% select(SampleID_short, Unknown), by = "SampleID_short") %>%
  select(-SampleID_short)

func_cols_excl <- setdiff(colnames(results_excl), c("SampleID","Distance_M","WaterPercent","Unknown"))
results_excl[func_cols_excl] <- lapply(results_excl[func_cols_excl], as.numeric)
nonzero_excl <- func_cols_excl[colSums(results_excl[func_cols_excl], na.rm = TRUE) > 0]

cat(sprintf("\nCave-exclusive FAPROTAX: %d non-zero functions\n", length(nonzero_excl)))

means_excl <- data.frame(
  Function       = nonzero_excl,
  Mean_Abundance = round(colMeans(results_excl[nonzero_excl], na.rm = TRUE), 5)
) %>% arrange(desc(Mean_Abundance))

run_cor <- function(predictor, label) {
  lapply(nonzero_excl, function(f) {
    test <- cor.test(predictor, as.numeric(results_excl[[f]]),
                     method = "spearman", exact = FALSE)
    data.frame(Function = f, rho = round(test$estimate, 3), P_value = round(test$p.value, 4))
  }) %>% bind_rows() %>%
    mutate(P_adjusted  = round(p.adjust(P_value, method = "fdr"), 4),
           Significant = sig(P_adjusted),
           Direction   = ifelse(rho > 0,
                                paste("Increases with", label),
                                paste("Decreases with", label))) %>%
    arrange(P_value)
}

dist_cor_excl <- run_cor(results_excl$Distance_M, "depth")
unk_cor_excl  <- run_cor(results_excl$Unknown,     "Unknown")

cat("\nDistance correlations (p < 0.10):\n"); print(dist_cor_excl %>% filter(P_value < 0.10))
cat("\nUnknown correlations (p < 0.10):\n");  print(unk_cor_excl  %>% filter(P_value < 0.10))

write_xlsx(
  list(
    FAPROTAX_Mean        = as.data.frame(means_excl),
    FAPROTAX_vs_Distance = as.data.frame(dist_cor_excl),
    FAPROTAX_vs_Unknown  = as.data.frame(unk_cor_excl),
    FAPROTAX_Raw         = as.data.frame(results_excl)
  ),
  file.path(out_dir, "FAPROTAX_CaveExclusive.xlsx")
)
cat("Saved: FAPROTAX_CaveExclusive.xlsx\n")

# =============================================================================
# SECTION 3 — PREDATORY ACTIVITY VS FEAST UNKNOWN (SUPPLEMENTARY TEST)
# =============================================================================

cat("\n=== SECTION 3: Predatory Activity vs FEAST Unknown ===\n")

faprotax_stats <- read_excel(file.path(out_dir, "FAPROTAX_Statistical_Summary.xlsx"),
                             sheet = "Raw_Data")

feast1 <- read_excel("FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx") %>%
  mutate(SampleID_short = sub("_Cave$", "", SampleID),
         SampleID_short = sub("-S.*", "",  SampleID_short))

faprotax_cave <- faprotax_stats %>%
  filter(Site_Type == "Cave") %>%
  mutate(SampleID_short = sub("-S.*", "", SampleID))

merged <- faprotax_cave %>%
  left_join(feast1 %>% select(SampleID_short, Unknown) %>%
              rename(FEAST_Unknown = Unknown), by = "SampleID_short")

cat(sprintf("Cave samples: %d | Matched with FEAST: %d\n",
            nrow(faprotax_cave), sum(!is.na(merged$FEAST_Unknown))))

if ("predatory_or_exoparasitic" %in% colnames(merged)) {
  test <- cor.test(merged$FEAST_Unknown, merged$predatory_or_exoparasitic,
                   method = "spearman", exact = FALSE)
  cat(sprintf("\nPredatory vs FEAST Unknown:\n  rho = %.3f, p = %.4f  %s\n",
              test$estimate, test$p.value, sig(test$p.value)))

  result <- data.frame(
    Test         = "Predatory_or_Exoparasitic vs FEAST_Unknown",
    n            = nrow(merged),
    Spearman_rho = round(test$estimate, 3),
    P_value      = round(test$p.value, 4),
    Significant  = sig(test$p.value),
    Interpretation = ifelse(test$p.value >= 0.05,
                            "Independent — both track depth but not each other",
                            "Direct relationship detected")
  )

  write_xlsx(
    list(
      Correlation_Result = result,
      Raw_Data           = as.data.frame(merged %>%
                             select(SampleID, Distance_M, FEAST_Unknown,
                                    predatory_or_exoparasitic))
    ),
    file.path(out_dir, "Predatory_vs_FEAST_Unknown.xlsx")
  )
  cat("Saved: Predatory_vs_FEAST_Unknown.xlsx\n")
} else {
  cat("Column 'predatory_or_exoparasitic' not found in FAPROTAX results — skipping test.\n")
}

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("============================\n")
