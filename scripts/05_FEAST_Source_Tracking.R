setwd("C:/Users/Diego/OneDrive/Desktop/Moonshiner_Microbial_Ecology/Moonshiner_16S_Code")
# =============================================================================
# 05_FEAST_Source_Tracking.R
# Purpose: Run FEAST microbial source-tracking on cave samples using surface
#          and water environments as potential sources. Runs two iterations:
#          (1) separate surface sub-types (Cave Entrance, Soil Transect, Water);
#          (2) collapsed surface sources. Runs statistical analysis on how source
#          contributions vary with cave depth and moisture. Produces Figures 8a,
#          8b, and 9.
#
# Requires:
#   Inputs/Karst_Input_Metadata.TSV
#   Inputs/feature-table.tsv
#   phyloseq_filtered_RA_DNA.rds    (for figure metadata; from 01_Phyloseq_Setup.R)
#
# Outputs (written to FEAST_output/ and Outputs/FEAST/):
#   FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx
#   FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx
#   FEAST_output/Karst_FEAST_iteration2_Combined.xlsx
#   FEAST_output/FEASTV2_Full_Statistical_Summary.xlsx
#   Outputs/FEAST/Figure8a_FEAST_Iteration1.png / .pdf
#   Outputs/FEAST/Figure8b_FEAST_Iteration2.png / .pdf
#   Outputs/FEAST/Figure9_FEAST_vs_Distance.png / .pdf
# =============================================================================

library(tidyverse)
library(writexl)
library(readxl)
library(FEAST)
library(phyloseq)
library(forcats)

# Save project root — FEAST changes the working directory internally when it
# runs. We restore to this path immediately after every FEAST() call so that
# all subsequent file paths resolve correctly from the project root.
project_root <- getwd()

dir.create("FEAST_output", showWarnings = FALSE)
out_dir <- "Outputs/FEAST"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1 — PREPARE METADATA AND ASV TABLE FOR FEAST
# =============================================================================

metadata <- read.delim("Inputs/Karst_Input_Metadata.TSV") %>%
  mutate(
    SourceSink = ifelse(Site_Type == "Cave", "Sink", "Source"),
    Env = case_when(
      SampleName %in% c("C0-1", "C0-2", "C0-3") ~ "Cave_Entrance",
      grepl("^T", SampleName)                    ~ "Soil_Transect",
      Site_Type == "Water"                        ~ "Water",
      Site_Type == "Cave"                         ~ "Cave"
    ),
    id = Site_Number
  )

rownames(metadata) <- metadata$Sampel.ID

otu_raw   <- read.delim("Inputs/feature-table.tsv", skip = 1, row.names = 1, check.names = FALSE)
otu_table <- t(otu_raw) %>% as.data.frame()

cat("Row check — OTU in metadata:", all(rownames(otu_table) %in% rownames(metadata)), "\n")
cat("Row check — metadata in OTU:", all(rownames(metadata) %in% rownames(otu_table)), "\n")

write.table(metadata,  "FEAST_output/Karst_metadata_FEAST.txt", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(otu_table, "FEAST_output/Karst_ASV_FEAST.txt",      sep = "\t", quote = FALSE, row.names = TRUE)

# =============================================================================
# STEP 2 — FEAST ITERATION 1 (Separate surface sub-types)
# =============================================================================

meta_feast <- read.delim("FEAST_output/Karst_metadata_FEAST.txt", row.names = 1)
otu_feast  <- read.delim("FEAST_output/Karst_ASV_FEAST.txt",      row.names = 1, check.names = FALSE)

meta_feast1        <- meta_feast[, c("SourceSink", "Env", "id")]
meta_feast1        <- meta_feast1[rownames(otu_feast), ]
meta_feast1$id     <- ave(seq_len(nrow(meta_feast1)), meta_feast1$Env, FUN = seq_along)

otu_matrix         <- as.matrix(otu_feast)
mode(otu_matrix)   <- "integer"

feast_results <- FEAST(
  C                      = otu_matrix,
  metadata               = meta_feast1,
  different_sources_flag = 0,
  dir_path               = "FEAST_output/",
  outfile                = "Karst_FEAST_iteration1"
)

# FEAST drops you into FEAST_output/ — restore to project root before reading
setwd(project_root)

results <- read.delim(
  "FEAST_output/Karst_FEAST_iteration1_source_contributions_matrix.txt",
  row.names = 1
)

results_collapsed <- data.frame(
  Cave_Entrance = rowSums(results[, grepl("Cave_Entrance", colnames(results)), drop = FALSE]),
  Soil_Transect = rowSums(results[, grepl("Soil_Transect", colnames(results)), drop = FALSE]),
  Water         = rowSums(results[, grepl("Water",         colnames(results)), drop = FALSE]),
  Unknown       = results$Unknown
)
rownames(results_collapsed) <- rownames(results)

cat("\nIteration 1 results preview:\n")
print(round(results_collapsed, 3))

results_export <- cbind(SampleID = rownames(results_collapsed), round(results_collapsed, 4))
write_xlsx(results_export, "FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx")
cat("Saved: FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx\n")

# =============================================================================
# STEP 3 — FEAST ITERATION 2 (Collapsed surface sources)
# =============================================================================

meta_feast2      <- meta_feast[, c("SourceSink", "Env", "id")]
meta_feast2$Env  <- case_when(
  meta_feast$Site_Type == "Surface" ~ "Surface_All",
  meta_feast$Site_Type == "Water"   ~ "Water",
  meta_feast$Site_Type == "Cave"    ~ "Cave"
)
meta_feast2      <- meta_feast2[rownames(otu_feast), ]
meta_feast2$id   <- ave(seq_len(nrow(meta_feast2)), meta_feast2$Env, FUN = seq_along)

feast_results2 <- FEAST(
  C                      = otu_matrix,
  metadata               = meta_feast2,
  different_sources_flag = 0,
  dir_path               = "FEAST_output/",
  outfile                = "Karst_FEAST_iteration2"
)

# Restore again after the second FEAST call
setwd(project_root)

results2 <- read.delim(
  "FEAST_output/Karst_FEAST_iteration2_source_contributions_matrix.txt",
  row.names = 1
)

results_collapsed2 <- data.frame(
  Surface_All = rowSums(results2[, grepl("Surface_All", colnames(results2)), drop = FALSE]),
  Water       = rowSums(results2[, grepl("Water",       colnames(results2)), drop = FALSE]),
  Unknown     = results2$Unknown
)
rownames(results_collapsed2) <- rownames(results2)

cat("\nIteration 2 results preview:\n")
print(round(results_collapsed2, 3))

results_export2 <- cbind(SampleID = rownames(results_collapsed2), round(results_collapsed2, 4))
write_xlsx(results_export2, "FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx")
cat("Saved: FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx\n")

# =============================================================================
# STEP 4 — MERGE ITERATION 2 WITH METADATA
# =============================================================================

meta_for_merge <- metadata %>%
  select(Sampel.ID, Site_Number, Distance_Datum_M, Water_Percantage) %>%
  rename(
    SampleID              = Sampel.ID,
    SiteNumber            = Site_Number,
    DistanceFromEntranceM = Distance_Datum_M,
    WaterPercent          = Water_Percantage
  )

iter2_with_meta <- results_export2 %>%
  as.data.frame() %>%
  mutate(SampleID = sub("_Cave$", "", SampleID)) %>%
  left_join(meta_for_merge, by = "SampleID") %>%
  mutate(across(c(Surface_All, Water, Unknown, DistanceFromEntranceM, WaterPercent), as.numeric))

write_xlsx(iter2_with_meta, "FEAST_output/Karst_FEAST_iteration2_Combined.xlsx")
cat("Saved: FEAST_output/Karst_FEAST_iteration2_Combined.xlsx\n")

df <- iter2_with_meta

# =============================================================================
# STEP 5 — STATISTICAL ANALYSIS
# =============================================================================

sources <- c("Surface_All", "Water", "Unknown")
sig     <- function(p) ifelse(p < 0.05, "YES ***", ifelse(p < 0.10, "TRENDING ~", "NO"))

cat("\n=== SECTION 1: Source Contributions ~ Distance from Entrance ===\n")
corr_distance <- lapply(sources, function(s) {
  test <- cor.test(df[[s]], df$DistanceFromEntranceM, method = "spearman")
  data.frame(
    Source      = s,
    Spearman_r  = round(test$estimate, 3),
    P_value     = round(test$p.value,  4),
    Significant = sig(test$p.value),
    Trend       = ifelse(test$estimate > 0, "Increases with distance", "Decreases with distance")
  )
}) %>% bind_rows()
print(corr_distance)

cat("\n=== SECTION 2: Source Contributions ~ Moisture Content ===\n")
corr_moisture <- lapply(sources, function(s) {
  test <- cor.test(df[[s]], df$WaterPercent, method = "spearman")
  data.frame(
    Source      = s,
    Spearman_r  = round(test$estimate, 3),
    P_value     = round(test$p.value,  4),
    Significant = sig(test$p.value),
    Trend       = ifelse(test$estimate > 0, "Increases with moisture", "Decreases with moisture")
  )
}) %>% bind_rows()
print(corr_moisture)

cat("\n=== SECTION 3: Near vs Far Zone Comparison (threshold 60 m) ===\n")
df$Zone <- ifelse(df$DistanceFromEntranceM < 60, "Near", "Far")
print(table(df$Zone))

zone_results <- lapply(sources, function(s) {
  test      <- wilcox.test(df[[s]] ~ df$Zone)
  near_mean <- mean(df[[s]][df$Zone == "Near"], na.rm = TRUE)
  far_mean  <- mean(df[[s]][df$Zone == "Far"],  na.rm = TRUE)
  data.frame(
    Source      = s,
    Near_mean   = round(near_mean,    3),
    Far_mean    = round(far_mean,     3),
    W_statistic = test$statistic,
    P_value     = round(test$p.value, 4),
    Significant = sig(test$p.value)
  )
}) %>% bind_rows()
print(zone_results)

replicate_summary <- df %>%
  group_by(SiteNumber) %>%
  summarise(
    Distance_mean = round(mean(DistanceFromEntranceM, na.rm = TRUE), 1),
    n_reps        = n(),
    Surface_mean  = round(mean(Surface_All, na.rm = TRUE), 3),
    Surface_sd    = round(sd(Surface_All,   na.rm = TRUE), 3),
    Surface_CV    = round(sd(Surface_All)   / mean(Surface_All)   * 100, 1),
    Water_mean    = round(mean(Water,        na.rm = TRUE), 3),
    Water_sd      = round(sd(Water,          na.rm = TRUE), 3),
    Water_CV      = round(sd(Water)          / mean(Water)         * 100, 1),
    Unknown_mean  = round(mean(Unknown,      na.rm = TRUE), 3),
    Unknown_sd    = round(sd(Unknown,        na.rm = TRUE), 3),
    Unknown_CV    = round(sd(Unknown)        / mean(Unknown)       * 100, 1)
  ) %>%
  arrange(Distance_mean)

cat("\nUnknown fraction summary:\n")
cat(sprintf("  Mean:   %.3f (%.1f%%)\n", mean(df$Unknown,   na.rm = TRUE), mean(df$Unknown,   na.rm = TRUE) * 100))
cat(sprintf("  Median: %.3f (%.1f%%)\n", median(df$Unknown, na.rm = TRUE), median(df$Unknown, na.rm = TRUE) * 100))
cat(sprintf("  Range:  %.3f - %.3f\n",   min(df$Unknown,    na.rm = TRUE), max(df$Unknown,    na.rm = TRUE)))

write_xlsx(
  list(
    Corr_Distance     = corr_distance,
    Corr_Moisture     = corr_moisture,
    Zone_Comparison   = zone_results,
    Replicate_Summary = as.data.frame(replicate_summary),
    Raw_Data          = as.data.frame(df)
  ),
  "FEAST_output/FEASTV2_Full_Statistical_Summary.xlsx"
)
cat("Saved: FEAST_output/FEASTV2_Full_Statistical_Summary.xlsx\n")

# =============================================================================
# FIGURES 8a, 8b, 9 — FEAST BARPLOTS AND DISTANCE SCATTERPLOT
# =============================================================================

iter1 <- read_excel("FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx")
iter2 <- read_excel("FEAST_output/Karst_FEAST_iteration2_collapsed.xlsx")

phy.filtered.DNA.rel <- readRDS("phyloseq_filtered_RA_DNA.rds")
meta_fig <- data.frame(sample_data(phy.filtered.DNA.rel)) %>%
  tibble::rownames_to_column("SampleID_raw") %>%
  mutate(SampleID_clean = sub("_Cave$", "", SampleID_raw))

prep_feast <- function(dat, source_cols, labels) {
  dat %>%
    mutate(SampleID_clean = sub("_Cave$", "", SampleID)) %>%
    left_join(
      meta_fig %>% select(SampleID_clean, SampleName, Site_Number, Distance_Datum_M),
      by = "SampleID_clean"
    ) %>%
    mutate(
      Distance_Datum_M = as.numeric(Distance_Datum_M),
      Site_label       = paste0("Site ", Site_Number),
      SampleName       = fct_reorder(SampleName, Distance_Datum_M)
    ) %>%
    pivot_longer(cols = all_of(source_cols), names_to = "Source", values_to = "Proportion") %>%
    mutate(Source = recode(Source, !!!labels))
}

# ── Figure 8a ────────────────────────────────────────────────────────────────

iter1_long <- prep_feast(
  iter1,
  source_cols = c("Cave_Entrance", "Soil_Transect", "Water", "Unknown"),
  labels      = c(
    Cave_Entrance = "Cave Entrance",
    Soil_Transect = "Soil Transect",
    Water         = "Water",
    Unknown       = "Resident Cave"
  )
)

iter1_long$Source <- factor(iter1_long$Source,
                            levels = c("Resident Cave", "Water", "Cave Entrance", "Soil Transect"))

feast1_colors <- c(
  "Resident Cave" = "#3B1F6B",
  "Water"         = "#4A90C4",
  "Cave Entrance" = "#4A8C3F",
  "Soil Transect" = "#A8C97A"
)

p8a <- ggplot(iter1_long, aes(x = SampleName, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(~ Site_label, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0), limits = c(0, 1.001)) +
  scale_fill_manual(values = feast1_colors, name = "Microbial source") +
  labs(
    title = "Figure 8a. FEAST Source Tracking — Iteration 1 (Separate Surface Sources)",
    x     = "Sample",
    y     = "Proportional contribution"
  ) +
  annotate("text", x = -Inf, y = -Inf,
           label    = expression(paste("\u2190", " Increasing distance from cave entrance ", "\u2192")),
           hjust    = -0.05, vjust = 2.8,
           size     = 3.2, fontface = "italic", color = "grey40") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    panel.spacing    = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(face = "bold", size = 11),
    axis.line.y      = element_line(color = "grey40", linewidth = 0.4),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin      = margin(10, 15, 20, 10)
  )

ggsave(file.path(out_dir, "Figure8a_FEAST_Iteration1.png"), p8a, width = 14, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure8a_FEAST_Iteration1.pdf"), p8a, width = 14, height = 6)
cat("Saved: Figure 8a\n")

# ── Figure 8b ────────────────────────────────────────────────────────────────

iter2_long <- prep_feast(
  iter2,
  source_cols = c("Surface_All", "Water", "Unknown"),
  labels      = c(
    Surface_All = "Surface",
    Water       = "Water",
    Unknown     = "Resident Cave"
  )
)

iter2_long$Source <- factor(iter2_long$Source,
                            levels = c("Resident Cave", "Water", "Surface"))

feast2_colors <- c(
  "Resident Cave" = "#3B1F6B",
  "Water"         = "#4A90C4",
  "Surface"       = "#4A8C3F"
)

p8b <- ggplot(iter2_long, aes(x = SampleName, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(~ Site_label, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0), limits = c(0, 1.001)) +
  scale_fill_manual(values = feast2_colors, name = "Microbial source") +
  labs(
    title = "Figure 8b. FEAST Source Tracking — Iteration 2 (Collapsed Surface Sources)",
    x     = "Sample",
    y     = "Proportional contribution"
  ) +
  annotate("text", x = -Inf, y = -Inf,
           label    = expression(paste("\u2190", " Increasing distance from cave entrance ", "\u2192")),
           hjust    = -0.05, vjust = 2.8,
           size     = 3.2, fontface = "italic", color = "grey40") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    panel.spacing    = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(face = "bold", size = 11),
    axis.line.y      = element_line(color = "grey40", linewidth = 0.4),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin      = margin(10, 15, 20, 10)
  )

ggsave(file.path(out_dir, "Figure8b_FEAST_Iteration2.png"), p8b, width = 14, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure8b_FEAST_Iteration2.pdf"), p8b, width = 14, height = 6)
cat("Saved: Figure 8b\n")

# ── Figure 9 — Source contributions vs distance ───────────────────────────────

df_iter2 <- read_excel("FEAST_output/Karst_FEAST_iteration2_Combined.xlsx") %>%
  mutate(across(c(Surface_All, Water, Unknown, DistanceFromEntranceM), as.numeric))

stat_surface <- cor.test(df_iter2$Surface_All, df_iter2$DistanceFromEntranceM, method = "spearman", exact = FALSE)
stat_water   <- cor.test(df_iter2$Water,       df_iter2$DistanceFromEntranceM, method = "spearman", exact = FALSE)
stat_unknown <- cor.test(df_iter2$Unknown,     df_iter2$DistanceFromEntranceM, method = "spearman", exact = FALSE)

fmt_stat <- function(test) {
  rho   <- round(test$estimate, 3)
  p     <- test$p.value
  p_lab <- ifelse(p < 0.001, "< 0.001", paste0("= ", round(p, 3)))
  paste0("rho = ", rho, ", p ", p_lab)
}

df9_long <- df_iter2 %>%
  select(SampleID, DistanceFromEntranceM, Surface_All, Water, Unknown) %>%
  pivot_longer(
    cols      = c(Surface_All, Water, Unknown),
    names_to  = "Source",
    values_to = "Proportion"
  ) %>%
  mutate(
    Source = recode(Source,
                    Surface_All = "Surface",
                    Water       = "Water",
                    Unknown     = "Resident Cave"),
    Source = factor(Source, levels = c("Surface", "Water", "Resident Cave"))
  )

annot_9 <- data.frame(
  Source = factor(c("Surface", "Water", "Resident Cave"),
                  levels = c("Surface", "Water", "Resident Cave")),
  label  = c(fmt_stat(stat_surface), fmt_stat(stat_water), fmt_stat(stat_unknown)),
  sig    = c(TRUE, FALSE, FALSE)
)

source_colors <- c(
  "Surface"       = "#4A8C3F",
  "Water"         = "#4A90C4",
  "Resident Cave" = "#3B1F6B"
)

p9 <- ggplot(df9_long, aes(x = DistanceFromEntranceM, y = Proportion, color = Source)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(
    data      = df9_long %>% filter(Source == "Surface"),
    method    = "lm", se = TRUE,
    color     = "#4A8C3F", fill = "#4A8C3F",
    alpha     = 0.15, linewidth = 0.8
  ) +
  geom_text(
    data        = annot_9 %>% filter(sig == TRUE),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.6,
    size = 3.2, fontface = "bold.italic", color = "#B2182B",
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = annot_9 %>% filter(sig == FALSE),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.05, vjust = 1.6,
    size = 3.2, fontface = "italic", color = "grey40",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Source, scales = "free_y", nrow = 1) +
  scale_color_manual(values = source_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(
    name   = "Distance from cave entrance (m)",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    title = "Figure 9. FEAST Source Contributions vs Distance from Cave Entrance",
    y     = "Proportional contribution"
  ) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10),
    axis.title       = element_text(face = "bold", size = 11),
    legend.position  = "none",
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin      = margin(10, 15, 10, 10)
  )

ggsave(file.path(out_dir, "Figure9_FEAST_vs_Distance.png"), p9, width = 12, height = 5, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure9_FEAST_vs_Distance.pdf"), p9, width = 12, height = 5)
cat("Saved: Figure 9\n")

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("FEAST output folder: FEAST_output/\n")
cat("Figure output folder:", out_dir, "\n")
cat("============================\n")