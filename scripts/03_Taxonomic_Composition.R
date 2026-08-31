# =============================================================================
# 03_Taxonomic_Composition.R
# Purpose: Characterise and visualise the taxonomic composition of the full
#          community (all site types) and identify ASVs exclusively detected
#          in cave samples. Produces publication Figures 5, 6, and 13.
#
# Requires: phyloseq_filtered_counts_DNA.rds (from 01_Phyloseq_Setup.R)
#
# Outputs (written to Outputs/Taxonomic_Composition/):
#   Figure5_Top15_Phyla_Barplot.png / .pdf
#   Figure6_Top15_Orders_Barplot.png / .pdf
#   Figure13_CaveExclusive_Waffle.png / .pdf
#   Taxonomic_Abundance_Summary.xlsx
#   CaveExclusive_Taxonomy.xlsx
# =============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(writexl)

# Output folder
out_dir <- "Outputs/Taxonomic_Composition"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
phy.filtered.DNA <- readRDS("phyloseq_filtered_counts_DNA.rds")

# Shared color palette (16 colors — 15 taxa + Other)
taxa_colors <- c(
  "#3B1F6B", "#4A90C4", "#4A8C3F", "#A8C97A", "#C0504D",
  "#F4A460", "#8B6914", "#5B9BD5", "#70AD47", "#9E480E",
  "#997300", "#43682B", "#264478", "#843C0C", "#636363",
  "#D9D9D9"
)

# Helper: clean SILVA prefixes from taxonomy
clean_tax <- function(df) {
  for (col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
    df[[col]] <- gsub("^[a-z]__", "", df[[col]])
    df[[col]][df[[col]] == ""] <- NA
  }
  df
}

# =============================================================================
# SECTION 1 — TAXONOMIC BARPLOTS (ALL SITE TYPES)
# =============================================================================

# ── Phylum ───────────────────────────────────────────────────────────────────

ps_phylum <- tax_glom(phy.filtered.DNA, taxrank = "Phylum")
ps_phylum <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

phylum_long <- psmelt(ps_phylum) %>%
  mutate(
    Phylum           = gsub("p__", "", Phylum),
    Phylum           = ifelse(is.na(Phylum) | Phylum == "", "Unclassified", Phylum),
    Site_Number      = as.numeric(as.character(Site_Number)),
    Distance_Datum_M = as.numeric(as.character(Distance_Datum_M))
  )

top15_phylum <- phylum_long %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Phylum)

phylum_long <- phylum_long %>%
  mutate(Phylum_plot = ifelse(Phylum %in% top15_phylum, Phylum, "Other")) %>%
  arrange(Site_Type, Site_Number, Distance_Datum_M)

sample_order_phylum <- phylum_long %>%
  distinct(Sample, SampleName, Site_Type, Site_Number, Distance_Datum_M) %>%
  arrange(Site_Type, Site_Number, Distance_Datum_M) %>%
  pull(SampleName)

phylum_long$SampleName <- factor(phylum_long$SampleName, levels = unique(sample_order_phylum))

phylum_level_order <- c(
  phylum_long %>%
    filter(Phylum_plot != "Other") %>%
    group_by(Phylum_plot) %>%
    summarise(m = mean(Abundance), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Phylum_plot),
  "Other"
)

phylum_long$Phylum_plot <- factor(
  phylum_long$Phylum_plot,
  levels = c("Other", rev(phylum_level_order[phylum_level_order != "Other"]))
)

phylum_colors <- taxa_colors
names(phylum_colors) <- c(phylum_level_order[1:15], "Other")

p5 <- ggplot(phylum_long, aes(x = SampleName, y = Abundance, fill = Phylum_plot)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(~ Site_Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 1.001)) +
  scale_fill_manual(values = phylum_colors, name = "Phylum",
                    guide  = guide_legend(reverse = TRUE)) +
  labs(title = "Figure 5. Relative Abundance of Top 15 Bacterial Phyla Across Site Types",
       x = NULL, y = "Relative Abundance") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.spacing = unit(0, "lines"),
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.line.y  = element_line(color = "grey40", linewidth = 0.4),
    legend.position = "right", legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "Figure5_Top15_Phyla_Barplot.png"), p5, width = 16, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure5_Top15_Phyla_Barplot.pdf"), p5, width = 16, height = 6)
cat("Saved: Figure 5 — Phylum\n")

# ── Order ─────────────────────────────────────────────────────────────────────

ps_order <- tax_glom(phy.filtered.DNA, taxrank = "Order")
ps_order <- transform_sample_counts(ps_order, function(x) x / sum(x))

order_long <- psmelt(ps_order) %>%
  mutate(
    Order            = gsub("o__", "", Order),
    Order            = ifelse(is.na(Order) | Order == "", "Unclassified", Order),
    Site_Number      = as.numeric(as.character(Site_Number)),
    Distance_Datum_M = as.numeric(as.character(Distance_Datum_M))
  )

top15_order <- order_long %>%
  group_by(Order) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 15) %>%
  pull(Order)

order_long <- order_long %>%
  mutate(Order_plot = ifelse(Order %in% top15_order, Order, "Other")) %>%
  arrange(Site_Type, Site_Number, Distance_Datum_M)

sample_order_order <- order_long %>%
  distinct(Sample, SampleName, Site_Type, Site_Number, Distance_Datum_M) %>%
  arrange(Site_Type, Site_Number, Distance_Datum_M) %>%
  pull(SampleName)

order_long$SampleName <- factor(order_long$SampleName, levels = unique(sample_order_order))

order_level_order <- c(
  order_long %>%
    filter(Order_plot != "Other") %>%
    group_by(Order_plot) %>%
    summarise(m = mean(Abundance), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Order_plot),
  "Other"
)

order_long$Order_plot <- factor(
  order_long$Order_plot,
  levels = c("Other", rev(order_level_order[order_level_order != "Other"]))
)

order_colors <- taxa_colors
names(order_colors) <- c(order_level_order[1:15], "Other")

p6 <- ggplot(order_long, aes(x = SampleName, y = Abundance, fill = Order_plot)) +
  geom_bar(stat = "identity", width = 0.85) +
  facet_grid(~ Site_Type, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, 1.001)) +
  scale_fill_manual(values = order_colors, name = "Order",
                    guide  = guide_legend(reverse = TRUE)) +
  labs(title = "Figure 6. Relative Abundance of Top 15 Bacterial Orders Across Site Types",
       x = NULL, y = "Relative Abundance") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), panel.spacing = unit(0, "lines"),
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.line.y  = element_line(color = "grey40", linewidth = 0.4),
    legend.position = "right", legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_dir, "Figure6_Top15_Orders_Barplot.png"), p6, width = 16, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure6_Top15_Orders_Barplot.pdf"), p6, width = 16, height = 6)
cat("Saved: Figure 6 — Order\n")

# ── Summary tables ────────────────────────────────────────────────────────────

phylum_summary <- phylum_long %>%
  group_by(Site_Type, Phylum_plot) %>%
  summarise(Mean_Abundance = round(mean(Abundance), 4),
            SD_Abundance   = round(sd(Abundance),   4), .groups = "drop") %>%
  arrange(Site_Type, desc(Mean_Abundance))

phylum_wide <- phylum_summary %>%
  pivot_wider(names_from = Site_Type, values_from = c(Mean_Abundance, SD_Abundance)) %>%
  arrange(desc(Mean_Abundance_Cave))

order_summary <- order_long %>%
  group_by(Site_Type, Order_plot) %>%
  summarise(Mean_Abundance = round(mean(Abundance), 4),
            SD_Abundance   = round(sd(Abundance),   4), .groups = "drop") %>%
  arrange(Site_Type, desc(Mean_Abundance))

order_wide <- order_summary %>%
  pivot_wider(names_from = Site_Type, values_from = c(Mean_Abundance, SD_Abundance)) %>%
  arrange(desc(Mean_Abundance_Cave))

write_xlsx(
  list(
    Phylum_by_SiteType = as.data.frame(phylum_wide),
    Phylum_long        = as.data.frame(phylum_summary),
    Phylum_per_sample  = as.data.frame(phylum_long),
    Order_by_SiteType  = as.data.frame(order_wide),
    Order_long         = as.data.frame(order_summary),
    Order_per_sample   = as.data.frame(order_long)
  ),
  file.path(out_dir, "Taxonomic_Abundance_Summary.xlsx")
)
cat("Saved: Taxonomic_Abundance_Summary.xlsx\n")

# =============================================================================
# SECTION 2 — CAVE-EXCLUSIVE ASV TAXONOMY
# =============================================================================

otu_raw <- as.data.frame(otu_table(phy.filtered.DNA))
if (!taxa_are_rows(phy.filtered.DNA)) otu_raw <- t(otu_raw)
otu_raw <- as.data.frame(otu_raw)

meta           <- as.data.frame(as.matrix(sample_data(phy.filtered.DNA)))
cave_samples   <- rownames(meta)[meta$Site_Type == "Cave"]
source_samples <- rownames(meta)[meta$Site_Type != "Cave"]

cave_present        <- rowSums(otu_raw[, cave_samples])   > 0
source_present      <- rowSums(otu_raw[, source_samples]) > 0
cave_exclusive_asvs <- rownames(otu_raw)[cave_present & !source_present]

cat(sprintf("\nCave-exclusive ASVs: %d\n", length(cave_exclusive_asvs)))

tax_full      <- as.data.frame(tax_table(phy.filtered.DNA))
tax_full      <- clean_tax(tax_full)
tax_exclusive <- tax_full[cave_exclusive_asvs, ]
tax_all_cave  <- tax_full[rownames(otu_raw)[cave_present], ]

# Enrichment helper
enrich <- function(tax_sub, tax_all, col) {
  e <- tax_sub %>%
    count(.data[[col]], name = "N_Exclusive") %>%
    mutate(Pct_Exclusive = round(100 * N_Exclusive / sum(N_Exclusive), 2))
  a <- tax_all %>%
    count(.data[[col]], name = "N_All_Cave") %>%
    mutate(Pct_All_Cave  = round(100 * N_All_Cave  / sum(N_All_Cave),  2))
  e %>%
    left_join(a, by = col) %>%
    mutate(Enrichment = round(Pct_Exclusive / Pct_All_Cave, 2)) %>%
    arrange(desc(N_Exclusive))
}

phylum_compare <- enrich(tax_exclusive, tax_all_cave, "Phylum")
class_compare  <- enrich(tax_exclusive, tax_all_cave, "Class")
order_compare  <- enrich(tax_exclusive, tax_all_cave, "Order")

cat("\n=== Phylum enrichment — Exclusive vs All Cave ===\n"); print(phylum_compare)
cat("\n=== Order enrichment (top 25) ===\n"); print(head(order_compare, 25))
cat("\n=== Classification rates ===\n")
cat(sprintf("Phylum: %d / %d (%.1f%%)\n",
            sum(!is.na(tax_exclusive$Phylum)), nrow(tax_exclusive),
            100 * sum(!is.na(tax_exclusive$Phylum)) / nrow(tax_exclusive)))
cat(sprintf("Genus:  %d / %d (%.1f%%)\n",
            sum(!is.na(tax_exclusive$Genus)),  nrow(tax_exclusive),
            100 * sum(!is.na(tax_exclusive$Genus))  / nrow(tax_exclusive)))

write_xlsx(
  list(
    Phylum_Comparison = as.data.frame(phylum_compare),
    Class_Comparison  = as.data.frame(class_compare),
    Order_Comparison  = as.data.frame(order_compare),
    Raw_Taxonomy      = tax_exclusive %>% rownames_to_column("ASV_ID")
  ),
  file.path(out_dir, "CaveExclusive_Taxonomy.xlsx")
)
cat("Saved: CaveExclusive_Taxonomy.xlsx\n")

# =============================================================================
# FIGURE 13 — WAFFLE CHART: CAVE ASV COMPOSITION
# =============================================================================

total_asvs     <- nrow(otu_raw)
in_cave        <- sum(cave_present)
excl           <- length(cave_exclusive_asvs)
shared         <- in_cave - excl
not_in_cave    <- total_asvs - in_cave

square_size <- 100
n_not_cave  <- ceiling(not_in_cave / square_size)
n_shared    <- ceiling(shared      / square_size)
n_exclusive <- ceiling(excl        / square_size)
n_total     <- n_not_cave + n_shared + n_exclusive

ncols <- 12
nrows <- ceiling(n_total / ncols)

waffle_df <- data.frame(
  id       = seq_len(n_total),
  Category = factor(
    c(rep("Not detected in cave",         n_not_cave),
      rep("Shared with surface or water", n_shared),
      rep("Cave-exclusive",               n_exclusive)),
    levels = c("Not detected in cave", "Shared with surface or water", "Cave-exclusive")
  )
) %>%
  mutate(col = (id - 1) %% ncols + 1,
         row = (id - 1) %/% ncols + 1)

waffle_colors <- c(
  "Not detected in cave"         = "#D9D9D9",
  "Shared with surface or water" = "#9090B8",
  "Cave-exclusive"               = "#3B1F6B"
)

legend_labels <- c(
  "Not detected in cave"         = paste0("Not detected in cave\n(n = ", not_in_cave,
                                           " | ", round(100 * not_in_cave / total_asvs, 1), "%)"),
  "Shared with surface or water" = paste0("Shared with surface or water\n(n = ", shared,
                                           " | ", round(100 * shared / total_asvs, 1), "%)"),
  "Cave-exclusive"               = paste0("Cave-exclusive\n(n = ", excl,
                                           " | ", round(100 * excl / total_asvs, 1), "%)")
)

p13 <- ggplot(waffle_df, aes(x = col, y = -row, fill = Category)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.9, height = 0.9) +
  scale_fill_manual(values = waffle_colors, labels = legend_labels, name = NULL) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  annotate("text",
           x = ncols / 2 + 0.5, y = -(nrows + 0.8),
           label    = paste0("Each square represents ~100 ASVs   |   Total = ", total_asvs, " ASVs"),
           size = 3.2, fontface = "italic", color = "grey40", hjust = 0.5) +
  coord_equal(clip = "off") +
  labs(title    = "Figure 13. Composition of Cave-Associated ASVs",
       subtitle = "Partitioned by presence across cave, surface, and water environments") +
  theme_void() +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13, margin = margin(b = 4)),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey40", margin = margin(b = 12)),
    legend.position  = "bottom", legend.direction = "vertical",
    legend.text      = element_text(size = 10, lineheight = 1.3),
    legend.key.size  = unit(0.9, "lines"), legend.spacing.y = unit(0.4, "cm"),
    plot.margin      = margin(15, 15, 30, 15)
  ) +
  guides(fill = guide_legend(byrow = TRUE))

ggsave(file.path(out_dir, "Figure13_CaveExclusive_Waffle.png"), p13, width = 7, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure13_CaveExclusive_Waffle.pdf"), p13, width = 7, height = 6)
cat("Saved: Figure 13 — Waffle Chart\n")

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("============================\n")
