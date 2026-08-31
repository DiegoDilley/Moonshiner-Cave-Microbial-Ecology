# =============================================================================
# 04_Alpha_Diversity.R
# Purpose: Calculate alpha diversity metrics (Shannon, Simpson, Richness) for
#          all samples, run statistical comparisons across site types, test
#          correlations with cave depth and FEAST Unknown fraction, and quantify
#          diversity of cave-exclusive ASVs. Produces Figures 7 and 10.
#
# Requires:
#   phyloseq_filtered_counts_DNA.rds        (from 01_Phyloseq_Setup.R)
#   FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx  (from 05_FEAST_Source_Tracking.R)
#
# Outputs (written to Outputs/Alpha_Diversity/):
#   Figure7_AlphaDiversity_Boxplots.png / .pdf
#   Figure10_AlphaVsFEAST_Unknown.png / .pdf
#   Karst_AlphaDiversity.xlsx
#   AlphaDiversity_Stats.xlsx
#   CaveExclusive_Diversity.xlsx
# =============================================================================

library(phyloseq)
library(vegan)
library(writexl)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# Output folder
out_dir <- "Outputs/Alpha_Diversity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
phy.filtered.DNA <- readRDS("phyloseq_filtered_counts_DNA.rds")

# Shared color palette
site_colors <- c(
  "Cave"    = "#3B1F6B",
  "Surface" = "#4A8C3F",
  "Water"   = "#4A90C4"
)

# =============================================================================
# STEP 1 — CALCULATE ALPHA DIVERSITY
# =============================================================================

otu <- as.data.frame(t(otu_table(phy.filtered.DNA)))

alpha_df <- data.frame(
  SampleID = rownames(otu),
  Shannon  = vegan::diversity(otu, index = "shannon"),
  Simpson  = vegan::diversity(otu, index = "simpson"),
  Richness = vegan::specnumber(otu)
)

meta <- as.data.frame(as.matrix(sample_data(phy.filtered.DNA)))
meta$SampleID <- rownames(meta)

alpha_df <- alpha_df %>%
  left_join(
    meta %>% select(SampleID, SampleName, Site_Type, Site_Number,
                    Distance_Datum_M, Water_Percantage),
    by = "SampleID"
  )

cat("Site_Type check:\n")
print(table(alpha_df$Site_Type, useNA = "always"))

# =============================================================================
# STEP 2 — APPEND FEAST UNKNOWN FRACTION
# =============================================================================

feast <- read_excel("FEAST_output/Karst_FEAST_iteration1_collapsed.xlsx")

feast <- feast %>%
  mutate(SampleID_short = sub("_Cave$", "", SampleID),
         SampleID_short = sub("-S.*",   "", SampleID_short))

alpha_df <- alpha_df %>%
  mutate(SampleID_short = sub("-S.*", "", SampleID)) %>%
  left_join(
    feast %>%
      select(SampleID_short, Unknown, Water, Cave_Entrance, Soil_Transect) %>%
      rename(FEAST_Unknown       = Unknown,
             FEAST_Water         = Water,
             FEAST_Cave_Entrance = Cave_Entrance,
             FEAST_Soil_Transect = Soil_Transect),
    by = "SampleID_short"
  ) %>%
  select(-SampleID_short)

cat("\nSamples with FEAST data:\n")
print(alpha_df %>% filter(!is.na(FEAST_Unknown)) %>%
        select(SampleID, Site_Type, FEAST_Unknown))

write_xlsx(alpha_df, file.path(out_dir, "Karst_AlphaDiversity.xlsx"))
cat("\nSaved: Karst_AlphaDiversity.xlsx\n")

# =============================================================================
# STEP 3 — STATISTICAL TESTS
# =============================================================================

df <- alpha_df
df$Shannon          <- as.numeric(df$Shannon)
df$Simpson          <- as.numeric(df$Simpson)
df$Richness         <- as.numeric(df$Richness)
df$Distance_Datum_M <- as.numeric(df$Distance_Datum_M)
df$FEAST_Unknown    <- as.numeric(df$FEAST_Unknown)

metrics <- c("Shannon", "Simpson", "Richness")
sig     <- function(p) ifelse(p < 0.05, "YES ***", ifelse(p < 0.10, "TRENDING ~", "NO"))

# Q1: Does site type influence alpha diversity?
cat("=================================================================\n")
cat("Q1: Alpha Diversity by Site Type (Kruskal-Wallis)\n")
cat("=================================================================\n")

kw_results <- lapply(metrics, function(m) {
  test <- kruskal.test(df[[m]] ~ df$Site_Type)
  data.frame(
    Metric      = m,
    H_statistic = round(test$statistic, 3),
    df          = test$parameter,
    P_value     = round(test$p.value, 4),
    Significant = sig(test$p.value)
  )
}) %>% bind_rows()
print(kw_results)

sig_metrics <- kw_results %>% filter(P_value < 0.05) %>% pull(Metric)

pairwise_results <- lapply(sig_metrics, function(m) {
  pairs <- combn(unique(df$Site_Type), 2, simplify = FALSE)
  lapply(pairs, function(p) {
    g1   <- df[[m]][df$Site_Type == p[1]]
    g2   <- df[[m]][df$Site_Type == p[2]]
    test <- wilcox.test(g1, g2)
    data.frame(
      Metric      = m, Group1 = p[1], Group2 = p[2],
      Mean_G1     = round(mean(g1, na.rm = TRUE), 3),
      Mean_G2     = round(mean(g2, na.rm = TRUE), 3),
      P_value     = round(test$p.value, 4),
      Significant = sig(test$p.value)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

print(pairwise_results)

means <- df %>%
  group_by(Site_Type) %>%
  summarise(
    N             = n(),
    Mean_Shannon  = round(mean(Shannon,  na.rm = TRUE), 3),
    SD_Shannon    = round(sd(Shannon,    na.rm = TRUE), 3),
    Mean_Simpson  = round(mean(Simpson,  na.rm = TRUE), 3),
    SD_Simpson    = round(sd(Simpson,    na.rm = TRUE), 3),
    Mean_Richness = round(mean(Richness, na.rm = TRUE), 1),
    SD_Richness   = round(sd(Richness,   na.rm = TRUE), 1)
  )
print(means)

# Q2: Does distance influence cave alpha diversity?
cat("\n=================================================================\n")
cat("Q2: Alpha Diversity vs Distance from Entrance (Cave only)\n")
cat("=================================================================\n")

cave <- df %>% filter(Site_Type == "Cave")

dist_results <- lapply(metrics, function(m) {
  test <- cor.test(cave$Distance_Datum_M, cave[[m]], method = "spearman", exact = FALSE)
  data.frame(
    Metric      = m,
    rho         = round(test$estimate, 3),
    P_value     = round(test$p.value, 4),
    Direction   = ifelse(test$estimate > 0, "Increases with depth", "Decreases with depth"),
    Significant = sig(test$p.value)
  )
}) %>% bind_rows()
print(dist_results)

# Q3: Does FEAST Unknown correlate with alpha diversity?
cat("\n=================================================================\n")
cat("Q3: FEAST Unknown Fraction vs Alpha Diversity (Cave only)\n")
cat("=================================================================\n")

unknown_results <- lapply(metrics, function(m) {
  test <- cor.test(cave$FEAST_Unknown, cave[[m]], method = "spearman", exact = FALSE)
  data.frame(
    Metric      = m,
    rho         = round(test$estimate, 3),
    P_value     = round(test$p.value, 4),
    Direction   = ifelse(test$estimate > 0, "Higher Unknown = Higher diversity",
                         "Higher Unknown = Lower diversity"),
    Significant = sig(test$p.value)
  )
}) %>% bind_rows()
print(unknown_results)

write_xlsx(
  list(
    Q1_KruskalWallis = as.data.frame(kw_results),
    Q1_Pairwise      = as.data.frame(pairwise_results),
    Q1_Means         = as.data.frame(means),
    Q2_Distance      = as.data.frame(dist_results),
    Q3_Unknown       = as.data.frame(unknown_results)
  ),
  file.path(out_dir, "AlphaDiversity_Stats.xlsx")
)
cat("\nSaved: AlphaDiversity_Stats.xlsx\n")

# =============================================================================
# STEP 4 — CAVE-EXCLUSIVE ASV DIVERSITY
# =============================================================================

otu_raw <- as.data.frame(otu_table(phy.filtered.DNA))
if (!taxa_are_rows(phy.filtered.DNA)) otu_raw <- t(otu_raw)
otu_raw <- as.data.frame(otu_raw)

meta2          <- as.data.frame(as.matrix(sample_data(phy.filtered.DNA)))
cave_samples   <- rownames(meta2)[meta2$Site_Type == "Cave"]
source_samples <- rownames(meta2)[meta2$Site_Type != "Cave"]

cave_present        <- rowSums(otu_raw[, cave_samples])   > 0
source_present      <- rowSums(otu_raw[, source_samples]) > 0
cave_exclusive_asvs <- rownames(otu_raw)[cave_present & !source_present]
shared_asvs         <- rownames(otu_raw)[cave_present &  source_present]

cat(sprintf("\nTotal ASVs:              %d\n", nrow(otu_raw)))
cat(sprintf("ASVs in cave:            %d\n", sum(cave_present)))
cat(sprintf("Cave-exclusive ASVs:     %d\n", length(cave_exclusive_asvs)))
cat(sprintf("Shared with sources:     %d\n", length(shared_asvs)))
cat(sprintf("Cave-exclusive %%:       %.1f%%\n",
            100 * length(cave_exclusive_asvs) / sum(cave_present)))

otu_exclusive_t <- as.data.frame(t(otu_raw[cave_exclusive_asvs, cave_samples]))

exclusive_alpha <- data.frame(
  SampleID           = rownames(otu_exclusive_t),
  Shannon_Exclusive  = vegan::diversity(otu_exclusive_t, index = "shannon"),
  Richness_Exclusive = vegan::specnumber(otu_exclusive_t)
)

otu_all_t <- as.data.frame(t(otu_raw))

all_alpha <- data.frame(
  SampleID = rownames(otu_all_t),
  Shannon  = vegan::diversity(otu_all_t, index = "shannon"),
  Richness = vegan::specnumber(otu_all_t)
) %>%
  left_join(meta2 %>% select(Site_Type) %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  left_join(exclusive_alpha, by = "SampleID")

summary_table <- all_alpha %>%
  group_by(Site_Type) %>%
  summarise(
    N                       = n(),
    Mean_Richness_Full      = round(mean(Richness,           na.rm = TRUE), 1),
    SD_Richness_Full        = round(sd(Richness,             na.rm = TRUE), 1),
    Mean_Shannon_Full       = round(mean(Shannon,            na.rm = TRUE), 3),
    Mean_Richness_Exclusive = round(mean(Richness_Exclusive, na.rm = TRUE), 1),
    SD_Richness_Exclusive   = round(sd(Richness_Exclusive,   na.rm = TRUE), 1),
    Mean_Shannon_Exclusive  = round(mean(Shannon_Exclusive,  na.rm = TRUE), 3)
  )

print(summary_table)

write_xlsx(
  list(
    Summary         = as.data.frame(summary_table),
    Per_Sample      = as.data.frame(all_alpha),
    Exclusive_ASV_N = data.frame(
      Category = c("Total ASVs", "Cave present", "Cave exclusive",
                   "Cave shared with sources", "Cave exclusive %"),
      Value    = c(nrow(otu_raw), sum(cave_present), length(cave_exclusive_asvs),
                   length(shared_asvs),
                   round(100 * length(cave_exclusive_asvs) / sum(cave_present), 1))
    )
  ),
  file.path(out_dir, "CaveExclusive_Diversity.xlsx")
)
cat("\nSaved: CaveExclusive_Diversity.xlsx\n")

# =============================================================================
# FIGURE 7 — ALPHA DIVERSITY BOXPLOTS BY SITE TYPE
# =============================================================================

alpha_long <- df %>%
  select(SampleID, Site_Type, Shannon, Simpson, Richness) %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness),
               names_to  = "Metric",
               values_to = "Value") %>%
  mutate(
    Metric    = factor(Metric, levels = c("Richness", "Shannon", "Simpson")),
    Site_Type = factor(Site_Type, levels = c("Cave", "Surface", "Water"))
  )

y_max       <- alpha_long %>% group_by(Metric) %>% summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
richness_ym <- y_max %>% filter(Metric == "Richness") %>% pull(ymax)
shannon_ym  <- y_max %>% filter(Metric == "Shannon")  %>% pull(ymax)
simpson_ym  <- y_max %>% filter(Metric == "Simpson")  %>% pull(ymax)

bracket_df <- data.frame(
  Metric  = factor(
    c("Richness","Richness","Richness","Shannon","Shannon","Shannon","Simpson","Simpson","Simpson"),
    levels = c("Richness","Shannon","Simpson")),
  group1  = rep(c("Cave",    "Water",    "Cave"),  3),
  group2  = rep(c("Surface", "Surface",  "Water"), 3),
  p_label = c(
    "p < 0.001", "p = 0.001", "p = 0.087~",
    "p = 0.113",  "p = 0.073~", "p = 0.491",
    "p = 0.317",  "p = 0.167",  "p = 0.491"
  ),
  y = c(
    richness_ym * 1.05, richness_ym * 1.15, richness_ym * 1.25,
    shannon_ym  * 1.03, shannon_ym  * 1.08, shannon_ym  * 1.13,
    simpson_ym  * 1.03, simpson_ym  * 1.08, simpson_ym  * 1.13
  )
)

x_pos      <- c("Cave" = 1, "Surface" = 2, "Water" = 3)
bracket_df <- bracket_df %>%
  mutate(x_mid = (x_pos[group1] + x_pos[group2]) / 2)

p7 <- ggplot(alpha_long, aes(x = Site_Type, y = Value, fill = Site_Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(color = Site_Type), width = 0.15, size = 2, alpha = 0.8, show.legend = FALSE) +
  geom_segment(data = bracket_df,
               aes(x = group1, xend = group2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_segment(data = bracket_df,
               aes(x = group1, xend = group1, y = y, yend = y * 0.99),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_segment(data = bracket_df,
               aes(x = group2, xend = group2, y = y, yend = y * 0.99),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = bracket_df,
            aes(x = x_mid, y = y * 1.02, label = p_label),
            inherit.aes = FALSE, size = 2.8, color = "grey20", fontface = "italic") +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values  = site_colors, name = "Site Type") +
  scale_color_manual(values = site_colors) +
  coord_cartesian(clip = "off") +
  labs(title = "Figure 7. Alpha Diversity Across Site Types", x = NULL, y = "Diversity value") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "none",
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin = margin(10, 10, 25, 10)
  )

ggsave(file.path(out_dir, "Figure7_AlphaDiversity_Boxplots.png"), p7, width = 10, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure7_AlphaDiversity_Boxplots.pdf"), p7, width = 10, height = 6)
cat("Saved: Figure 7\n")

# =============================================================================
# FIGURE 10 — ALPHA DIVERSITY VS FEAST UNKNOWN FRACTION (CAVE ONLY)
# =============================================================================

cave_df <- df %>% filter(Site_Type == "Cave", !is.na(FEAST_Unknown))

stat_shannon  <- cor.test(cave_df$Shannon,  cave_df$FEAST_Unknown, method = "spearman", exact = FALSE)
stat_richness <- cor.test(cave_df$Richness, cave_df$FEAST_Unknown, method = "spearman", exact = FALSE)

fmt_stat <- function(test) {
  rho   <- round(test$estimate, 3)
  p     <- test$p.value
  p_lab <- ifelse(p < 0.001, "< 0.001", paste0("= ", round(p, 3)))
  paste0("rho = ", rho, ", p ", p_lab)
}

cave_long <- cave_df %>%
  select(SampleID, FEAST_Unknown, Shannon, Richness) %>%
  pivot_longer(cols = c(Shannon, Richness), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Richness", "Shannon")))

annot_10 <- data.frame(
  Metric = factor(c("Richness", "Shannon"), levels = c("Richness", "Shannon")),
  label  = c(fmt_stat(stat_richness), fmt_stat(stat_shannon))
)

p10 <- ggplot(cave_long, aes(x = FEAST_Unknown, y = Value)) +
  geom_point(color = "#3B1F6B", size = 3.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#3B1F6B", fill = "#3B1F6B", alpha = 0.15, linewidth = 0.8) +
  geom_text(data    = annot_10,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.6,
            size = 3.2, fontface = "bold.italic", color = "#B2182B", inherit.aes = FALSE) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     name   = "Resident cave fraction (FEAST Unknown)") +
  labs(title = "Figure 10. Alpha Diversity vs Resident Cave Fraction (Cave Samples Only)",
       y = "Diversity value") +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11),
    axis.text  = element_text(size = 10), axis.title = element_text(face = "bold", size = 11),
    legend.position = "none",
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave(file.path(out_dir, "Figure10_AlphaVsFEAST_Unknown.png"), p10, width = 10, height = 5, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure10_AlphaVsFEAST_Unknown.pdf"), p10, width = 10, height = 5)
cat("Saved: Figure 10\n")

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("============================\n")
