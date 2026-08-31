# =============================================================================
# 02_Beta_Diversity.R
# Purpose: Beta-diversity ordination (NMDS) and associated statistical tests
#          (PERMANOVA, PERMDISP, envfit) for all sites and cave-only samples.
#          Also produces publication-ready Figures 3 and 4.
#
# Requires: phyloseq_filtered_RA_DNA.rds (from 01_Phyloseq_Setup.R)
#
# Outputs (written to Outputs/Beta_Diversity/):
#   Figure3_NMDS_AllSites.png / .pdf
#   Figure4_NMDS_Cave_Envfit.png / .pdf
#   Cave_Metadata_Corrplot.png
#   NMDS_Statistical_Summary.xlsx
# =============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(pairwiseAdonis)
library(corrplot)
library(writexl)

# Output folder
out_dir <- "Outputs/Beta_Diversity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
phy.filtered.DNA.rel <- readRDS("phyloseq_filtered_RA_DNA.rds")

# Shared color palette
site_colors <- c(
  "Cave"    = "#3B1F6B",
  "Surface" = "#4A8C3F",
  "Water"   = "#4A90C4"
)

# =============================================================================
# STEP 1 — ORDINATION (FULL DATASET)
# =============================================================================

phy.ord <- ordinate(phy.filtered.DNA.rel, method = "NMDS", distance = "bray")

# =============================================================================
# STEP 2 — PERMANOVA + BETADISPERSION (ALL SAMPLES)
# =============================================================================

otu_mat   <- as.data.frame(t(otu_table(phy.filtered.DNA.rel)))
meta_df   <- data.frame(sample_data(phy.filtered.DNA.rel))
meta_df   <- meta_df[rownames(otu_mat), ]
bray_dist <- vegdist(otu_mat, method = "bray")

set.seed(123)
permanova_result <- adonis2(bray_dist ~ Site_Type, data = meta_df, permutations = 999)
pairwise_result  <- pairwise.adonis2(bray_dist ~ Site_Type, data = meta_df, permutations = 999)

cat("\n=== PERMANOVA — All Samples ===\n")
print(permanova_result)
cat("\n=== Pairwise PERMANOVA ===\n")
print(pairwise_result)

# Betadispersion
set.seed(123)
betadisp       <- betadisper(bray_dist, meta_df$Site_Type)
betadisp_test  <- permutest(betadisp, permutations = 999)
betadisp_tukey <- TukeyHSD(betadisp)

cat("\n=== BETADISPERSION (PERMDISP) ===\n")
print(betadisp_test)
print(betadisp_tukey)

betadisp_summary <- data.frame(
  Site_Type     = names(betadisp$group.distances),
  Mean_Distance = round(betadisp$group.distances, 4)
)

betadisp_overall <- data.frame(
  Test        = "PERMDISP overall",
  F_statistic = round(betadisp_test$tab$F[1], 3),
  P_value     = round(betadisp_test$tab$`Pr(>F)`[1], 4),
  Significant = ifelse(betadisp_test$tab$`Pr(>F)`[1] < 0.05, "YES ***",
                       ifelse(betadisp_test$tab$`Pr(>F)`[1] < 0.10, "TRENDING ~", "NO")),
  Interpretation = ifelse(betadisp_test$tab$`Pr(>F)`[1] < 0.05,
                          "Groups differ in dispersion — PERMANOVA may partly reflect variance",
                          "Groups do not differ in dispersion — PERMANOVA reflects centroid differences")
)

# =============================================================================
# STEP 3 — BUILD ANNOTATION STRINGS
# =============================================================================

r2_overall      <- round(permanova_result$R2[1], 3)
f_overall       <- round(permanova_result$F[1], 2)
p_overall       <- permanova_result$`Pr(>F)`[1]
p_overall_label <- ifelse(p_overall < 0.001, "< 0.001", round(p_overall, 3))

overall_label <- paste0(
  "PERMANOVA  R² = ", r2_overall,
  ",  F = ", f_overall,
  ",  p ", p_overall_label
)

pw_cv <- pairwise_result$Cave_vs_Surface
pw_cw <- pairwise_result$Cave_vs_Water
pw_sw <- pairwise_result$Surface_vs_Water

fmt_pw <- function(pw, label) {
  r2    <- round(pw$R2[1], 3)
  p     <- pw$`Pr(>F)`[1]
  p_lab <- ifelse(p < 0.001, "< 0.001", round(p, 3))
  paste0(label, ":  R² = ", r2, ", p ", p_lab)
}

pairwise_label <- paste0(
  "Pairwise PERMANOVA\n",
  fmt_pw(pw_cv, "Cave vs Surface"), "\n",
  fmt_pw(pw_cw, "Cave vs Water"),   "\n",
  fmt_pw(pw_sw, "Surface vs Water")
)

# =============================================================================
# STEP 4 — CAVE-ONLY NMDS + ENVFIT
# =============================================================================

phy.cave     <- subset_samples(phy.filtered.DNA.rel, Site_Type == "Cave")
phy.cave.ord <- ordinate(phy.cave, method = "NMDS", distance = "bray")
cave.meta    <- data.frame(sample_data(phy.cave))

set.seed(123)
env_fit <- envfit(
  phy.cave.ord$points,
  cave.meta[, "Distance_Datum_M", drop = FALSE],
  permutations = 999,
  na.rm = TRUE
)

cat("\n=== envfit — Distance from Entrance (Cave only) ===\n")
print(env_fit)

# =============================================================================
# STEP 5 — CAVE-ONLY PERMANOVA (DISTANCE AS PREDICTOR)
# =============================================================================

otu_cave  <- as.data.frame(t(otu_table(phy.cave)))
meta_cave <- data.frame(sample_data(phy.cave))
meta_cave <- meta_cave[rownames(otu_cave), ]
meta_cave$Distance_Datum_M <- as.numeric(meta_cave$Distance_Datum_M)
bray_cave <- vegdist(otu_cave, method = "bray")

set.seed(123)
permanova_cave <- adonis2(bray_cave ~ Distance_Datum_M,
                          data         = meta_cave,
                          permutations = 999)

cat("\n=== PERMANOVA — Cave only, Distance as predictor ===\n")
print(permanova_cave)

cave_permanova_df <- data.frame(
  Predictor   = "Distance_Datum_M",
  R2          = round(permanova_cave$R2[1], 3),
  F_statistic = round(permanova_cave$F[1], 3),
  P_value     = permanova_cave$`Pr(>F)`[1],
  Significant = ifelse(permanova_cave$`Pr(>F)`[1] < 0.05, "YES ***",
                       ifelse(permanova_cave$`Pr(>F)`[1] < 0.10, "TRENDING ~", "NO"))
)

# =============================================================================
# STEP 6 — METADATA CORRELATION (CAVE)
# =============================================================================

cor_mat <- cor(cave.meta[, c("Distance_Datum_M", "Depth_Surface_M", "Water_Percantage")],
               method = "spearman")
cat("\n=== Metadata correlations (Cave only) ===\n")
print(round(cor_mat, 3))

png(file.path(out_dir, "Cave_Metadata_Corrplot.png"), width = 600, height = 600)
corrplot(cor_mat, method = "color", addCoef.col = "black", tl.col = "black")
dev.off()

# =============================================================================
# FIGURE 3 — FULL NMDS (publication-ready)
# =============================================================================

p3 <- plot_ordination(phy.filtered.DNA.rel, phy.ord, color = "Site_Type") +
  stat_ellipse(
    aes(fill = Site_Type),
    geom      = "polygon",
    alpha     = 0.12,
    level     = 0.95,
    linetype  = "dashed",
    linewidth = 0.5
  ) +
  geom_point(aes(color = Site_Type), size = 3.5, alpha = 0.9) +
  scale_color_manual(values = site_colors, name = "Site Type") +
  scale_fill_manual(values  = site_colors, name = "Site Type") +
  annotate("text",
           x = -Inf, y = Inf,
           label = overall_label,
           hjust = -0.04, vjust = 1.6,
           size = 3.2, fontface = "italic", color = "grey20") +
  annotate("text",
           x = -Inf, y = -Inf,
           label    = pairwise_label,
           hjust    = -0.05, vjust = -0.3,
           size     = 2.9, fontface = "italic", color = "grey20") +
  annotate("text",
           x = Inf, y = -Inf,
           label = paste0("Stress = ", round(phy.ord$stress, 3)),
           hjust = 1.08, vjust = -0.6,
           size = 3.2, fontface = "italic", color = "grey30") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10),
    panel.border     = element_rect(color = "grey40", linewidth = 0.6),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "Figure 3. Microbial Community Composition Across Site Types",
    x = "NMDS1", y = "NMDS2"
  )

ggsave(file.path(out_dir, "Figure3_NMDS_AllSites.png"), p3, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure3_NMDS_AllSites.pdf"), p3, width = 8, height = 6)
cat("Saved: Figure 3\n")

# =============================================================================
# FIGURE 4 — CAVE NMDS + ENVFIT (publication-ready)
# =============================================================================

arrow_scale <- 0.38
env_r2    <- round(env_fit$vectors$r[1], 3)
env_p     <- env_fit$vectors$pvals[1]
env_p_lab <- ifelse(env_p < 0.001, "< 0.001", paste0("= ", round(env_p, 3)))

env_vec           <- as.data.frame(scores(env_fit, display = "vectors"))
colnames(env_vec) <- c("NMDS1", "NMDS2")
env_vec$label     <- paste0("Distance from entrance\nr² = ", env_r2, ", p ", env_p_lab)

ord_scores <- as.data.frame(scores(phy.cave.ord, display = "sites"))
ord_scores$Distance_Datum_M <- cave.meta[rownames(ord_scores), "Distance_Datum_M"]

p4 <- ggplot(ord_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(fill = Distance_Datum_M),
    shape = 21, size = 4.5, color = "white", stroke = 0.4, alpha = 0.95
  ) +
  scale_fill_gradientn(
    colors = c("#D8D8E8", "#9090B8", "#5A4A8A", "#3B1F6B"),
    name   = "Distance from\nEntrance (m)"
  ) +
  geom_segment(
    data    = env_vec,
    aes(x = 0, y = 0,
        xend = NMDS1 * arrow_scale,
        yend = NMDS2 * arrow_scale),
    arrow       = arrow(length = unit(0.22, "cm"), type = "closed"),
    color       = "black", linewidth = 0.9, inherit.aes = FALSE
  ) +
  geom_text(
    data    = env_vec,
    aes(x     = NMDS1 * (arrow_scale + 0.10),
        y     = NMDS2 * (arrow_scale + 0.10),
        label = label),
    size = 3.3, fontface = "italic", color = "grey20", inherit.aes = FALSE
  ) +
  annotate("text",
           x = Inf, y = -Inf,
           label = paste0("Stress = ", round(phy.cave.ord$stress, 3)),
           hjust = 1.08, vjust = -0.6,
           size = 3.2, fontface = "italic", color = "grey30") +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title       = element_text(face = "bold", size = 11),
    axis.text        = element_text(size = 10),
    panel.border     = element_rect(color = "grey40", linewidth = 0.6),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "Figure 4. Cave Microbial Community Composition Along Depth Gradient",
    x = "NMDS1", y = "NMDS2"
  )

ggsave(file.path(out_dir, "Figure4_NMDS_Cave_Envfit.png"), p4, width = 8, height = 6, dpi = 500, bg = "white")
ggsave(file.path(out_dir, "Figure4_NMDS_Cave_Envfit.pdf"), p4, width = 8, height = 6)
cat("Saved: Figure 4\n")

# =============================================================================
# STEP 7 — EXPORT ALL STATISTICS
# =============================================================================

pw_df <- bind_rows(
  data.frame(Comparison = "Cave_vs_Surface",
             R2 = round(pw_cv$R2[1], 3), F = round(pw_cv$F[1], 3), P_value = pw_cv$`Pr(>F)`[1]),
  data.frame(Comparison = "Cave_vs_Water",
             R2 = round(pw_cw$R2[1], 3), F = round(pw_cw$F[1], 3), P_value = pw_cw$`Pr(>F)`[1]),
  data.frame(Comparison = "Surface_vs_Water",
             R2 = round(pw_sw$R2[1], 3), F = round(pw_sw$F[1], 3), P_value = pw_sw$`Pr(>F)`[1])
) %>%
  mutate(Significant = ifelse(P_value < 0.05, "YES ***",
                              ifelse(P_value < 0.10, "TRENDING ~", "NO")))

overall_df <- data.frame(
  Predictor   = "Site_Type",
  R2          = round(permanova_result$R2[1], 3),
  F_statistic = round(permanova_result$F[1], 3),
  P_value     = permanova_result$`Pr(>F)`[1],
  Significant = "YES ***"
)

envfit_df <- data.frame(
  Variable    = "Distance_Datum_M",
  r2          = round(env_fit$vectors$r[1], 3),
  MDS1        = round(scores(env_fit, display = "vectors")[1, 1], 3),
  MDS2        = round(scores(env_fit, display = "vectors")[1, 2], 3),
  P_value     = env_fit$vectors$pvals[1],
  Significant = ifelse(env_fit$vectors$pvals[1] < 0.05, "YES ***", "NO")
)

write_xlsx(
  list(
    PERMANOVA_Overall       = overall_df,
    PERMANOVA_Pairwise      = pw_df,
    PERMDISP_Overall        = betadisp_overall,
    PERMDISP_GroupDistances = betadisp_summary,
    Cave_PERMANOVA_Distance = cave_permanova_df,
    Envfit_Distance         = envfit_df,
    Metadata_Correlations   = as.data.frame(round(cor_mat, 3))
  ),
  file.path(out_dir, "NMDS_Statistical_Summary.xlsx")
)

cat("\n============================\n")
cat("ALL ANALYSES COMPLETE\n")
cat("Output folder:", out_dir, "\n")
cat("============================\n")
