# Moonshiner Cave 16S Analysis Scripts

R scripts for the 16S amplicon analysis of cave, surface, and water samples. They start from QIIME2 outputs and cover phyloseq setup, beta and alpha diversity, taxonomic composition, FEAST source tracking, and FAPROTAX functional prediction.

Sequence processing happens upstream of this folder (QIIME2 with Deblur, SILVA 138.1 taxonomy, following github.com/QLFaber/16Samplicons).

## Inputs

Put these in `Inputs/`:

- `feature-table.tsv` - QIIME2 ASV count table
- `taxonomy.tsv` - QIIME2 taxonomy with 7 rank columns
- `Karst_Input_Metadata.tsv` - sample metadata

## Scripts

**01_Phyloseq_Setup.R**
Builds the phyloseq object. Drops samples under 1,000 reads, removes chloroplast and mitochondria, and removes ASVs with 10 or fewer total reads. Saves `phyloseq_filtered_counts_DNA.rds` and `phyloseq_filtered_RA_DNA.rds`, which every other script reads.

**02_Beta_Diversity.R**
Bray-Curtis NMDS for all samples and for cave samples only. PERMANOVA, pairwise PERMANOVA, PERMDISP, and envfit for distance from the entrance. Makes Figures 3 and 4.

**03_Taxonomic_Composition.R**
Top 15 phyla and orders across site types. Finds ASVs detected only in cave samples and compares their taxonomy to the full cave community. Makes Figures 5, 6, and 13.

**04_Alpha_Diversity.R**
Shannon, Simpson, and richness. Kruskal-Wallis and pairwise Wilcoxon by site type, Spearman correlations with distance and with the FEAST Unknown fraction, and diversity of cave-exclusive ASVs. Makes Figures 7 and 10.

**05_FEAST_Source_Tracking.R**
FEAST with cave samples as sinks. Iteration 1 keeps Cave Entrance, Soil Transect, and Water as separate sources. Iteration 2 collapses the surface sources into one. Tests how source contributions change with distance and moisture. Makes Figures 8a, 8b, and 9.

**06_FEAST_Resident_Analysis.R**
Correlates each cave ASV with the FEAST Unknown fraction to find ASVs associated with it, reports their taxonomy, and compares Bray-Curtis dissimilarity between resident and non-resident ASV sets.

**07_Functional_Annotation.R**
FAPROTAX predictions through microeco for the full community and for cave-exclusive ASVs. Kruskal-Wallis by site type, Spearman correlations with distance, and a test of predatory function against the FEAST Unknown fraction.

## Run order

01, 02, 03, **05**, 04, 06, 07

04 reads the FEAST results from 05, so 05 has to run before it even though it is numbered after.

## Before running

- Set the working directory to this folder. 05 has a hardcoded `setwd()` on line 1 that needs to be changed to your path.
- Create `FEAST_output/` and `Outputs/FEAST/` by hand. 05 writes to both but does not create them. The other scripts make their own output folders.
- The metadata column names `Sampel.ID` and `Water_Percantage` are misspelled in the metadata file and the scripts use those spellings. Don't fix one without fixing the other.
- 05 reads the metadata as `Karst_Input_Metadata.TSV` and 01 reads it as `.tsv`. This works on Windows but will break on Mac or Linux.

## Outputs

- `.rds` phyloseq objects are saved in the main folder
- FEAST tables go to `FEAST_output/`
- Everything else goes to `Outputs/<analysis name>/` as figures (PNG and PDF) and Excel tables

In the figures, the FEAST "Unknown" fraction is labeled "Resident Cave".

## Packages

CRAN: `tidyverse`, `data.table`, `vegan`, `corrplot`, `writexl`, `readxl`, `microeco`, `scales`

Bioconductor: `phyloseq`

GitHub: `FEAST` (cozygene/FEAST), `pairwiseAdonis` (pmartinezarbizu/pairwiseAdonis)
