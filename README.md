# Microbial Ecology of Moonshiner Cave

Undergraduate thesis research characterizing the structure, function, and assembly of microbial communities in a karst cave system, with a manuscript in preparation for the *Journal of Geomicrobiology*.



## Authorship:

Diego Dilley, Quincy Faber, Alexander B. Michaud 



Affiliations: Byrd Polar and Climate Research Center, Ohio State University, Columbus, OH, USA 

School of Earth Sciences, Ohio State University, Columbus, OH, USA 



## Overview

Caves are simultaneously connected to and isolated from the surface — water percolates in, debris falls through the entrance, but the interior stays dark, low-carbon, and thermally stable. This project asks whether that connection is enough to explain what lives underground, or whether the cave is actively selecting its own community.

I sampled 31 sites across Moonshiner Cave (Rockcastle County, Kentucky) — cave sediment along a depth transect, a paired surface soil transect directly above it, and cave groundwater sources — and used 16S rRNA sequencing, microbial source tracking (FEAST), functional annotation (FAPROTAX), enzyme assays, and phylogenetic null modeling (βNTI/Stegen framework) to figure out where the cave community comes from and what's shaping it.



## Key findings

* **Cave, surface, and groundwater are three distinct communities**, not a cave community that's just a diluted version of the surface (environment type explains \~31% of compositional variance, PERMANOVA p = 0.001).
* **Most of the cave community is unexplained by any sampled source.** FEAST source tracking attributed a mean of 68% of cave sediment composition to no detectable surface or groundwater origin — i.e., a genuinely resident cave community.
* **Distance from the entrance matters, but less than it looks like at first.** An NMDS/envfit analysis suggested entrance distance explained \~78% of intra-cave variation; a more conservative PERMANOVA on the same data put the real effect at \~12%. Distance is a useful proxy for the environmental gradient (light, temperature, carbon) radiating from the entrance, but it's not the whole story.
* **Karst hydrology breaks the distance-decay pattern.** Sites next to sinking streams and drip pools showed elevated surface- or water-derived signal regardless of how deep into the cave they were — physical distance from the entrance didn't predict phylogenetic community similarity (bMTI vs. distance: no relationship).
* **Function shifts sharply the moment you enter the cave.** Surface communities are dominated by heterotrophic carbon degradation; cave sediment shifts toward aerobic ammonia oxidation and nitrification — consistent with a chemolithoautotrophic cave model. Enzyme assays back this up directly: β-galactosidase activity (an indicator of organic matter breakdown) drops \~10-fold from surface to cave sediment.
* **Assembly is driven mostly by selection, not chance.** Using the Stegen et al. null-modeling framework, both cave and surface communities showed >94% homogeneous selection within their own habitat type, meaning the cave isn't just accumulating whoever wanders in — it's filtering for taxa that can survive there.



## Tools \& methods

* **Bioinformatics:** cutadapt, fastp, Deblur (ASV denoising), SILVA v138.1 taxonomy, phyloseq
* **Phylogenetics:** MUSCLE (Super5), ModelTest-NG, IQ-TREE (GTR+I+G4, 1000 bootstrap)
* **Statistics (R):** vegan (NMDS, PERMANOVA, betadisper), pairwiseAdonis, FEAST (source tracking), microeco/FAPROTAX (functional annotation), picante (βMNTD/βNTI null modeling)
* **Wet lab:** DNA extraction, Qubit quantification, enzyme activity assays (β-galactosidase)
* **Sequencing:** Illumina MiSeq, 16S rRNA V4–V5 region

## Repository structure

```
├── scripts/          # R pipeline: QC, ASV calling, phylogenetics, statistics, FEAST, βNTI
├── figures/          # Site map, ordinations, diversity plots, source-tracking and assembly-process figures
├── data/             # \[see note below]

├── Presentations/    # Undergrad thesis presentation material
└── README.md
```

## Data availability

Raw sequencing reads have been deposited in the NCBI GenBank Sequence Read Archive under the Bioproject PRJNA1517589. Any additional data is available upon request

## Status

Thesis defended, April 2026. Manuscript in preparation.
