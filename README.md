# Microbial Ecology of Moonshiner Cave


Undergraduate thesis research characterizing the structure, function, and assembly of microbial communities in a moonshiner cave.


## Authorship:

Diego Dilley, Quincy Faber, Alexander B. Michaud 



Affiliations: Byrd Polar and Climate Research Center, Ohio State University, Columbus, OH, USA 

School of Earth Sciences, Ohio State University, Columbus, OH, USA 



## Overview

Caves are simultaneously connected to and isolated from the surface — water percolates in, debris falls through the entrance, but the interior stays dark, low-carbon, and thermally stable. This project asks whether that connection is enough to explain what lives underground, or whether the cave is actively selecting its own community.

I sampled 31 sites across Moonshiner Cave (Rockcastle County, Kentucky) — cave sediment along a depth transect, a paired surface soil transect directly above it, and cave groundwater sources — and used 16S rRNA sequencing, microbial source tracking (FEAST), functional annotation (FAPROTAX), enzyme assays, and phylogenetic null modeling (βNTI/Stegen framework) to figure out where the cave community comes from and what's shaping it.

Code used for this project has been uploaded to this GitHub repo


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

