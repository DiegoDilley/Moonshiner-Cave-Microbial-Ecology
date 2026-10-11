# Moonshiner βNTI / RCbray pipeline (scripts 01–06)

This is the code I used to build a phylogenetic tree from 16S ASVs and then run the Stegen et al. (2012, 2013) null model framework on it. Every pair of samples gets a βNTI value and an RCbray value, and those two numbers sort the pair into one assembly process: selection, dispersal, or drift.

The scripts run in order. Some run locally in R and some run on OSC (Pitzer) through SLURM. The scripts call βNTI `bMTI` and βMNTD `bMNTD`. They mean the same thing.

| Step | File | Runs on | What it does |
|---|---|---|---|
| 01 | `01_DataPrep.R` | Local R | Pulls out the sequences for ASVs that are in the filtered phyloseq object |
| 02 | `02_CondaSetUp.txt` | OSC, once | Builds the `phylo_tree` conda environment |
| 03 | `03_phylo_tree.txt` | OSC, SLURM | MUSCLE alignment and ModelTest-NG model selection |
| 04 | `04_iq_tree.txt` | OSC, SLURM | Maximum likelihood tree in IQ-TREE |
| 05 | `05_iqTree.R` | Local R | Roots the tree and adds it to the phyloseq object |
| 06 | `06_bNTI_analysis_parallel_v3.R` | OSC, SLURM, parallel | βMNTD, βNTI, RCbray, and process classification |

```
phyloseq_filtered_counts_DNA.rds + sequences.fasta
  └─ 01 ──> asvs_for_tree.fasta
             └─ 03 ──> asvs_aligned.fasta + ModelTest-NG results
                        └─ 04 ──> final_tree.treefile
                                   └─ 05 ──> phyloseq_filtered_counts_DNA_tree.rds
                                              └─ 06 ──> bMTI_RCbray_results.csv + matrices
```

## Inputs

You need two files in `Inputs/` before starting:

- `phyloseq_filtered_counts_DNA.rds`: the filtered phyloseq object (ASV counts, taxonomy, sample data) from the upstream 16S processing (Deblur/QIIME2).
- `sequences.fasta`: the representative sequence for every ASV. The sequence names have to match `taxa_names()` in the phyloseq object.

The data files are not included in this repo.

---

## 01_DataPrep.R (local R)

Packages: `phyloseq`, `Biostrings`

`sequences.fasta` has every ASV from the original run, but the phyloseq object has already been filtered. This script keeps only the sequences for ASVs that are still in the phyloseq object, so the tree and the count table have the same ASVs.

It prints two numbers: ASVs in the phyloseq object and ASVs matched in the fasta. They should be the same. If they aren't, the ASV names in the fasta don't match the phyloseq names, and step 06 will fail at the tree/OTU check.

Outputs:

- `Inputs/asvs_for_tree.fasta`: the sequences used to build the tree.
- `Inputs/asvs_for_tree_TEST.fasta`: the first 250 sequences. Use this for a test SLURM run so you find problems before spending OSC hours on the full set.

Upload `asvs_for_tree.fasta` to your `PhyloAnalysis` folder on OSC.

## 02_CondaSetUp.txt (OSC, run once)

Run these commands in an OSC terminal. They load miniconda and create an environment called `phylo_tree` with MUSCLE, ModelTest-NG, and PhyML from bioconda, then add IQ-TREE. You only do this once. Steps 03 and 04 activate the environment themselves.

## 03_phylo_tree.txt (OSC, SLURM)

This file has the job script plus notes on making and submitting it. Copy everything from `#!/bin/bash` to `echo "All steps finished!"` into `phylo_tree.sh` (`nano phylo_tree.sh`), then submit with `sbatch phylo_tree.sh`.

Before running, change `--account`, `--mail-user`, and `WORKDIR` to your own.

Resources: 1 node, 16 tasks, 64 GB, 24 hours.

What the job does:

1. **MUSCLE alignment.** Aligns `asvs_for_tree.fasta` with MUSCLE v5 in `-super5` mode (made for large sets of sequences). Writes `asvs_aligned.fasta`.
2. **ModelTest-NG.** Tests nucleotide substitution models on the alignment and writes results to `Tree_Output/modeltest_results.*`. GTR+I+G4 came out best by both BIC and AIC (weight = 1.0).
3. **PhyML.** This step fails. By the time it runs, the alignment and the ModelTest-NG results are already written and fine, so don't rerun the whole job. Go to step 04, which builds the tree in IQ-TREE instead.

Check on it while it runs:

```bash
squeue -u $USER            # job status
tail -f phylo_tree_*.log   # live log
```

## 04_iq_tree.txt (OSC, SLURM)

This replaces the PhyML step. Same setup as 03: copy the job part into `run_iqtree.sh`, edit `--account`, `--mail-user`, and `WORKDIR`, then `sbatch run_iqtree.sh`.

It builds a maximum likelihood tree from `asvs_aligned.fasta` using the model ModelTest-NG picked:

- `-m GTR+I+G4`: substitution model from step 03
- `-B 1000`: 1000 ultrafast bootstrap replicates
- `-T 16`: 16 threads

Outputs (in `Tree_Output/`):

- `final_tree.treefile`: the tree used in the rest of the pipeline
- `final_tree.iqtree` and `final_tree.log`: IQ-TREE's report and log
- `best_model_summary.txt`: records the model, bootstrap settings, input file, and date

The top of the file has a Biopython snippet that converts the alignment from FASTA to PHYLIP. You can skip it. IQ-TREE reads FASTA directly, and the job uses `asvs_aligned.fasta` anyway.

Download `final_tree.treefile` into your local `Inputs/` folder.

## 05_iqTree.R (local R)

Packages: `phyloseq`, `ape`, `phangorn`

1. Reads `Inputs/final_tree.treefile`.
2. Midpoint roots it with `phangorn::midpoint()`. IQ-TREE writes an unrooted tree, and phyloseq tools like UniFrac expect a rooted one.
3. Merges the rooted tree into the filtered phyloseq object.
4. Saves `phyloseq_filtered_counts_DNA_tree.rds`.

Check that the number of tips it prints matches the number of ASVs from step 01.

Step 06 reads this file from `Inputs/` in the BNTI folder on OSC, so upload it there.

## 06_bNTI_analysis_parallel_v3.R (OSC, SLURM, parallel)

Packages (installed in your OSC R library): `phyloseq`, `picante`, `vegan`, `dplyr`, `foreach`, `doParallel`

This is the main analysis. It runs 999 randomizations for βMNTD and another 999 for RCbray, so it runs on OSC with multiple cores.

Before running:

- The script sets `.libPaths()` and `setwd()` to my OSC folders. Change both to yours.
- The number of cores comes from `SLURM_CPUS_PER_TASK`. Submit with `#SBATCH --cpus-per-task=N`. If you only set `--ntasks` (like in 03 and 04), that variable is empty and the script runs on 1 core.

### What it does

**Step 1. Load data.** Reads `Inputs/phyloseq_filtered_counts_DNA_tree.rds`, turns the OTU table into samples × ASVs, and stops if any ASV is missing from the tree.

**Step 2. Observed βMNTD.** For each pair of samples, `picante::comdistnt(abundance.weighted = TRUE)` finds the phylogenetic distance from each ASV in one community to its closest relative in the other community, weighted by abundance. Low values mean the two communities are made of close relatives.

**Step 3. βMNTD null model.** Shuffles the ASV labels across the tips of the tree 999 times and recalculates βMNTD each time. Which ASVs are in each sample and how abundant they are stays the same. Only their positions on the tree get randomized. This gives the βMNTD you'd expect if phylogeny didn't matter. Saved as `Outputs/null_bMNTD_array.rds`.

**Step 4. βNTI.** How far the observed βMNTD is from the null, in standard deviations:

```
βNTI = (observed βMNTD − mean of null βMNTD) / SD of null βMNTD
```

**Step 5. RCbray null model.** Builds 999 null communities for every sample, following the abundance-weighted Raup-Crick method in the Stegen et al. 2013 supplement:

- A. Draw the same number of ASVs as the sample's observed richness, without replacement. ASVs are weighted by how common they are in the whole dataset (occurrence frequency × mean relative abundance).
- B. Spread the sample's real total read count across those ASVs with a multinomial draw, using the same weights rescaled to the chosen ASVs. This is the change in v3. v2 gave every chosen ASV an equal chance at each read, which doesn't match Stegen 2013.

Then it calculates Bray-Curtis on the null communities and compares the observed Bray-Curtis for each pair to its 999 null values. RCbray is that comparison scaled from −1 to +1. Saved as `Outputs/null_BC_array.rds`.

**Step 6. Classify each pair.**

| Condition | Process |
|---|---|
| βNTI > +2 | Variable selection |
| βNTI < −2 | Homogeneous selection |
| \|βNTI\| ≤ 2 and RCbray > +0.95 | Dispersal limitation + drift |
| \|βNTI\| ≤ 2 and RCbray < −0.95 | Homogenizing dispersal |
| \|βNTI\| ≤ 2 and \|RCbray\| ≤ 0.95 | Drift (undominated) |

βNTI is checked first. RCbray only matters for pairs where selection isn't significant.

**Step 7. Save outputs** (in `Outputs/`):

- `bMTI_RCbray_results.csv`: one row per sample pair with `Sample1`, `Sample2`, `SiteType1`, `SiteType2`, `Comparison`, `bMNTD_obs`, `bMTI`, `RCbray`, `Process`
- `bMNTD_observed_matrix.csv`, `bMTI_matrix.csv`, `RCbray_matrix.csv`: full sample × sample matrices
- `Process_Summary.csv`: count and percent of each process for each comparison type
- `null_bMNTD_array.rds`, `null_BC_array.rds`: the raw null distributions. These files are big, but keep them. They let you recalculate βNTI or RCbray without rerunning the randomizations.

**Site type note:** the `Comparison` column is built from `Site_Type` in the phyloseq sample data, and that column labels C001, C002, and C003 wrong (they're surface soils). The figure scripts that come after step 06 rebuild the comparison groups from `metadata.xlsx`, so use those groups, not the `Comparison` column in this file. The bMTI, RCbray, and Process values aren't affected.

---

## Software

| Tool | Version |
|---|---|
| R (local) | 4.5.1 |
| phyloseq | 1.52.0 |
| Biostrings | 2.76.0 |
| ape | 5.8-1 |
| phangorn | 2.12.1 |
| picante | 1.8.2 |
| vegan | 2.7-3 |
| dplyr | 1.2.0 |
| foreach | 1.5.2 |
| doParallel | 1.0.17 |
| miniconda (OSC module) | miniconda3/24.1.2-py310 |
| MUSCLE | v5 |
| ModelTest-NG | from bioconda |
| IQ-TREE | 3.1.1 |

R package versions come from my local library. The versions in the OSC library may be different.

## References

Stegen JC, Lin X, Konopka AE, Fredrickson JK. 2012. Stochastic and deterministic assembly processes in subsurface microbial communities. *ISME J* 6:1653–1664.

Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ, Rockhold ML, Konopka A. 2013. Quantifying community assembly processes and identifying features that impose them. *ISME J* 7:2069–2079.

Chase JM, Kraft NJB, Smith KG, Vellend M, Inouye BD. 2011. Using null models to disentangle variation in community dissimilarity from variation in α-diversity. *Ecosphere* 2:art24.
