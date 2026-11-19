# Plant telomere length and warm-season temperature

Data and R code for the phylogenetic comparative analyses in:

> **Nguyen L, Choi JY.** Plant telomere length is associated with warm-season temperature rather than life history. *Submitted.* Preprint: [DOI to be added]

## Quick start

Requires R (tested with 4.3.3) and six packages:

```r
install.packages(c("ape", "phylolm", "geiger", "MASS", "MCMCglmm", "coda"))
```

Then, from the repository root:

```bash
Rscript run_all.R
```

This runs every analysis in order (about 30 minutes on a laptop), writes the results to `results/` and records the R and package versions in `results/session_info.txt`. Each step can also be run alone, for example:

```bash
Rscript R/05_pgls_34_predictors.R
```

The Bayesian models (step 08) take hours to refit, so `run_all.R` summarizes the MCMC chains used in the paper, which are included in `results/mrpmm_chains/`. The header of `R/08_mrpmm.R` explains how to refit them.

## Repository contents

```
data/
  TableS15_species_data.csv    one row per species: telomere length, sequencing, traits, climate (Table S15)
  tree_plantTL.tre             time-calibrated phylogeny of the 257 land plants (Newick, millions of years)
R/
  00_setup.R                   packages, data loading, variable definitions, shared functions
  01_sequencing_covariates.R   read length, sequencing platform and assembly quality
  02_major_clades.R            telomere length across major plant lineages
  03_phylogenetic_signal.R     models of trait evolution; Pagel's lambda with confidence intervals
  04_life_history.R            lifespan, growth form and ploidy
  05_pgls_34_predictors.R      telomere length against 34 geographic, climatic and trait predictors
  06_robustness.R              clade removal, one species per genus, sequencing platform
  07_predictor_comparison.R    model selection, multi-predictor models, slopes within groups
  08_mrpmm.R                   Bayesian multi-response phylogenetic mixed models
  09_phylogenetic_timing.R     disparity through time, node-height tests, contrasts by node age
run_all.R                      runs steps 01-09
docs/analysis_walkthrough.md   what each analysis does, with code excerpts and the numbers it produces
```

## Data

`data/TableS15_species_data.csv` is Table S15 of the manuscript. Main columns:

| Column | Description |
|---|---|
| Telomere length (bp) | Median telomere length across telomeric reads (Topsicle) |
| SRA spot length (bp), Sequencing platform | Median read length of the libraries; PacBio or Oxford Nanopore. Included as covariates. |
| Latitude, Longitude, Median altitude (m) | Representative location of each species |
| bio1–bio19, gst, gsl, gdd0, ngd10, scd | CHELSA v2.1 bioclimatic variables at that location |
| Lifespan, Growth form, Ploidy, Plant height, ... | Life-history and genome traits |
| In comparative analyses | `yes` for the 257 land plants in the phylogeny; `no` for the 7 green algae |


## Contact

Jae Young Choi, Department of Ecology and Evolutionary Biology, University of Kansas. Lab website: [jychoilab.github.io](https://jychoilab.github.io). 
