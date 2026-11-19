## run_all.R — reproduce every comparative analysis in the manuscript.
## Run from the repository root:  Rscript run_all.R
Sys.setenv(PREDS = "all")          # robustness analyses for all 34 predictors (Table S6)
scripts <- c("01_sequencing_covariates.R", "02_major_clades.R", "03_phylogenetic_signal.R", "04_life_history.R",
             "05_pgls_34_predictors.R", "06_robustness.R", "07_predictor_comparison.R", "09_phylogenetic_timing.R")
for (s in scripts) {
  cat("\n==========", s, "==========\n")
  status <- system2("Rscript", file.path("R", s))
  if (status != 0) stop("R/", s, " failed (exit status ", status, "); later scripts were not run")   # do not carry on after a failure
}
## MR-PMM: the MCMC chains used in the manuscript are in results/mrpmm_chains; this summarizes them (Tables S11, S12).
## To refit the chains yourself, see the header of R/08_mrpmm.R.
system2("Rscript", c("R/08_mrpmm.R", "summary"))
## Record the R and package versions used for this run
writeLines(capture.output({ suppressPackageStartupMessages({ library(ape); library(phylolm); library(geiger); library(MASS); library(MCMCglmm); library(coda); library(caper) }); print(sessionInfo()) }),
           "results/session_info.txt")
