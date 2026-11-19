## 08_mrpmm.R
## Bivariate multi-response phylogenetic mixed models (MR-PMM; Halliwell et al. 2025) with MCMCglmm
## (Fig 3C, Fig S5, Tables S11 and S12). Telomere length and one predictor are joint responses; the
## covariance between them is split into a phylogenetic part (random effect with the phylogenetic
## relatedness matrix) and a residual part. log10 spot length is a fixed effect on telomere length only.
##
## Why MR-PMM in addition to PGLS: PGLS treats the predictor as fixed and asks about one overall slope. MR-PMM treats
## both traits as outcomes and splits their association into
##   r_phy: correlation of the parts of each trait shared through common ancestry (co-evolution along the tree), and
##   r_res: correlation of the species-specific parts not explained by phylogeny (recent or within-lineage association).
## Result: bio5 r_phy = 0.56 (95% HPD 0.26-0.84) with r_res = 0.15 (-0.02-0.32), i.e. the temperature-telomere
## association is mainly phylogenetic: lineages that moved into warmer climates evolved longer telomeres.
##
## Each model is one MCMC chain (110,000 iterations, burn-in 10,000, thinning 100) and takes several
## minutes, so the script runs one job per call and a summary step at the end:
##   Rscript R/08_mrpmm.R main <predictor> <chain>            # e.g. main bio5 1   (Table S11: 13 predictors x 2 chains)
##   Rscript R/08_mrpmm.R jackknife <predictor> <drop_set>     # e.g. jackknife bio5 Rosaceae   (Table S12: single chain)
##   Rscript R/08_mrpmm.R summary
## Shorter chains for testing: NITT=11000 BURN=1000 THIN=10 Rscript R/08_mrpmm.R main bio5 1
## CAUTION: a test run writes to results/mrpmm_chains/ and overwrites the stored chain for that model; run tests in a copy.
## The 13 jackknife models that were below ESS 200 at the default length (log_alt without Asteraceae, Brassicaceae, Fabaceae,
## Lamiaceae, Malvaceae, Myrtaceae, Poaceae, Rosaceae, Rosaceae+Brassicaceae, Salicaceae or Solanaceae; bio5 without
## Bryophytes+lycophyte or Malvaceae) were run with NITT=1100000 BURN=100000 THIN=1000 (Myrtaceae with SEED_OFFSET=1).
## Stored results: results/mrpmm_chains/ holds 101 chain files (26 main + 75 jackknife), 1,000 posterior samples each;
## the summary step reads them, so Tables S11-S12 can be rebuilt without refitting. results/mrpmm_rerun_check/ holds
## independent refits of a few main models showing that refitting gives the same posteriors up to Monte Carlo error.
source("R/00_setup.R")
suppressPackageStartupMessages({ library(MCMCglmm); library(coda) })
dir.create("results/mrpmm_chains", showWarnings = FALSE)

## Table S11: 13 predictors, 2 chains each (26 chains)
main_predictors <- c("abs_lat", "log_alt", "genome_size_log10", "phenotype_num_chr", "bio1", "bio5", "bio7",
                     "bio8", "bio10", "bio15", "bio17", "gdd0", "gst")
## Table S12: 5 predictors x 15 drop sets, 1 chain each (75 chains)
jk_predictors   <- c("abs_lat", "bio5", "bio10", "log_alt", "phenotype_num_chr")
## Drop sets: each family with >= 5 species (13 families), Rosaceae + Brassicaceae together, and the non-seed plants.
## Each entry is a vector of family names to remove (bryophyte/lycophyte families are looked up from bryo_lyco).
fam5 <- names(which(table(d$Family) >= 5))
jk_sets <- c(setNames(as.list(fam5), fam5), list("Rosaceae+Brassicaceae" = c("Rosaceae", "Brassicaceae"),
             "Bryophytes+lycophyte" = unique(d$Family[d$bryo_lyco])))

## spot length standardized over all 264 species (as in the manuscript)
## (mean and SD from all_species, the full Table S15, not from the 257 species in d)
d$log_spot_z <- (d$log_spot - mean(all_species$log_spot)) / sd(all_species$log_spot)

## Fit one bivariate MR-PMM and save its posterior samples.
## x = predictor, chain = chain number (sets the seed), drop = family names to remove (NULL = none),
## tag = "main" or "jk", label = name of the drop set for the output file.
run_chain <- function(x, chain, drop = NULL, tag, label = NULL) {
  ## Chain length from environment variables (defaults: 110,000 iterations, first 10,000 discarded, every 100th kept
  ## = 1,000 posterior samples)
  NITT <- as.integer(Sys.getenv("NITT", 110000))
  BURN <- as.integer(Sys.getenv("BURN", 10000))
  THIN <- as.integer(Sys.getenv("THIN", 100))
  dd <- d[!d$Family %in% drop, ]                                                 # remove the dropped families (none if drop = NULL)
  dd <- dd[complete.cases(dd[, c("tel_log10", x, "log_spot_z")]), ]              # species with both traits and spot length
  Y  <- scale(as.matrix(dd[, c("tel_log10", x)])); colnames(Y) <- c("y1", "y2")    # both traits standardized
                                                                                 # y1 = telomere length, y2 = predictor; mean 0, SD 1,
                                                                                 # so the two traits' variances are on the same scale
  df <- data.frame(Y, animal = dd$tip, log_spot_z = dd$log_spot_z)              # "animal" = species name, linking rows to tree tips
  tr <- keep.tip(tree, dd$tip)
  Ainv <- inverseA(tr, nodes = "ALL", scale = TRUE)$Ainv                          # phylogenetic relatedness
  ## inverseA returns the inverse of the phylogenetic covariance matrix (covariance = shared branch length from the
  ## root, scaled so root-to-tip = 1; the tree is ultrametric, 490 Myr deep). With nodes = "ALL" the matrix covers the
  ## internal nodes as well as the species (with all 257 species: 257 tips + 255 internal nodes = 512 rows; the root is
  ## left out). Its species block reproduces ape::vcv(tree) / 490, so the model for the species is the same.
  ## Why "ALL": under Brownian motion a node's value depends only on its parent, so the inverse is non-zero only between
  ## a node and its parent/children (about 6 entries per species; 1,508 non-zero entries for the 253-species bio5 model)
  ## and is built directly from the tree. With nodes = "TIPS" the species are all linked through ancestors that are not
  ## in the model, so the inverse is dense (61,045 non-zero entries) and must be obtained by inverting a dense matrix.
  ## MCMCglmm uses sparse linear algebra, so its cost per iteration scales with the non-zero entries (Hadfield &
  ## Nakagawa 2010, J Evol Biol 23:494-508). Measured here (bio5 model, 11,000 iterations): 4.8 s with "ALL" vs 431 s
  ## with "TIPS", about 90x faster.
  ##
  ## Priors for the two 2 x 2 covariance matrices:
  ##   G (phylogenetic): parameter-expanded prior (V = identity, nu = 2, working parameters alpha with mean 0 and
  ##     variance 1000), weakly informative and well behaved when a variance is near zero.
  ##   R (residual): inverse-Wishart with V = identity, nu = 2.
  ## Note: with nu = 2 for 2 traits, the implied prior on each correlation is (1 - r^2)^(-1/2), mildly favoring values
  ## near +/-1 (nu = 3 would give a flat prior). This matters most when a trait has little phylogenetic variance
  ## (e.g. altitude, phylo_share_x = 0.21), where the data say little about r_phy.
  prior <- list(G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000)),  # parameter-expanded
                R = list(V = diag(2), nu = 2))                                                                # inverse-Wishart
  ## Seed = 1000 x chain + number of characters in the predictor name (+ optional offset), so chains 1 and 2 differ.
  ## (Predictors with names of equal length, e.g. bio1 and bio5, share a seed; their data differ, so this is harmless.)
  set.seed(1000 * chain + nchar(x) + as.integer(Sys.getenv("SEED_OFFSET", 0)))   # SEED_OFFSET: restart a chain that stops on an ill-conditioned G matrix
  ## The model:
  ##   cbind(y1, y2) ~ trait - 1          a separate intercept for each trait (no shared intercept)
  ##     + at.level(trait, 1):log_spot_z  spot length affects telomere length (trait 1) only, not the predictor
  ##   random = ~ us(trait):animal         unstructured 2 x 2 phylogenetic covariance (2 variances + 1 covariance),
  ##                                       with relatedness given by Ainv (ginverse)
  ##   rcov   = ~ us(trait):units          unstructured 2 x 2 residual covariance (species-specific parts)
  ##   pr = FALSE                          do not store the species-level random effects (saves memory)
  m <- MCMCglmm(cbind(y1, y2) ~ trait - 1 + at.level(trait, 1):log_spot_z,
                random = ~ us(trait):animal, rcov = ~ us(trait):units, family = rep("gaussian", 2),
                ginverse = list(animal = Ainv), prior = prior, data = df,
                nitt = NITT, burnin = BURN, thin = THIN, verbose = FALSE, pr = FALSE)
  ## m$VCV: one row per posterior sample, one column per (co)variance, e.g. "traity1:traity2.animal" is the
  ## phylogenetic covariance between telomere length and the predictor.
  V <- m$VCV
  ## Correlations computed for every posterior sample: covariance / sqrt(variance 1 x variance 2)
  r_phy <- V[, "traity1:traity2.animal"] / sqrt(V[, "traity1:traity1.animal"] * V[, "traity2:traity2.animal"])   # phylogenetic correlation
  r_res <- V[, "traity1:traity2.units"]  / sqrt(V[, "traity1:traity1.units"]  * V[, "traity2:traity2.units"])    # residual correlation
  ## Share of the predictor's variance that is phylogenetic (phylogenetic heritability, comparable to Pagel's lambda).
  ## Near 1: closely related species have similar values (e.g. bio5 ~0.55); near 0: little phylogenetic pattern.
  share <- V[, "traity2:traity2.animal"] / (V[, "traity2:traity2.animal"] + V[, "traity2:traity2.units"])
  out <- data.frame(tag = tag, var = x, dropped = if (!is.null(label)) label else if (is.null(drop)) "(none)" else paste(drop, collapse = "+"), chain = chain,
                    n = nrow(df), r_phy = as.numeric(r_phy), r_res = as.numeric(r_res), phylo_share_x = as.numeric(share))
  ## One CSV per chain, named tag__predictor__dropset__cN.csv (spaces and "+" replaced by "-"), 1,000 rows
  f <- sprintf("results/mrpmm_chains/%s__%s__%s__c%d.csv", tag, x, gsub("[+ ]", "-", out$dropped[1]), chain)
  write.csv(out, f, row.names = FALSE)
  cat(tag, x, out$dropped[1], "chain", chain, "n =", nrow(df), "median r_phy =", round(median(r_phy), 2), "\n")
}

## Combine all chain files into Tables S11 (main) and S12 (jackknife).
summarize <- function() {
  x <- do.call(rbind, lapply(list.files("results/mrpmm_chains", full.names = TRUE), read.csv))   # all posterior samples
  ## One summary row per model (tag + predictor + drop set); chains of the same model are pooled.
  res <- do.call(rbind, lapply(split(x, paste(x$tag, x$var, x$dropped)), function(s) {
    ch <- split(s, s$chain)                                  # list of chains for this model (2 for main, 1 for jackknife)
    ## Posterior median and 95% highest posterior density (HPD) interval: the narrowest interval containing 95% of
    ## the samples. An interval that excludes 0 is the Bayesian analogue of a significant correlation.
    q <- function(v) { h <- HPDinterval(mcmc(v)); c(median(v), h[1], h[2]) }
    ## Effective sample size, summed over chains: the number of independent samples the autocorrelated chain is worth
    ## (maximum 1,000 per chain here). >= 200 was the target.
    ess  <- function(col) sum(sapply(ch, function(c) effectiveSize(mcmc(c[[col]]))))
    ## Gelman-Rubin R-hat (potential scale reduction factor, point estimate): compares within- and between-chain
    ## variance; ~1.00 means the chains agree (< 1.1, ideally < 1.01). Needs >= 2 chains, so NA for jackknife models.
    rhat <- function(col) if (length(ch) > 1) gelman.diag(mcmc.list(lapply(ch, function(c) mcmc(c[[col]]))), autoburnin = FALSE)$psrf[1] else NA
    a <- q(s$r_phy)
    b <- q(s$r_res)
    data.frame(tag = s$tag[1], var = s$var[1], dropped = s$dropped[1], n = s$n[1], chains = length(ch),
               r_phy = a[1], r_phy_lo = a[2], r_phy_hi = a[3], r_res = b[1], r_res_lo = b[2], r_res_hi = b[3],
               phylo_share_x = median(s$phylo_share_x), ESS_r_phy = ess("r_phy"), ESS_r_res = ess("r_res"),
               Rhat_r_phy = rhat("r_phy"), Rhat_r_res = rhat("r_res"))
  }))
  ## Expected (main, r_phy with 95% HPD): bio5 0.56 (0.26-0.84); bio10 0.52 (0.24-0.81); bio1 0.41; gdd0 0.37;
  ## latitude -0.29 (-0.56 to -0.01); altitude -0.71 (-1.00 to -0.34; phylo_share_x 0.21, ESS 422, R-hat 1.014, so the
  ## least precise); chromosome number 0.48 (0.13-0.81). r_res intervals include 0 except bio17 (-0.21) and bio7 (0.20).
  ## Jackknife: bio5 and bio10 r_phy stay positive with HPD excluding 0 for every drop set except Rosaceae + Brassicaceae
  ## (bio5 0.42, -0.04 to 0.74); altitude and chromosome number exclude 0 for every drop set; latitude includes 0 for 7 of 15.
  print(res[res$tag == "main", ], digits = 2, row.names = FALSE)
  save_csv(res[res$tag == "main", ], "08_mrpmm_main.csv")
  save_csv(res[res$tag == "jk", ], "08_mrpmm_jackknife.csv")
}

## Command-line dispatch: the first argument chooses the job
a <- commandArgs(TRUE)
if (length(a) == 0) stop("usage: main <predictor> <chain> | jackknife <predictor> <drop_set> | summary")
if (a[1] == "main") run_chain(a[2], as.integer(a[3]), NULL, "main")
if (a[1] == "jackknife") run_chain(a[2], 1L, jk_sets[[a[3]]], "jk", label = a[3])
if (a[1] == "summary") summarize()
## To run everything (e.g. on a cluster, one job per line):
##   for (x in main_predictors) for (ch in 1:2) system(paste("Rscript R/08_mrpmm.R main", x, ch))
##   for (x in jk_predictors) for (s in names(jk_sets)) system(paste("Rscript R/08_mrpmm.R jackknife", x, shQuote(s)))
