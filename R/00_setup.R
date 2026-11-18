## 00_setup.R
## Packages, data loading and helper functions shared by every analysis script.
## Run all scripts from the repository root, e.g. Rscript R/05_pgls_34_predictors.R

suppressPackageStartupMessages({
  library(ape)        # tree handling, independent contrasts
  library(phylolm)    # PGLS with Pagel's lambda
  library(geiger)     # fitContinuous, disparity through time
  library(MASS)       # robust regression (node-height tests)
})
dir.create("results", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Data: Table S15 of the manuscript exported to CSV, and the time-calibrated tree
## ---------------------------------------------------------------------------
raw <- read.csv("data/TableS15_species_data.csv", check.names = FALSE, na.strings = c("", "NA"),
                fileEncoding = "UTF-8")

## Short names used in the code (left) for the Table S15 column headers (right)
rename <- c(
  Species = "Species", APG.IV = "Major clade (APG IV)", Order = "Order", Family = "Family",
  in_comparative = "In comparative analyses", TL_bp = "Telomere length (bp)", tel_log10 = "log10 telomere length",
  Instrument = "Sequencing instrument", platform = "Sequencing platform", spot_bp = "SRA spot length (bp)",
  N50_kb = "Contig N50 (kb)", genome_size_log10 = "log10 genome size", phenotype_num_chr = "Chromosome number (2n)",
  ploidy = "Ploidy (diploid/polyploid)", lifespan = "Lifespan", growthform = "Growth form",
  phenotype_height = "Plant height (m)", phenotype_root_structure = "Root structure",
  phenotype_shoot_structure = "Shoot structure", phenotype_genome_diversity = "Genome diversity",
  rangesize_log10 = "Range size (log10 km²)", native_range = "Native range size (km²)",
  latitude = "Latitude", longitude = "Longitude", abs_lat = "Absolute latitude", median_altitude = "Median altitude (m)")
bioclim <- c(paste0("bio", 1:19), "gst", "gsl", "gdd0", "ngd10", "scd")
all_species <- raw[, c(rename, bioclim)]
names(all_species) <- c(names(rename), bioclim)

## Derived variables used in the analyses
all_species <- within(all_species, {
  tip       <- gsub(" ", "_", trimws(Species))           # tree tip labels
  genus     <- sub(" .*", "", trimws(Species))
  log_spot  <- log10(spot_bp)                             # technical covariate in every comparative model
  log_alt   <- log10(pmax(median_altitude, 0) + 1)
  log_native_range <- log10(native_range)
  ONT       <- as.numeric(platform == "ONT")              # NA where platform is unknown
})

## Clade groupings derived from APG IV clade and order
rosid_orders   <- c("Brassicales", "Cucurbitales", "Fabales", "Fagales", "Malpighiales", "Malvales",
                    "Myrtales", "Rosales", "Sapindales", "Vitales")
asterid_orders <- c("Apiales", "Asterales", "Boraginales", "Cornales", "Dipsacales", "Ericales",
                    "Garryales", "Lamiales", "Solanales")
fabid_orders   <- c("Rosales", "Fabales", "Cucurbitales", "Fagales", "Malpighiales")
malvid_orders  <- c("Brassicales", "Malvales", "Sapindales", "Myrtales")
all_species$clade <- with(all_species, ifelse(Order %in% rosid_orders, "Rosids",
  ifelse(Order %in% asterid_orders, "Asterids", ifelse(APG.IV == "Monocots", "Monocots",
  ifelse(Order %in% c("Caryophyllales", "Santalales", "Proteales"), "Other eudicots", "Other lineages")))))
all_species$bryo_lyco <- all_species$APG.IV %in% c("Anthocerotophyta", "Bryophytes", "Marchantiophyta", "Lycophytes")

## Tree and land-plant data set used in all comparative analyses (n = 257)
tree <- read.tree("data/tree_plantTL.tre"); tree$node.label <- NULL
d <- all_species[all_species$in_comparative == "yes" & all_species$tip %in% tree$tip.label, ]
rownames(d) <- d$tip
stopifnot(nrow(d) == 257, all(tree$tip.label %in% d$tip))

## The 34 predictors tested against telomere length (Table 1 and Table S5)
predictors <- c(
  abs_lat = "Latitude (absolute)", log_alt = "Altitude (log10)", genome_size_log10 = "Genome size (log10 Mb)",
  phenotype_num_chr = "Chromosome number (2n)", bio1 = "Mean annual temperature (bio1)",
  bio2 = "Mean diurnal temperature range (bio2)", bio3 = "Isothermality (bio3)",
  bio4 = "Temperature seasonality (bio4)", bio5 = "Max temperature of warmest month (bio5)",
  bio6 = "Min temperature of coldest month (bio6)", bio7 = "Annual temperature range (bio7)",
  bio8 = "Mean temperature of wettest quarter (bio8)", bio9 = "Mean temperature of driest quarter (bio9)",
  bio10 = "Mean temperature of warmest quarter (bio10)", bio11 = "Mean temperature of coldest quarter (bio11)",
  bio12 = "Annual precipitation (bio12)", bio13 = "Precipitation of wettest month (bio13)",
  bio14 = "Precipitation of driest month (bio14)", bio15 = "Precipitation seasonality (bio15)",
  bio16 = "Precipitation of wettest quarter (bio16)", bio17 = "Precipitation of driest quarter (bio17)",
  bio18 = "Precipitation of warmest quarter (bio18)", bio19 = "Precipitation of coldest quarter (bio19)",
  gst = "Growing season temperature (gst)", gsl = "Growing season length (gsl)",
  gdd0 = "Growing degree days > 0 °C (gdd0)", ngd10 = "Growing degree days > 10 °C (ngd10)",
  scd = "Snow cover days (scd)", phenotype_shoot_structure = "Shoot structure",
  phenotype_root_structure = "Root structure", phenotype_height = "Plant height",
  log_native_range = "Native range size (log10)", phenotype_genome_diversity = "Genome diversity",
  rangesize_log10 = "Range size (log10)")

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
## Fit a PGLS (Pagel's lambda) to the complete cases of the variables in formula f.
## Variables named in `standardize` are scaled to unit variance first, so slopes are per SD.
## Example: pgls(tel_log10 ~ bio5 + log_spot, standardize = "bio5")
## Returns a list: model (phylolm fit; lambda in $optpar), data (species used; nrow = n),
## tree (pruned tree), coef (estimate, SE, t, p for each term).
pgls <- function(f, data = d, standardize = NULL) {
  v  <- all.vars(f)                                              # variable names in the formula
  dd <- data[complete.cases(data[, v]), v, drop = FALSE]         # species with no missing value for these variables,
                                                                 #   so n can differ between models
  for (x in standardize) dd[[x]] <- as.numeric(scale(dd[[x]]))   # mean 0, SD 1 among the species in this model
  tr <- keep.tip(tree, rownames(dd))                             # prune the tree to those species (rownames are tip labels)
  dd <- dd[tr$tip.label, , drop = FALSE]                         # put rows in tree tip order
  m  <- suppressWarnings(phylolm(f, dd, tr, model = "lambda"))   # GLS with lambda estimated by maximum likelihood, separately
                                                                 #   for each model; warnings when lambda is at 0 or 1 are hidden
  list(model = m, data = dd, tree = tr, coef = summary(m)$coefficients)
}
## Small-sample corrected AIC (AICc) for a phylolm fit; lower is better. Used to rank predictors (Table S8).
## m = phylolm fit (r$model from pgls()), n = number of species in the model (nrow(r$data)).
## Compare AICc only between models fitted to the same species and the same response.
## Example (230 species): bio5 + spot length, k = 5, AICc = -147.0; spot length only, k = 4, AICc = -135.5
aicc <- function(m, n) {
  k <- length(coef(m)) + 2                                   # parameters: regression coefficients (intercept + slopes)
                                                             #   + residual variance (sigma2) + Pagel's lambda
  -2 * m$logLik + 2 * k +                                    # AIC: fit (-2 log-likelihood) + 2 per parameter
    2 * k * (k + 1) / (n - k - 1)                            # small-sample correction; near 0 when n is large
}
## Spot-length adjusted telomere length: subtract the spot-length term of the PGLS of telomere length
## on bio5 + log10 spot length, centred on the mean spot length (used in figures and single-trait analyses)
## The value is the log10 telomere length each species would have at the average spot length.
## Model: log10 TL = a + b1*bio5 + b2*log_spot + e. Replacing a species' log_spot with the mean and
## subtracting leaves a, b1*bio5 and e unchanged, so
##   adjusted = log10 TL - b2 * (log_spot - mean log_spot)
## Only the spot-length effect is removed. This is not a residual: the bio5 effect is kept.
## bio5 is in the model so that b2 measures spot length holding temperature constant.
## Example: Arabidopsis thaliana, 6,300 bp reads (shorter than the mean, about 12,300 bp):
##   2,445 bp observed -> 2,760 bp adjusted (b2 = 0.181, 253 species).
## Returns a vector named by tip label, for the species with bio5 and spot-length data.
adjust_for_spot <- function(data = d) {
  r <- pgls(tel_log10 ~ bio5 + log_spot, data)                    # fit the main model; b2 = spot-length slope
  dd <- r$data                                                    # species used in the fit
  setNames(data[rownames(dd), "tel_log10"] -                      # observed log10 telomere length
             coef(r$model)["log_spot"] * (dd$log_spot - mean(dd$log_spot)),   # minus b2 x deviation from mean spot length
           rownames(dd))
}
save_csv <- function(x, file) write.csv(x, file.path("results", file), row.names = FALSE)
