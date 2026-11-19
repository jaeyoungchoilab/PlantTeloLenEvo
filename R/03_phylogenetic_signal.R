## 03_phylogenetic_signal.R
## Evolutionary model comparison for telomere length (Table S3) and Pagel's lambda with
## profile-likelihood 95% CIs for telomere length subsets and other traits (Fig 2, Table S4).
source("R/00_setup.R")

## ---- 1. Model comparison with fitContinuous (Table S3) ----
## Eight models of how telomere length evolved along the tree, fitted by maximum likelihood and ranked by AICc.
##   BM      Brownian motion: random drift at a constant rate
##   lambda  Pagel's lambda: BM with shared history down-weighted (lambda < 1 = extra variation at the tips)
##   kappa   Pagel's kappa: branch lengths raised to the power kappa (change tracks number of splits more than time)
##   OU      Ornstein-Uhlenbeck: drift pulled toward an optimum (alpha = strength of pull)
##   EB      early burst: rate of change slows through time (a < 0)
##   trend   rate trend: rate of change speeds up or slows down steadily (geiger's "rate_trend")
##   mean    mean trend: trait drifts in one direction over time (geiger's "mean_trend")
##   white   white noise: no phylogenetic signal, species independent
## Raw log10 telomere length is used (not adjusted for spot length): the question is how the trait itself is
## distributed on the tree.
y <- setNames(d$tel_log10, d$tip)                                   # trait values named by tree tip label
mods <- c("BM", "lambda", "kappa", "OU", "EB", "trend", "mean", "white")
set.seed(1)                                                         # fitContinuous draws random starting values
fits <- lapply(mods, function(m) suppressWarnings(fitContinuous(tree, y, model = m, ncores = 1)))   # warnings: parameters at bounds (see below)
## ncores = 1: by default the random starts run in parallel, where set.seed() does not reach, so results vary between runs
## Table: log-likelihood, AICc and each model's own parameter. f$opt holds z0 (root value), sigsq (rate), lnL, method,
## k, aic, aicc, plus the model-specific parameter, which is what is left after removing the shared entries.
fc <- data.frame(model = mods, lnL = sapply(fits, function(f) f$opt$lnL), AICc = sapply(fits, function(f) f$opt$aicc),
                 parameter = sapply(fits, function(f) { p <- f$opt[setdiff(names(f$opt), c("z0", "sigsq", "lnL", "method", "k", "aic", "aicc"))]
                   if (length(p)) paste(names(p), signif(unlist(p), 3), sep = " = ", collapse = "; ") else "" }))
fc$dAICc <- fc$AICc - min(fc$AICc); fc <- fc[order(fc$AICc), ]       # difference from the best model; sort best first
fc$lnL <- round(fc$lnL, 2); fc$AICc <- round(fc$AICc, 2); fc$dAICc <- round(fc$dAICc, 2)
print(fc, row.names = FALSE)
## Result: lambda (0.90) is best, 17.8 AICc units ahead of kappa; BM and white noise fit poorly.
## Three models reach the limits of what can be estimated, which is why warnings are suppressed above:
##   EB: a sits at its upper bound (~0), so the fit collapses to Brownian motion (same log-likelihood): no early burst.
##   trend: the slope is at the optimizer's limit, so its log-likelihood is a lower bound. The reported slope
##          (100, the limit) depends on the random starting values (other seeds give ~97-100); its log-likelihood does not.
##   mean: not identifiable on an ultrametric tree (drift cannot be separated from the root value); the parameter
##         estimates are meaningless and only the log-likelihood (equal to BM) is usable.
## None of these affects the ranking.

## ---- 2. Pagel's lambda with a profile-likelihood CI (caper) ----
## lambda is estimated by maximum likelihood with caper::pgls on an intercept-only model (v ~ 1), i.e. the
## phylogenetic signal of the trait itself. caper reports a profile-likelihood 95% CI: all values of lambda
## whose log-likelihood is within qchisq(0.95, 1) / 2 = 1.92 of the maximum.
## caper is called as caper::pgls because 00_setup.R defines its own pgls() (phylolm wrapper).
## Example, annual plants (n = 72): lambda = 0.57, 95% CI 0.28-0.83.
## Safeguard: caper's optimizer starts at lambda = 0.5 and searches locally, so when the likelihood has two
## peaks it can stop on the lower one. We therefore also evaluate the likelihood on a grid of fixed lambda
## (0, 0.05, ..., 1); if a grid point beats caper's maximum, lambda is re-estimated near that point and caper's
## profile CI is recomputed over the full 0-1 range. In this dataset this happens only for SRA spot length
## (peaks at lambda ~0 and 0.41, differing by 0.03 log-likelihood units; either way, no significant signal).
lambda_ci <- function(v, name) {
  v  <- v[!is.na(v)]                                             # species with data for this trait
  tr <- keep.tip(tree, names(v))                                 # tree pruned to those species
  cd <- caper::comparative.data(tr, data.frame(tip = names(v), v = v), names.col = "tip", vcv = TRUE)
  m  <- caper::pgls(v ~ 1, cd, lambda = "ML")                    # lambda by maximum likelihood
  grid <- c(1e-6, seq(0.05, 1, 0.05))                            # caper's lower bound is 1e-6, not 0
  ll   <- sapply(grid, function(l) as.numeric(logLik(caper::pgls(v ~ 1, cd, lambda = l))))
  if (max(ll) > as.numeric(logLik(m)) + 1e-6) {                  # caper missed the global maximum
    g <- grid[which.max(ll)]
    m <- caper::pgls(v ~ 1, cd, lambda = "ML", param.CI = NULL,  # re-estimate within +/- 0.05 of the best grid point
                     bounds = list(lambda = c(max(1e-6, g - 0.05), min(1, g + 0.05))))
    m$bounds$lambda <- c(1e-6, 1)                                # restore full range, then compute the CI with caper
    m$param.CI <- list(lambda = caper::pgls.confint(m, "lambda", 0.95))
  }
  ci <- m$param.CI$lambda$ci.val                                 # profile-likelihood 95% CI
  ci[is.na(ci)] <- c(0, 1)[is.na(ci)]                            # NA = the limit reaches the bound (0 or 1)
  ## likelihood-ratio test against lambda = 0 (no phylogenetic signal): lambda = 0 is an ordinary mean, fitted by lm
  ll0 <- as.numeric(logLik(lm(v ~ 1)))
  data.frame(trait = name, n = nrow(cd$data), lambda = unname(m$param["lambda"]), CI_lower = ci[1], CI_upper = ci[2],
             p_vs_lambda0 = pchisq(2 * (as.numeric(logLik(m)) - ll0), 1, lower.tail = FALSE))
}
tl  <- setNames(d$tel_log10, d$tip)
sub <- function(keep) tl[which(keep)]
subsets <- rbind(
  lambda_ci(tl, "Telomere length, all land plants"),
  lambda_ci(sub(d$lifespan == "annual"), "Annual"),     lambda_ci(sub(d$lifespan == "perennial"), "Perennial"),
  lambda_ci(sub(d$growthform == "herb"), "Herb"),       lambda_ci(sub(d$growthform == "woody"), "Woody"),
  lambda_ci(sub(d$APG.IV %in% c("Core eudicots", "Eudicots")), "Eudicots"),
  lambda_ci(sub(d$APG.IV == "Monocots"), "Monocots"),
  lambda_ci(sub(d$ploidy == "diploid"), "Diploid"),     lambda_ci(sub(d$ploidy == "polyploid"), "Polyploid"))
traits <- c(genome_size_log10 = "Genome size", phenotype_num_chr = "Chromosome number",
            phenotype_genome_diversity = "Genome diversity", rangesize_log10 = "Range size",
            log_native_range = "Native range size", phenotype_height = "Plant height",
            phenotype_shoot_structure = "Shoot structure", phenotype_root_structure = "Root structure",
            abs_lat = "Latitude (absolute)", bio5 = "Peak temperature (bio5)",
            bio10 = "Warmest-quarter temperature (bio10)", gst = "Growing season temperature (gst)",
            bio17 = "Precipitation of driest quarter (bio17)", log_alt = "Altitude", log_spot = "SRA spot length")
others <- do.call(rbind, lapply(names(traits), function(v) lambda_ci(setNames(d[[v]], d$tip), traits[[v]])))
lam <- rbind(cbind(panel = "Telomere length subsets", subsets), cbind(panel = "Other traits", others))
print(lam, digits = 3, row.names = FALSE)

## lambda = 1 test for telomere length (lambda vs Brownian motion)
m_l <- pgls(tel_log10 ~ 1)$model; m_bm <- phylolm(tel_log10 ~ 1, d, tree, model = "BM")
cat("Telomere length lambda vs 1: p =", signif(pchisq(2 * (m_l$logLik - m_bm$logLik), 1, lower.tail = FALSE), 3), "\n")

save_csv(fc, "03_model_comparison.csv"); save_csv(lam, "03_lambda.csv")
