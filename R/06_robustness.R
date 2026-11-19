## 06_robustness.R
## Robustness of the temperature associations (Fig 3D, Tables S6 and S7):
##   1. dropping each family or order with >= 3 species, Rosaceae + Brassicaceae, and bryophytes + lycophyte
##   2. 200 random subsamples with one species per genus
##   3. separate fits for PacBio and ONT species, a predictor x platform interaction, and platform as an
##      additional covariate, for all 34 predictors
## By default steps 1-2 use bio5 and bio10 (fast). Set PREDS=all to run all 34 predictors, as in Table S6:
##   PREDS=all Rscript R/06_robustness.R
## run_all.R sets PREDS=all. CAUTION: running this script on its own WITHOUT PREDS=all overwrites the 34-predictor
## CSVs in results/ with bio5/bio10-only versions. (Step 3 always runs all 34 predictors.)
## Result: bio5 stays significant in all 45 drop-clade data sets (worst p = 0.006, dropping Brassicaceae), in 97.5% of
## the 200 one-species-per-genus subsamples, and within PacBio and ONT species separately.
source("R/00_setup.R"); set.seed(1)                          # seed fixes the random genus draws in step 2
preds <- if (Sys.getenv("PREDS") == "all") names(predictors) else c("bio5", "bio10")

## Fit the main-analysis model, tel_log10 ~ x + log_spot (x standardized, Pagel's lambda), to any subset of species
## and return n, slope per SD, SE and p. tryCatch turns a failed fit (e.g. too few species) into NAs, so one failure
## does not stop the hundreds of fits below.
slope <- function(data, x, covar = "log_spot") {
  r <- tryCatch(pgls(reformulate(c(x, covar), "tel_log10"), data, standardize = x), error = function(e) NULL)
  if (is.null(r)) return(c(n = NA, slope = NA, SE = NA, p = NA))
  c(n = nrow(r$data), slope = r$coef[x, 1], SE = r$coef[x, 2], p = r$coef[x, 4])
}

## ---- 1. Drop-clade jackknife ----
## If a single clade drives the association, removing it should weaken or remove it.
## `drops` is a named list of TRUE/FALSE vectors over the species in d; TRUE = species removed in that data set.
## 23 families + 20 orders with >= 3 species, plus Rosaceae + Brassicaceae together and the non-seed plants
## (bryophytes + lycophyte; bryo_lyco defined in 00_setup.R) = 45 reduced data sets.
fam <- names(which(table(d$Family) >= 3))                    # families with at least 3 species
ord <- names(which(table(d$Order) >= 3))                     # orders with at least 3 species
drops <- c(setNames(lapply(fam, function(f) d$Family == f), paste("Family:", fam)),
           setNames(lapply(ord, function(o) d$Order == o), paste("Order:", ord)),
           list("Rosaceae + Brassicaceae" = d$Family %in% c("Rosaceae", "Brassicaceae"),
                "Bryophytes + lycophyte" = d$bryo_lyco))
cat("Number of reduced data sets:", length(drops), "\n")
## For each predictor and each reduced data set: refit on d[!drops[[nm]], ] (everything except the dropped group).
## t() turns the named vector from slope() into a one-row data frame. One row per predictor x dropped group.
jk <- do.call(rbind, lapply(preds, function(x) do.call(rbind, lapply(names(drops), function(nm)
  data.frame(predictor = x, dropped = nm, t(slope(d[!drops[[nm]], ], x)))))))
## Summary per predictor: how many of the 45 fits are still p < 0.05, the range of slopes, and the worst p.
## bio5: 45/45 significant, slopes 0.031-0.043 per SD, max p = 0.006. bio10: 45/45, max p = 0.035. Latitude: 0/45.
jk_sum <- do.call(rbind, lapply(split(jk, jk$predictor), function(s)
  data.frame(predictor = s$predictor[1], n_sets = nrow(s), n_significant = sum(s$p < 0.05),
             slope_min = min(s$slope), slope_max = max(s$slope), max_p = max(s$p))))
print(jk_sum, digits = 3, row.names = FALSE)

## ---- 2. One species per genus ----
## Densely sampled genera could dominate the fit; thinning to one species per genus tests whether the result needs them.
## split() groups row numbers by genus (161 genera, 34 with more than one species); each genus contributes its only
## species, or one chosen at random. (if/else is needed because sample(i, 1) with a single number i would sample 1:i.)
genus_draw <- function() unlist(lapply(split(seq_len(nrow(d)), d$genus), function(i) if (length(i) == 1) i else sample(i, 1)))
## 200 draws; each draw is refitted for every predictor (157-159 species per fit after missing values).
## The draws do not depend on `preds`, so bio5 results are the same with or without PREDS=all.
gs <- do.call(rbind, lapply(1:200, function(rep) { dd <- d[genus_draw(), ]
  do.call(rbind, lapply(preds, function(x) data.frame(rep = rep, predictor = x, t(slope(dd, x))))) }))
## Fraction of the 200 draws with p < 0.05 (mean of TRUE/FALSE). bio5 = 0.975, bio10 = 0.885, latitude = 0.34.
gs_sum <- aggregate(cbind(fraction_significant = p < 0.05) ~ predictor, gs, mean)
print(gs_sum, digits = 3, row.names = FALSE)

## ---- 3. Sequencing platform (Table S7): within-platform slopes, predictor x platform interaction,
##         and platform as an additional covariate, for every predictor ----
## Concern: ONT species come from warmer sites than PacBio species (about 2 °C higher bio5, PGLS p = 0.019). If ONT
## reads also gave longer telomere estimates, platform alone could create a bio5 association. Three checks:
##   (a) fit PacBio and ONT species separately: does the association hold within each platform?
##   (b) x * ONT interaction: does the slope differ between platforms?
##   (c) ONT as an extra covariate: is x still associated after adjusting for platform?
dp <- d[!is.na(d$platform), ]                                # species with known platform (192 PacBio, 65 ONT)
cat("Species per platform:", paste(names(table(dp$platform)), table(dp$platform), collapse = ", "), "\n")
cat("Do ONT species occupy warmer sites? bio5 ~ ONT (PGLS) p =", signif(pgls(bio5 ~ ONT, dp)$coef["ONT", 4], 2), "\n")
platform_one <- function(x) {
  sdx <- sd(d[[x]], na.rm = TRUE)                 # full-data SD, so per-SD slopes are comparable between platforms
  ## (a) Within-platform fits. x is NOT standardized within each subset (each platform's own SD would make "1 SD" a
  ##     different amount for PacBio and ONT); instead the raw slope is multiplied by the full-data SD.
  ##     Returns c(n, slope per SD, p); NAs if the fit fails or x has no estimate in that subset.
  w <- lapply(c("PacBio", "ONT"), function(pl) {
    r <- tryCatch(pgls(reformulate(c(x, "log_spot"), "tel_log10"), dp[dp$platform == pl, ]), error = function(e) NULL)
    if (is.null(r) || !x %in% rownames(r$coef)) return(c(NA, NA, NA))
    c(nrow(r$data), r$coef[x, 1] * sdx, r$coef[x, 4]) })
  ## (b) Interaction: x * ONT expands to x + ONT + x:ONT; the p of x:ONT tests whether the slope differs by platform.
  ri <- pgls(as.formula(paste("tel_log10 ~", x, "* ONT + log_spot")), dp, standardize = x)
  ## (c) Platform as an additional covariate (ONT = 1, PacBio = 0). Note: because ONT species really do live in warmer
  ##     places, this adjustment can remove some genuine climate signal along with any platform artifact, so (a) is
  ##     the cleaner test.
  rc <- pgls(reformulate(c(x, "log_spot", "ONT"), "tel_log10"), dp, standardize = x)
  data.frame(variable = x, label = predictors[[x]],
             PacBio_n = w[[1]][1], PacBio_slope_per_SD = w[[1]][2], PacBio_p = w[[1]][3],
             ONT_n = w[[2]][1], ONT_slope_per_SD = w[[2]][2], ONT_p = w[[2]][3],
             interaction_p = ri$coef[paste0(x, ":ONT"), 4],
             platform_cov_n = nrow(rc$data), platform_cov_slope_per_SD = rc$coef[x, 1], platform_cov_p = rc$coef[x, 4])
}
plat <- do.call(rbind, lapply(names(predictors), platform_one))   # one row per predictor (all 34)
## BH across the 34 predictors for the platform-covariate p-values. For bio5, p = 0.002 but BH p = 0.075 (no predictor
## passes BH here); the robustness question is whether bio5 survives adjustment, not a new 34-predictor screen.
plat$platform_cov_p_BH <- p.adjust(plat$platform_cov_p, "BH")
## Expected for bio5: PacBio 0.028 per SD (p = 0.018, n = 188); ONT 0.071 (p = 0.019, n = 65); interaction p = 0.10;
## with platform covariate 0.033 (p = 0.002). bio10 is similar.
print(plat[plat$variable %in% c("bio5", "bio10"), ], digits = 3, row.names = FALSE)

save_csv(jk, "06_jackknife.csv"); save_csv(jk_sum, "06_jackknife_summary.csv"); save_csv(gs, "06_genus_subsampling.csv")
save_csv(plat, "06_platform.csv")
