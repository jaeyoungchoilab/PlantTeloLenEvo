## 01_sequencing_covariates.R
## Technical covariates of long-read telomere length estimates (Results, Fig S1, Table S2).
source("R/00_setup.R")

## PGLS models: telomere length on spot length and platform, and spot length on the main predictors
models <- list(
  "tel_log10 ~ log_spot"       = tel_log10 ~ log_spot,
  "tel_log10 ~ ONT"            = tel_log10 ~ ONT,
  "tel_log10 ~ log_spot + ONT" = tel_log10 ~ log_spot + ONT,
  "log_spot ~ ONT"             = log_spot ~ ONT,
  "log_spot ~ bio5"            = log_spot ~ bio5,
  "log_spot ~ abs_lat"         = log_spot ~ abs_lat,
  "log_spot ~ bio10"           = log_spot ~ bio10)
## Fit each model and stack the results into one table, one row per predictor
## (the model with two predictors gives two rows). Slopes are in raw units.
tech <- do.call(rbind, lapply(names(models), function(nm) {        # loop over the 7 models; rbind stacks the outputs
  r <- pgls(models[[nm]])                                           # fit the PGLS (helper in 00_setup.R)
  s <- r$coef[-1, , drop = FALSE]                                   # coefficient table without the intercept row
  data.frame(model = nm, n = nrow(r$data), term = rownames(s), estimate = s[, 1], SE = s[, 2], p = s[, 4],
             lambda = r$model$optpar, row.names = NULL)
}))
print(tech, digits = 3)
cat("ONT / PacBio fold difference in telomere length:", round(10^tech$estimate[tech$model == "tel_log10 ~ ONT"], 2), "\n")
cat("Species per platform:\n"); print(table(d$platform, useNA = "ifany"))

## Contig N50 of the genome assembly vs telomere length (Spearman)
n50 <- cor.test(d$N50_kb, d$tel_log10, method = "spearman", exact = FALSE)
cat(sprintf("Contig N50: Spearman rho = %.3f, p = %.2f, n = %d\n", n50$estimate, n50$p.value, sum(complete.cases(d$N50_kb, d$tel_log10))))

tech <- rbind(tech, data.frame(model = "Spearman: tel_log10 vs contig N50", n = sum(complete.cases(d$N50_kb, d$tel_log10)),
                               term = "N50_kb", estimate = unname(n50$estimate), SE = NA, p = n50$p.value, lambda = NA))
save_csv(tech, "01_technical_covariates.csv")
