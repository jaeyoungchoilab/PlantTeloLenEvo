## 05_pgls_34_predictors.R
## Telomere length regressed on each of 34 predictors (Table 1, Table S5, Fig 3A):
## OLS and PGLS, with and without log10 spot length as a covariate; Benjamini-Hochberg and
## Li & Ji effective-number-of-tests corrections.
## Result: 17 of 34 predictors are significant by OLS (e.g. latitude, genome size), but under PGLS only the warm-season
## temperatures remain after correction: bio5 (BH p = 0.018, Li & Ji p = 0.009) and, marginally, bio10 (0.050 / 0.050).
source("R/00_setup.R")

## One predictor at a time. Predictors are standardized (slope = change in log10 telomere length per SD).
## x = predictor name, covar = covariate name ("log_spot") or NULL for no covariate.
fit_one <- function(x, covar = NULL) {
  f <- reformulate(c(x, covar), "tel_log10")                 # builds the formula, e.g. tel_log10 ~ bio5 + log_spot
  r <- pgls(f, standardize = x)                              # PGLS helper (00_setup.R): complete cases, x scaled to SD 1, lambda by ML
  o <- summary(lm(f, r$data))$coefficients                   # OLS on the same species and the same standardized x,
                                                             # so OLS and PGLS differ only in the phylogeny
  data.frame(variable = x, label = predictors[[x]], N = nrow(r$data),   # N differs by predictor (253 for most climate
                                                             # variables; 36-171 for phenotype and range variables)
             OLS_slope = o[x, 1], OLS_p = o[x, 4],           # column 1 = estimate, 4 = p
             PGLS_slope = r$coef[x, 1], PGLS_SE = r$coef[x, 2], PGLS_p = r$coef[x, 4], PGLS_lambda = r$model$optpar)
}
## All 34 predictors (list `predictors` in 00_setup.R) for one specification, then Benjamini-Hochberg (false discovery
## rate) correction across the 34 tests, separately for PGLS and OLS p-values.
run <- function(covar = NULL, tag) {
  r <- do.call(rbind, lapply(names(predictors), fit_one, covar = covar))   # one row per predictor
  r$PGLS_p_BH <- p.adjust(r$PGLS_p, "BH")
  r$OLS_p_BH  <- p.adjust(r$OLS_p, "BH")
  cbind(spec = tag, r)
}
main  <- run("log_spot", "log10 spot length covariate (main analysis)")
nocov <- run(NULL, "no covariate")                           # check that the result does not depend on the covariate

## Effective number of independent tests (Li and Ji 2005)
## From the eigenvalues of the predictors' Spearman correlation matrix: each eigenvalue contributes
## 1 if it is >= 1, plus its fractional part. Correlated predictors give fewer effective tests.
## Uncorrelated predictors would give 34 eigenvalues of 1 (Meff = 34); correlated groups (e.g. bio1, bio5, bio10, gdd0,
## gst: all "how warm") pile into a few large eigenvalues. Here: 8 eigenvalues >= 1 (contributing 8 + 4.24) and 23
## between 0 and 1 (contributing 4.98), so Meff = 17.2: the 34 predictors behave like about 17 independent tests.
## Each correlation uses the species that have both variables (pairwise.complete.obs), so the matrix is not exactly a
## valid correlation matrix and has 3 small negative eigenvalues (down to -0.15); pmax(ev, 0) sets them to 0.
C  <- cor(d[, names(predictors)], use = "pairwise.complete.obs", method = "spearman")
ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
Meff <- sum((ev >= 1) + (pmax(ev, 0) - floor(pmax(ev, 0))))
cat("Effective number of independent tests, Meff =", round(Meff, 1), "\n")
## Sidak correction for Meff tests: chance of at least one p this small among Meff independent tests.
## e.g. bio5: 1 - (1 - 0.00053)^17.2 = 0.0091
sidak <- function(p) pmin(1, 1 - (1 - p)^Meff)
main$PGLS_p_LiJi <- sidak(main$PGLS_p); nocov$PGLS_p_LiJi <- sidak(nocov$PGLS_p)

## Strongest PGLS associations in the main analysis
top <- main[order(main$PGLS_p), c("variable", "N", "PGLS_slope", "PGLS_SE", "PGLS_p", "PGLS_p_BH", "PGLS_p_LiJi", "PGLS_lambda")]
print(head(top, 8), digits = 3, row.names = FALSE)

## Effect size for bio5 on the original scale: % change in telomere length per 10 °C
## Refit without standardizing, so the slope is per °C (0.00642 log10 bp per °C; x SD of bio5, 5.95 °C, = 0.0382 per SD).
## Over 10 °C: 10^(10 x 0.00642) - 1 = 0.16, i.e. about 16% longer telomeres.
b5 <- pgls(tel_log10 ~ bio5 + log_spot)
cat(sprintf("bio5 slope = %.6f log10 bp per °C (lambda = %.3f): %.0f%% longer telomeres per 10 °C\n",
            b5$coef["bio5", 1], b5$model$optpar, 100 * (10^(10 * b5$coef["bio5", 1]) - 1)))

save_csv(rbind(main, nocov), "05_pgls_34_predictors.csv")   # 68 rows: 34 predictors x 2 specifications
save_csv(data.frame(component = seq_along(ev), eigenvalue = ev), "05_spearman_eigenvalues.csv")
