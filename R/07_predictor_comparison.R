## 07_predictor_comparison.R
## Which of the correlated climate variables best explains telomere length? (Results, Figs S3 and S4, Tables S8-S10)
##   1. Spearman correlations among predictors
##   2. AICc comparison of single-predictor PGLS models on a common set of species
##   3. Multi-predictor PGLS models (bio5 with latitude/altitude, or with another temperature variable)
##   4. bio5 and bio10 slopes within clades and life-history groups
##   5. Temperature x growth form and temperature x lifespan interactions
## Result: bio5 has the best AICc (Akaike weight 0.81; bio10 next at dAICc = 4.0); latitude and altitude do worse than
## the spot-length-only null model and add nothing once bio5 is in the model; the bio5 slope is positive in 13 of 14
## subgroups and does not differ significantly between herbs and woody plants or between annuals and perennials.
source("R/00_setup.R")

## ---- 1. Correlations among predictors ----
## Same 34 x 34 Spearman matrix as in script 05. Expected: bio5-bio10 = 0.95, bio5-bio1 = 0.82, bio5-|latitude| = -0.59.
## bio5 and bio10 are nearly interchangeable, which is why the comparisons below are needed.
C <- cor(d[, names(predictors)], use = "pairwise.complete.obs", method = "spearman")
cat(sprintf("Spearman rho: bio5-bio10 = %.2f, bio5-bio1 = %.2f, bio5-|latitude| = %.2f\n",
            C["bio5", "bio10"], C["bio5", "bio1"], C["bio5", "abs_lat"]))

## ---- 2. AICc of single climate/geographic predictors, same species for every model ----
## Predictors are standardized within the common set; every model includes log10 spot length.
## AICc can only compare models fitted to exactly the same species, so `common` keeps the species that have all 26
## climate/geographic variables plus telomere length and spot length (230 species).
clim <- c(paste0("bio", 1:19), "gst", "gsl", "gdd0", "ngd10", "scd", "abs_lat", "log_alt")   # 26 variables
common <- d[complete.cases(d[, c("tel_log10", "log_spot", clim)]), ]
## One model per variable, tel_log10 ~ x + log_spot, plus a null model with spot length only ("none").
## setdiff(x, "none") returns x for a real variable (so it is standardized) and nothing for the null model.
ai <- do.call(rbind, lapply(c(clim, "none"), function(x) {
  f <- if (x == "none") tel_log10 ~ log_spot else reformulate(c(x, "log_spot"), "tel_log10")
  r <- pgls(f, common, standardize = setdiff(x, "none"))
  data.frame(variable = ifelse(x == "none", "spot length only (null)", x), N = nrow(r$data), AICc = aicc(r$model, nrow(r$data)),
             p = if (x == "none") NA else r$coef[x, 4])
}))
## dAICc = how much worse than the best model (rule of thumb: < 2 similar support, 4-7 considerably less, > 10 essentially
## none). Akaike weight = relative support, summing to 1 over all 27 models; read as the probability that a model is the
## best in the set. Expected: bio5 dAICc 0, weight 0.81; bio10 4.0, 0.11; null 11.5; latitude 11.9; altitude 13.5
## (latitude and altitude are worse than the null: they cost a parameter and explain nothing).
ai$dAICc <- ai$AICc - min(ai$AICc)
ai$weight <- exp(-ai$dAICc / 2) / sum(exp(-ai$dAICc / 2))
ai <- ai[order(ai$AICc), ]                                   # best model first
print(head(ai, 6), digits = 3, row.names = FALSE)
cat("Null model (spot length only): dAICc =", round(ai$dAICc[ai$variable == "spot length only (null)"], 1),
    "| latitude:", round(ai$dAICc[ai$variable == "abs_lat"], 1), "| altitude:", round(ai$dAICc[ai$variable == "log_alt"], 1), "\n")

## ---- 3. Multi-predictor PGLS (predictors standardized) ----
## bio5 in the same model as a competitor: does bio5 still matter once the other variable is accounted for?
## Each model uses the species with complete data for its own variables (234 with latitude/altitude, 253 otherwise).
multi <- list(
  "bio5 + |lat| + altitude + spot"       = tel_log10 ~ bio5 + abs_lat + log_alt + log_spot,
  "bio5 + |lat| + altitude + spot + ONT" = tel_log10 ~ bio5 + abs_lat + log_alt + log_spot + ONT,
  "bio5 + bio1 + spot"                   = tel_log10 ~ bio5 + bio1 + log_spot,
  "bio5 + bio6 + spot"                   = tel_log10 ~ bio5 + bio6 + log_spot,
  "bio5 + bio7 + spot"                   = tel_log10 ~ bio5 + bio7 + log_spot,
  "bio5 + bio10 + spot"                  = tel_log10 ~ bio5 + bio10 + log_spot)
## Note: bio7 (annual temperature range) is defined as bio5 - bio6, so "bio5 + bio6" and "bio5 + bio7" are the same
## model reparameterized (identical fit; the bio6 and bio7 p-values are both 0.277). They are not independent tests.
## variance inflation factor of each predictor (spot length included as a regressor)
## VIF = 1 / (1 - R^2) from an OLS regression of predictor x on the other predictors + spot length.
## VIF = 1: no overlap; VIF = 10: 90% of x is predictable from the others, so its coefficient is unstable (> 5-10 is
## usually a problem). Plain OLS among predictors, i.e. phylogeny is ignored here (standard practice).
vif <- function(dd, xs) sapply(xs, function(x) 1 / (1 - summary(lm(reformulate(setdiff(c(xs, "log_spot"), x), x), dd))$r.squared))
mm <- do.call(rbind, lapply(names(multi), function(nm) {
  xs <- setdiff(all.vars(multi[[nm]])[-1], c("log_spot", "ONT"))   # climate/geographic predictors ([-1] drops the response)
  r <- pgls(multi[[nm]], standardize = c(xs, "log_spot"))    # predictors and spot length per SD; ONT stays 0/1
  s <- r$coef[-1, , drop = FALSE]                            # drop the intercept row
  v <- vif(r$data, xs)
  data.frame(model = nm, term = rownames(s), slope_per_SD = s[, 1], p = s[, 4], VIF = unname(v[rownames(s)]),   # VIF is NA for spot length and ONT
             N = nrow(r$data), row.names = NULL)
}))
## Expected bio5 p: with latitude + altitude 0.0005 (latitude p = 0.66, altitude p = 0.49; VIFs 1.1-1.6, so not a
## collinearity artifact); + ONT 0.002 (ONT = 0.076, i.e. 10^0.076 - 1 = 19% longer estimates from ONT, p = 0.003);
## with bio1 0.002; with bio6 0.0004; with bio7 0.001; with bio10 0.046 but VIF = 11.8, so bio5 vs bio10 cannot really
## be separated in one model (the AICc comparison in part 2 is the better evidence).
print(mm[mm$term != "log_spot", ], digits = 3, row.names = FALSE)

## ---- 4. Slopes within clades and life-history groups ----
## Each group is a TRUE/FALSE vector over species. fabid_orders and malvid_orders are defined in 00_setup.R.
## %in% (not ==) for growth form and lifespan, so missing values become FALSE instead of NA.
groups <- list(
  "All species" = rep(TRUE, nrow(d)), Monocots = d$clade == "Monocots", Poales = d$Order == "Poales",
  Rosids = d$clade == "Rosids", Fabids = d$Order %in% fabid_orders, Malvids = d$Order %in% malvid_orders,
  Rosales = d$Order == "Rosales", Fabales = d$Order == "Fabales", Brassicales = d$Order == "Brassicales",
  Malvales = d$Order == "Malvales", Asterids = d$clade == "Asterids",
  Herbs = d$growthform %in% "herb", Woody = d$growthform %in% "woody",
  Annual = d$lifespan %in% "annual", Perennial = d$lifespan %in% "perennial")
## Temperature is NOT standardized here: each group has a different temperature range, so per-°C slopes (converted to
## % per 10 °C as in script 05) are what is comparable across groups.
wg <- do.call(rbind, lapply(names(groups), function(g) do.call(rbind, lapply(c("bio5", "bio10"), function(x) {
  r <- pgls(reformulate(c(x, "log_spot"), "tel_log10"), d[groups[[g]], ])   # fit within this group only
  b <- r$coef[x, 1]                                   # log10 bp per °C
  data.frame(group = g, variable = x, N = nrow(r$data), slope_per_degree = b, p = r$coef[x, 4],
             pct_per_10C = 100 * (10^(10 * b) - 1), lambda = r$model$optpar)
}))))
## Expected for bio5 (% per 10 °C): all species +16% (p = 0.0005); Rosids +20%; herbs +27%; woody +7% (p = 0.24);
## annual +33%; perennial +15%. Positive in 13 of 14 subgroups (only Malvales, n = 16, is negative).
## Small clades give noisy slopes (e.g. Brassicales +70% from 18 species), and lambda at exactly 0 or 1 in small
## clades means there is too little data to estimate it there.
print(wg[wg$variable == "bio5", ], digits = 3, row.names = FALSE)

## ---- 5. Temperature x life-history interactions ----
## Does the temperature slope differ between woody plants and herbs, or between perennials and annuals?
## 0/1 coding as in script 04 (NA where growth form or lifespan is unknown; those species are dropped by pgls()).
d$woody <- as.numeric(d$growthform == "woody")              # 1 = woody, 0 = herb
d$perennial <- as.numeric(d$lifespan == "perennial")        # 1 = perennial, 0 = annual
ix <- do.call(rbind, lapply(c("bio5", "bio10"), function(x) do.call(rbind, lapply(c("woody", "perennial"), function(h) {
  ## x * h expands to x + h + x:h; the p of x:h tests whether the slope differs between the two groups
  r <- pgls(as.formula(paste("tel_log10 ~", x, "*", h, "+ log_spot")))
  data.frame(variable = x, life_history = h, N = nrow(r$data), interaction_p = r$coef[paste0(x, ":", h), 4])
}))))
## Expected: bio5 x woody p = 0.069, bio5 x perennial p = 0.091, bio10 0.14 and 0.27. None significant, so the herb vs
## woody difference in part 4 (+27% vs +7%) is not statistically supported.
print(ix, digits = 3, row.names = FALSE)

save_csv(data.frame(variable = rownames(C), C, check.names = FALSE), "07_spearman_matrix.csv")
save_csv(ai, "07_AICc_single_predictors.csv"); save_csv(mm, "07_multi_predictor_pgls.csv")
save_csv(wg, "07_within_group_slopes.csv"); save_csv(ix, "07_interactions.csv")
