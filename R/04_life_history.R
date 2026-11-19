## 04_life_history.R
## Telomere length and life history: lifespan, growth form and ploidy (Results, Fig S2).
## Mann-Whitney tests ignore phylogeny; PGLS tests include phylogeny and log10 spot length.
## Result: no lifespan difference by either test. Herb vs woody and diploid vs polyploid differ by Mann-Whitney
## but not by PGLS, because woody species and polyploids are clustered in particular clades of the tree.
source("R/00_setup.R")

## The comparisons: name = c(column in d, group A, group B)
lh <- list(Lifespan = c("lifespan", "annual", "perennial"), "Growth form" = c("growthform", "herb", "woody"),
           Ploidy = c("ploidy", "diploid", "polyploid"))
## Run the same tests on each comparison; lapply gives one row per comparison and rbind stacks them
res <- do.call(rbind, lapply(names(lh), function(nm) {
  z <- lh[[nm]]                                              # e.g. c("lifespan", "annual", "perennial")
  dd <- d[d[[z[1]]] %in% z[2:3], ]                           # species in group A or B; NA (no data) is dropped,
                                                             # so n differs: 231 lifespan, 240 growth form, 182 ploidy
  dd$groupB <- as.numeric(dd[[z[1]]] == z[3])               # 0 = first group, 1 = second group
                                                             # e.g. annual = 0, perennial = 1; herb = 0, woody = 1; diploid = 0, polyploid = 1
  ## Test 1, Mann-Whitney (Wilcoxon rank-sum): ranks all species by telomere length and asks whether one group ranks
  ## higher. Ignores phylogeny (every species counted as independent) and spot length. Raw bp is used; ranks are the
  ## same on the log scale. exact = FALSE: normal approximation, needed because of tied values.
  mw <- wilcox.test(TL_bp ~ groupB, dd, exact = FALSE)
  ## Test 2, PGLS (helper in 00_setup.R): log10 telomere length with spot length as covariate and Pagel's lambda.
  ## Because groupB is 0/1, its coefficient is the difference in mean log10 TL, group B minus group A, at the same
  ## spot length. Percent difference = 10^coefficient - 1, e.g. woody vs herb: 10^0.065 - 1 = +16% (p = 0.088).
  r  <- pgls(tel_log10 ~ groupB + log_spot, dd)
  data.frame(trait = nm, groupA = z[2], groupB = z[3], nA = sum(dd$groupB == 0), nB = sum(dd$groupB == 1),
             median_A_bp = median(dd$TL_bp[dd$groupB == 0]), median_B_bp = median(dd$TL_bp[dd$groupB == 1]),
             MannWhitney_p = mw$p.value,
             PGLS_effect_log10 = r$coef["groupB", 1], PGLS_p = r$coef["groupB", 4],   # column 1 = estimate, 4 = p
             PGLS_lambda = r$model$optpar)                   # lambda of the model residuals
}))
print(res, digits = 3, row.names = FALSE)
## Expected (Mann-Whitney p / PGLS p): lifespan 0.27 / 0.96; growth form 0.0008 / 0.088; ploidy 0.023 / 0.77
save_csv(res, "04_life_history.csv")
