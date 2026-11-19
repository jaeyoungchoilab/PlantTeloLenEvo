## 02_major_clades.R
## Telomere length across major plant lineages (Results, Fig 1A). These are descriptive comparisons:
## each contrasts the two sides of a single evolutionary transition, so they are not phylogenetically replicated.
source("R/00_setup.R")

x <- all_species                      # all 264 species, including the 7 green algae
algae <- x$in_comparative == "no"         # the 7 green algae (not in the land-plant phylogeny)
bryo  <- x$APG.IV %in% c("Anthocerotophyta", "Bryophytes", "Marchantiophyta")
lyco  <- x$APG.IV %in% "Lycophytes"
eudicot <- x$APG.IV %in% c("Core eudicots", "Eudicots")

compare <- function(name, a, b) {
  data.frame(contrast = name, nA = sum(a), median_A_bp = median(x$TL_bp[a]), nB = sum(b), median_B_bp = median(x$TL_bp[b]),
             MannWhitney_p = wilcox.test(x$TL_bp[a], x$TL_bp[b], exact = FALSE)$p.value)
}
res <- rbind(
  compare("Green algae (A) vs land plants (B)", algae, !algae),
  compare("Bryophytes (A) vs vascular plants (B)", bryo, !algae & !bryo),
  compare("Bryophytes + lycophyte (A) vs seed plants (B)", bryo | lyco, !algae & !bryo & !lyco),
  compare("Monocots (A) vs eudicots (B)", x$APG.IV %in% "Monocots", eudicot))
cat("Median telomere length, all 264 species:", median(x$TL_bp), "bp\n")
print(res, digits = 3)
save_csv(res, "02_major_clades.csv")
