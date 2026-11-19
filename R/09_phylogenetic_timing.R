## 09_phylogenetic_timing.R
## Is the temperature association driven by a few ancient divergences? (Results, Fig 4, Fig S6, Table S13)
##   1. Disparity through time (Harmon et al. 2003) with the morphological disparity index (MDI)
##   2. Node-height tests (Freckleton and Harvey 2006; Slater and Pennell 2014)
##   3. Regression of telomere length contrasts on bio5 contrasts through the origin (Garland et al. 1992),
##      fitted separately to younger and older nodes
##   4. Family means
## The DTT simulations take ~4 minutes per trait.
##
## Why this matters: a correlation across species can come from one or a few events deep in the tree, e.g. one old
## lineage that happened to move into warm climates and also had long telomeres (Uyeda et al. 2018). Such a pattern is
## really a single unreplicated event, however many species descend from it. These analyses ask whether trait change,
## and the telomere-temperature association, are spread across the tree, including among recent divergences.
## Result: neither trait shows an early burst (DTT); change is not concentrated at old nodes (node-height tests); and
## the contrast slope is positive and significant in both younger (< 25 Ma, p = 0.016) and older (> 25 Ma, p = 0.010)
## nodes, so the association is replicated among recent divergences, not driven by a few ancient ones.
source("R/00_setup.R")

## ---- Data and trees ----
## Species with telomere length, spot length and bio5 (253), with the tree pruned to them and rows in tip order.
dd <- d[complete.cases(d[, c("tel_log10", "log_spot", "bio5")]), ]
t2 <- keep.tip(tree, rownames(dd)); dd <- dd[t2$tip.label, ]
## Main PGLS model; its lambda (0.90) is used for the lambda-transformed tree below, and its bio5 slope (0.00642 log10
## bp per °C) is the value the all-node contrast regression in part 3 should reproduce.
m  <- phylolm(tel_log10 ~ bio5 + log_spot, dd, t2, model = "lambda"); lam <- m$optpar     # lambda = 0.90
## Contrast methods take one trait at a time, so spot length cannot be a covariate. Instead telomere length is adjusted
## for spot length first (adjust_for_spot in 00_setup.R: log10 TL - b_spot x (log_spot - mean)). Because b_spot comes
## from the same model as the bio5 slope, regressing the adjusted values on bio5 gives back exactly that bio5 slope.
tl <- adjust_for_spot(dd)[t2$tip.label]            # telomere length adjusted for spot length
b5 <- setNames(dd$bio5, rownames(dd))              # bio5 named by tip label
## Two versions of the tree (branch lengths in Myr; root at 490 Ma; 253 tips, 252 internal nodes, fully bifurcating):
## time-calibrated tree (very short branches set to 0.01 Myr so that contrasts are defined)
## A standardized contrast divides a difference by the square root of the branch lengths below it, so near-zero
## branches give huge, unstable contrasts. 16 internal branches are < 0.01 Myr; no terminal branch is.
tT <- t2; tT$edge.length[tT$edge.length < 0.01] <- 0.01
## lambda-transformed tree: internal branches multiplied by lambda, terminal branches extended to keep it ultrametric
## This is Pagel's lambda transformation (as in geiger::rescale(tree, "lambda", lam)). Shared history is shrunk by 10%
## and the difference is given to the terminal branches, so each species keeps total depth H = 490 Myr. Contrasts on
## this tree match the covariance structure of the PGLS model; the 0.01 floor adds at most 0.03 Myr to the depth.
H  <- max(node.depth.edgelength(t2)); tL <- t2; tL$edge.length <- tL$edge.length * lam
tip_edge <- tL$edge[, 2] <= Ntip(tL)               # TRUE for terminal branches (child node is a tip)
tL$edge.length[tip_edge] <- tL$edge.length[tip_edge] + (1 - lam) * H
tL$edge.length[tL$edge.length < 0.01] <- 0.01

## ---- 1. Disparity through time ----
## Relative disparity at each branching time = mean disparity (average squared pairwise difference) of the
## clades alive at that time / total disparity. MDI = area between the observed curve and the median of
## 10,000 Brownian-motion simulations on the same tree (negative = less disparity within clades than expected,
## the signature of an early burst). geiger's p-value is the share of simulations below the observed curve;
## we report it two-sided.
## Reading the curve: at the root (relative time 0) disparity is 1 by definition. Moving toward the present, each
## branching time splits the tree into more clades, and the curve shows how much of the total variation lies WITHIN
## those clades. Early burst: most variation arose early, between old clades, so within-clade disparity drops quickly
## and the curve falls below the Brownian expectation (MDI < 0). Late change: the curve stays above it (MDI > 0).
## The MDI is integrated over the whole tree (limits c(0, 1)). geiger's .area.between.curves is an internal function
## (hence :::) used by dtt() itself; abc(x, y1, y2) is the area of y2 minus y1.
abc <- geiger:::.area.between.curves
run_dtt <- function(x, seed) {
  set.seed(seed)                                                   # simulations are random
  g   <- dtt(tT, x, index = "avg.sq", nsim = 10000, plot = FALSE)  # observed curve (g$dtt) and 10,000 BM curves (g$sim)
  med <- apply(g$sim, 1, median)                                   # median simulated curve at each time
  MDI <- abc(g$times, med, g$dtt, c(0, 1))                         # area of observed minus median
  ## p1 = share of simulations whose curve lies below the observed one (area of observed minus simulation >= 0);
  ## two-sided p = 2 x the smaller tail.
  p1  <- mean(apply(g$sim, 2, function(s) abc(g$times, s, g$dtt, c(0, 1))) >= 0)
  ## g$times is relative time (0 = root, 1 = present); converted to age in Ma for plotting (Fig 4A), with the
  ## 2.5% and 97.5% quantiles of the simulations as the 95% envelope.
  list(curve = data.frame(age = (1 - g$times) * max(branching.times(tT)), dtt = g$dtt, median = med,
                          lo = apply(g$sim, 1, quantile, 0.025), hi = apply(g$sim, 1, quantile, 0.975)),
       test = data.frame(MDI = MDI, p_two_sided = 2 * min(p1, 1 - p1)))
}
dtt_tl <- run_dtt(tl, 1)
dtt_b5 <- run_dtt(b5, 2)
## Expected: telomere length MDI = -0.19, p = 0.30; bio5 MDI = 0.21, p = 0.71. Neither departs from Brownian motion.
mdi <- rbind(cbind(trait = "Telomere length", dtt_tl$test), cbind(trait = "bio5", dtt_b5$test))
print(mdi, digits = 3, row.names = FALSE)

## ---- 2. Node-height tests ----
## log |standardized contrast| regressed on node height above the root (OLS and robust regression),
## as in geiger::nh.test; zero contrasts (|u| < 1e-8; three for bio5) are set to 1e-5.
## Each internal node gives one standardized contrast: the difference between its two daughter lineages divided by
## its expected standard deviation under Brownian motion. Under Brownian motion their size should not depend on when
## the node occurred. Negative slope (big contrasts near the root, small near the tips): early burst / slowdown.
## Positive slope: more change toward the present than the tree's branch lengths predict.
## Absolute values because the sign of a contrast is arbitrary; log because |u| is strongly right-skewed (zeros would
## give log(0) = -Inf, hence the 1e-5 floor). Robust regression (MASS::rlm) down-weights outlying contrasts; rlm gives
## no p-value, so it is computed from its t value with n_nodes - 2 df.
nh <- function(x, tr, tree_label, trait) {
  bt <- branching.times(tr); ht <- max(bt) - bt                    # node age (Ma) -> height above the root
  u  <- abs(pic(x, tr)); u[u < 1e-8] <- 1e-5                       # |standardized contrasts|, one per node, same order as bt
  ols <- summary(lm(log(u) ~ ht))$coefficients
  rob <- summary(rlm(log(u) ~ ht, maxit = 100))$coefficients
  data.frame(trait = trait, tree = tree_label, n_nodes = Nnode(tr), OLS_slope = ols[2, 1], OLS_p = ols[2, 4],
             robust_slope = rob[2, 1], robust_p = 2 * pt(-abs(rob[2, 3]), Nnode(tr) - 2))
}
## Expected (robust): time-calibrated tree, telomere length slope > 0, p = 3.3e-5; bio5 > 0, p = 1.8e-12, i.e. no sign
## of early change, if anything more recent change. On the lambda-transformed tree telomere length shows no trend
## (p = 0.53): shrinking shared history by lambda accounts for the extra recent variation, so the contrasts used in
## part 3 meet the Brownian-motion assumption. (bio5 keeps a positive trend, p = 0.005; that tree uses the lambda of
## the telomere-length model, and in part 3 bio5 is the predictor, so its contrasts need not be standardized.)
NH <- rbind(nh(tl, tT, "time-calibrated", "Telomere length"), nh(b5, tT, "time-calibrated", "bio5"),
            nh(tl, tL, "lambda-transformed", "Telomere length"), nh(b5, tL, "lambda-transformed", "bio5"))
print(NH, digits = 3, row.names = FALSE)

## ---- 3. Contrast regression through the origin, by node age (Table S13) ----
## On the lambda-transformed tree the contrast slope equals the PGLS slope (Garland and Ives 2000; Blomberg et al. 2012).
## Telomere-length contrasts (ux) are regressed on bio5 contrasts (uy), one pair per node. Through the origin (- 1)
## because the sign of each contrast depends only on which daughter is subtracted from which, so the line must pass
## through (0, 0). Each node's age comes from the untransformed tree (names(ux) are node numbers).
ux <- pic(tl, tL); uy <- pic(b5, tL); age <- branching.times(t2)[names(ux)]
## Fit the through-origin regression to the nodes selected by k (TRUE/FALSE). With no intercept there is 1 parameter,
## so the 95% CI uses a t distribution with n - 1 df.
fit_nodes <- function(k, label) {
  f <- summary(lm(ux[k] ~ uy[k] - 1))$coefficients
  n <- sum(k)
  data.frame(nodes = label, n_nodes = n, slope = f[1, 1], SE = f[1, 2],
             CI_lower = f[1, 1] - qt(0.975, n - 1) * f[1, 2], CI_upper = f[1, 1] + qt(0.975, n - 1) * f[1, 2], p = f[1, 4])
}
med_age <- median(age)                                       # 24.9 Ma
## All 252 nodes: slope 0.0064220, matching the PGLS slope 0.0064219 (tiny difference from the 0.01 branch floor)
CR <- fit_nodes(rep(TRUE, length(ux)), "All nodes")
cc <- 25   # median node age (rounded), which splits the 252 nodes 126/126
## If the association came only from a few ancient divergences, the younger half would show no slope.
## Expected: younger than 25 Ma slope 0.0056 (95% CI 0.0011-0.0102), p = 0.016; older 0.0079 (0.0019-0.0138), p = 0.010.
## The two slopes overlap broadly; a formal test (lm(ux ~ uy + uy:older - 1)) gives p = 0.56 for the difference.
CR <- rbind(CR, fit_nodes(age <= cc, sprintf("Younger than %d Ma", cc)), fit_nodes(age > cc, sprintf("Older than %d Ma", cc)))
cat(sprintf("Median node age = %.0f Ma; PGLS slope = %.6f\n", med_age, coef(m)["bio5"]))
print(CR, digits = 3, row.names = FALSE)

## ---- 4. Family means (families with >= 3 species) ----
## A simple, non-phylogenetic illustration (Fig 4C): mean spot-adjusted telomere length (mean of log10 values, i.e. the
## geometric mean on the bp scale) and mean bio5 per family, for the 22 families with at least 3 species.
fam <- do.call(rbind, lapply(split(rownames(dd), dd$Family), function(k)      # k = species in one family
  data.frame(Family = dd[k[1], "Family"], n = length(k), tl = mean(tl[k]), bio5 = mean(dd[k, "bio5"]))))
fam3 <- fam[fam$n >= 3, ]
## Spearman rank correlation; exact = FALSE because of tied ranks. Expected: rho = 0.46, p = 0.031.
rho <- cor.test(fam3$bio5, fam3$tl, method = "spearman", exact = FALSE)
cat(sprintf("Family means (n = %d families): Spearman rho = %.2f, p = %.3f\n", nrow(fam3), rho$estimate, rho$p.value))

save_csv(mdi, "09_dtt_mdi.csv"); save_csv(rbind(cbind(trait = "Telomere length", dtt_tl$curve), cbind(trait = "bio5", dtt_b5$curve)), "09_dtt_curves.csv")
save_csv(NH, "09_node_height_tests.csv"); save_csv(CR, "09_contrast_regression_by_age.csv"); save_csv(fam3, "09_family_means.csv")
