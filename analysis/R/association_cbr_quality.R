#!/usr/bin/env Rscript
# CBR quality model (13.10; renamed from explore_cbr_quality.R 27.10): does the CBR "quality" of the DBXN copies (similarity to MperDBXN) predict green larvae
# beyond the number of copies? Quadrifid species WITH DBXN and known colour. Per species: max and mean CBR
# similarity over intact DBXN copies, all CBRs / bilin CBRs / lutein CBRs. Raw medians by colour, then phyloglm
# (logistic_MPLE) green ~ n_DBXN + z(score), one model per score, on the main placeholder tree (Grafen; NOT for the paper).
# Control: the same score for the canonical Vg of the species (vg_all, vg_bilin) should not predict colour.
# Usage: Rscript analysis/R/association_cbr_quality.R [TREE]  -> printed only
suppressMessages({library(ape); library(phylolm)})
args <- commandArgs(trailingOnly = TRUE)
tf <- if (length(args)) args[1] else "analysis/results/placeholder_trees/c4_b.grafen.nwk"
d <- read.delim("revision/data/lep_species_copies_2026-10-10.tsv", stringsAsFactors = FALSE)
d <- d[d$family %in% c("Noctuidae", "Erebidae", "Nolidae") & d$pheno_simple %in% c("green", "not_green") & d$n_DBXN > 0, ]
pc <- read.delim("analysis/data/cbr/cbr.lep_per_copy.tsv", stringsAsFactors = FALSE)
vg <- pc[pc$class_tree == "Vg" & pc$intact == 1, ]       # negative control: CBR similarity of the canonical Vg
pc <- pc[pc$class_tree == "DBXN" & pc$intact == 1, ]
for (k in c("all", "bilin", "lutein")) {
  v <- pc[[paste0(k, "_similarity_pct")]]
  d[[paste0(k, "_max")]]  <- sapply(d$tree_label, function(s) max(c(v[pc$species == s], NA), na.rm = TRUE))
  d[[paste0(k, "_mean")]] <- sapply(d$tree_label, function(s) mean(v[pc$species == s]))
}
d$vg_all  <- sapply(d$tree_label, function(s) max(c(vg$all_similarity_pct[vg$species == s], NA), na.rm = TRUE))
d$vg_bilin <- sapply(d$tree_label, function(s) max(c(vg$bilin_similarity_pct[vg$species == s], NA), na.rm = TRUE))
d <- d[is.finite(d$all_max) & is.finite(d$vg_all), ]
d$green <- as.integer(d$pheno_simple == "green")
rownames(d) <- d$tree_label
tr <- read.tree(tf); tr <- keep.tip(tr, intersect(tr$tip.label, d$tree_label)); d <- d[tr$tip.label, ]
cat(sprintf("species with DBXN and colour: %d (%d green)\n\n", nrow(d), sum(d$green)))
sc <- c("all_max", "all_mean", "bilin_max", "bilin_mean", "lutein_max", "lutein_mean", "vg_all", "vg_bilin")  # vg_*: control
cat("score\tmedian_green\tmedian_not_green\tOR_per_SD(adj. n_DBXN)\t95%CI\tP\tOR_copies\n")
for (s in sc) {
  d$z <- as.numeric(scale(d[[s]]))
  fits <- lapply(c(0.1, 1, 10, 100), function(a0) tryCatch(suppressWarnings(
    phyloglm(green ~ n_DBXN + z, data = d, phy = tr, method = "logistic_MPLE", btol = 30, start.alpha = a0)), error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits); m <- fits[[which.max(sapply(fits, function(x) x$penlogLik))]]
  co <- summary(m)$coefficients
  cat(sprintf("%s\t%.1f\t%.1f\t%.2f\t%.2f-%.2f\t%.3g\t%.2f\n", s, median(d[[s]][d$green == 1]), median(d[[s]][d$green == 0]),
              exp(co["z", 1]), exp(co["z", 1] - 1.96 * co["z", 2]), exp(co["z", 1] + 1.96 * co["z", 2]), co["z", 4], exp(co["n_DBXN", 1])))
}
cat(sprintf("\ncor(n_DBXN, all_max) = %.2f; cor(n_DBXN, all_mean) = %.2f\n", cor(d$n_DBXN, d$all_max), cor(d$n_DBXN, d$all_mean)))
