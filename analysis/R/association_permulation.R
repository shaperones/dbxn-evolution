#!/usr/bin/env Rscript
# Analysis 4 of revision/05_association_methods.md: permulation test for the primary phyloglm (quadrifids,
# green ~ intact DBXN copies) and the classification-free CBR predictor (copies of any class with CBR >= 70 %).
# Null: a liability evolves by Brownian motion on the tree independently of the copy number; the species with the
# highest liability get "green", as many as observed — keeps the phylogenetic clustering and the frequency of
# green coloration (permulation for a binary trait; Saputra et al. 2021). The regression is refitted on every
# simulated trait; empirical P = share of simulations with beta >= observed (one-sided, the hypothesis is beta > 0).
# Trees: simulations cycle over the trees of the main set.
#
# Usage: association_permulation.R TREESETS.tsv OUT.tsv [NPERM=1000]
suppressMessages({library(ape); library(phylolm)})
args  <- commandArgs(trailingOnly = TRUE)
sets  <- read.delim(args[1], stringsAsFactors = FALSE)
out   <- args[2]
nperm <- if (length(args) > 2) as.integer(args[3]) else 1000L
set.seed(1)

d <- read.delim("revision/data/lep_species_copies_2026-10-10.tsv", stringsAsFactors = FALSE)
cb <- read.delim("analysis/data/cbr/cbr.lep_species.tsv", stringsAsFactors = FALSE)
d$n_cbr_70 <- cb$n_cbr_70[match(d$tree_label, cb$tree_label)]
q <- d[d$family %in% c("Noctuidae", "Erebidae", "Nolidae") & d$pheno_simple %in% c("green", "not_green"), ]
q$green <- as.integer(q$pheno_simple == "green"); rownames(q) <- q$tree_label
trees <- read.tree(sets$file[sets$main == 1][1]); if (inherits(trees, "phylo")) trees <- list(trees)
trees <- lapply(trees, function(t) keep.tip(t, intersect(t$tip.label, q$tree_label)))

fit_beta <- function(y, x, tr, a0) {
  dd <- data.frame(y = y, x = x, row.names = tr$tip.label)
  m <- tryCatch(suppressWarnings(phyloglm(y ~ x, data = dd, phy = tr, method = "logistic_MPLE", btol = 30, start.alpha = a0)),
                error = function(e) NULL)
  if (is.null(m)) NA else coef(m)[2]
}
res <- list()
for (pred in c("n_DBXN", "n_cbr_70")) {
  tr <- trees[[1]]; H <- mean(node.depth.edgelength(tr)[seq_len(Ntip(tr))])
  y <- q[tr$tip.label, "green"]; x <- q[tr$tip.label, pred]
  obs <- fit_beta(y, x, tr, 1 / H)
  ngreen <- sum(y)
  null <- sapply(seq_len(nperm), function(i) {
    t <- trees[[(i - 1) %% length(trees) + 1]]
    liab <- rTraitCont(t, model = "BM")[t$tip.label]
    ys <- as.integer(rank(-liab, ties.method = "random") <= ngreen)
    fit_beta(ys, q[t$tip.label, pred], t, 1 / H)
  })
  null <- null[!is.na(null)]
  res[[pred]] <- data.frame(predictor = pred, n = length(y), n_green = ngreen, beta_obs = obs,
                            null_mean = mean(null), null_q95 = quantile(null, 0.95), null_q99 = quantile(null, 0.99),
                            p_empirical = (sum(null >= obs) + 1) / (length(null) + 1), n_perm = length(null), trees = length(trees))
  cat(sprintf("%s: beta %.3f; null mean %.3f, 95%% %.3f; empirical P = %.4f (%d permulations)\n",
              pred, obs, mean(null), quantile(null, 0.95), res[[pred]]$p_empirical, length(null)))
}
res <- do.call(rbind, res); num <- sapply(res, is.numeric); res[num] <- lapply(res[num], signif, 4)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
