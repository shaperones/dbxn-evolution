#!/usr/bin/env Rscript
# Analysis 3 of revision/05_association_methods.md: events behind the association (quadrifids).
#  3a. Correlated evolution of green larval coloration and "many DBXN copies" (>= K intact copies; K = 3, the
#      median among species with DBXN; variants 2 and 4): corHMM::fitCorrelationTest — independent and dependent
#      models, each with one and two rate classes. Correlation is supported only if a dependent model has a lower
#      AICc than BOTH independent models (incl. the hidden-rates one).
#  3b. Stochastic maps (NSIM) on the correlated model (K = 3, first tree of the main set): gains and losses of
#      green coloration on branches with few / many copies. corHMM states 1-4 = (not green, few), (not green, many),
#      (green, few), (green, many) — checked against the tip data below; in the maps they are named "green_many"
#      ("0_0", "0_1", "1_0", "1_1"). Counts are also divided by the time the maps spend in the source state
#      (rate per unit of tree depth): branches with few copies are much longer in total, so raw counts of gains in
#      the two backgrounds are not comparable.
#  3c. Independent losses of DBXN: Dollo parsimony on the species tree (DBXN gained once at the quadrifid crown):
#      every maximal clade without intact DBXN = one loss; with colour and pseudogene status.
#
# Usage: association_corhmm.R TREESETS.tsv OUT_PREFIX [NSIM=500] [NSTARTS=10] [SEEDS=1,2,3] [CORES=8]
#   writes OUT_PREFIX.corr_models.tsv, OUT_PREFIX.simmap_counts.tsv, OUT_PREFIX.dbxn_losses.tsv
# 25.10: the four models of fitCorrelationTest are fitted here with corHMM(nstarts = NSTARTS) - several random starts,
# best likelihood kept, and this is repeated for every seed in SEEDS; the best fit over the seeds is used and the
# number of seeds that reached it (AICc within 0.1) is written as a convergence check (column seeds_at_best). With the
# single default start the optimum depended on the random stream (AICc of the same model on the same tree differed by
# up to 13 between runs; with 10 starts still by 3.9 on the main tree; revision/09 §8).
suppressMessages({library(ape); library(corHMM)})
args <- commandArgs(trailingOnly = TRUE)
sets <- read.delim(args[1], stringsAsFactors = FALSE)
pre  <- args[2]
nsim <- if (length(args) > 2) as.integer(args[3]) else 500L
nstarts <- if (length(args) > 3) as.integer(args[4]) else 10L
seeds <- if (length(args) > 4) as.integer(strsplit(args[5], ",")[[1]]) else 1:3
cores <- if (length(args) > 5) as.integer(args[6]) else 8L
set.seed(seeds[1])

d <- read.delim("revision/data/lep_species_copies_2026-10-10.tsv", stringsAsFactors = FALSE)
quad <- c("Noctuidae", "Erebidae", "Nolidae")
d$green <- c(green = 1L, not_green = 0L)[d$pheno_simple]
q <- d[d$family %in% quad & !is.na(d$green), ]
rownames(q) <- q$tree_label

# ---- 3a -------------------------------------------------------------------------------------------------------
corr_test <- function(tr, k) {
  tr <- keep.tip(tr, intersect(tr$tip.label, q$tree_label))
  dat <- data.frame(sp = tr$tip.label, green = q[tr$tip.label, "green"], many = as.integer(q[tr$tip.label, "n_DBXN"] >= k))
  rc <- getRateCatMat(2)                             # the models of corHMM::fitCorrelationTest (v2.8), with nstarts
  ind1 <- getStateMat4Dat(dat, "ARD", collapse = FALSE, indep = TRUE)$rate.mat
  cor1 <- getStateMat4Dat(dat, "ARD", collapse = FALSE, indep = FALSE)$rate.mat
  fit <- function(rcat, mat) {                       # best over seeds; attribute seeds_at_best
    f <- lapply(seeds, function(sd) { set.seed(sd); suppressMessages(suppressWarnings(
      corHMM(phy = tr, data = dat, rate.cat = rcat, rate.mat = mat, nstarts = nstarts, n.cores = cores))) })
    a <- sapply(f, function(m) m$AICc)
    best <- f[[which.min(a)]]; best$seeds_at_best <- sum(a - min(a) < 0.1); best
  }
  sink(nullfile()); on.exit(sink())
  fits <- list(independent_model_fit = fit(1, ind1), hidden_Markov_independent_model_fit = fit(2, getFullMat(list(ind1, ind1), RateClassMat = rc)),
               correlated_model_fit = fit(1, cor1), hidden_Markov_correlated_model_fit = fit(2, getFullMat(list(cor1, cor1), RateClassMat = rc)))
  sink(); on.exit()
  a <- sapply(fits, function(m) m$AICc)
  w <- exp(-0.5 * (a - min(a))); w <- w / sum(w)
  list(tab = data.frame(K = k, model = c("independent", "independent, hidden rates", "correlated", "correlated, hidden rates"),
                        AICc = round(a, 2), dAICc = round(a - min(a), 2), weight = round(w, 3),
                        supported = min(a[3:4]) < min(a[1:2]),
                        seeds_at_best = sapply(fits, function(m) m$seeds_at_best), seeds = length(seeds)),
       fits = fits, tree = tr, dat = dat)
}
rows <- list(); keep <- NULL
for (i in seq_len(nrow(sets))) {
  tr <- read.tree(sets$file[i]); if (!inherits(tr, "phylo")) tr <- tr[[1]]
  for (k in if (sets$main[i] == 1) c(3, 2, 4) else 3) {
    r <- corr_test(tr, k)
    rows[[length(rows) + 1]] <- cbind(tree_set = sets$set[i], r$tab)
    if (sets$main[i] == 1 && k == 3) keep <- r
    cat(sprintf("%s K=%d: best %s (dAICc to best independent %.1f)\n", sets$set[i], k, r$tab$model[which.min(r$tab$AICc)],
                min(r$tab$AICc[1:2]) - min(r$tab$AICc[3:4])))
  }
}
write.table(do.call(rbind, rows), paste0(pre, ".corr_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 3b -------------------------------------------------------------------------------------------------------
m <- keep$fits$correlated_model_fit
tip_state <- 1 + 2 * keep$dat$green + keep$dat$many                    # 1 = 00, 2 = 01, 3 = 10, 4 = 11
obs <- apply(m$tip.states, 1, which.max)
stopifnot(all(obs[keep$dat$sp] == tip_state))                         # state coding as assumed
maps <- makeSimmap(keep$tree, keep$dat, m$solution, rate.cat = 1, nSim = nsim)
ev <- c(green_gain_few = "0_0>1_0", green_gain_many = "0_1>1_1", green_loss_few = "1_0>0_0", green_loss_many = "1_1>0_1",
        copies_up_notgreen = "0_0>0_1", copies_up_green = "1_0>1_1", copies_down_notgreen = "0_1>0_0", copies_down_green = "1_1>1_0")
stopifnot(all(unique(unlist(lapply(maps[[1]]$maps, names))) %in% c("0_0", "0_1", "1_0", "1_1")))
cnt <- t(sapply(maps, function(s) {
  tr <- unlist(lapply(s$maps, function(x) { st <- names(x); if (length(st) > 1) paste(head(st, -1), tail(st, -1), sep = ">") }))
  sapply(ev, function(e) sum(tr == e))
}))
tim <- t(sapply(maps, function(s) {                                   # time in each state, summed over branches
  x <- unlist(s$maps); sapply(c("0_0", "0_1", "1_0", "1_1"), function(st) sum(x[names(x) == st]))
}))
rate <- cnt / tim[, sub(">.*", "", ev), drop = FALSE]
sm <- data.frame(event = names(ev), transition = ev, mean = round(colMeans(cnt), 2),
                 lo95 = apply(cnt, 2, quantile, 0.025), hi95 = apply(cnt, 2, quantile, 0.975),
                 time_in_source = round(colMeans(tim[, sub(">.*", "", ev), drop = FALSE]), 4),
                 rate_mean = signif(colMeans(rate), 4), rate_lo95 = signif(apply(rate, 2, quantile, 0.025), 4),
                 rate_hi95 = signif(apply(rate, 2, quantile, 0.975), 4))
write.table(sm, paste0(pre, ".simmap_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sm, row.names = FALSE)

# ---- 3c -------------------------------------------------------------------------------------------------------
tr <- read.tree(sets$file[sets$main == 1][1]); if (!inherits(tr, "phylo")) tr <- tr[[1]]
qa <- d[d$family %in% quad, ]; rownames(qa) <- qa$tree_label           # all quadrifids incl. phenotype no_data
tr <- keep.tip(tr, intersect(tr$tip.label, qa$tree_label))
absent <- qa[tr$tip.label, "n_DBXN"] == 0; names(absent) <- tr$tip.label
desc <- prop.part(tr)                                                   # tips under each internal node
clades <- c(lapply(seq_len(Ntip(tr)), function(i) tr$tip.label[i]), lapply(desc, function(x) attr(desc, "labels")[x]))
allabs <- which(sapply(clades, function(x) all(absent[x])))
maximal <- allabs[!sapply(allabs, function(i) any(sapply(setdiff(allabs, i), function(j) all(clades[[i]] %in% clades[[j]]) && length(clades[[j]]) > length(clades[[i]]))))]
loss <- do.call(rbind, lapply(maximal, function(i) {
  x <- qa[clades[[i]], ]
  data.frame(n_species = nrow(x), subfamilies = paste(unique(x$subfamily), collapse = ", "),
             green = sum(x$pheno_simple == "green"), not_green = sum(x$pheno_simple == "not_green"),
             no_data = sum(!x$pheno_simple %in% c("green", "not_green")),
             with_DBXN_pseudogene = sum(x$n_DBXN_pseudo > 0), species = paste(x$species, collapse = "; "))
}))
loss <- loss[order(-loss$n_species), ]
write.table(loss, paste0(pre, ".dbxn_losses.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("DBXN losses (Dollo): %d events, %d species; green among them: %d\n", nrow(loss), sum(loss$n_species), sum(loss$green)))
print(loss[, 1:6], row.names = FALSE)
