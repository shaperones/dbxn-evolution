#!/usr/bin/env Rscript
# Specification curve: odds ratio of green larvae per DBXN copy for every combination of sample (5), color coding (2),
# copy definition (7), dated tree (first tree of each set) and estimator (3).
#
# Usage: association_speccurve.R TREESETS.tsv OUT.tsv [CORES]
suppressMessages({library(ape); library(phylolm); library(parallel)})
args  <- commandArgs(trailingOnly = TRUE)
sets  <- read.delim(args[1], stringsAsFactors = FALSE)
out   <- args[2]
cores <- if (length(args) > 2) as.integer(args[3]) else 4L

d <- read.delim("revision/data/lep_species_copies_2026-10-10.tsv", stringsAsFactors = FALSE)
rownames(d) <- d$tree_label
quad <- c("Noctuidae", "Erebidae", "Nolidae")
cb <- read.delim("analysis/data/cbr/cbr.lep_species.tsv", stringsAsFactors = FALSE)
d <- cbind(d, cb[match(d$tree_label, cb$tree_label), c("n_cbr_65", "n_cbr_70", "n_cbr_75")])
d$indel <- grepl("assembly_indel", d$flags)

samples <- list(   # name -> subset; the first is the main choice in every dimension
  "quadrifids"             = function(x) x[x$family %in% quad, ],
  "all species"            = function(x) x,
  "species with DBXN"      = function(x) x[x$family %in% quad & x$n_DBXN > 0, ],
  "Noctuidae"              = function(x) x[x$family == "Noctuidae", ],
  "without assembly indels" = function(x) x[x$family %in% quad & !x$indel, ])
colours <- c("all descriptions" = "pheno_simple", "unambiguous only" = "pheno_strict")
copies  <- c("DBXN clade, intact"        = "n_DBXN",
             "narrow DBXN clade"         = "n_DBXN_strict",
             "near-identical merged"     = "n_DBXN_dedup",
             "with disrupted full-size"  = "n_DBXN_fullsize",
             "any Vg, CBR >= 65 %"       = "n_cbr_65",
             "any Vg, CBR >= 70 %"       = "n_cbr_70",
             "any Vg, CBR >= 75 %"       = "n_cbr_75")
estimators <- c("phyloglm MPLE", "phyloglm IG10", "binaryPGLMM")
sets <- sets[order(-sets$main), ]
trees <- setNames(lapply(sets$file, function(f) { t <- read.tree(f); if (inherits(t, "phylo")) t else t[[1]] }), sets$set)

grid <- expand.grid(estimator = estimators, tree_set = sets$set, copies = names(copies), colour = names(colours),
                    sample = names(samples), stringsAsFactors = FALSE)[, 5:1]

fit <- function(i) {
  g <- grid[i, ]
  dd <- samples[[g$sample]](d)
  dd$green <- c(green = 1L, not_green = 0L)[dd[[colours[[g$colour]]]]]
  dd$x <- dd[[copies[[g$copies]]]]
  dd <- dd[!is.na(dd$green) & !is.na(dd$x), ]
  tr <- keep.tip(trees[[g$tree_set]], intersect(trees[[g$tree_set]]$tip.label, dd$tree_label)); dd <- dd[tr$tip.label, ]
  h <- mean(node.depth.edgelength(tr)[seq_len(Ntip(tr))])
  r <- tryCatch(suppressWarnings(switch(g$estimator,
    "phyloglm MPLE" = {   # several starting alphas, keep the best penalised likelihood (as in association_phyloglm.R)
      ms <- Filter(Negate(is.null), lapply(c(0.1, 1, 10, 100) / h, function(a0) tryCatch(
        phyloglm(green ~ x, data = dd, phy = tr, method = "logistic_MPLE", btol = 30, start.alpha = a0), error = function(e) NULL)))
      m <- ms[[which.max(sapply(ms, function(z) z$penlogLik))]]; summary(m)$coefficients[2, 1:2] },
    "phyloglm IG10" = summary(phyloglm(green ~ x, data = dd, phy = tr, method = "logistic_IG10", btol = 30))$coefficients[2, 1:2],
    "binaryPGLMM"   = { m <- binaryPGLMM(green ~ x, data = dd, phy = tr); c(m$B[2], m$B.se[2]) })),
    error = function(e) c(NA, NA))
  b <- unname(r[1]); se <- unname(r[2])
  data.frame(g, n = nrow(dd), n_green = sum(dd$green), beta = b, se = se, OR = exp(b), lo = exp(b - 1.96 * se),
             hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), fit, mc.cores = cores))
res$main <- as.integer(res$sample == names(samples)[1] & res$colour == names(colours)[1] & res$copies == names(copies)[1] &
                       res$tree_set == sets$set[1] & res$estimator == estimators[1])
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], signif, 4)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
ok <- !is.na(res$beta)
cat(sprintf("%d specifications, %d fitted; OR > 1: %d; 95%% CI above 1: %d; median OR %.2f (range %.2f-%.2f); main OR %.2f\n",
            nrow(res), sum(ok), sum(res$OR[ok] > 1), sum(res$lo[ok] > 1), median(res$OR[ok]), min(res$OR[ok]), max(res$OR[ok]),
            res$OR[res$main == 1]))
