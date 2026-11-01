#!/usr/bin/env Rscript
# Phylogenetic logistic regression of green larval coloration on DBXN copy number (phylolm, logistic_MPLE) for several
# samples, color codings, predictors and tree sets; estimates pooled over the trees of a set (Rubin's rules).
#
# Usage: association_phyloglm.R TREESETS.tsv OUT.tsv [BOOT]
#   TREESETS.tsv: set, file (newick; tips = tree_label), main (1 = primary set)
#   BOOT: parametric bootstrap replicates for the primary model (first tree of the main set; default 0)
suppressMessages({library(ape); library(phylolm)})
args <- commandArgs(trailingOnly = TRUE)
sets <- read.delim(args[1], stringsAsFactors = FALSE)
out  <- args[2]
boot <- if (length(args) > 2) as.integer(args[3]) else 0L

d <- read.delim("revision/data/lep_species_copies_2026-10-10.tsv", stringsAsFactors = FALSE)
rownames(d) <- d$tree_label
quad <- c("Noctuidae", "Erebidae", "Nolidae")
d$green        <- c(green = 1L, not_green = 0L)[d$pheno_simple]
d$green_strict <- c(green = 1L, not_green = 0L)[d$pheno_strict]
# CCR predictors (cbr_summary.py)
cb <- read.delim("analysis/data/cbr/cbr.lep_species.tsv", stringsAsFactors = FALSE)
d <- cbind(d, cb[match(d$tree_label, cb$tree_label), c("cbr_sum_DBXN", "n_cbr_65", "n_cbr_70", "n_cbr_75")])
d$indel        <- grepl("assembly_indel", d$flags)
# N_DBXN_FROM=<tsv: species, panel_intact_DBXN>: copy counts of the panel-only miniprot run instead of the main ones
if (nzchar(Sys.getenv("N_DBXN_FROM"))) {
  po <- read.delim(Sys.getenv("N_DBXN_FROM"), stringsAsFactors = FALSE)
  k <- match(d$tree_label, po$species); stopifnot(!any(is.na(k)))
  cat(sprintf("n_DBXN from %s: changed in %d species\n", Sys.getenv("N_DBXN_FROM"), sum(d$n_DBXN != po$panel_intact_DBXN[k])))
  d$n_DBXN <- po$panel_intact_DBXN[k]
}

fit_one <- function(dd, tr, f, boot = 0) {
  tr <- keep.tip(tr, intersect(tr$tip.label, dd$tree_label)); dd <- dd[tr$tip.label, ]
  # the penalised likelihood can have several optima in alpha: start from several values, keep the best
  fits <- lapply(c(0.1, 1, 10, 100) / mean(node.depth.edgelength(tr)[seq_len(Ntip(tr))]), function(a0)
    tryCatch(suppressWarnings(phyloglm(f, data = dd, phy = tr, method = "logistic_MPLE", btol = 30, start.alpha = a0)), error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  m <- fits[[which.max(sapply(fits, function(x) x$penlogLik))]]
  if (boot > 0) m <- phyloglm(f, data = dd, phy = tr, method = "logistic_MPLE", btol = 30, start.alpha = m$alpha, boot = boot)
  s <- summary(m)$coefficients
  c(beta = s[2, 1], se = s[2, 2], alpha = m$alpha * mean(node.depth.edgelength(tr)[seq_len(Ntip(tr))]), n = nrow(dd), n_green = sum(m$y),
    boot_lo = if (boot > 0) s[2, "lowerbootCI"] else NA, boot_hi = if (boot > 0) s[2, "upperbootCI"] else NA)
}
pool <- function(dd, trees, f) {           # Rubin's rules over the trees of one set
  r <- t(sapply(trees, function(tr) tryCatch(fit_one(dd, tr, f), error = function(e) rep(NA, 7))))
  r <- r[!is.na(r[, 1]), , drop = FALSE]; m <- nrow(r)
  if (m == 0) return(NULL)
  q <- mean(r[, "beta"]); W <- mean(r[, "se"]^2); B <- if (m > 1) var(r[, "beta"]) else 0
  se <- sqrt(W + (1 + 1 / m) * B)
  data.frame(n = r[1, "n"], n_green = r[1, "n_green"], trees = m, beta = q, se = se, OR = exp(q), lo = exp(q - 1.96 * se),
             hi = exp(q + 1.96 * se), p = 2 * pnorm(-abs(q / se)), alphaT_median = median(r[, "alpha"]),
             beta_min = min(r[, "beta"]), beta_max = max(r[, "beta"]))
}

specs <- list(   # group, label, data subset, formula
  list("primary",   "quadrifids, intact DBXN copies",           function(x) x[x$family %in% quad, ],                 green ~ n_DBXN),
  list("sample",    "all species",                              function(x) x,                                      green ~ n_DBXN),
  list("sample",    "species with DBXN only",                   function(x) x[x$family %in% quad & x$n_DBXN > 0, ],  green ~ n_DBXN),
  list("sample",    "Noctuidae only",                           function(x) x[x$family == "Noctuidae", ],           green ~ n_DBXN),
  list("sample",    "without assembly_indel species",           function(x) x[x$family %in% quad & !x$indel, ],      green ~ n_DBXN),
  list("phenotype", "strict (unambiguous descriptions)",        function(x) { x <- x[x$family %in% quad, ]; x$green <- x$green_strict; x }, green ~ n_DBXN),
  list("predictor", "narrow DBXN clade (Noctuidae subclade)",   function(x) x[x$family %in% quad, ],                 green ~ n_DBXN_strict),
  list("predictor", "near-identical copies collapsed",          function(x) x[x$family %in% quad, ],                 green ~ n_DBXN_dedup),
  list("predictor", "full-size copies incl. disrupted",         function(x) x[x$family %in% quad, ],                 green ~ n_DBXN_fullsize),
  list("predictor", "log2(copies + 1)",                         function(x) x[x$family %in% quad, ],                 green ~ log2(n_DBXN + 1)),
  list("predictor", "three or more copies (yes/no)",            function(x) x[x$family %in% quad, ],                 green ~ I(n_DBXN >= 3)),
  list("predictor", "CBR similarity summed over DBXN copies",   function(x) x[x$family %in% quad, ],                 green ~ cbr_sum_DBXN),
  list("predictor", "copies of any class with CBR >= 70 %",     function(x) x[x$family %in% quad, ],                 green ~ n_cbr_70),
  list("predictor", "copies of any class with CBR >= 65 %",     function(x) x[x$family %in% quad, ],                 green ~ n_cbr_65),
  list("predictor", "copies of any class with CBR >= 75 %",     function(x) x[x$family %in% quad, ],                 green ~ n_cbr_75),
  list("sample",    "all species, copies with CBR >= 70 %",     function(x) x,                                      green ~ n_cbr_70),
  list("control",   "negative control: Vg copies",              function(x) x[x$family %in% quad, ],                 green ~ n_Vg)
)
for (sf in names(which(table(d$subfamily[d$family %in% quad]) >= 5)))
  specs[[length(specs) + 1]] <- local({ s <- sf; list("leave-one-out", paste("without", s), function(x) x[x$family %in% quad & x$subfamily != s, ], green ~ n_DBXN) })

res <- list()
for (i in seq_len(nrow(sets))) {
  trees <- read.tree(sets$file[i]); if (inherits(trees, "phylo")) trees <- list(trees)
  use <- if (sets$main[i] == 1) specs else specs[1:2]      # full grid on the primary tree set only
  for (sp in use) {
    dd <- sp[[3]](d); dd <- dd[!is.na(dd$green), ]
    r <- pool(dd, trees, sp[[4]])
    if (!is.null(r)) res[[length(res) + 1]] <- cbind(tree_set = sets$set[i], group = sp[[1]], model = sp[[2]], r)
  }
  # other estimators for the primary model, first tree of every set: GEE-type phyloglm and a phylogenetic GLMM
  dd <- d[d$family %in% quad & !is.na(d$green), ]; tr <- keep.tip(trees[[1]], dd$tree_label); dd <- dd[tr$tip.label, ]
  ig <- tryCatch(suppressWarnings(phyloglm(green ~ n_DBXN, data = dd, phy = tr, method = "logistic_IG10", btol = 30)), error = function(e) NULL)
  pg <- tryCatch(suppressWarnings(binaryPGLMM(green ~ n_DBXN, data = dd, phy = tr)), error = function(e) NULL)
  est <- function(lab, b, se, extra) cbind(tree_set = sets$set[i], group = "estimator", model = lab,
    data.frame(n = nrow(dd), n_green = sum(dd$green), trees = 1, beta = b, se = se, OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
               p = 2 * pnorm(-abs(b / se)), alphaT_median = extra, beta_min = b, beta_max = b))
  if (!is.null(ig)) res[[length(res) + 1]] <- est("phyloglm IG10 (tree 1)", coef(ig)[2], summary(ig)$coefficients[2, 2], NA)
  if (!is.null(pg)) { res[[length(res) + 1]] <- est("binaryPGLMM (tree 1)", pg$B[2], pg$B.se[2], NA)
    cat(sprintf("%s: binaryPGLMM phylogenetic variance s2 = %.2f (P for s2 = 0: %.2g)\n", sets$set[i], pg$s2, pg$P.H0.s2)) }
  cat(sprintf("%s: shortest branch %.2g of tree height\n", sets$set[i], min(trees[[1]]$edge.length) / max(node.depth.edgelength(trees[[1]]))))
  if (sets$main[i] == 1 && boot > 0) {
    dd <- d[d$family %in% quad & !is.na(d$green), ]
    b <- fit_one(dd, trees[[1]], green ~ n_DBXN, boot = boot)
    cat(sprintf("bootstrap (%d) on %s, tree 1: OR %.2f (%.2f-%.2f)\n", boot, sets$set[i], exp(b["beta"]), exp(b["boot_lo"]), exp(b["boot_hi"])))
  }
}
res <- do.call(rbind, res); num <- sapply(res, is.numeric); res[num] <- lapply(res[num], signif, 4)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res[, c("tree_set", "group", "model", "n", "n_green", "trees", "OR", "lo", "hi", "p", "alphaT_median")], row.names = FALSE, right = FALSE)
