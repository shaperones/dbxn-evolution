# Evolution of DBXN vitellogenins and green larval coloration — analysis code

Custom code for the evolutionary analyses of *Cryo-EM structures of spider and insect vitellogenin suggest their independent neofunctionalization in green camouflage evolution* (Egorkin, Burtseva, Latypova et al., Nature Communications, in revision). Runs of standard
programs (BUSCO, IQ-TREE, ASTRAL-IV, MAFFT, ClipKIT, MCMCTree, GeneRax, Cactus, halSynteny, InterProScan) are described in
Methods with their versions and options and are not repeated here.

The data — genome-derived sequences and annotations (GFF3), alignments, gene and species trees, MCMCTree, GeneRax and
synteny outputs, template sets and intermediate tables — are in [Zenodo DOI]. Unpacked into this directory, the archive
provides `analysis/data/`, `analysis/results/` and `revision/data/`; all paths in the scripts are relative to the
repository root.

## Environment

```bash
micromamba create -f environment/vg-local.yml       # R 4.4 (phylolm, ape, phytools), Python 3.12 (dendropy)
Rscript environment/install_cran.R                  # corHMM and other CRAN packages not on conda for osx-arm64
```
The annotation scripts (`analysis/annotation/`) use only the Python standard library; the programs they call
(miniprot 0.18, BLAST+ 2.17.0, MAFFT 7.526) are listed in `environment/annotation.yml`.

## Scripts

Every script starts with a description of its inputs, outputs and options.

| Script | What it does | Paper |
|---|---|---|
| `analysis/annotation/miniprot_copies.py` | Locus-level calling of Vg-family copies from miniprot protein-to-genome alignments (`miniprot -I -j2 --outs=0.4 -N 50 --outc=0.2 --gff --trans`): splits alignments at long introns, clusters them into copies, assigns the provisional class by the best template, flags disrupted copies; writes the copy table, proteins and one GFF3 model per copy | Methods: annotation; Supplementary Data 5 |
| `analysis/annotation/assign_dbxn_clade.py` | DBXN / Vg class of every copy from the gene tree alone: the largest supported clade with MperDBXN and the DBXN-like reference proteins that excludes all canonical Vgs | Fig. 7a |
| `analysis/annotation/ctck_genomic.py` | Presence of the C-terminal cystine-knot exon downstream of each copy, searched in the genome | Figs. 7b, 7e |
| `analysis/annotation/define_vg_loci.py` | Vg locus of each genome: the contig with most Vg-family copies, between two flanking BUSCO genes where both are found on it, extended to all Vg-family copies of the contig (input of the Cactus alignment) | Fig. 7e, Supplementary Fig. 11 |
| `analysis/annotation/synteny/make_synteny_inputs.py` | Copies in locus coordinates (input of the whole-locus alignment and of the synteny figures) | Fig. 7e, Supplementary Fig. 11 |
| `analysis/py/mcmctree_chronograms.py` (+ `mcmctree_nodes.py`) | 100 chronograms sampled from the MCMCTree posterior (`mcmc.txt`) on the fixed topology, used to propagate dating uncertainty | Methods: association |
| `analysis/py/make_copy_table.py` | Species × copies table: intact and disrupted DBXN and Vg, near-identical copies, full-length copies, assembly flags | Figs. 8a–b |
| `analysis/py/vg_family_order.py` | Genomic order of the copies and their position in or outside the Vg locus | Fig. 8a |
| `analysis/py/generax_summary.py` | Branches of the species tree to which the key duplications and losses map in the GeneRax reconciliations | Fig. 7d |
| `analysis/py/cbr_matrix.py`, `cbr_summary.py` | Chromophore-contact residues (CCRs; within 5 Å of a chromophore) of MperDBXN / MvirDBXN mapped onto every protein of the alignment; CCR identity and CCR similarity (identical or similar residues) per copy and per species (scripts keep the internal prefix `cbr`) | Figs. 7a, 7c, 7f, 8c |
| `analysis/py/asr_cbr.py` | Ancestral states of the CCRs at the key nodes of the gene tree (from IQ-TREE ancestral reconstruction) | Results, Supplementary Data 2 |
| `analysis/R/association_phyloglm.R` | Phylogenetic logistic regression of green coloration on DBXN copy number over posterior chronograms (Rubin's rules), robustness checks, parametric bootstrap | Fig. 8b |
| `analysis/R/association_permulation.R` | Phylogenetic permutation test | Results |
| `analysis/R/association_speccurve.R` | All combinations of analysis choices (specification curve) | Supplementary Fig. 12a |
| `analysis/R/association_corhmm.R` | Correlated evolution of coloration and copy number (corHMM, the four models of fitCorrelationTest, several random starts and seeds) and stochastic maps | Results, Supplementary Fig. 12b,c |
| `analysis/R/association_cbr_quality.R` | CCR similarity of the best DBXN copy as a predictor of green coloration, with copy number in the model; the Vg of the same species as control | Fig. 8c |

Example (main association analysis):
```bash
Rscript analysis/R/association_phyloglm.R analysis/results/treesets_mcmctree.tsv analysis/results/association_phyloglm.MCMCTREE.tsv 1000
Rscript analysis/R/association_corhmm.R analysis/results/treesets_mcmctree_mean.tsv analysis/results/association_corhmm.MCMCTREE 500 10 1,2,3 8
```

## License

MIT (see `LICENSE`).
