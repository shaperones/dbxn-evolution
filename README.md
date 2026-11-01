# DBXN evolution — analysis code

Code for the evolutionary analyses of vitellogenin-derived pigments (dibilinoxanthinins, DBXNs) in Lepidoptera and spiders.
Data: [Zenodo DOI].

## Environment

```bash
micromamba create -f environment/vg-local.yml
Rscript environment/install_cran.R
```

## Scripts

| Script | Purpose |
|---|---|
| `analysis/annotation/miniprot_copies.py` | Vg-family copies from miniprot alignments |
| `analysis/annotation/assign_dbxn_clade.py` | DBXN / Vg class from the gene tree |
| `analysis/annotation/ctck_genomic.py` | CTCK exon search in the genome |
| `analysis/annotation/define_vg_loci.py` | Vg loci for the whole-locus alignment |
| `analysis/annotation/synteny/make_synteny_inputs.py` | Copies in locus coordinates |
| `analysis/py/mcmctree_chronograms.py` | Posterior chronograms from MCMCTree |
| `analysis/py/make_copy_table.py`, `vg_family_order.py` | Copy counts and genomic order |
| `analysis/py/generax_summary.py` | Duplications and losses from GeneRax |
| `analysis/py/cbr_matrix.py`, `cbr_summary.py` | CCR identity and similarity |
| `analysis/py/asr_cbr.py` | Ancestral states of the CCRs |
| `analysis/R/association_*.R` | Association of coloration with DBXN copies |

## License

MIT
