#!/usr/bin/env python3
"""Species x Vg/DBXN copies table for the association analysis (Lepidoptera, 194 genomes).

Sources: miniprot copy tables (all loci, intact or not), gene-tree classes of intact copies
(broad and strict DBXN clade), genomic CTCK, phenotypes, Vg locus flags, HALPER components.
Counts per species:
  n_DBXN / n_DBXN_strict / n_Vg   intact copies by gene-tree class (the copies in the gene tree)
  n_DBXN_dedup                    n_DBXN after collapsing near-identical intact copies (see below)
  n_DBXN_pseudo / n_Vg_pseudo     loci that are not intact (class by best template)
  n_DBXN_fullsize                 DBXN loci covering >= 80 % of the template, intact or not
  n_DBXN_ctck                     intact DBXN copies with a genomic CTCK (expected 0)
  n_neg / n_neg_broken            negative-control genes (apoLp, MTP, cv-d) found / not intact
Flags:
  identical_pair    >= 1 pair of intact same-class copies with >= --ident protein identity over
                    >= 90 % of the shorter copy (assembly duplication or very recent duplication)
  assembly_indel    >= 2 broken negative-control genes or a broken canonical Vg with no intact Vg:
                    frameshifts are likely assembly errors, so intactness is unreliable
  locus:<flag>      flag of the Vg locus definition (copies outside the BUSCO anchors, etc.)

Usage: make_copy_table.py --out revision/data/lep_species_copies.tsv   (needs blastp, makeblastdb)
       make_copy_table.py --out revision/data/lep_species_copies_2026-10-10.tsv --refresh-pheno
           only copies pheno_simple / pheno_3level / pheno_strict from the color table into the existing table (27.10)
"""
import argparse
import csv
import os
import subprocess
import tempfile
from collections import Counter, defaultdict

P = dict(
    species="revision/data/lepidoptera_species_v2.tsv",
    pheno="revision/data/lepidoptera_larval_color.tsv",   # 27.10: the single source of larval color (build_color_table.py)
    copies="analysis/data/miniprot",
    classes="analysis/data/gene_tree/lep/lep_classes_broad.tsv",
    ctck="analysis/data/ctck/all.ctck.tsv",
    loci="analysis/data/synteny/loci/vg_loci.tsv",
    halper="analysis/data/synteny/halper_vg_orthologs_all_vs_all/summary/halper_components.tsv",
    tree="analysis/data/species_trees/lep_02_noEut/c4_b.nwk",
)


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def read_fasta(path):
    name, seq = None, {}
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seq[name] = ""
        else:
            seq[name] += line.strip()
    return seq


def identical_pairs(seqs, cls, ident, work):
    """Pairs of same-species, same-class copies with >= ident % identity over >= 90 % of the shorter."""
    fa = os.path.join(work, "intact.faa")
    with open(fa, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "prot", "-out", fa], check=True, stdout=subprocess.DEVNULL)
    out = subprocess.run(["blastp", "-query", fa, "-db", fa, "-evalue", "1e-50", "-num_threads", "4", "-max_target_seqs", "50",
                          "-outfmt", "6 qseqid sseqid pident length qlen slen"], check=True, capture_output=True, text=True).stdout
    pairs = set()
    for line in out.splitlines():
        q, s, pid, ln, ql, sl = line.split("\t")
        if q < s and q.split("_mp")[0] == s.split("_mp")[0] and cls[q] == cls[s] \
                and float(pid) >= ident and int(ln) >= 0.9 * min(int(ql), int(sl)):
            pairs.add((q, s))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ident", type=float, default=99.0)
    ap.add_argument("--refresh-pheno", action="store_true", help="update only the color columns of the existing --out table")
    a = ap.parse_args()
    if a.refresh_pheno:
        pheno = {r["tree_label"]: r for r in tsv(P["pheno"])}
        rows = tsv(a.out); cols = list(rows[0])
        changed = 0
        for r in rows:
            new = {k: pheno[r["tree_label"]][k] for k in ("pheno_simple", "pheno_3level", "pheno_strict")}
            changed += any(r[k] != v for k, v in new.items())
            r.update(new)
        with open(a.out, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        print(f"{a.out}: color columns from {P['pheno']}; {changed} species changed")
        return

    import re
    tips = set(re.findall(r"[(,]([^(),:;]+)", open(P["tree"]).read()))
    sp_rows = [r for r in tsv(P["species"]) if r["tree_label"] in tips]
    pheno = {r["tree_label"]: r for r in tsv(P["pheno"])}
    tree_cls = {r["name"].split("|")[1]: r for r in tsv(P["classes"])}
    ctck = {r["copy_id"]: r["ctck"] for r in tsv(P["ctck"])}
    locus = defaultdict(set)
    for r in tsv(P["loci"]):
        locus[r["species"]].add(r["flags"])
    comp = {r["gene_id"].split("__")[1]: r["halper_component"] for r in tsv(P["halper"])}

    copies, seqs, cls = defaultdict(list), {}, {}
    for r in sp_rows:
        sp = r["tree_label"]
        copies[sp] = tsv(f"{P['copies']}/{sp}.miniprot.copies.tsv")
        fa = read_fasta(f"{P['copies']}/{sp}.miniprot.copies.faa")
        for c in copies[sp]:
            if c["copy_id"] in tree_cls:
                seqs[c["copy_id"]] = fa[c["copy_id"]]
                cls[c["copy_id"]] = tree_cls[c["copy_id"]]["class_tree"]
    with tempfile.TemporaryDirectory() as work:
        pairs = identical_pairs(seqs, cls, a.ident, work)
    pair_n, dup = Counter(), defaultdict(set)
    for q, s in pairs:
        pair_n[q.split("_mp")[0]] += 1
        dup[q.split("_mp")[0]].add(s)  # the second member of each pair is dropped in the dedup count

    cols = ("tree_label species family subfamily assembly pheno_simple pheno_3level pheno_strict hairy n_DBXN n_DBXN_strict n_Vg n_DBXN_dedup "
            "n_DBXN_pseudo n_Vg_pseudo n_DBXN_fullsize n_DBXN_ctck n_Vg_ctck n_DBXN_main_component n_neg n_neg_broken "
            "n_identical_pairs flags").split()
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sp_rows:
            sp = r["tree_label"]
            c = Counter()
            for x in copies[sp]:
                cid, t = x["copy_id"], tree_cls.get(x["copy_id"])
                if t:
                    k = t["class_tree"]
                    c["n_" + k] += 1
                    c["n_DBXN_strict"] += t["class_tree_strict"] == "DBXN"
                    c[f"n_{k}_ctck"] += ctck.get(cid) == "yes"
                    c["n_DBXN_dedup"] += k == "DBXN" and cid not in dup[sp]
                    c["n_DBXN_main_component"] += comp.get(cid) == "halper_component_002"
                elif x["class"] in ("DBXN", "Vg"):
                    c[f"n_{x['class']}_pseudo"] += 1
                elif x["class"] == "negative":
                    c["n_neg"] += 1
                    c["n_neg_broken"] += x["intact"] == "0"
                if x["class"] == "DBXN" and float(x["best_cov"]) >= 0.8:
                    c["n_DBXN_fullsize"] += 1
            flags = []
            if pair_n[sp]:
                flags.append("identical_pair")
            if c["n_neg_broken"] >= 2 or (c["n_Vg_pseudo"] and not c["n_Vg"]):
                flags.append("assembly_indel")
            flags += [f"locus:{f}" for f in sorted(locus[sp]) if f not in ("ok", ".")]
            if sp not in locus:
                flags.append("locus:none")
            row = {k: c.get(k, 0) for k in cols if k.startswith("n_")}
            row.update(tree_label=sp, species=r["species"], family=r["family"], subfamily=r["subfamily"], assembly=r["assembly"],
                       pheno_simple=pheno[sp]["pheno_simple"], pheno_3level=pheno[sp]["pheno_3level"],
                       pheno_strict=pheno[sp]["pheno_strict"], hairy=pheno[sp]["hairy"],
                       n_identical_pairs=pair_n[sp], flags=";".join(flags) or ".")
            w.writerow(row)
    print(f"{len(sp_rows)} species -> {a.out}")


if __name__ == "__main__":
    main()
