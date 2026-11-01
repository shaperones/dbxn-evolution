#!/usr/bin/env python3
"""Genomic order of the Vg-family copies of each lepidopteran genome and their position in or outside the Vg locus.

Vg locus: the chain of copies on the contig of the canonical Vg with neighbours <= --gap bp apart (default 100 kb);
oriented so that Vg is transcribed left to right.

Usage: vg_family_order.py [--gap 100000]  -> analysis/results/lep_vgfamily_order.tsv
"""
import argparse
import csv
from collections import defaultdict

SPECIES = "revision/data/lep_species_copies_2026-10-10.tsv"
CLASSES = "analysis/data/gene_tree/lep/lep_classes_broad.tsv"
COPIES = "analysis/data/miniprot/{}.miniprot.copies.tsv"
CBR = "analysis/data/cbr/cbr.lep_per_copy.tsv"
OUT = "analysis/results/lep_vgfamily_order.tsv"


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gap", type=int, default=100000, help="largest gap between neighbouring copies of one locus, bp")
    a = ap.parse_args()
    cls = {r["name"].split("|")[1]: r["class_tree"] for r in tsv(CLASSES)}
    sim = {r["name"].split("|")[1]: r["all_similarity_pct"] for r in tsv(CBR) if r["intact"] == "1"}
    out, n_out = [], defaultdict(int)
    for sp, d in ((r["tree_label"], r) for r in tsv(SPECIES)):
        rows = [r for r in tsv(COPIES.format(sp)) if r["copy_id"] in cls or r["class"] in ("DBXN", "Vg")]
        for r in rows:
            r["gene"] = cls.get(r["copy_id"], r["class"])
            r["intact"] = int(r["copy_id"] in cls)
            r["s"], r["e"] = int(r["start"]), int(r["end"])
        rows = [r for r in rows if r["gene"] in ("DBXN", "Vg")]
        for g in ("DBXN", "Vg"):                       # same counts as the species table
            assert sum(r["gene"] == g and r["intact"] for r in rows) == int(d[f"n_{g}"]), (sp, g)
            assert sum(r["gene"] == g and not r["intact"] for r in rows) == int(d[f"n_{g}_pseudo"]), (sp, g, "pseudo")
        vg = [r for r in rows if r["gene"] == "Vg"]
        anchor = sorted(vg, key=lambda r: -r["intact"])[0] if vg else None
        locus = set()
        if anchor:
            same = sorted((r for r in rows if r["contig"] == anchor["contig"]), key=lambda r: r["s"])
            i = j = same.index(anchor)
            while i > 0 and same[i]["s"] - same[i - 1]["e"] <= a.gap:
                i -= 1
            while j < len(same) - 1 and same[j + 1]["s"] - same[j]["e"] <= a.gap:
                j += 1
            seg = same[i:j + 1]
            if anchor["strand"] == "-":
                seg = seg[::-1]
            locus = {r["copy_id"] for r in seg}
        else:
            seg = []
        rest = sorted((r for r in rows if r["copy_id"] not in locus), key=lambda r: (r["contig"], r["s"]))
        for k, r in enumerate(seg + rest, 1):
            inl = r["copy_id"] in locus
            n_out[sp] += not inl
            out.append(dict(species=sp, family=d["family"], order=k, copy_id=r["copy_id"], gene=r["gene"], intact=r["intact"],
                            in_vg_locus=int(inl), contig=r["contig"], start=r["start"], end=r["end"], strand=r["strand"],
                            same_strand_as_vg=int(bool(anchor) and r["strand"] == anchor["strand"]),
                            cbr_similarity_pct=sim.get(r["copy_id"], "NA") if r["gene"] == "DBXN" and r["intact"] else "NA"))
    with open(OUT, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} copies of {len({r['species'] for r in out})} species -> {OUT}; outside the Vg locus: "
          + ", ".join(f"{s} {n}" for s, n in n_out.items() if n))


if __name__ == "__main__":
    main()
