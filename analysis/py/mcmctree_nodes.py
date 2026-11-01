#!/usr/bin/env python3
"""Posterior mean ages and 95% HPD (Ma) of named nodes (MRCA of two species) from MCMCTree FigTree.tre files.

Usage: mcmctree_nodes.py --dir analysis/data/mcmctree/araneae --schemes A B1 C1
"""
import argparse
import csv
import re


class N:
    def __init__(s, p=None):
        s.ch, s.name, s.bl, s.hpd, s.p = [], None, 0.0, None, p

    def leaves(s):
        return [s] if not s.ch else [l for c in s.ch for l in c.leaves()]

    def nodes(s):
        yield s
        for c in s.ch:
            yield from c.nodes()


def parse_figtree(txt):
    start = txt.index("(")
    t = txt[start:txt.index(";", start)]
    toks = re.findall(r"\[&95%HPD=\{[^}]*\}\]|[(),:]|[^(),:\[\s]+", t)
    root = cur = N()
    i = 0
    while i < len(toks):
        x = toks[i]
        if x == "(":
            n = N(cur); cur.ch.append(n); cur = n
        elif x == ",":
            n = N(cur.p); cur.p.ch.append(n); cur = n
        elif x == ")":
            cur = cur.p
        elif x == ":":
            i += 1; cur.bl = float(toks[i])
        elif x.startswith("[&95%HPD"):
            cur.hpd = [float(v) for v in re.findall(r"[\d.]+", x[x.index("{"):])]
        else:
            cur.name = x
        i += 1
    return root.ch[0] if len(root.ch) == 1 else root


def height(n):
    return 0.0 if not n.ch else n.ch[0].bl + height(n.ch[0])


def mrca(root, names):
    best = None
    for n in root.nodes():
        s = {l.name for l in n.leaves()}
        if names <= s and (best is None or len(s) < len(best[1])):
            best = (n, s)
    return best[0]


SPIDER_NODES = [
    ("R root: OCC vs Sparassidae+marronoids", "spParPse", "spHetDav"),
    ("SM Sparassidae vs marronoids", "spEraAtri", "spHetDav"),
    ("S Sparassidae crown", "spTheSev", "spHetDav"),
    ("HD Heteropodinae vs Micrommata+Deleninae", "spMicVir", "spHetDav"),
    ("MD Micrommata vs Deleninae", "spMicVir", "spDelCan"),
    ("Heteropodinae crown", "spPseKas", "spHetDav"),
    ("HS Sinopoda vs Heteropoda", "spSinOki", "spHetDav"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--schemes", nargs="+", default=["A", "B1", "C1"])
    ap.add_argument("--unit", type=float, default=100.0, help="Ma per MCMCTree time unit")
    a = ap.parse_args()
    mp = {r["new_name"]: r["tree_label"] for r in csv.DictReader(open(f"{a.dir}/mapping_species.tsv"), delimiter="\t")}
    trees = {}
    for s in a.schemes:
        t = parse_figtree(open(f"{a.dir}/{s}/FigTree.tre").read())
        for l in t.leaves():
            l.name = mp[l.name]
        trees[s] = t
    print("node\t" + "\t".join(a.schemes))
    for label, x, y in SPIDER_NODES:
        cells = []
        for s in a.schemes:
            n = mrca(trees[s], {x, y})
            lo, hi = (n.hpd[0] * a.unit, n.hpd[1] * a.unit) if n.hpd else (float("nan"), float("nan"))
            cells.append(f"{height(n) * a.unit:.1f} ({lo:.1f}-{hi:.1f})")
        print(label + "\t" + "\t".join(cells))


if __name__ == "__main__":
    main()
