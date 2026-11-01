#!/usr/bin/env python3
"""Placement of the key Vg duplications and losses on the species tree, from GeneRax reconciliations.

Usage: generax_summary.py --nhx lep_b_reconciliated.nhx --species-tree starting_species_tree.newick \
           --classes lep_classes_broad.tsv --species-table lepidoptera_species_v2.tsv
"""
import argparse
import csv
import re
from collections import Counter

TOK = re.compile(r"\[&&NHX[^\]]*\]|[(),:;]|[^(),:;\[]+")


class N:
    def __init__(self, parent=None):
        self.ch, self.name, self.parent, self.tags = [], None, parent, {}

    def leaves(self):
        return [self] if not self.ch else [l for c in self.ch for l in c.leaves()]

    def nodes(self):
        yield self
        for c in self.ch:
            yield from c.nodes()


def parse(txt):
    root = cur = N()
    toks = TOK.findall(txt.strip())
    i = 0
    while i < len(toks):
        t = toks[i]
        if t == "(":
            n = N(cur); cur.ch.append(n); cur = n
        elif t == ",":
            n = N(cur.parent); cur.parent.ch.append(n); cur = n
        elif t == ")":
            cur = cur.parent
        elif t == ":":
            i += 1  # branch length
        elif t.startswith("[&&NHX"):
            cur.tags = dict(x.split("=", 1) for x in t[6:-1].split(":") if "=" in x)
        elif t != ";":
            cur.name = t.strip()
        i += 1
    return root.ch[0] if len(root.ch) == 1 and not root.name else root


def mrca(root, names):
    best = None
    for n in root.nodes():
        s = {l.name for l in n.leaves()}
        if names <= s and (best is None or len(s) < len(best[1])):
            best = (n, s)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nhx", required=True)
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--classes", help="lep_classes_broad.tsv: groups = class_tree (DBXN / Vg)")
    ap.add_argument("--groups", help="TSV gene<TAB>group (alternative to --classes; gene names as in the nhx)")
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--rank", default="family", help="column used to describe species-tree nodes")
    a = ap.parse_args()

    gt = parse(open(a.nhx).read())
    st = parse(open(a.species_tree).read())
    info = {r["tree_label"]: r for r in csv.DictReader(open(a.species_table), delimiter="\t")}
    st_nodes = {n.name: n for n in st.nodes()}

    def describe(sp_node):
        n = st_nodes.get(sp_node)
        if n is None:
            return sp_node
        c = Counter(info.get(l.name, {}).get(a.rank, l.name) for l in n.leaves())
        return f"{sp_node} [{len(n.leaves())} spp: " + ", ".join(f"{k} {v}" for k, v in c.most_common()) + "]"

    group = {}
    if a.classes:
        for r in csv.DictReader(open(a.classes), delimiter="\t"):
            group[r["name"].replace("|", "__")] = r["class_tree"]
    if a.groups:
        for line in open(a.groups):
            g, k = line.rstrip("\n").split("\t")[:2]
            group[g] = k
    tips = {l.name for l in gt.leaves()}
    total = Counter(n.tags.get("D") for n in gt.nodes() if n.ch)
    print(f"gene tree: {len(tips)} genes; internal nodes D=Y {total.get('Y', 0)}, D=N {total.get('N', 0)}")
    for g in sorted(set(group.values())):
        genes = {x for x, k in group.items() if k == g and x in tips}
        if not genes:
            continue
        node, s = mrca(gt, genes)
        extra = s - genes
        print(f"\n== {g}: {len(genes)} genes; MRCA has {len(s)} tips ({len(extra)} from other groups)")
        print(f"   MRCA event D={node.tags.get('D')} -> species node {describe(node.tags.get('S'))}")
        p = node.parent
        if p is not None:
            print(f"   split from sister lineage: event D={p.tags.get('D')} -> species node {describe(p.tags.get('S'))}")
        dup = Counter(n.tags.get("S") for n in node.nodes() if n.ch and n.tags.get("D") == "Y")
        print(f"   duplications inside the group: {sum(dup.values())}; at the 8 most frequent species nodes:")
        for sp, k in dup.most_common(8):
            print(f"      {k:3d}  {describe(sp)}")


if __name__ == "__main__":
    main()
