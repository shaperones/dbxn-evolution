#!/usr/bin/env python3
"""Final DBXN / Vg class of every copy from the gene tree (no CBR involved).

Rule (plan §3.1): clades are defined relative to the outgroup Vg (other insect orders); the DBXN clade
is the largest clade with UFBoot >= --min-support that contains MperDBXN and all RefSeq DBXN-like
references and no Lepidoptera canonical Vg reference (--definition broad, default). The smallest such
clade (--definition strict) is also reported: because every DBXN reference is from Noctuidae, it stops
at the Noctuidae subclade by construction and must not be read as "DBXN only in Noctuidae".
DBXN arose from Vg within Lepidoptera, so the canonical Vg are paraphyletic around it:
copies inside the smallest clade holding all Lepidoptera canonical Vg references but outside
the DBXN clade are Vg; anything else is 'unassigned' and is reported with its neighbours.
If no clade satisfies the rule, the script says why (which canonical Vg sits inside, or
which reference falls outside) instead of guessing.

Usage: assign_dbxn_clade.py --tree lep.treefile --set lep_tree_set.tsv --out lep_classes
IQ-TREE node labels 'SH-aLRT/UFBoot' or 'UFBoot' are both read.
Outputs: <out>.tsv (per copy), <out>.summary.txt
"""
import argparse
import csv
import re
from collections import Counter, defaultdict

LABEL = re.compile(r"\s*('(?:[^']|'')*'|[(),:;]|[^(),:;']+)")


class Node:
    def __init__(self, parent=None):
        self.ch, self.name, self.parent, self.ufb = [], None, parent, None

    def leaves(self):
        return [self] if not self.ch else [l for c in self.ch for l in c.leaves()]

    def nodes(self):
        yield self
        for c in self.ch:
            yield from c.nodes()


def parse(txt):
    root = cur = Node()
    toks = LABEL.findall(txt.strip())
    i = 0
    while i < len(toks):
        t = toks[i].strip()
        if t == "(":
            n = Node(cur); cur.ch.append(n); cur = n
        elif t == ",":
            n = Node(cur.parent); cur.parent.ch.append(n); cur = n
        elif t == ")":
            cur = cur.parent
        elif t == ":":
            i += 1
        elif t == ";":
            break
        elif t:
            lab = t.strip("'")
            if cur.ch:
                try:
                    cur.ufb = float(lab.split("/")[-1])
                except ValueError:
                    pass
            else:
                cur.name = lab
        i += 1
    return root


def clades(root, outgroup):
    """Rooted clades implied by the outgroup on the unrooted tree: for every edge take the side
    without any outgroup sequence. Returns [(leaf_set, support)]; edges splitting the outgroup
    are skipped. Independent of how the tree file happens to be rooted."""
    allt = frozenset(l.name for l in root.leaves())
    out = []
    for n in root.nodes():
        if n is root:
            continue
        s = frozenset(l.name for l in n.leaves())
        c = allt - s
        if outgroup <= c:
            side = s
        elif outgroup <= s:
            side = c
        else:
            continue
        if len(side) > 1:
            out.append((side, n.ufb))
    return out


def largest_clade(cl, must, forbid, min_support):
    best = None
    for side, ufb in cl:
        if must <= side and not (side & forbid) and (ufb is None or ufb >= min_support):
            if best is None or len(side) > len(best[0]):
                best = (side, ufb)
    return best


def smallest_clade(cl, must, forbid, min_support):
    best = None
    for side, ufb in cl:
        if must <= side and not (side & forbid) and (ufb is None or ufb >= min_support):
            if best is None or len(side) < len(best[0]):
                best = (side, ufb)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--set", required=True, help="select_tree_set.py TSV")
    ap.add_argument("--min-support", type=float, default=95)
    ap.add_argument("--definition", choices=["broad", "strict"], default="broad",
                    help="broad: largest supported clade with all DBXN references and no Lepidoptera canonical "
                         "Vg (the DBXN lineage as a whole); strict: smallest such clade. All DBXN references "
                         "come from Noctuidae, so 'strict' stops at the Noctuidae subclade by construction.")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = {r["name"]: r for r in csv.DictReader(open(a.set), delimiter="\t")}
    root = parse(open(a.tree).read())
    names = {l.name for l in root.leaves()}
    missing = set(rows) - names
    role = {n: rows[n]["role"] for n in names if n in rows}
    og = {n for n, r in role.items() if r == "outgroup_Vg"}
    dbxn_ref = {n for n, r in role.items() if r == "DBXN_reference"}
    vg_ref = {n for n, r in role.items() if r == "lep_canonical_Vg"}
    cl = clades(root, og)
    og_mono = any(side == frozenset(names) - og for side, _ in cl) or len(og) <= 1

    out = []
    msg = [f"tips {len(names)}; sequences from the set missing in the tree (should be 0): {len(missing)}",
           f"outgroup Vg ({len(og)} seqs) monophyletic: {'yes' if og_mono else 'NO - check the tree'}"]
    d_strict = smallest_clade(cl, dbxn_ref, vg_ref, a.min_support)
    d_broad = largest_clade(cl, dbxn_ref, vg_ref, a.min_support)
    d = d_broad if a.definition == "broad" else d_strict
    v = smallest_clade(cl, vg_ref, set(), 0)  # Lepidoptera Vg clade; contains the DBXN clade
    if d is None:
        loose = smallest_clade(cl, dbxn_ref, set(), 0)
        bad = sorted(loose[0] & vg_ref) if loose else []
        msg.append("NO DBXN clade satisfies the rule. Smallest clade with all DBXN references also contains "
                   f"canonical Vg: {bad[:10]}")
    else:
        msg.append(f"DBXN clade ({a.definition}): {len(d[0])} tips, UFBoot {d[1]}")
        if d_strict and d_broad and d_strict[0] != d_broad[0]:
            msg.append(f"  strict clade {len(d_strict[0])} tips (UFBoot {d_strict[1]}) is nested in broad "
                       f"{len(d_broad[0])} tips (UFBoot {d_broad[1]}); column class_tree_strict marks the difference")
    if v is None:
        msg.append("NO clade holding all Lepidoptera canonical Vg (outgroup placement problem?)")
    else:
        inside = "contains" if d and d[0] <= v[0] else "does NOT contain"
        msg.append(f"Lepidoptera Vg clade: {len(v[0])} tips, UFBoot {v[1]}; {inside} the DBXN clade")

    by_sp = defaultdict(Counter)
    for n, r in rows.items():
        if r["role"] != "copy":
            continue
        if n not in names:
            cls = "not_in_tree"
        elif d and n in d[0]:
            cls = "DBXN"
        elif v and n in v[0]:
            cls = "Vg"
        else:
            cls = "unassigned"
        sister = "."
        if cls == "unassigned" and n in names:
            # smallest clade holding this copy and at least one reference: what references are in it
            near = min((side for side, _ in cl if n in side and side & (dbxn_ref | vg_ref)), key=len, default=frozenset())
            sister = ",".join(f"{k}:{c}" for k, c in sorted(Counter(role[x] for x in near if role.get(x) != "copy").items()))
        strict = "DBXN" if (d_strict and n in d_strict[0]) else ("Vg" if cls in ("DBXN", "Vg") else cls)
        out.append([r["species"], n, r["class_besthit"], cls, strict, r["intact"], sister])
        by_sp[r["species"]][cls] += 1

    with open(a.out + ".tsv", "w") as fh:
        fh.write("species\tname\tclass_besthit\tclass_tree\tclass_tree_strict\tintact\tsister_if_unassigned\n")
        for r in out:
            fh.write("\t".join(r) + "\n")
    agree = Counter((r[2], r[3]) for r in out)
    msg.append("best-hit class -> tree class: " + ", ".join(f"{k[0]}->{k[1]}: {c}" for k, c in sorted(agree.items())))
    with open(a.out + ".summary.txt", "w") as fh:
        fh.write("\n".join(msg) + "\n")
    print("\n".join(msg))


if __name__ == "__main__":
    main()
