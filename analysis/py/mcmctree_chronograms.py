#!/usr/bin/env python3
"""Sample chronograms from an MCMCTree posterior (mcmc.txt) on the fixed topology of FigTree.tre.

MCMCTree stores one row per MCMC sample with the ages of all internal nodes (t_n<k>) but no trees.
Each t_n column is matched to a FigTree.tre node by its posterior mean (FigTree branch lengths are
posterior mean times); the match must be one-to-one within --tol, otherwise the script stops. Then
--n rows are drawn at random, ages are put on the nodes and branch length = parent age - child age.

Usage: mcmctree_chronograms.py --dir analysis/data/mcmctree/araneae/A \
           --mapping analysis/data/mcmctree/araneae/mapping_species.tsv --n 100 --seed 1
Writes <dir>/chronograms_<n>.nwk (tips = tree_label codes, branch lengths in Ma) and
<dir>/chronogram_mean.nwk (posterior mean ages).
"""
import argparse
import csv
import random
import sys

from mcmctree_nodes import parse_figtree, height


def to_newick(n, ages, unit):
    if not n.ch:
        return n.name
    return "(" + ",".join(f"{to_newick(c, ages, unit)}:{(ages[id(n)] - ages.get(id(c), 0.0)) * unit:.6f}"
                          for c in n.ch) + ")"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--unit", type=float, default=100.0, help="Ma per MCMCTree time unit")
    ap.add_argument("--tol", type=float, default=0.002, help="max |mean(t_n) - FigTree age| in MCMCTree units")
    a = ap.parse_args()

    mp = {r["new_name"]: r["tree_label"] for r in csv.DictReader(open(a.mapping), delimiter="\t")}
    tree = parse_figtree(open(f"{a.dir}/FigTree.tre").read())
    for l in tree.leaves():
        l.name = mp[l.name]
    internal = [n for n in tree.nodes() if n.ch]

    with open(f"{a.dir}/mcmc.txt") as fh:
        header = fh.readline().split()
        tcols = [i for i, c in enumerate(header) if c.startswith("t_n")]
        rows = [[float(v) for v in (line.split()[i] for i in tcols)] for line in fh if line.strip()]
    means = [sum(r[j] for r in rows) / len(rows) for j in range(len(tcols))]

    if len(tcols) != len(internal):
        sys.exit(f"{len(internal)} internal nodes but {len(tcols)} t_n columns")
    # MCMCTree numbers internal nodes in pre-order of the input tree (root = t_n<ntips+1>), and FigTree.tre keeps
    # that order; checked against the posterior means. Fallback for other layouts: nearest posterior mean.
    if all(abs(height(n) - m) <= a.tol for n, m in zip(internal, means)):
        match = dict(enumerate(internal))
    else:
        match, used = {}, set()
        for j, m in enumerate(means):
            best = min(internal, key=lambda n: abs(height(n) - m))
            if abs(height(best) - m) > a.tol or id(best) in used:
                sys.exit(f"cannot match {header[tcols[j]]} (mean {m:.5f}) uniquely to a FigTree node")
            match[j] = best
            used.add(id(best))
    worst = max(abs(height(match[j]) - means[j]) for j in match)

    random.seed(a.seed)
    pick = random.sample(range(len(rows)), a.n)
    with open(f"{a.dir}/chronograms_{a.n}.nwk", "w") as fh:
        for k in pick:
            ages = {id(match[j]): rows[k][j] for j in match}
            fh.write(to_newick(tree, ages, a.unit) + ";\n")
    ages = {id(match[j]): means[j] for j in match}
    with open(f"{a.dir}/chronogram_mean.nwk", "w") as fh:
        fh.write(to_newick(tree, ages, a.unit) + ";\n")
    print(f"{a.dir}: {len(rows)} samples, {len(tcols)} nodes matched (max |mean - FigTree| {worst:.2g}); "
          f"wrote {a.n} chronograms (root age mean {max(means) * a.unit:.1f} Ma)")


if __name__ == "__main__":
    main()
