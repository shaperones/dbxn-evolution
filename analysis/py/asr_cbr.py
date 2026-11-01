#!/usr/bin/env python3
"""Ancestral states of the chromophore-binding residues (CBRs) at key nodes of the Vg/DBXN gene trees.

Question: are the CBRs that DBXN shares with Notodontidae Vg (but not with quadrifid Vg) ancestral for
noctuoid Vg (kept in DBXN, lost in quadrifid Vg after the duplication), or derived in DBXN?
Marginal reconstruction (IQ-TREE --ancestral) on a fixed topology; model as in the gene tree.

Runs (Lepidoptera on the alignment the trees were built from, columns after ClipKIT; spiders on the untrimmed MAFFT alignment):
  lep_ml      614 genes (GeneRax / AU set), ML tree of the full-set AU test (ml614.treefile)
  lep_constr  same genes, ML tree under the constraint "DBXN inside quadrifid Vg" (constr614.treefile)
  spider      75 sequences of tree v4 (spider.contree; 27.10), untrimmed MAFFT alignment (31.10: ClipKIT kpic-gappy dropped
              2 of the 57 CBRs of MvirDBXN, C652 and L1053 - nearly invariant columns with singleton changes; the run on the
              ClipKIT alignment is kept in analysis/data/asr/spider_clip/). Lepidoptera: all 37 CBRs are in the ClipKIT alignment
Nodes (MRCA after rooting: Lepidoptera — on the non-noctuoid outgroup Vg (Pyraloidea, Tineidae);
spiders — midpoint, as in the CBR matrices of Fig. 7):
  Lepidoptera: DBXN; Vg_quadrifid; Vg_Notodontidae; Vg_Noctuoidea (quadrifid + Notodontidae Vg; in the
  constrained tree it contains DBXN); DBXN_split (MRCA of DBXN and quadrifid Vg = the duplication);
  Vg_Noctuoidea_Geometridae
  Spiders: MvirDBXN lineage (MvirDBXN + Neosparassus ortholog); its parent (split from other Sparassidae Vg);
  clade A
Outputs in analysis/data/asr/<run>/: IQ-TREE files, <run>.cbr_nodes.tsv (per node x CBR: MAP residue, its
posterior, posterior of the reference residue); summary printed and written to <run>.summary.tsv.

Usage (repository root):
  asr_cbr.py prepare|run|summarise lep_ml|lep_constr|spider [--iqtree PATH] [--threads 8]
  (`run` = prepare + IQ-TREE + summarise)
"""
import argparse
import csv
import os
import subprocess

import dendropy

from cbr_matrix import read_fasta

G = "analysis/data/gene_tree"
RUNS = {
    "lep_ml": dict(aln=f"{G}/lep/lep_set.clip", linsi=f"{G}/lep/lep_set.linsi.aln", tree=f"{G}/lep_au/full614/ml614.treefile",
                   ref="REF|DBXN_reference|spMelPer_pVg1", pos="analysis/data/cbr/positions/MperDBXN_cbr.tsv", col="coordinate",
                   model="JTT+F+I+R10", rename=lambda n: n.replace("|", "__", 1)),
    "lep_constr": dict(aln=f"{G}/lep/lep_set.clip", linsi=f"{G}/lep/lep_set.linsi.aln", tree=f"{G}/lep_au/full614/constr614.treefile",
                       ref="REF|DBXN_reference|spMelPer_pVg1", pos="analysis/data/cbr/positions/MperDBXN_cbr.tsv", col="coordinate",
                       model="JTT+F+I+R10", rename=lambda n: n.replace("|", "__", 1)),
    # spiders: tree v4 (27.10; v3 without the apolipophorin-like clade C; run on v3 11.10 kept in analysis/data/asr/spider_v3/)
    "spider": dict(aln=f"{G}/spider_v4/spider_set.linsi.aln", linsi=f"{G}/spider_v4/spider_set.linsi.aln", tree=f"{G}/spider_v4/spider.contree",
                   ref="REF|spider_reference|spMicVir_pVg1", pos="analysis/data/cbr/positions/MvirDBXN_cbr.tsv",
                   col="real_coordinate_plus_10aa", model="JTT+F+I+R4", rename=lambda n: n),
}
QUAD = {"Erebidae", "Nolidae", "Noctuidae"}
EDGES = {  # ancestor -> descendant comparisons reported in <run>.changes.tsv
    "lep": [("Vg_Noctuoidea_Geometridae", "Vg_Noctuoidea"), ("Vg_Noctuoidea", "Vg_Notodontidae"),
            ("DBXN_split", "DBXN"), ("DBXN_split", "Vg_quadrifid")],
    "spider": [("clade_A", "MvirDBXN_parent"), ("MvirDBXN_parent", "MvirDBXN_lineage")],
}
AA = "ARNDCQEGHILKMFPSTWYV"


def out_dir(run):
    return f"analysis/data/asr/{run}"


def cbr_columns(R):
    """CBR -> 1-based column of the ClipKIT alignment (via the reference in the MAFFT alignment)."""
    L, C = read_fasta(R["linsi"]), read_fasta(R["aln"])
    names = list(L)
    lcols, ccols = list(zip(*[L[n] for n in names])), list(zip(*[C[n] for n in names]))
    l2c, j = {}, 0
    for i, c in enumerate(lcols):
        if j < len(ccols) and c == ccols[j]:
            l2c[i], j = j, j + 1
    assert j == len(ccols), "ClipKIT columns are not an ordered subset of the MAFFT columns"
    p2l, k = {}, 0
    for i, ch in enumerate(L[R["ref"]]):
        if ch != "-":
            k += 1
            p2l[k] = i
    out = []
    for r in csv.DictReader(open(R["pos"]), delimiter="\t"):
        p = int(r[R["col"]])
        out.append((f"{r['ligand_binding_residue']}{p}", r["type"], l2c[p2l[p]] + 1 if p2l[p] in l2c else None))
    return out


def prepare(run):
    R, d = RUNS[run], out_dir(run)
    os.makedirs(d, exist_ok=True)
    tips = {l.taxon.label for l in dendropy.Tree.get(path=R["tree"], schema="newick", preserve_underscores=True).leaf_node_iter()}
    C = read_fasta(R["aln"])
    seqs = {R["rename"](n): s for n, s in C.items() if R["rename"](n) in tips}
    assert len(seqs) == len(tips), f"{len(tips) - len(seqs)} tree tips not in the alignment"
    with open(f"{d}/{run}.aln.fasta", "w") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in seqs.items())
    with open(f"{d}/{run}.cbr_columns.tsv", "w") as fh:
        fh.write("cbr\ttype\tclip_column\n")
        fh.writelines(f"{c}\t{t}\t{x if x else 'NA'}\n" for c, t, x in cbr_columns(R))
    print(f"{run}: {len(seqs)} sequences -> {d}/")


def run_iqtree(run, iqtree, threads):
    R, d = RUNS[run], out_dir(run)
    subprocess.run([iqtree, "-s", f"{d}/{run}.aln.fasta", "-te", R["tree"], "-m", R["model"], "--ancestral",
                    "-T", str(threads), "--prefix", f"{d}/{run}", "-redo", "--quiet"], check=True)


def node_sets(run, tree):
    """Named nodes -> MRCA, on the tree rooted as described in the module docstring."""
    labs = [l.taxon.label for l in tree.leaf_node_iter()]
    if run.startswith("lep"):
        meta = {r["name"].replace("|", "__", 1): r for r in csv.DictReader(open("analysis/data/cbr/cbr.lep_per_copy.tsv"), delimiter="\t")}
        fam = {n: meta[n]["family"] for n in labs}
        cls = {n: meta[n]["class_tree"] for n in labs}
        og = [n for n in labs if cls[n] == "Vg" and fam[n] not in QUAD | {"Notodontidae", "Geometridae"}]
        dbxn = [n for n in labs if cls[n] == "DBXN"]
        qvg = [n for n in labs if cls[n] == "Vg" and fam[n] in QUAD]
        nvg = [n for n in labs if cls[n] == "Vg" and fam[n] == "Notodontidae"]
        gvg = [n for n in labs if cls[n] == "Vg" and fam[n] == "Geometridae"]
        root_on = og
        sets = {"DBXN": dbxn, "Vg_quadrifid": qvg, "Vg_Notodontidae": nvg, "Vg_Noctuoidea": qvg + nvg,
                "DBXN_split": dbxn + qvg, "Vg_Noctuoidea_Geometridae": qvg + nvg + gvg + dbxn}
        return root_on, sets
    mv = ["REF|spider_reference|spMicVir_pVg1", "spMicVir|spMicVir_asBTCV1_sq000533t1.p1", "spNeoDia|spNeoDia_asBTCV1_sq1933ct1.p4"]
    clade = {r["name"]: r["clade"] for r in csv.DictReader(open("analysis/data/cbr/cbr.spider_per_copy.tsv"), delimiter="\t")}
    return None, {"MvirDBXN_lineage": mv, "MvirDBXN_parent": mv, "clade_A": [n for n in labs if clade.get(n) == "A"]}


def summarise(run):
    R, d = RUNS[run], out_dir(run)
    cols = list(csv.DictReader(open(f"{d}/{run}.cbr_columns.tsv"), delimiter="\t"))
    tree = dendropy.Tree.get(path=f"{d}/{run}.treefile", schema="newick", preserve_underscores=True, rooting="force-unrooted")
    root_on, sets = node_sets(run, tree)
    if root_on:
        labs = [l.taxon.label for l in tree.leaf_node_iter()]
        inner = next(x for x in labs if x not in root_on)
        tree.reroot_at_edge(tree.find_node_with_taxon_label(inner).edge, update_bipartitions=True)
        tree.is_rooted = True
        tree.reroot_at_edge(tree.mrca(taxon_labels=root_on).edge, update_bipartitions=True)
    else:
        tree.reroot_at_midpoint(update_bipartitions=True)
    tree.is_rooted = True
    node_of = {}
    for name, tips in sets.items():
        n = tree.mrca(taxon_labels=tips)
        if name.endswith("_parent"):
            n = n.parent_node
        node_of[name] = (n.label or "").split("/")[0]       # IQ-TREE labels internal nodes NodeN (or NodeN/support)
        node_of[name + "_size"] = len(n.leaf_nodes())
    want = {v: k for k, v in node_of.items() if not k.endswith("_size")}
    site2cbr = {int(c["clip_column"]): c for c in cols if c["clip_column"] != "NA"}
    state = {}
    with open(f"{d}/{run}.state") as fh:
        rows = (l for l in fh if not l.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            if r["Node"] in want and int(r["Site"]) in site2cbr:
                state[(want[r["Node"]], int(r["Site"]))] = r
    out = []
    for name in [k for k in node_of if not k.endswith("_size")]:
        for site, c in site2cbr.items():
            r = state[(name, site)]
            ref_aa = c["cbr"][0]
            out.append(dict(node=name, n_tips=node_of[name + "_size"], cbr=c["cbr"], type=c["type"], map=r["State"],
                            pp_map=f"{float(r['p_' + r['State']]) if r['State'] in AA else 0:.2f}",
                            pp_ref=f"{float(r['p_' + ref_aa]):.2f}"))
    with open(f"{d}/{run}.cbr_nodes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    names = [k for k in node_of if not k.endswith("_size")]
    by_node = {(r["node"], r["cbr"]): r for r in out}
    with open(f"{d}/{run}.cbr_nodes_wide.tsv", "w") as fh:      # CBR x node: "MAP (posterior)"
        fh.write("cbr\ttype\t" + "\t".join(f"{n}({node_of[n + '_size']})" for n in names) + "\n")
        for site, c in site2cbr.items():
            fh.write(f"{c['cbr']}\t{c['type']}\t" + "\t".join(
                f"{by_node[(n, c['cbr'])]['map']} ({by_node[(n, c['cbr'])]['pp_map']})" for n in names) + "\n")
    summ = []
    for name in names:
        rs = [r for r in out if r["node"] == name]
        summ.append(dict(node=name, n_tips=rs[0]["n_tips"], n_cbr=len(rs),
                         map_is_ref=sum(r["map"] == r["cbr"][0] for r in rs),
                         map_is_ref_pp95=sum(r["map"] == r["cbr"][0] and float(r["pp_map"]) >= 0.95 for r in rs),
                         uncertain_pp_lt_0_8=sum(float(r["pp_map"]) < 0.8 for r in rs)))
    with open(f"{d}/{run}.summary.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summ)
    for s in summ:
        print("\t".join(map(str, s.values())))

    # CBR changes along the paths of interest: ancestor -> descendant, reference residue gained / lost
    edges = EDGES["lep" if run.startswith("lep") else "spider"]
    by = {(r["node"], r["cbr"]): r for r in out}
    ch = []
    for anc, desc in edges:
        for c in [r["cbr"] for r in out if r["node"] == anc]:
            a_, d_ = by[(anc, c)], by[(desc, c)]
            if a_["map"] != d_["map"]:
                ref_aa = c[0]
                kind = "gain_ref" if d_["map"] == ref_aa else ("loss_ref" if a_["map"] == ref_aa else "other")
                ch.append(dict(ancestor=anc, descendant=desc, cbr=c, type=a_["type"], from_=f"{a_['map']} ({a_['pp_map']})",
                               to=f"{d_['map']} ({d_['pp_map']})", change=kind))
    with open(f"{d}/{run}.changes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ancestor", "descendant", "cbr", "type", "from_", "to", "change"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(ch)
    for anc, desc in edges:
        sub = [x for x in ch if x["ancestor"] == anc and x["descendant"] == desc]
        print(f"{anc} -> {desc}: " + ", ".join(f"{x['cbr']} {x['from_'].split()[0]}>{x['to'].split()[0]} [{x['change']}]" for x in sub))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "run", "summarise"])
    ap.add_argument("run", choices=RUNS)
    ap.add_argument("--iqtree", default=os.path.expanduser("~/micromamba/envs/vg-annot/bin/iqtree3"))
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    if a.step in ("prepare", "run"):
        prepare(a.run)
    if a.step == "run":
        run_iqtree(a.run, a.iqtree, a.threads)
    if a.step in ("run", "summarise"):
        summarise(a.run)


if __name__ == "__main__":
    main()
