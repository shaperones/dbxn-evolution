#!/usr/bin/env python3
"""CBR similarity by gene-tree class: per-copy tables, group summaries and species-level CBR predictors.

Inputs are the outputs of cbr_matrix.py on the alignments before ClipKIT (lep_set.linsi.aln vs MperDBXN,
spider_set.linsi.aln vs MvirDBXN). Classes come from the gene trees only (CBRs are not used to define them):
Lepidoptera — lep_classes_broad.tsv (DBXN / Vg, broad and strict DBXN clade); spiders — Vg clades A / B of tree v4
(27.10: v3 without its clade C, which groups with insect apolipophorins in a tree with the negative references -
analysis/data/gene_tree/spider_v3_with_negatives/; A holds PpVg3 and MvirDBXN, B PpVg1, PpVg2, PtVg4, PtVg6).

Species-level predictors (Lepidoptera, intact copies only), written to <out>.lep_species.tsv:
  cbr_sum_DBXN     sum of CBR similarity (fraction) over DBXN copies — class from the tree, score from CBRs
  n_cbr_high       copies of ANY class with CBR similarity >= --threshold (70 %) — no classification at all;
                   n_cbr_65, n_cbr_70, n_cbr_75 — the same for each threshold (sensitivity)
  cbr_max          highest CBR similarity among all intact copies

Usage (from the repository root):
  cbr_summary.py --out analysis/data/cbr/cbr
"""
import argparse
import csv
import statistics as st
from collections import defaultdict

from generax_summary import parse

P = dict(
    lep_scores="analysis/data/cbr/lep_cbr_MperDBXN.scores.tsv",
    lep_classes="analysis/data/gene_tree/lep/lep_classes_broad.tsv",
    species="revision/data/lep_species_copies_2026-10-10.tsv",
    sp_scores="analysis/data/cbr/spider_cbr_MvirDBXN.scores.tsv",
    sp_tree="analysis/data/gene_tree/spider_v4/spider.contree",   # 27.10: v4 = v3 without clade C (apolipophorin-like)
)
OUTGROUPS = {"spParPse", "spXysAce", "spAmauFer", "spEraAtri"}
MVIR = {"REF|spider_reference|spMicVir_pVg1": "MvirDBXN (structure)",
        "spMicVir|spMicVir_asBTCV1_sq000533t1.p1": "MvirDBXN (transcript)",
        "spNeoDia|spNeoDia_asBTCV1_sq1933ct1.p4": "Neosparassus ortholog of MvirDBXN"}
SETS = ("all", "bilin", "lutein")


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


CLADE_ANCHORS = {  # clade = smallest clade (either side of a branch; the tree is unrooted) holding all anchors:
    # the outgroup-genome copies of the clade + one panel or Sparassidae sequence; all three are UFBoot 100
    "A": ["spParPse|spParPse_mp1", "spXysAce|spXysAce_mp4", "spAmauFer|spAmauFer_mp4", "spEraAtri|spEraAtri_mp5",
          "REF|spider_reference|XP_087020500.1", "REF|spider_reference|spMicVir_pVg1"],
    "B": ["spParPse|spParPse_mp2", "spParPse|spParPse_mp3", "spXysAce|spXysAce_mp1", "spXysAce|spXysAce_mp2",
          "spAmauFer|spAmauFer_mp2", "spAmauFer|spAmauFer_mp3", "spEraAtri|spEraAtri_mp2", "spEraAtri|spEraAtri_mp3",
          "REF|spider_reference|AXN69712.1", "REF|spider_reference|AXN69713.1", "REF|spider_reference|XP_015921257.2",
          "REF|spider_reference|XP_015930207.2", "REF|spider_reference|XP_015930209.2"],
}
CLADE_SIZE = {"A": 41, "B": 31}  # tree v4 (as in v3)


def spider_clades(path):
    """Clades A / B of the spider Vg tree from CLADE_ANCHORS; stops if a size differs from CLADE_SIZE.
    A rule such as "the smallest clade with all four outgroups" fails here: in clade C a 5-tip subclade
    already holds all four, and in clade B every outgroup genome has two copies."""
    tree = parse(open(path).read())
    tips = {l.name for l in tree.leaves()}
    sides = []
    for n in tree.nodes():
        if n.ch:
            s = frozenset(l.name for l in n.leaves())
            sides += [s, frozenset(tips - s)]
    clade = {}
    for c, anchors in CLADE_ANCHORS.items():
        best = min((s for s in sides if set(anchors) <= s), key=len)
        if len(best) != CLADE_SIZE[c]:
            raise SystemExit(f"clade {c}: {len(best)} tips, expected {CLADE_SIZE[c]}")
        clade.update({x: c for x in best})
    return clade


def summarise(rows, keys):
    groups = defaultdict(list)
    for r in rows:
        groups[tuple(r[k] for k in keys)].append(r)
    out = []
    for g, rs in sorted(groups.items()):
        rec = dict(zip(keys, g), n=len(rs))
        for s in SETS:
            v = [float(r[f"{s}_similarity_pct"]) for r in rs]
            rec.update({f"{s}_median": f"{st.median(v):.1f}", f"{s}_min": f"{min(v):.1f}", f"{s}_max": f"{max(v):.1f}"})
        # 30.10: similarity is already over the CBRs a sequence covers (cbr_matrix.py: missing ends excluded, gaps inside
        # the span = differences); the former all_nongap_median (every gap excluded) is dropped
        rec.update(all_gap_median=st.median(int(r["all_gap"]) for r in rs),
                   all_missing_median=st.median(int(r["all_missing"]) for r in rs))
        out.append(rec)
    return out


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threshold", type=float, default=70.0,
                    help="CBR similarity, %%, for n_cbr_high (70: no intact quadrifid Vg reaches it; decided 11.10)")
    ap.add_argument("--alt-thresholds", type=float, nargs="*", default=[65.0, 75.0], help="sensitivity variants")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # Lepidoptera
    cls = {r["name"]: r for r in tsv(P["lep_classes"])}
    sp = {r["tree_label"]: r for r in tsv(P["species"])}
    lep = []
    for r in tsv(P["lep_scores"]):
        c = cls.get(r["name"])
        if c is None:
            group, fam, intact, species = r["name"].split("|")[1], "reference", ".", "."
        else:
            group, species, intact = c["class_tree"], c["species"], c["intact"]
            fam = sp[species]["family"] if species in sp else "?"
        lep.append(dict(name=r["name"], species=species, family=fam, class_tree=group, intact=intact,
                        **{k: r[k] for k in r if k != "name"}))
    write(a.out + ".lep_per_copy.tsv", lep)
    copies = [r for r in lep if r["family"] != "reference"]
    summ = summarise(copies, ["class_tree", "family"]) + summarise(copies, ["class_tree"]) \
        + summarise([r for r in lep if r["family"] == "reference"], ["class_tree"])
    for s in summ:
        s.setdefault("family", "all" if s["class_tree"] in ("DBXN", "Vg") else "panel")
    write(a.out + ".lep_summary.tsv", [{k: s[k] for k in ["class_tree", "family"] + [k for k in s if k not in ("class_tree", "family")]} for s in summ])

    thr = [a.threshold] + [t for t in a.alt_thresholds if t != a.threshold]
    per_sp = defaultdict(lambda: dict(cbr_sum_DBXN=0.0, cbr_max=0.0, **{f"n_cbr_{t:g}": 0 for t in thr}))
    for r in copies:
        if r["intact"] != "1":
            continue
        v = float(r["all_similarity_pct"])
        d = per_sp[r["species"]]
        if r["class_tree"] == "DBXN":
            d["cbr_sum_DBXN"] += v / 100
        for t in thr:
            d[f"n_cbr_{t:g}"] += v >= t
        d["cbr_max"] = max(d["cbr_max"], v)
    write(a.out + ".lep_species.tsv", [dict(tree_label=s, cbr_sum_DBXN=f"{per_sp[s]['cbr_sum_DBXN']:.3f}",
                                            n_cbr_high=per_sp[s][f"n_cbr_{a.threshold:g}"],
                                            **{f"n_cbr_{t:g}": per_sp[s][f"n_cbr_{t:g}"] for t in thr},
                                            cbr_max=f"{per_sp[s]['cbr_max']:.1f}") for s in sorted(sp)])

    # spiders
    clade = spider_clades(P["sp_tree"])
    spd = []
    for r in tsv(P["sp_scores"]):
        spd.append(dict(name=r["name"], species=r["name"].split("|")[0], clade=clade.get(r["name"], "."),
                        note=MVIR.get(r["name"], "reference" if r["name"].startswith("REF|") else "."),
                        **{k: r[k] for k in r if k != "name"}))
    write(a.out + ".spider_per_copy.tsv", spd)
    write(a.out + ".spider_summary.tsv", summarise([r for r in spd if r["note"] == "."], ["clade"]))

    dbxn = [float(r["all_similarity_pct"]) for r in copies if r["class_tree"] == "DBXN"]
    vg = [float(r["all_similarity_pct"]) for r in copies if r["class_tree"] == "Vg"]
    print(f"Lepidoptera vs MperDBXN: DBXN n={len(dbxn)} median {st.median(dbxn):.1f} ({min(dbxn):.1f}-{max(dbxn):.1f}); "
          f"Vg n={len(vg)} median {st.median(vg):.1f} ({min(vg):.1f}-{max(vg):.1f}); "
          f"DBXN below {a.threshold}: {sum(x < a.threshold for x in dbxn)}, Vg at/above: {sum(x >= a.threshold for x in vg)}")


if __name__ == "__main__":
    main()
