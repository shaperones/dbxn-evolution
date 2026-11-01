#!/usr/bin/env python3
"""Vg loci for synteny (Cactus): BUSCO anchors extended to every Vg/DBXN copy on the same contig.

The v1 locus was the interval between two M. persicariae BUSCO anchors (9200at7088, 839at7088).
With the new annotation 10 of 502 tree-classified copies in 6 species lie outside that interval
on the same chromosome (whole Vg+DBXN arrays shifted by 1-4 Mb), and one DBXN lies on another
contig. Here the locus = union(anchor interval, all copies on the anchor contig) +- --flank.
Copies on other contigs get their own small locus (copy span +- --flank) and a flag.

Anchors: BUSCO full_table.tsv per species (Complete or Duplicated; with Duplicated, the copy on
the contig with most Vg copies is used). If an anchor is missing, the locus is built from the
copies alone and flagged 'no_anchor' (the v1 workflow rescued 3 such genomes with nested BUSCOs).

Usage:
  define_vg_loci.py --samples samples.tsv --busco-pattern '/path/{species}/run_lepidoptera_odb10/full_table.tsv' \
      --copies-dir main/miniprot --classes lep_classes_broad.tsv --out vg_loci [--extract]
--extract writes <out>/<species>.locus.fasta (genome paths from samples.tsv) for Cactus.
Outputs: <out>/vg_loci.bed, <out>/vg_loci.tsv
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

from gff_utils import read_fasta, write_fasta

ANCHORS = ("9200at7088", "839at7088")


def read_busco(path):
    hits = defaultdict(list)
    if not os.path.exists(path):
        return hits
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and f[0] in ANCHORS and f[1] in ("Complete", "Duplicated"):
            hits[f[0]].append((f[2], min(int(f[3]), int(f[4])), max(int(f[3]), int(f[4]))))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--busco-pattern", required=True, help="path with {species} placeholder")
    ap.add_argument("--copies-dir", required=True)
    ap.add_argument("--classes", required=True, help="lep_classes_broad.tsv (name = species|copy_id)")
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--extract", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    cls = {r["name"].split("|")[-1]: r["class_tree"] for r in csv.DictReader(open(a.classes), delimiter="\t")}
    genomes = {r["tree_label"]: r["input"] for r in csv.DictReader(open(a.samples), delimiter="\t")
               if r["mode"] == "genome" and r["set"] == "main"}
    bed, rows = [], []
    for sp, gpath in sorted(genomes.items()):
        tsv = os.path.join(a.copies_dir, f"{sp}.miniprot.copies.tsv")
        if not os.path.exists(tsv):
            continue
        copies = [r for r in csv.DictReader(open(tsv), delimiter="\t") if cls.get(r["copy_id"]) in ("DBXN", "Vg")]
        if not copies:
            continue
        by_ctg = defaultdict(list)
        for r in copies:
            by_ctg[r["contig"]].append((int(r["start"]), int(r["end"])))
        main_ctg = Counter(r["contig"] for r in copies).most_common(1)[0][0]
        busco = read_busco(a.busco_pattern.format(species=sp))
        anchor = {}
        for aid in ANCHORS:
            on_main = [h for h in busco.get(aid, []) if h[0] == main_ctg]
            if on_main:
                anchor[aid] = on_main[0]
        flags = []
        spans = by_ctg[main_ctg][:]
        if len(anchor) == 2:
            spans += [(anchor[x][1], anchor[x][2]) for x in ANCHORS]
            lo_a = min(anchor[x][1] for x in ANCHORS); hi_a = max(anchor[x][2] for x in ANCHORS)
            if any(s < lo_a or e > hi_a for s, e in by_ctg[main_ctg]):
                flags.append("copies_outside_anchors")
        else:
            flags.append("no_anchor" if not anchor else "one_anchor")
            spans += [(anchor[x][1], anchor[x][2]) for x in anchor]
        lo = max(1, min(s for s, _ in spans) - a.flank)
        hi = max(e for _, e in spans) + a.flank
        loci = [(main_ctg, lo, hi, "main")]
        for ctg, ss in by_ctg.items():
            if ctg != main_ctg:
                flags.append(f"copies_on_other_contig:{ctg}")
                loci.append((ctg, max(1, min(s for s, _ in ss) - a.flank), max(e for _, e in ss) + a.flank, "extra"))
        for ctg, s, e, kind in loci:
            bed.append((ctg, s - 1, e, f"{sp}|{kind}"))
        rows.append([sp, main_ctg, lo, hi, hi - lo + 1, len(copies), ",".join(flags) or "ok", gpath])

        if a.extract:
            want = {c for c, *_ in loci}
            with open(os.path.join(a.out, f"{sp}.locus.fasta"), "w") as fh:
                for name, seq in read_fasta(gpath):
                    if name in want:
                        for ctg, s, e, kind in loci:
                            if ctg == name:
                                write_fasta(fh, f"{sp}.{kind}.{ctg}:{s}-{min(e, len(seq))}", seq[s - 1:e])

    with open(os.path.join(a.out, "vg_loci.bed"), "w") as fh:
        for b in bed:
            fh.write("\t".join(map(str, b)) + "\n")
    with open(os.path.join(a.out, "vg_loci.tsv"), "w") as fh:
        fh.write("species\tcontig\tstart\tend\tlength\tn_copies\tflags\tgenome\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    f = Counter(x for r in rows for x in r[6].split(",") if not x.startswith("copies_on_other_contig"))
    print(f"{len(rows)} species; flags: {dict(f)}; extra contigs: "
          f"{sum(1 for r in rows if 'other_contig' in r[6])} species")


if __name__ == "__main__":
    main()
