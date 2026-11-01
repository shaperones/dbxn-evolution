#!/usr/bin/env python3
"""Chromophore-contact residues (CCRs) of a reference DBXN mapped onto every sequence of an alignment.

Each CCR is identical, similar (same Clustal strong or weak group), different, gap (inside the aligned span) or missing
(outside it). Identity and similarity are % of the CCRs a sequence covers; a gap counts as a difference.
--label-col names the CCRs by another numbering (MvirDBXN: structure = transcript - 10).

Usage:
  cbr_matrix.py --alignment lep_set.linsi.aln --reference spMelPer_pVg1 \
      --residues ligand_binding_residues1.spMelPer_g4_clipkit_coords.tsv --coord-col coordinate --out lep_cbr
  cbr_matrix.py --alignment spider_set.linsi.aln --reference spMicVir_pVg1 --residues MvirDBXN_cbr.tsv \
      --coord-col real_coordinate_plus_10aa --label-col coordinate_from_doc --out spider_cbr_MvirDBXN
Outputs: <out>.matrix.tsv, <out>.scores.tsv
"""
import argparse
import csv
import sys

STRONG = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
WEAK = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"]


def read_fasta(path):
    name, seqs = None, {}
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def state(ref, res):
    if res == ".":
        return "missing"
    if res in "-X*":
        return "gap"
    if res == ref:
        return "identical"
    if any(ref in g and res in g for g in STRONG + WEAK):
        return "similar"
    return "different"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--reference", required=True, help="substring that identifies the reference sequence name")
    ap.add_argument("--residues", required=True)
    ap.add_argument("--coord-col", default="coordinate")
    ap.add_argument("--label-col", help="column with the residue numbers used in the labels (default: --coord-col)")
    ap.add_argument("--allow-mismatch", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    aln = read_fasta(a.alignment)
    refs = [k for k in aln if a.reference in k]
    if len(refs) != 1:
        sys.exit(f"reference '{a.reference}' matches {len(refs)} sequences: {refs}")
    ref = aln[refs[0]].upper()
    col_of, pos = {}, 0
    for i, c in enumerate(ref):
        if c != "-":
            pos += 1
            col_of[pos] = i
    rows = list(csv.DictReader(open(a.residues), delimiter="\t"))
    cbr = [(r["ligand_binding_residue"], r["type"], int(r[a.coord_col])) for r in rows]
    num = [int(r[a.label_col or a.coord_col]) for r in rows]
    bad = [f"{aa}{p} (reference has {ref[col_of[p]]})" for aa, _, p in cbr if ref[col_of[p]] != aa]
    if bad:
        print(f"{len(bad)} of {len(cbr)} CBRs do not match {refs[0]}: " + ", ".join(bad), file=sys.stderr)
        if not a.allow_mismatch:
            sys.exit("stop (use --allow-mismatch to keep the stated coordinates)")
    labels = [f"{ref[col_of[p]]}{k}" for (_, _, p), k in zip(cbr, num)]
    sets = {"all": [True] * len(cbr), "bilin": [t in ("bilin", "both") for _, t, _ in cbr],
            "lutein": [t in ("lutein", "both") for _, t, _ in cbr]}

    with open(a.out + ".matrix.tsv", "w") as fm, open(a.out + ".scores.tsv", "w") as fs:
        fm.write("#type\t" + "\t".join(t for _, t, _ in cbr) + "\n")
        fm.write("name\t" + "\t".join(labels) + "\n")
        fs.write("name\t" + "\t".join(f"{s}_{k}" for s in sets for k in ("n", "identical", "similar", "different", "gap", "missing",
                                                                         "identity_pct", "similarity_pct")) + "\n")
        for name, seq in aln.items():
            seq = seq.upper()
            span = [i for i, ch in enumerate(seq) if ch not in "-."]
            lo, hi = (span[0], span[-1]) if span else (1, 0)
            res = [seq[col_of[p]] if lo <= col_of[p] <= hi else "." for _, _, p in cbr]     # "." = missing (matrix.tsv)
            st = [state(ref[col_of[p]], x) for (_, _, p), x in zip(cbr, res)]
            fm.write(name + "\t" + "\t".join(res) + "\n")
            row = [name]
            for mask in sets.values():
                sub = [s for s, m in zip(st, mask) if m]
                n = len(sub)
                c = {k: sub.count(k) for k in ("identical", "similar", "different", "gap", "missing")}
                cov = n - c["missing"]
                pct = (lambda k: f"{100 * k / cov:.1f}") if cov else (lambda k: "NA")
                row += [n, c["identical"], c["similar"], c["different"], c["gap"], c["missing"], pct(c["identical"]),
                        pct(c["identical"] + c["similar"])]
            fs.write("\t".join(map(str, row)) + "\n")
    print(f"{len(aln)} sequences, {len(cbr)} CBRs -> {a.out}.matrix.tsv, {a.out}.scores.tsv")


if __name__ == "__main__":
    main()
