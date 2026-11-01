#!/usr/bin/env python3
"""Inputs for halSynteny / HALPER on the revision Cactus alignment.

Locus sequences in cactus.hal are named <species>.<main|extra>.<contig>:<s>-<e> and hold genome[s-1:e]
(define_vg_loci.py), so a copy at genomic start..end (1-based, inclusive) sits at locus 0-based
[start - s, end - s + 1). Six species have more than one locus sequence; every copy keeps its own.

Writes to OUT:
  mapped_vgs.tsv        new_id (= GeneRax id species__copy), species, seq, genomic coords, locus coords,
                        classification (broad tree class), class_strict
  dbxns.txt             new_id of broad-DBXN copies
  species_order.tsv     HAL leaves in seqfile tree order (C4-b)
  adjacent_pairs.tsv    tree-neighbour pairs for halSynteny
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter

SEQ_RE = re.compile(r"^(?P<sp>[^.]+)\.(?P<kind>main|extra)\.(?P<ctg>.+):(?P<s>\d+)-(?P<e>\d+)$")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seqfile", required=True, help="cactus_input.txt (Newick line + species fasta lines)")
    ap.add_argument("--loci-dir", required=True, help="synteny/loci with <species>.locus.fasta")
    ap.add_argument("--copies-dir", required=True, help="main/miniprot with <species>.miniprot.copies.tsv")
    ap.add_argument("--classes", required=True, help="gene_tree/lep_classes_broad.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    lines = [l.strip() for l in open(a.seqfile) if l.strip()]
    newick = lines[0]
    hal_species = [l.split()[0] for l in lines[1:]]
    order = [t for t in re.findall(r"[(,]\s*([A-Za-z0-9_]+)\s*:", newick) if t in set(hal_species)]
    if sorted(order) != sorted(hal_species):
        sys.exit("seqfile tree tips and species lines differ")

    cls = {}
    for r in csv.DictReader(open(a.classes), delimiter="\t"):
        cls[r["name"]] = (r["class_tree"], r["class_tree_strict"])

    rows, missing = [], []
    for sp in order:
        seqs = []
        for l in open(os.path.join(a.loci_dir, f"{sp}.locus.fasta")):
            if l.startswith(">"):
                m = SEQ_RE.match(l[1:].strip())
                if not m:
                    sys.exit(f"unexpected locus header: {l.strip()}")
                seqs.append((l[1:].strip(), m["ctg"], int(m["s"]), int(m["e"])))
        for r in csv.DictReader(open(os.path.join(a.copies_dir, f"{sp}.miniprot.copies.tsv")), delimiter="\t"):
            c = cls.get(f"{sp}|{r['copy_id']}")
            if not c or c[0] not in ("DBXN", "Vg"):
                continue
            st, en = int(r["start"]), int(r["end"])
            hit = [x for x in seqs if x[1] == r["contig"] and x[2] <= st and en <= x[3]]
            if len(hit) != 1:
                missing.append(f"{sp}|{r['copy_id']}")
                continue
            name, _, s, _ = hit[0]
            rows.append({
                "new_id": f"{sp}__{r['copy_id']}", "species": sp, "seq": name,
                "chromosome": r["contig"], "genomic_start": st, "genomic_end": en, "strand": r["strand"],
                "locus_start": st - s, "locus_end": en - s + 1,
                "classification": c[0], "class_strict": c[1],
            })
    if missing:
        sys.exit(f"{len(missing)} classified copies not inside a HAL locus: {missing[:5]}")

    with open(os.path.join(a.out, "mapped_vgs.tsv"), "w") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(a.out, "dbxns.txt"), "w") as fh:
        fh.writelines(r["new_id"] + "\n" for r in rows if r["classification"] == "DBXN")
    with open(os.path.join(a.out, "species_order.tsv"), "w") as fh:
        fh.writelines(sp + "\n" for sp in order)
    with open(os.path.join(a.out, "adjacent_pairs.tsv"), "w") as fh:
        fh.writelines(f"{x}\t{y}\n" for x, y in zip(order, order[1:]))

    per_sp = Counter(r["species"] for r in rows)
    print(f"species {len(order)}; copies {len(rows)} ({Counter(r['classification'] for r in rows)}); "
          f"species with copies {len(per_sp)}; multi-sequence species "
          f"{sorted({r['species'] for r in rows if '.extra.' in r['seq']})}; adjacent pairs {len(order) - 1}")


if __name__ == "__main__":
    main()
