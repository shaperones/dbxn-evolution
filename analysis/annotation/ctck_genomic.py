#!/usr/bin/env python3
"""Does each Vg/DBXN copy carry a C-terminal cystine knot (CTCK) exon? Checked in the genome.

miniprot proteins end where their template ends. DBXN templates lack CTCK, so a protein-level
domain scan (InterProScan) would report "no CTCK" for every DBXN copy by construction. Here the
genome itself is searched: CTCK queries (the C-terminal --tail aa of panel proteins that carry it:
Lepidoptera canonical Vg, MvirDBXN, spider Vg) are aligned with tblastn to a window from the copy
start to --downstream bp past its end, stopping before the next copy on the same strand. A hit on
the copy's strand inside the window = CTCK present.

Usage (per genome):
  ctck_genomic.py --copies SP.miniprot.copies.tsv --genome GENOME.fa --panel seed_panel_v2.faa \
      --out SP.ctck.tsv
Needs tblastn in PATH. Output: one row per DBXN / Vg / Vg_spider copy.
"""
import argparse
import csv
import os
import subprocess
import tempfile

from gff_utils import read_fasta, write_fasta

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LEP = {"Noctuidae", "Erebidae", "Geometridae", "Crambidae", "Pyralidae", "Plutellidae", "Saturniidae",
       "Tortricidae", "Psychidae", "Sphingidae", "Nymphalidae", "Papilionidae", "Bombycidae"}


def ctck_queries(panel, tail):
    q = []
    for header, s in read_fasta(panel):
        pid, pclass, taxon, group = (header.split("|") + [""] * 4)[:4]
        group = group.split()[0] if group else ""
        if (pclass == "Vg_positive" and group in LEP) or pclass in ("Vg_positive_spider",
                                                                      "DBXN_spider_reference_structure"):
            q.append((f"{pid}_Cterm", s.rstrip("*")[-tail:]))
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--tail", type=int, default=110, help="C-terminal length used as CTCK query (aa)")
    ap.add_argument("--downstream", type=int, default=10000)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t")
              if r["class"] in ("DBXN", "Vg", "Vg_spider")]
    if not copies:
        open(a.out, "w").write("species\tcopy_id\tclass\tintact\tctck\tbest_query\tevalue\tpident\n")
        return
    contigs = {r["contig"] for r in copies}
    genome = {n: s for n, s in read_fasta(a.genome) if n in contigs}
    allc = list(csv.DictReader(open(a.copies), delimiter="\t"))  # neighbours of any class bound windows

    tmp = tempfile.mkdtemp()
    qf = os.path.join(tmp, "q.faa")
    with open(qf, "w") as fh:
        for n, s in ctck_queries(a.panel, a.tail):
            write_fasta(fh, n, s)
    wf = os.path.join(tmp, "w.fa")
    with open(wf, "w") as fh:
        for r in copies:
            s, e, strand, ctg = int(r["start"]), int(r["end"]), r["strand"], r["contig"]
            same = [(int(x["start"]), int(x["end"])) for x in allc
                    if x["contig"] == ctg and x["strand"] == strand and x["copy_id"] != r["copy_id"]]
            if strand == "+":
                nxt = min([x[0] for x in same if x[0] > e], default=len(genome[ctg]) + 1)
                ws, we = s, min(e + a.downstream, nxt - 1)
            else:
                prv = max([x[1] for x in same if x[1] < s], default=0)
                ws, we = max(s - a.downstream, prv + 1), e
            seq = genome[ctg][ws - 1:we]
            if strand == "-":
                seq = seq.translate(COMP)[::-1]  # window in the gene's orientation
            write_fasta(fh, r["copy_id"], seq)
    res = subprocess.run(["tblastn", "-query", qf, "-subject", wf, "-evalue", str(a.evalue),
                          "-outfmt", "6 qseqid sseqid pident evalue sstart send", "-seg", "no"],
                         capture_output=True, text=True, check=True).stdout
    best = {}
    for line in res.splitlines():
        q, sid, pid, ev, ss, se = line.split("\t")
        if int(se) < int(ss):  # wrong strand relative to the copy
            continue
        if sid not in best or float(ev) < best[sid][1]:
            best[sid] = (q, float(ev), float(pid))
    with open(a.out, "w") as fh:
        fh.write("species\tcopy_id\tclass\tintact\tctck\tbest_query\tevalue\tpident\n")
        for r in copies:
            b = best.get(r["copy_id"])
            fh.write("\t".join([r["species"], r["copy_id"], r["class"], r["intact"], "yes" if b else "no",
                                b[0] if b else ".", f"{b[1]:.1e}" if b else ".", f"{b[2]:.1f}" if b else "."]) + "\n")
    n = sum(1 for r in copies if r["copy_id"] in best)
    print(f"{os.path.basename(a.copies)}: CTCK in {n}/{len(copies)} copies")


if __name__ == "__main__":
    main()
