#!/usr/bin/env python3
"""Vg-family copies of one genome from miniprot alignments.

  miniprot -t16 -I -j2 --outs=0.4 -N 50 --outc=0.2 --gff --trans genome.fa templates.faa > SPECIES.mp.gff

Alignments are split at long introns and clustered into copies by CDS overlap (best score first); partial copies with
complementary template ranges are joined. A copy is classed by its best template and is intact if one alignment covers
>= --intact-cov of its template without frameshifts or stop codons.
Outputs: <prefix>.copies.tsv, .copies.faa, .counts.tsv, .copies.gff3
"""
import argparse
import csv
import re
from collections import defaultdict

from gff_utils import parse_attrs, read_fasta, write_fasta


class Aln:
    __slots__ = ("id", "query", "qlen", "qs", "qe", "strand", "contig", "exons", "score",
                 "fs", "st", "prot", "cls", "note")

    def cds_len(self):
        return sum(e - s + 1 for s, e, _, _ in self.exons)

    def cov(self):
        return (self.qe - self.qs) / self.qlen

    def span(self):
        return min(x[0] for x in self.exons), max(x[1] for x in self.exons)


def parse_miniprot_gff(path, classes):
    alns, cur, paf, sta = {}, None, None, None
    for line in open(path):
        if line.startswith("##PAF"):
            paf, sta = line.rstrip("\n").split("\t")[1:], None
            continue
        if line.startswith("##STA"):
            sta = line.rstrip("\n").split("\t")[1]
            continue
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        a = parse_attrs(f[8])
        if f[2] == "mRNA":
            tags = dict(t.split(":", 2)[0::2] for t in paf[12:] if t.count(":") >= 2)
            x = Aln()
            x.id, x.query, x.qlen = a["ID"], paf[0], int(paf[1])
            x.qs, x.qe = int(paf[2]), int(paf[3])
            x.strand, x.contig = f[6], f[0]
            x.score = int(tags.get("AS", f[5] if f[5] != "." else 0))
            x.fs, x.st = int(tags.get("fs", 0)), int(tags.get("st", 0))
            x.prot, x.exons, x.note = sta or "", [], ""
            x.cls = classes.get(x.query, "unknown")
            alns[x.id] = cur = x
        elif f[2] == "CDS" and cur is not None and a.get("Parent") == cur.id:
            t = a.get("Target", "").split()
            qs, qe = (int(t[1]) - 1, int(t[2])) if len(t) == 3 else (None, None)
            cur.exons.append((int(f[3]), int(f[4]), qs, qe))
    return list(alns.values())


def split_long_introns(x, max_intron):
    ex = sorted(x.exons)
    groups, g = [], [ex[0]]
    for prev, nxt in zip(ex, ex[1:]):
        if nxt[0] - prev[1] - 1 > max_intron:
            groups.append(g)
            g = []
        g.append(nxt)
    groups.append(g)
    if len(groups) == 1:
        return [x]
    out = []
    total = x.cds_len()
    for i, g in enumerate(groups, 1):
        y = Aln()
        for k in Aln.__slots__:
            setattr(y, k, getattr(x, k))
        y.id, y.exons, y.note = f"{x.id}.s{i}", g, "split_long_intron"
        qs = [e[2] for e in g if e[2] is not None]
        qe = [e[3] for e in g if e[3] is not None]
        if qs:
            y.qs, y.qe = min(qs), max(qe)
        frac = sum(e[1] - e[0] + 1 for e in g) / total
        y.score = int(x.score * frac)
        y.prot = ""  # translation of the unsplit alignment no longer applies
        out.append(y)
    return out


def overlap(a_exons, b_exons):
    ov = 0
    for s1, e1, *_ in a_exons:
        for s2, e2, *_ in b_exons:
            ov += max(0, min(e1, e2) - max(s1, s2) + 1)
    return ov


class Copy:
    def __init__(self, founder):
        self.founder, self.members, self.bridges, self.joined = founder, [founder], 0, False
        self.founders = [founder]
        self.neighbours = []  # founders of other copies on the same contig and strand

    def best_by_class(self):
        best = {}
        for m in self.members:
            if m.cls not in best or m.score > best[m.cls].score:
                best[m.cls] = m
        return best

    @property
    def cls(self):
        return self.founder.cls

    def core(self, tol=1000, max_foreign_bp=30):
        """Members lying within a founder's span (+-tol) and not reaching into another copy's
        CDS (> max_foreign_bp): excludes alignments that jump into a neighbour, including those
        that borrow only its short first exon.
        If the founder itself is such a bridge, members free of foreign CDS are used instead.
        Joined copies keep one founder per piece."""
        def clean(m):
            return all(overlap(m.exons, n.exons) <= max_foreign_bp for n in self.neighbours)
        spans = [f.span() for f in self.founders]
        inside = [m for m in self.members
                  if any(fs - tol <= m.span()[0] and m.span()[1] <= fe + tol for fs, fe in spans)]
        ok = [m for m in inside if clean(m)] or [m for m in self.members if clean(m)]
        return ok or inside

    def span(self):
        s = [m.span() for m in self.core()]
        return min(x[0] for x in s), max(x[1] for x in s)

    def intact_aln(self, min_cov):
        ok = [m for m in self.core()
              if m.cls == self.cls and m.cov() >= min_cov and m.fs == 0 and m.st == 0 and m.note == ""]
        return max(ok, key=lambda m: m.score) if ok else None

    def best_cov(self):
        return max((m.cov() for m in self.core() if m.cls == self.cls), default=0.0)


def cluster(alns, min_overlap):
    copies = []
    by_key = defaultdict(list)
    for x in sorted(alns, key=lambda x: -x.score):
        hits = []
        for c in by_key[(x.contig, x.strand)]:
            ov = overlap(x.exons, c.founder.exons)
            if ov and ov / min(x.cds_len(), c.founder.cds_len()) >= min_overlap:
                hits.append(c)
        if len(hits) == 1:
            hits[0].members.append(x)
        elif len(hits) > 1:
            for c in hits:
                c.bridges += 1
        else:
            c = Copy(x)
            copies.append(c)
            by_key[(x.contig, x.strand)].append(c)
    return copies


def join_fragments(copies, intact_cov, max_join):
    """Join neighbouring partial copies whose template ranges are complementary."""
    out = []
    by_key = defaultdict(list)
    for c in copies:
        by_key[(c.founder.contig, c.founder.strand)].append(c)
    for (contig, strand), cs in by_key.items():
        cs.sort(key=lambda c: c.span()[0], reverse=(strand == "-"))  # transcription order
        merged = [cs[0]]
        for c in cs[1:]:
            p = merged[-1]
            dist = (c.span()[0] - p.span()[1]) if strand == "+" else (p.span()[0] - c.span()[1])
            if (p.cls == c.cls and p.best_cov() < intact_cov and c.best_cov() < intact_cov
                    and 0 <= dist <= max_join and complementary(p, c)):
                p.members += c.members
                p.founders += c.founders
                p.joined = True
                continue
            merged.append(c)
        out += merged
    return out


def complementary(p, c, tol=30):
    """Same template aligned to the upstream piece's 5' part and the downstream piece's 3' part."""
    q1 = {m.query: m for m in p.members}
    for m in c.members:
        u = q1.get(m.query)
        if u is not None and u.qe <= m.qs + tol and u.qs < m.qs:
            return True
    return False


def n_runs(genome_path, contigs, min_n=10):
    runs = defaultdict(list)
    for name, seq in read_fasta(genome_path):
        if name not in contigs:
            continue
        for m in re.finditer(r"[Nn]{%d,}" % min_n, seq):
            runs[name].append((m.start() + 1, m.end()))
    return runs


def gff_records(cid, cls, x, intact):
    """GFF3 lines of one copy: gene, mRNA and CDS from alignment x (CDS phase from the transcription order)."""
    ex = sorted((e[0], e[1]) for e in x.exons)
    s, e = ex[0][0], ex[-1][1]
    attrs = f"class={cls};template={x.query};intact={int(intact)};score={x.score};frameshifts={x.fs};stops={x.st}"
    out = [f"{x.contig}\tminiprot_copies\tgene\t{s}\t{e}\t.\t{x.strand}\t.\tID={cid};{attrs}",
           f"{x.contig}\tminiprot_copies\tmRNA\t{s}\t{e}\t.\t{x.strand}\t.\tID={cid}.t1;Parent={cid};{attrs}"]
    done = 0
    for a0, a1 in (ex if x.strand == "+" else ex[::-1]):
        out.append(f"{x.contig}\tminiprot_copies\tCDS\t{a0}\t{a1}\t.\t{x.strand}\t{(3 - done % 3) % 3}\tParent={cid}.t1")
        done += a1 - a0 + 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="miniprot --gff --trans output for one genome")
    ap.add_argument("--classes", required=True, help="query_classes.tsv (id, class)")
    ap.add_argument("--species", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--genome", help="genome FASTA, to flag assembly gaps near copies")
    ap.add_argument("--max-intron", type=int, default=20000)
    ap.add_argument("--min-overlap", type=float, default=0.3)
    ap.add_argument("--intact-cov", type=float, default=0.8)
    ap.add_argument("--min-cov", type=float, default=0.2, help="drop alignments covering less of their template")
    ap.add_argument("--max-join", type=int, default=50000)
    ap.add_argument("--gap-window", type=int, default=5000)
    a = ap.parse_args()

    classes = {r["id"]: r["class"] for r in csv.DictReader(open(a.classes), delimiter="\t")}
    alns = []
    for x in parse_miniprot_gff(a.gff, classes):
        if x.exons:
            alns += split_long_introns(x, a.max_intron)
    alns = [x for x in alns if x.cov() >= a.min_cov]
    copies = join_fragments(cluster(alns, a.min_overlap), a.intact_cov, a.max_join)
    for c in copies:
        c.neighbours = [f for d in copies if d is not c
                        and (d.founder.contig, d.founder.strand) == (c.founder.contig, c.founder.strand)
                        for f in d.founders]
    copies.sort(key=lambda c: (c.founder.contig, c.span()[0]))
    for c in copies:
        if c.founder not in c.core():
            c.bridge_founder = True

    gaps = n_runs(a.genome, {c.founder.contig for c in copies}) if a.genome else {}
    counts = defaultdict(lambda: [0, 0])
    with open(a.out_prefix + ".copies.tsv", "w") as tsv, open(a.out_prefix + ".copies.faa", "w") as fa, \
            open(a.out_prefix + ".copies.gff3", "w") as gff:
        gff.write("##gff-version 3\n")
        tsv.write("species\tcopy_id\tcontig\tstart\tend\tstrand\tclass\tbest_template\tbest_score\t"
                  "margin_vs_other_class\tsecond_class\tbest_cov\tintact\tintact_template\tn_alignments\t"
                  "n_bridges\tflags\n")
        for i, c in enumerate(copies, 1):
            cid = f"{a.species}_mp{i}"
            s, e = c.span()
            best = c.best_by_class()
            others = sorted((m for k, m in best.items() if k != c.cls), key=lambda m: -m.score)
            margin = c.founder.score - others[0].score if others else c.founder.score
            ia = c.intact_aln(a.intact_cov)
            flags = []
            if c.joined:
                flags.append("joined_fragments")
            if getattr(c, "bridge_founder", False):
                flags.append("bridge_founder")
            if any(m.note for m in c.members):
                flags.append("has_split_alignments")
            if not ia and any(m.cls == c.cls and m.cov() >= a.intact_cov for m in c.core()):
                flags.append("pseudo_or_frameshift")
            g = [r for r in gaps.get(c.founder.contig, []) if r[1] >= s - a.gap_window and r[0] <= e + a.gap_window]
            if g:
                flags.append("gap_near")
            tsv.write("\t".join(map(str, [
                a.species, cid, c.founder.contig, s, e, c.founder.strand, c.cls, c.founder.query,
                c.founder.score, margin, others[0].cls if others else ".", f"{c.best_cov():.3f}",
                int(ia is not None), ia.query if ia else ".", len(c.members), c.bridges,
                ",".join(flags) or "."])) + "\n")
            rep = ia or max((m for m in c.core() if m.cls == c.cls and m.prot), key=lambda m: m.score, default=None)
            if rep is not None and rep.prot:
                write_fasta(fa, cid, rep.prot)
            gff.write("\n".join(gff_records(cid, c.cls, rep or c.founder, ia is not None)) + "\n")
            counts[c.cls][0] += 1
            counts[c.cls][1] += int(ia is not None)
    with open(a.out_prefix + ".counts.tsv", "w") as fh:
        fh.write("species\tmethod\tclass\tn_copies\tn_intact\n")
        for k in sorted(counts):
            fh.write(f"{a.species}\tminiprot\t{k}\t{counts[k][0]}\t{counts[k][1]}\n")
    print(f"{a.species}: {len(alns)} alignments -> {len(copies)} copies; " +
          ", ".join(f"{k} {v[0]} ({v[1]} intact)" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
