"""Minimal GFF3 / FASTA helpers shared by the annotation scripts (stdlib only)."""
from collections import defaultdict


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(fh, name, seq):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i:i + 80] + "\n")


def parse_attrs(field):
    out = {}
    for item in field.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def read_cds(gff_path):
    """Return {transcript_id: Transcript} built from CDS features (Parent=...)."""
    tx = {}
    gene_of = {}
    for line in open(gff_path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = parse_attrs(f[8])
        if f[2] in ("mRNA", "transcript"):
            gene_of[a.get("ID")] = a.get("Parent")
        if f[2] != "CDS":
            continue
        for parent in a.get("Parent", "").split(","):
            t = tx.setdefault(parent, Transcript(parent, f[0], f[6]))
            phase = int(f[7]) if f[7] in "012" else 0
            t.cds.append((int(f[3]), int(f[4]), phase))
    for tid, t in tx.items():
        t.gene = gene_of.get(tid)
        t.sort()
    return tx


class Transcript:
    def __init__(self, tid, contig, strand):
        self.id, self.contig, self.strand, self.cds, self.gene = tid, contig, strand, [], None

    def sort(self):
        # transcription order
        self.cds.sort(key=lambda c: c[0], reverse=(self.strand == "-"))

    @property
    def start(self):
        return min(c[0] for c in self.cds)

    @property
    def end(self):
        return max(c[1] for c in self.cds)

    def introns(self):
        ex = sorted((s, e) for s, e, _ in self.cds)
        return [b[0] - a[1] - 1 for a, b in zip(ex, ex[1:])]

    def aa_to_genome(self, aa_start, aa_end):
        """Genomic (start, end) of codons aa_start..aa_end (1-based, inclusive)."""
        first_phase = self.cds[0][2]
        nt_from = first_phase + 3 * (aa_start - 1)          # 0-based offset in CDS
        nt_to = first_phase + 3 * aa_end - 1
        pos, coords, offset = None, [], 0
        for s, e, _ in self.cds:
            length = e - s + 1
            for target in (nt_from, nt_to):
                if offset <= target < offset + length:
                    k = target - offset
                    coords.append(s + k if self.strand == "+" else e - k)
            offset += length
        if len(coords) < 2:  # position beyond annotated CDS: clamp to model ends
            coords = [self.start, self.end] if not coords else coords + [self.end if self.strand == "+" else self.start]
        return min(coords), max(coords)


def group_by(items, key):
    d = defaultdict(list)
    for it in items:
        d[key(it)].append(it)
    return d
