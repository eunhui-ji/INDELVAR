#!/usr/bin/env python3
"""The MANE Select transcript model: coordinates, CDS, and translation."""
from __future__ import annotations

import bisect
import collections
import gzip
import pathlib
import re

from . import refs
from .variant import Rejected

ATTR = re.compile(r"(\w+)=([^;]*)")

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

CODON_TABLE = {}
for _i, _aa in enumerate(
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _b = "TCAG"
    CODON_TABLE[_b[_i // 16] + _b[(_i // 4) % 4] + _b[_i % 4]] = _aa
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def translate(nt: str, initiator_met: bool = False) -> str:
    """Translate the CDS with the standard genetic code, stopping at the first stop."""
    out = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        codon = nt[i:i + 3]
        try:
            aa = CODON_TABLE[codon]
        except KeyError:
            raise Rejected("UNTRANSLATABLE_CODON", codon) from None
        if i == 0 and initiator_met:
            out.append("M")
            continue
        if aa == "*":
            return "".join(out)
        out.append(aa)
    return "".join(out)


class Cds:
    """CDS blocks of one transcript, ordered 5'->3' in transcript space."""

    def __init__(self, chrom_ucsc: str, strand: str, blocks: list[tuple[int, int]]):
        self.chrom_ucsc, self.strand = chrom_ucsc, strand
        self.blocks = sorted(blocks, reverse=(strand == "-"))
        self.cum: list[int] = []
        n = 0
        for a, b in self.blocks:
            self.cum.append(n)
            n += b - a + 1
        self.n = n

    def genomic(self, i: int) -> int:
        """Genomic position of CDS base i (1-based)."""
        if not 1 <= i <= self.n:
            raise Rejected("CDS_INDEX_OUT_OF_RANGE", f"{i} outside 1..{self.n}")
        k = bisect.bisect_right(self.cum, i - 1) - 1
        a, b = self.blocks[k]
        off = i - 1 - self.cum[k]
        return a + off if self.strand == "+" else b - off

    def contains(self, g: int) -> bool:
        return any(a <= g <= b for a, b in self.blocks)

    def span(self) -> tuple[int, int]:
        lo = min(a for a, _ in self.blocks)
        hi = max(b for _, b in self.blocks)
        return lo, hi


class Transcript:
    def __init__(self, tx_versioned: str, gene: str, strand: str,
                 chrom_ucsc: str, cds_blocks, exon_blocks, has_stop_codon: bool):
        self.id = tx_versioned
        self.base_id = tx_versioned.split(".")[0]
        self.gene = gene
        self.strand = strand
        self.chrom_ucsc = chrom_ucsc
        self.cds = Cds(chrom_ucsc, strand, cds_blocks)
        self.exons = sorted(exon_blocks)
        self.has_stop_codon = has_stop_codon

    @property
    def cds_length(self) -> int | None:
        """Summed CDS length in nucleotides, excluding the annotated stop codon."""
        if not self.has_stop_codon:
            return None
        return self.cds.n - 3

    def internal_exon_boundaries(self) -> list[int]:
        """Every exon start and end except the two outermost transcript ends."""
        if len(self.exons) < 2:
            return []
        sites: list[int] = []
        lo = self.exons[0][0]
        hi = self.exons[-1][1]
        for a, b in self.exons:
            if a != lo:
                sites.append(a)
            if b != hi:
                sites.append(b)
        return sorted(set(sites))


class ManeIndex:
    """Every MANE Select transcript, indexed by id and by genomic interval."""

    def __init__(self, gff: pathlib.Path, protein_faa: pathlib.Path):
        cds = collections.defaultdict(list)
        exon = collections.defaultdict(list)
        meta: dict[str, tuple[str, str, str]] = {}
        select: set[str] = set()
        stop: set[str] = set()
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                kind = f[2]
                if kind not in ("transcript", "CDS", "exon", "stop_codon"):
                    continue
                a = dict(ATTR.findall(f[8]))
                tx = a.get("transcript_id", "")
                if not tx:
                    continue
                base = tx.split(".")[0]
                if kind == "transcript":
                    if "MANE_Select" in a.get("tag", ""):
                        select.add(base)
                        meta[base] = (tx, a.get("gene_name", ""), f[6])
                    continue
                if kind == "CDS":
                    cds[base].append((int(f[3]), int(f[4]), f[0], f[6]))
                elif kind == "stop_codon":
                    stop.add(base)
                else:
                    exon[base].append((int(f[3]), int(f[4]), f[0]))

        self.by_id: dict[str, Transcript] = {}
        for base in select:
            if base not in cds:
                continue
            tx_v, gene, strand = meta[base]
            chrom = cds[base][0][2]
            self.by_id[base] = Transcript(
                tx_v, gene, strand, chrom,
                [(a, b) for a, b, _, _ in cds[base]],
                [(a, b) for a, b, _ in exon.get(base, [])],
                has_stop_codon=base in stop,
            )
        self.by_contig: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
        for base, t in self.by_id.items():
            lo, hi = t.cds.span()
            self.by_contig[t.chrom_ucsc].append((lo, hi, base))
        for v in self.by_contig.values():
            v.sort()

        self.protein: dict[str, str] = {}
        cur, buf = None, []
        with gzip.open(protein_faa, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur:
                        self.protein[cur] = "".join(buf)
                    m = re.search(r"(ENST\d+)", line)
                    cur, buf = (m.group(1) if m else None), []
                elif cur:
                    buf.append(line.strip())
        if cur:
            self.protein[cur] = "".join(buf)

    def containing(self, chrom_ucsc: str, lo: int, hi: int) -> list[Transcript]:
        """MANE Select transcripts whose CDS covers any base of [lo, hi]."""
        out = []
        for a, b, base in self.by_contig.get(chrom_ucsc, ()):
            if a > hi:
                break
            if b < lo:
                continue
            t = self.by_id[base]
            if any(t.cds.contains(g) for g in range(lo, hi + 1)):
                out.append(t)
        return out


_INDEX: ManeIndex | None = None


def index() -> ManeIndex:
    global _INDEX
    if _INDEX is None:
        _INDEX = ManeIndex(refs.path("mane_gff"), refs.path("mane_protein"))
    return _INDEX
