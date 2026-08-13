#!/usr/bin/env python3
"""Validate, normalize, and scope GRCh38 in-frame indels."""
from __future__ import annotations

import dataclasses
import pathlib
import re

from . import refs

VALID_BASES = re.compile(r"^[ACGTN]+$")


class Rejected(Exception):
    """Rejection carrying a machine-readable reason code."""

    def __init__(self, code: str, detail: str = ""):
        super().__init__(f"{code}{': ' + detail if detail else ''}")
        self.code, self.detail = code, detail


class Faidx:
    """Random-access FASTA reader, 1-based inclusive."""

    def __init__(self, fasta: pathlib.Path, fai: pathlib.Path):
        self.fh = fasta.open("rb")
        self.idx: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text().splitlines():
            name, length, offset, line_bases, line_width = line.split("\t")[:5]
            self.idx[name] = (int(length), int(offset), int(line_bases), int(line_width))

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig not in self.idx:
            raise Rejected("CHROM_NOT_IN_REFERENCE", contig)
        length, offset, lb, lw = self.idx[contig]
        if start < 1 or end > length or end < start:
            raise Rejected("POS_OUT_OF_RANGE", f"{contig}:{start}-{end} outside 1..{length}")
        n = end - start + 1
        s0 = start - 1
        self.fh.seek(offset + (s0 // lb) * lw + (s0 % lb))
        got = b""
        while len(got.replace(b"\n", b"")) < n:
            chunk = self.fh.read(n * 2 + 64)
            if not chunk:
                break
            got += chunk
        return got.replace(b"\n", b"").decode()[:n].upper()


class ContigMap:
    """Translates contig names between the Ensembl and UCSC spellings."""

    def __init__(self, tsv: pathlib.Path):
        self.to_ensembl: dict[str, str] = {}
        self.to_ucsc: dict[str, str] = {}
        for line in tsv.read_text().splitlines():
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.split("\t")
            if f[0] == "ensembl_name":
                continue
            ens, ucsc = f[0], f[1]
            self.to_ensembl[ens] = ens
            self.to_ensembl[ucsc] = ens
            self.to_ucsc[ens] = ucsc
            self.to_ucsc[ucsc] = ucsc

    def ensembl(self, name: str) -> str:
        try:
            return self.to_ensembl[name]
        except KeyError:
            raise Rejected("CHROM_NOT_CANONICAL", name) from None

    def ucsc(self, name: str) -> str:
        try:
            return self.to_ucsc[name]
        except KeyError:
            raise Rejected("CHROM_NOT_CANONICAL", name) from None


@dataclasses.dataclass(frozen=True)
class NormalisedVariant:
    """A normalised variant; `chrom` is the Ensembl name and `submitted_*` keep
    the input."""
    chrom: str
    pos: int
    ref: str
    alt: str
    indel_type: str          # "deletion" | "insertion"
    size_aa: int
    submitted_chrom: str
    submitted_pos: int
    submitted_ref: str
    submitted_alt: str

    @property
    def variant_id(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def submitted_id(self) -> str:
        return (f"{self.submitted_chrom}:{self.submitted_pos}:"
                f"{self.submitted_ref}:{self.submitted_alt}")

    @property
    def genomic_span(self) -> tuple[int, int]:
        """Inclusive genomic span the variant touches, anchor base included."""
        if self.indel_type == "deletion":
            return self.pos, self.pos + len(self.ref) - 1
        return self.pos, self.pos + 1


def _validate(chrom, pos, ref, alt) -> tuple[str, int, str, str]:
    if not isinstance(chrom, str) or not chrom.strip():
        raise Rejected("INPUT_SCHEMA_INVALID", "CHROM empty")
    try:
        pos = int(pos)
    except (TypeError, ValueError):
        raise Rejected("INPUT_SCHEMA_INVALID", f"POS not an integer: {pos!r}") from None
    if pos < 1:
        raise Rejected("INPUT_SCHEMA_INVALID", f"POS {pos} < 1")
    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()
    for name, seq in (("REF", ref), ("ALT", alt)):
        if not seq:
            raise Rejected("INPUT_SCHEMA_INVALID", f"{name} empty")
        if not VALID_BASES.match(seq):
            raise Rejected("INPUT_SCHEMA_INVALID", f"{name} not ACGTN: {seq!r}")
    if ref == alt:
        raise Rejected("INPUT_SCHEMA_INVALID", "REF equals ALT")
    return chrom.strip(), pos, ref, alt


def validate_schema(chrom, pos, ref, alt) -> None:
    """Raise unless CHROM/POS/REF/ALT is well formed."""
    _validate(chrom, pos, ref, alt)


def left_normalise(chrom: str, pos: int, ref: str, alt: str,
                   fa: Faidx) -> tuple[int, str, str]:
    """Left-align and trim REF/ALT against the reference."""
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    # The indel can be written one base earlier while the trailing bases agree.
    while (len(ref) == 1 or len(alt) == 1) and ref[-1] == alt[-1] and pos > 1:
        prev = fa.fetch(chrom, pos - 1, pos - 1)
        ref, alt, pos = prev + ref[:-1], prev + alt[:-1], pos - 1
    return pos, ref, alt


def parse(chrom, pos, ref, alt, fa: Faidx, cmap: ContigMap) -> NormalisedVariant:
    """Normalise a variant, or raise Rejected with a reason code."""
    s_chrom, s_pos, s_ref, s_alt = _validate(chrom, pos, ref, alt)
    ens = cmap.ensembl(s_chrom)

    # REF is checked before normalisation, which moves the coordinates.
    observed = fa.fetch(ens, s_pos, s_pos + len(s_ref) - 1)
    if observed != s_ref:
        raise Rejected("REF_MISMATCH",
                       f"{ens}:{s_pos} submitted REF {s_ref} but reference has {observed}")

    n_pos, n_ref, n_alt = left_normalise(ens, s_pos, s_ref, s_alt, fa)

    if len(n_ref) > 1 and len(n_alt) > 1:
        raise Rejected("COMPLEX_INDEL",
                       f"REF and ALT both longer than one base after normalisation: "
                       f"{n_ref}>{n_alt}")
    if len(n_ref) == 1 and len(n_alt) == 1:
        raise Rejected("NOT_AN_INDEL", f"{n_ref}>{n_alt} is a substitution")
    if n_ref[0] != n_alt[0]:
        raise Rejected("ANCHOR_BASE_SUBSTITUTED",
                       f"{n_ref}>{n_alt} substitutes the anchor base as well "
                       f"as changing length; this is a delins")
    delta = len(n_alt) - len(n_ref)
    if delta % 3 != 0:
        raise Rejected("FRAMESHIFT", f"length change {delta} is not a multiple of three")
    size_aa = abs(delta) // 3
    itype = "insertion" if delta > 0 else "deletion"
    if not 1 <= size_aa <= 10:
        raise Rejected("SIZE_OUT_OF_RANGE", f"{size_aa} aa outside 1..10")

    return NormalisedVariant(
        chrom=ens, pos=n_pos, ref=n_ref, alt=n_alt,
        indel_type=itype, size_aa=size_aa,
        submitted_chrom=s_chrom, submitted_pos=s_pos,
        submitted_ref=s_ref, submitted_alt=s_alt,
    )


def open_reference() -> tuple[Faidx, ContigMap]:
    return (Faidx(refs.path("grch38_fasta"), refs.path("grch38_fai")),
            ContigMap(refs.path("contig_map")))
