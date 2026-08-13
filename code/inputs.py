#!/usr/bin/env python3
"""Read variants from VCF, TSV, CSV, or a single CHROM:POS:REF:ALT string."""
from __future__ import annotations

import csv
import dataclasses
import pathlib
import re

ASSEMBLY_OK = "GRCh38"
# Accepted GRCh38 assembly names.
ASSEMBLY_ALIASES = {"grch38", "hg38", "grch38.p14", "grch38.p13"}

_REQUIRED = ("chrom", "pos", "ref", "alt")
_ALIASES = {
    "chrom": {"chrom", "chr", "#chrom", "chromosome", "seqname"},
    "pos": {"pos", "position", "start"},
    "ref": {"ref", "reference", "ref_allele"},
    "alt": {"alt", "alternate", "alt_allele"},
    "assembly": {"assembly", "build", "genome_build"},
    "submitted_id": {"submitted_id", "id", "variant_id", "name"},
}


class InputError(ValueError):
    """Raised when an input record cannot be parsed."""


@dataclasses.dataclass
class Submission:
    submitted_id: str
    chrom: str
    pos: str
    ref: str
    alt: str
    assembly: str
    source_line: int
    reject: tuple[str, str] | None = None  # (code, detail) set at read time


def _norm_assembly(value: str | None, default: str) -> tuple[str, tuple[str, str] | None]:
    v = (value or default or "").strip()
    if not v:
        return "", ("ASSEMBLY_NOT_GRCH38", "no assembly given")
    if v.lower() not in ASSEMBLY_ALIASES:
        return v, ("ASSEMBLY_NOT_GRCH38",
                   f"assembly {v!r} is not GRCh38; this tool scores GRCh38 only")
    return ASSEMBLY_OK, None


def _make(sid: str, chrom: str, pos: str, ref: str, alt: str,
          assembly: str | None, default_assembly: str, line: int) -> Submission:
    asm, bad = _norm_assembly(assembly, default_assembly)
    s = Submission(submitted_id=sid, chrom=chrom.strip(), pos=str(pos).strip(),
                   ref=(ref or "").strip(), alt=(alt or "").strip(),
                   assembly=asm, source_line=line, reject=bad)
    if s.reject is None and (not s.ref or not s.alt or s.ref == "." or s.alt == "."):
        s.reject = ("INCOMPLETE_VARIANT", "REF and ALT are both required")
    return s


def read_single(spec: str, assembly: str) -> list[Submission]:
    """`CHROM:POS:REF:ALT`, or the same four separated by whitespace."""
    parts = spec.split(":") if ":" in spec else spec.split()
    if len(parts) != 4:
        raise InputError(
            f"cannot read {spec!r}: expected CHROM:POS:REF:ALT (four fields), "
            f"got {len(parts)}")
    return [_make(spec, parts[0], parts[1], parts[2], parts[3], None, assembly, 1)]


def _sniff(path: pathlib.Path) -> str:
    head = path.read_text(errors="replace")[:8192]
    if path.suffix.lower() == ".csv":
        return ","
    if path.suffix.lower() in (".tsv", ".txt"):
        return "\t"
    first = next((l for l in head.splitlines() if l.strip()), "")
    return "," if first.count(",") > first.count("\t") else "\t"


def read_table(path: pathlib.Path, assembly: str) -> list[Submission]:
    delim = _sniff(path)
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            raise InputError(f"{path} is empty") from None
        lower = [h.strip().lstrip("﻿").lower() for h in header]
        col: dict[str, int] = {}
        for key, names in _ALIASES.items():
            for i, h in enumerate(lower):
                if h in names and key not in col:
                    col[key] = i
        missing = [k for k in _REQUIRED if k not in col]
        if missing:
            raise InputError(
                f"{path}: no column found for {', '.join(missing)}. "
                f"Header was: {', '.join(header)}")
        out = []
        for n, row in enumerate(reader, start=2):
            if not any(f.strip() for f in row):
                continue
            def get(key, default=""):
                i = col.get(key)
                return row[i] if i is not None and i < len(row) else default
            if len(row) <= max(col[k] for k in _REQUIRED):
                sid = get("submitted_id") or f"{path.name}:line{n}"
                s = Submission(sid, "", "", "", "", "", n,
                               reject=("MALFORMED_RECORD",
                                       f"row has {len(row)} fields, fewer than "
                                       f"the header's {len(header)}"))
                out.append(s)
                continue
            sid = get("submitted_id") or ":".join(
                [get("chrom"), get("pos"), get("ref"), get("alt")])
            out.append(_make(sid, get("chrom"), get("pos"), get("ref"),
                             get("alt"), get("assembly") or None, assembly, n))
    if not out:
        raise InputError(f"{path} has a header but no data rows")
    return out


_VCF_ASSEMBLY = re.compile(r"(?:assembly|reference)=\S*?(GRCh3[78]|hg1[89])",
                           re.IGNORECASE)


def read_vcf(path: pathlib.Path, assembly: str) -> list[Submission]:
    out: list[Submission] = []
    declared: str | None = None
    seen_header = False
    with path.open() as fh:
        for n, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if line.startswith("##"):
                m = _VCF_ASSEMBLY.search(line)
                if m and declared is None:
                    declared = m.group(1)
                continue
            if line.startswith("#CHROM"):
                seen_header = True
                continue
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) < 5:
                out.append(Submission(f"{path.name}:line{n}", "", "", "", "", "",
                                      n, reject=("MALFORMED_RECORD",
                                                 f"VCF data line has {len(f)} "
                                                 "columns, fewer than the 5 "
                                                 "required")))
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            # Each alternate is normalised on its own.
            alts = [a for a in alt.split(",") if a != ""]
            if not alts:
                alts = [""]
            for k, a in enumerate(alts):
                base = vid if vid and vid != "." else f"{chrom}:{pos}:{ref}:{alt}"
                sid = base if len(alts) == 1 else f"{base}#alt{k + 1}"
                if a.startswith("<") or "[" in a or "]" in a or "*" in a:
                    out.append(Submission(sid, chrom, pos, ref, a, ASSEMBLY_OK, n,
                                          reject=("UNSUPPORTED_ALT_REPRESENTATION",
                                                  f"ALT {a!r} is a symbolic, "
                                                  "breakend or overlapping-deletion "
                                                  "variant, not a literal sequence")))
                    continue
                out.append(_make(sid, chrom, pos, ref, a,
                                 declared, assembly, n))
    if not seen_header:
        raise InputError(f"{path} has no #CHROM header line; it is not a VCF")
    if not out:
        raise InputError(f"{path} has a VCF header but no data lines")
    return out


def read_any(path: pathlib.Path | None, single: str | None,
             assembly: str) -> list[Submission]:
    if single is not None:
        return read_single(single, assembly)
    if path is None:
        raise InputError("give --input or --variant")
    if not path.exists():
        raise InputError(f"input file not found: {path}")
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".vcf") or suffix.endswith(".vcf.gz"):
        if suffix.endswith(".gz"):
            raise InputError("gzipped VCF is not read directly; decompress first")
        return read_vcf(path, assembly)
    return read_table(path, assembly)


def count_duplicate_ids(subs: list[Submission]) -> int:
    """Count repeated submitted ids."""
    counts: dict[str, int] = {}
    for s in subs:
        counts[s.submitted_id] = counts.get(s.submitted_id, 0) + 1
    return sum(1 for v in counts.values() if v > 1)
