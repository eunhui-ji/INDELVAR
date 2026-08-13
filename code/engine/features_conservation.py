#!/usr/bin/env python3
"""phylop100_mean and gerp_mean over the affected bases."""
from __future__ import annotations

import dataclasses

import pyBigWig

from . import refs
from .protein_edit import ProteinEdit
from .transcripts import Transcript
from .variant import ContigMap, NormalisedVariant


@dataclasses.dataclass
class ConservationResult:
    phylop_mean: float | None
    gerp_mean: float | None
    qc: str


_BW: dict[str, object] = {}


def _bw(key: str):
    if key not in _BW:
        _BW[key] = pyBigWig.open(str(refs.path(key)))
    return _BW[key]


def deleted_reference_bases(variant: NormalisedVariant) -> list[int]:
    """The deleted bases; the anchor at POS is excluded."""
    return list(range(variant.pos + 1, variant.pos + len(variant.ref)))


def insertion_context_bases(edit: ProteinEdit, t: Transcript) -> tuple[list[int], str]:
    """The twelve CDS bases of the four-residue context, in transcript order."""
    residues = edit.structure_context_residues()
    if len(residues) != 4:
        return [], "INCOMPLETE_CONTEXT_PROTEIN_TERMINAL"
    out: list[int] = []
    for r in residues:
        first = 3 * r - 2
        if first < 1 or first + 2 > t.cds.n:
            return [], "INCOMPLETE_CONTEXT_TRANSCRIPT_TERMINAL"
        out.extend(t.cds.genomic(first + k) for k in range(3))
    return out, ""


def _values(bw, contig: str, positions: list[int]) -> list[float]:
    out = []
    for p in positions:
        try:
            v = bw.values(contig, p - 1, p)
        except (RuntimeError, TypeError):
            continue
        if v and v[0] == v[0]:
            out.append(v[0])
    return out


def compute(variant: NormalisedVariant, edit: ProteinEdit, t: Transcript,
            cmap: ContigMap) -> ConservationResult:
    if variant.indel_type == "deletion":
        positions, qc = deleted_reference_bases(variant), ""
    else:
        positions, qc = insertion_context_bases(edit, t)
    if not positions:
        return ConservationResult(None, None, qc or "NO_BASES_RESOLVED")

    ucsc = cmap.ucsc(variant.chrom)
    ens = cmap.ensembl(variant.chrom)
    pv = _values(_bw("phylop"), ucsc, positions)
    gv = _values(_bw("gerp"), ens, positions)
    return ConservationResult(
        phylop_mean=(sum(pv) / len(pv)) if pv else None,
        gerp_mean=(sum(gv) / len(gv)) if gv else None,
        qc=qc or ("" if pv and gv else "TRACK_LOOKUP_EMPTY"),
    )
