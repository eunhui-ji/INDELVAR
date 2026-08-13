#!/usr/bin/env python3
"""The eight sequence-context features, including the three deltas."""
from __future__ import annotations

import dataclasses

from . import scales
from .protein_edit import ProteinEdit
from .variant import NormalisedVariant, Rejected

FLANK = 10


@dataclasses.dataclass(frozen=True)
class LocalSequences:
    wt: str
    mut: str
    left_flank: str
    right_flank: str
    changed: str


def local_sequences(edit: ProteinEdit, flank: int = FLANK) -> LocalSequences:
    wt = edit.wt_protein
    i = edit.offset
    if edit.kind == "deletion":
        j = i + edit.size_aa
        left = wt[max(0, i - flank):i]
        right = wt[j:j + flank]
        return LocalSequences(wt=left + edit.changed_residues + right,
                              mut=left + right,
                              left_flank=left, right_flank=right,
                              changed=edit.changed_residues)
    left = wt[max(0, i - flank):i]
    right = wt[i:i + flank]
    return LocalSequences(wt=left + right,
                          mut=left + edit.changed_residues + right,
                          left_flank=left, right_flank=right,
                          changed=edit.changed_residues)


def delta_local_hydropathy(loc: LocalSequences) -> float:
    kd = scales.scale(scales.HYDROPATHY_SCALE)
    return scales.mean_over(loc.mut, kd) - scales.mean_over(loc.wt, kd)


def delta_local_entropy(loc: LocalSequences) -> float:
    return scales.shannon_bits(loc.mut) - scales.shannon_bits(loc.wt)


def delta_local_charge(loc: LocalSequences) -> float:
    q = scales.charge_scale()
    return scales.mean_over(loc.mut, q) - scales.mean_over(loc.wt, q)


def local_entropy(edit: ProteinEdit, uniprot_sequence: str | None) -> float | None:
    """Shannon entropy in bits over the wild-type window; NA below three residues."""
    if uniprot_sequence is None:
        return None
    if edit.kind == "deletion":
        a, b = edit.deleted_residue_span
    else:
        ctx = edit.structure_context_residues()
        a, b = min(ctx), max(ctx)
    lo = max(1, a - FLANK)
    hi = min(len(uniprot_sequence), b + FLANK)
    window = uniprot_sequence[lo - 1:hi]
    kept = "".join(c for c in window if c in scales.scale(scales.HYDROPATHY_SCALE))
    if len(kept) < 3:
        return None
    return scales.shannon_bits(kept)


def indel_size_relative(variant: NormalisedVariant, edit: ProteinEdit) -> float:
    n = len(edit.wt_protein)
    if n <= 0:
        raise Rejected("PROTEIN_LENGTH_UNRESOLVED", edit.transcript.id)
    return variant.size_aa / n


def dist_to_splice(variant: NormalisedVariant, transcript) -> int | None:
    """Return the minimum distance, in bases, from the variant span to an
    internal boundary."""
    sites = transcript.internal_exon_boundaries()
    if not sites:
        return None
    lo, hi = variant.genomic_span
    return min(min(abs(b - s) for s in sites) for b in range(lo, hi + 1))
