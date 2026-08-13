#!/usr/bin/env python3
"""Apply a variant to the MANE Select CDS and reconstruct both proteins."""
from __future__ import annotations

import dataclasses
import re

from . import transcripts as tx_mod
from .transcripts import Transcript, revcomp, translate, STANDARD_AA
from . import scope as SCOPE
from .variant import Faidx, NormalisedVariant, Rejected

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
       "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
       "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V", "Ter": "*"}
AA1 = {v: k for k, v in AA3.items()}

HGVSP_DEL = re.compile(r"p\.([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?del$")
HGVSP_INS = re.compile(r"p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins([A-Za-z]+)$")
HGVSP_DUP = re.compile(r"p\.([A-Z][a-z]{2})(\d+)(?:_([A-Z][a-z]{2})(\d+))?dup$")


@dataclasses.dataclass(frozen=True)
class ProteinEdit:
    transcript: Transcript
    wt_protein: str
    mut_protein: str
    kind: str                      # deletion | insertion
    # 0-based: first removed residue, or residues before an insertion.
    offset: int
    changed_residues: str          # deleted residues, or inserted residues
    hgvsp_canonical: str           # rendered under the HGVS 3' rule

    @property
    def size_aa(self) -> int:
        return len(self.changed_residues)

    @property
    def deleted_residue_span(self) -> tuple[int, int]:
        """1-based inclusive [start, end] of the deleted residues."""
        if self.kind != "deletion":
            raise Rejected("NOT_A_DELETION", self.kind)
        return self.offset + 1, self.offset + self.size_aa

    @property
    def left_junction_residue(self) -> int:
        """1-based index of the wild-type residue immediately 5' of an insertion."""
        if self.kind != "insertion":
            raise Rejected("NOT_AN_INSERTION", self.kind)
        return self.offset

    def structure_context_residues(self) -> list[int]:
        """The window the structure and annotation features aggregate over."""
        if self.kind == "deletion":
            a, b = self.deleted_residue_span
            return list(range(a, b + 1))
        j = self.left_junction_residue
        n = len(self.wt_protein)
        return [r for r in (j - 1, j, j + 1, j + 2) if 1 <= r <= n]


def _edited_cds(variant: NormalisedVariant, t: Transcript, fa: Faidx) -> str:
    """Edit the CDS in genomic space, then reverse for transcript orientation."""
    inserted_at = inserted_before = None
    if variant.indel_type == "deletion":
        removed = set(range(variant.pos + 1, variant.pos + len(variant.ref)))
        inserted = ""
        missing = [g for g in sorted(removed) if not t.cds.contains(g)]
        if missing:
            raise Rejected("DELETION_OUTSIDE_CDS",
                           f"{len(missing)} deleted base(s) are not in the CDS of "
                           f"{t.id}, first at {missing[0]}")
    else:
        removed = set()
        inserted = variant.alt[1:]
        left, right = t.cds.contains(variant.pos), t.cds.contains(variant.pos + 1)
        if not left and not right:
            raise Rejected("INSERTION_AT_CDS_EDGE",
                           f"insertion point {variant.pos}/{variant.pos + 1} has "
                           f"neither anchor in a CDS block of {t.id}")
        if left and right and not any(a <= variant.pos and variant.pos + 1 <= b
                                      for a, b in t.cds.blocks):
            raise Rejected("INSERTION_AT_EXON_BOUNDARY",
                           f"insertion point {variant.pos}/{variant.pos + 1} straddles "
                           f"two CDS blocks of {t.id}")
        # One coding anchor puts the point on an exon edge.
        if left:
            inserted_at = variant.pos
        else:
            inserted_before = variant.pos + 1

    pieces = []
    for a, b in sorted(t.cds.blocks):
        block = fa.fetch(variant.chrom, a, b)
        buf = []
        for k, g in enumerate(range(a, b + 1)):
            if g == inserted_before:
                buf.append(inserted)
            if g in removed:
                continue
            buf.append(block[k])
            if g == inserted_at:
                buf.append(inserted)
        pieces.append("".join(buf))
    seq = "".join(pieces)
    return seq if t.strand == "+" else revcomp(seq)


def junction_projected(variant: NormalisedVariant, t: Transcript) -> bool:
    """True when an insertion anchor sits on a CDS block edge rather than inside."""
    if variant.indel_type != "insertion":
        return False
    return t.cds.contains(variant.pos) != t.cds.contains(variant.pos + 1)


def _wt_cds(t: Transcript, chrom: str, fa: Faidx) -> str:
    seq = "".join(fa.fetch(chrom, a, b) for a, b in sorted(t.cds.blocks))
    return seq if t.strand == "+" else revcomp(seq)


def _shift_3prime_del(wt: str, i: int, length: int) -> int:
    while i + length < len(wt) and wt[i] == wt[i + length]:
        i += 1
    return i


def _shift_3prime_ins(wt: str, i: int, ins: str) -> tuple[int, str]:
    while i < len(wt) and ins and ins[0] == wt[i]:
        ins = ins[1:] + wt[i]
        i += 1
    return i, ins


def _render_hgvsp(wt: str, kind: str, i: int, changed: str) -> str:
    """HGVS protein string under the most-3' rule."""
    if kind == "deletion":
        j = _shift_3prime_del(wt, i, len(changed))
        seg = wt[j:j + len(changed)]
        if len(seg) == 1:
            return f"p.{AA1[seg[0]]}{j + 1}del"
        return f"p.{AA1[seg[0]]}{j + 1}_{AA1[seg[-1]]}{j + len(seg)}del"
    j, ins = _shift_3prime_ins(wt, i, changed)
    # A pure duplication is written dup by HGVS, not ins.
    if j - len(ins) >= 0 and wt[j - len(ins):j] == ins:
        a = j - len(ins)
        if len(ins) == 1:
            return f"p.{AA1[ins[0]]}{a + 1}dup"
        return f"p.{AA1[ins[0]]}{a + 1}_{AA1[ins[-1]]}{j}dup"
    # HGVS writes an insertion after the final residue against the stop codon.
    left = wt[j - 1]
    right_aa = AA1[wt[j]] if j < len(wt) else "Ter"
    return (f"p.{AA1[left]}{j}_{right_aa}{j + 1}ins"
            + "".join(AA1[c] for c in ins))


def build(variant: NormalisedVariant, t: Transcript, fa: Faidx) -> ProteinEdit:
    """Reconstruct both proteins and locate the edit."""
    wt_nt = _wt_cds(t, variant.chrom, fa)
    mut_nt = _edited_cds(variant, t, fa)
    # The annotated start translates as Met when the initiator codon is intact.
    init = wt_nt[:3] == mut_nt[:3]
    wt = translate(wt_nt, initiator_met=True)
    mut = translate(mut_nt, initiator_met=init)
    if not init:
        raise Rejected("INITIATOR_CODON_ALTERED",
                       f"{t.id}: the edit changes the start codon "
                       f"{wt_nt[:3]} -> {mut_nt[:3]}")

    reference = tx_mod.index().protein.get(t.base_id)
    if reference is None:
        raise Rejected("MANE_PROTEIN_ABSENT", t.id)
    if wt != reference:
        raise Rejected("CDS_PROTEIN_MISMATCH",
                       f"{t.id}: CDS translation length {len(wt)} vs MANE protein "
                       f"length {len(reference)}")

    kind, p, changed, detail = SCOPE.classify_protein_edit(
        wt, mut, variant.size_aa, variant.indel_type)
    if kind == SCOPE.DELINS:
        raise Rejected(SCOPE.DELINS_REASON, f"{t.id}: {detail}")

    bad = sorted(set(changed) - STANDARD_AA)
    if bad:
        raise Rejected("NON_STANDARD_RESIDUE", f"{t.id}: {bad}")

    return ProteinEdit(
        transcript=t, wt_protein=wt, mut_protein=mut,
        kind=variant.indel_type, offset=p, changed_residues=changed,
        hgvsp_canonical=_render_hgvsp(wt, variant.indel_type, p, changed),
    )


def cross_check_hgvsp(edit: ProteinEdit, vep_hgvsp: str) -> tuple[bool, str]:
    """Validate the reconstructed protein against the VEP HGVSp."""
    if not vep_hgvsp or vep_hgvsp == "-":
        return False, "VEP returned no HGVSp"
    got = vep_hgvsp.split(":")[-1]
    ours = edit.hgvsp_canonical
    if got == ours:
        return True, got

    # VEP shifts 3' along the genome, which is 5' along the protein.
    theirs = apply_hgvsp(edit.wt_protein, got)
    if theirs is None:
        return False, f"VEP {got} could not be applied to the reference protein"
    if theirs == edit.mut_protein:
        return True, (f"VEP {got} and ours {ours} are the same edit written at "
                      f"different positions of a repeat; the mutant proteins are "
                      f"identical")
    return False, f"VEP {got} vs ours {ours}"


def apply_hgvsp(wt: str, s: str) -> str | None:
    """The mutant protein a well-formed HGVSp describes, or None."""
    def at(three: str, one_based: str) -> int | None:
        i = int(one_based) - 1
        if not 0 <= i < len(wt) or AA3.get(three) != wt[i]:
            return None
        return i

    m = HGVSP_DEL.match(s)
    if m:
        a = at(m.group(1), m.group(2))
        b = a if m.group(3) is None else at(m.group(3), m.group(4))
        if a is None or b is None or b < a:
            return None
        return wt[:a] + wt[b + 1:]

    m = HGVSP_DUP.match(s)
    if m:
        a = at(m.group(1), m.group(2))
        b = a if m.group(3) is None else at(m.group(3), m.group(4))
        if a is None or b is None or b < a:
            return None
        return wt[:b + 1] + wt[a:b + 1] + wt[b + 1:]

    m = HGVSP_INS.match(s)
    if m:
        a = at(m.group(1), m.group(2))
        # the right flank may be the stop codon, which is not in wt
        b = int(m.group(4)) - 1
        if a is None or b != a + 1 or b > len(wt):
            return None
        if b < len(wt) and AA3.get(m.group(3)) != wt[b]:
            return None
        if b == len(wt) and m.group(3) != "Ter":
            return None
        res = m.group(5)
        ins = "".join(AA3.get(res[i:i + 3], "?") for i in range(0, len(res), 3))
        if len(res) % 3 or "?" in ins:
            return None
        return wt[:a + 1] + ins + wt[a + 1:]
    return None
