#!/usr/bin/env python3
"""Decide whether a variant is an in-frame indel this tool scores."""
from __future__ import annotations

PURE_DELETION = "pure_deletion"
PURE_INSERTION = "pure_insertion"
DELINS = "delins"

DELINS_REASON = "PROTEIN_LEVEL_DELINS_OUT_OF_SCOPE"


def classify_protein_edit(wt: str, mut: str, size_aa: int,
                          indel_type: str) -> tuple[str, int, str, str]:
    """Classify the protein consequence as an in-frame deletion, insertion, or out of scope."""
    expected = len(wt) - size_aa if indel_type == "deletion" else len(wt) + size_aa
    if len(mut) != expected:
        return (DELINS, -1, "",
                f"mutant protein length {len(mut)}, a pure {indel_type} of "
                f"{size_aa} aa would give {expected}; the pure equation does not "
                f"explain this protein")

    p = 0
    while p < min(len(wt), len(mut)) and wt[p] == mut[p]:
        p += 1

    if indel_type == "deletion":
        changed = wt[p:p + size_aa]
        if wt[p + size_aa:] == mut[p:]:
            return PURE_DELETION, p, changed, ""
        return (DELINS, p, changed,
                f"residues after the deleted span differ: WT "
                f"{wt[p + size_aa:p + size_aa + 4]!r} vs mutant {mut[p:p + 4]!r}; "
                f"the junction codon substitutes a residue, so the protein edit "
                f"is a delins rather than a {size_aa} aa deletion")
    changed = mut[p:p + size_aa]
    if mut[p + size_aa:] == wt[p:]:
        return PURE_INSERTION, p, changed, ""
    return (DELINS, p, changed,
            f"residues after the inserted span differ: mutant "
            f"{mut[p + size_aa:p + size_aa + 4]!r} vs WT {wt[p:p + 4]!r}; "
            f"the junction codon substitutes a residue, so the protein edit is a "
            f"delins rather than a {size_aa} aa insertion")


CROSSCHECK_OK = "OK"
CROSSCHECK_ALT_OK = "VEP_HGVSP_UNAVAILABLE_ALTERNATE_CROSSCHECK_PASS"
CROSSCHECK_FAILED = "HGVSP_CROSS_CHECK_FAILED"
IDENTITY_NOT_ESTABLISHED = "PROTEIN_EDIT_IDENTITY_NOT_ESTABLISHED"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def vep_mutant_protein(wt: str, protein_position: str,
                       amino_acids: str) -> tuple[str | None, str]:
    """The mutant protein VEP's Protein_position and Amino_acids describe, or None."""
    pp = (protein_position or "").strip()
    aa = (amino_acids or "").strip()
    if not pp or pp == "-" or "/" not in aa:
        return None, (f"VEP gave no HGVSp and no usable Protein_position/"
                      f"Amino_acids (pp={protein_position!r}, aa={amino_acids!r})")
    ref, alt = (s.strip() for s in aa.split("/", 1))
    ref = "" if ref == "-" else ref
    alt = "" if alt == "-" else alt
    if set(ref) - STANDARD_AA or set(alt) - STANDARD_AA:
        return None, (f"Amino_acids {aa!r} contains a residue outside the standard "
                      f"twenty, so the mutant sequence is not established")

    try:
        if "-" in pp:
            a_s, b_s = pp.split("-", 1)
            start, end = int(a_s), int(b_s)
        else:
            start = end = int(pp)
    except ValueError:
        return None, f"Protein_position {pp!r} is not a residue range"

    if not ref:
        # a pure insertion: the range names the two residues it falls between
        if end != start + 1 or not 1 <= start <= len(wt) or not alt:
            return None, (f"Protein_position {pp!r} does not name two adjacent "
                          f"residues of a {len(wt)} aa protein for an insertion")
        return wt[:start] + alt + wt[start:], f"insert {alt} between {start} and {end}"

    if end - start + 1 != len(ref) or start < 1 or end > len(wt):
        return None, (f"Protein_position {pp!r} spans {end - start + 1} residues "
                      f"but Amino_acids names {len(ref)}")
    if wt[start - 1:end] != ref:
        return None, (f"Amino_acids reference {ref!r} is not what the protein has "
                      f"at {start}-{end} ({wt[start - 1:end]!r}), so the string "
                      f"describes a different sequence")
    return wt[:start - 1] + alt + wt[end:], f"replace {start}-{end} {ref} with {alt}"


def crosscheck_protein(edit_hgvsp: str, wt_protein: str, mut_protein: str,
                       vep_hgvsp: str, vep_protein_position: str,
                       vep_amino_acids: str,
                       hgvsp_comparator) -> tuple[str, str]:
    """Compare the reconstructed protein with the VEP annotation."""
    if vep_hgvsp and vep_hgvsp != "-":
        agrees, detail = hgvsp_comparator(vep_hgvsp)
        return (CROSSCHECK_OK if agrees else CROSSCHECK_FAILED), detail

    theirs, detail = vep_mutant_protein(wt_protein, vep_protein_position,
                                        vep_amino_acids)
    if theirs is None:
        return IDENTITY_NOT_ESTABLISHED, detail
    if theirs != mut_protein:
        i = next((k for k in range(min(len(theirs), len(mut_protein)))
                  if theirs[k] != mut_protein[k]), min(len(theirs), len(mut_protein)))
        return (CROSSCHECK_FAILED,
                f"VEP HGVSp absent; its Protein_position/Amino_acids ({detail}) "
                f"give a mutant protein of {len(theirs)} aa that first differs "
                f"from our {len(mut_protein)} aa reconstruction at residue {i + 1} "
                f"({theirs[i:i + 6]!r} vs {mut_protein[i:i + 6]!r}); ours is "
                f"{edit_hgvsp}")
    return (CROSSCHECK_ALT_OK,
            f"VEP HGVSp absent; Protein_position/Amino_acids ({detail}) rebuild "
            f"the identical {len(theirs)} aa mutant protein")


# VEP emitting two rows for one feature with different protein edits.
VEP_TRANSCRIPT_AMBIGUOUS = "VEP_TRANSCRIPT_ANNOTATION_AMBIGUOUS"

JUNCTION_UNRESOLVED = "INSERTION_AT_EXON_EXON_JUNCTION_UNRESOLVED"
JUNCTION_PROJECTED = "INSERTION_PROJECTED_ACROSS_EXON_JUNCTION"


def _same_edit(a, b) -> bool:
    return (a.consequence == b.consequence and a.hgvsp == b.hgvsp
            and a.protein_position == b.protein_position
            and a.amino_acids == b.amino_acids)


def select_inframe_mane(candidates, vep_rows, indel_type, rows_for):
    """The MANE Select transcripts carrying an in-frame consequence."""
    want = "inframe_deletion" if indel_type == "deletion" else "inframe_insertion"
    hits, ambiguous = [], []
    for t in candidates:
        rows = [r for r in rows_for(vep_rows, t.id) if want in r.consequence]
        if not rows:
            continue
        first = rows[0]
        odd = [r for r in rows[1:] if not _same_edit(first, r)]
        if odd:
            ambiguous.append((t.id, (
                f"VEP returned {len(rows)} rows for {t.id} that describe "
                f"different edits: {first.hgvsp or first.amino_acids!r} vs "
                f"{odd[0].hgvsp or odd[0].amino_acids!r}")))
            continue
        hits.append((t, first))
    hits.sort(key=lambda p: p[0].id)
    ambiguous.sort()
    return hits, ambiguous
