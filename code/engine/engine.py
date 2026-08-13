#!/usr/bin/env python3
"""Turn a GRCh38 variant into the 37-feature table the model scores."""
from __future__ import annotations

import dataclasses

from . import (features_annotation as FA, features_conservation as FC,
               features_constraint as FK, features_sequence as FS,
               features_structure as FT, protein_edit as PE, refs,
               transcripts as TX, vep as VEP)
from . import scope as SCOPE
from .variant import Rejected, open_reference, parse as parse_variant

@dataclasses.dataclass
class Row:
    identity: dict[str, object]
    features: dict[str, object]


def _reject(vid: str, submitted: str, code: str, detail: str,
            feature_order: list[str]) -> Row:
    return Row({"variant_id": vid, "submitted_id": submitted, "assembly": "GRCh38",
                "chrom": None, "pos": None, "ref": None, "alt": None,
                "gene_symbol": None, "transcript_id": None, "uniprot_acc": None,
                "hgvs_p": None, "status": code, "qc": detail},
               {f: None for f in feature_order})


def run_grouped(inputs: list[tuple[str, object, str, str]],
                use_vep: bool = True) -> list[list[Row]]:
    """Feature rows per input; one row per MANE Select transcript."""
    feature_order = refs.contract_features()
    fa, cmap = open_reference()
    idx = TX.index()

    staged: list[tuple[int, object, object]] = []
    # One row per MANE Select transcript.
    out: list[list[Row]] = [[] for _ in inputs]

    for k, (chrom, pos, ref, alt) in enumerate(inputs):
        submitted = f"{chrom}:{pos}:{ref}:{alt}"
        try:
            v = parse_variant(chrom, pos, ref, alt, fa, cmap)
        except Rejected as exc:
            out[k].append(_reject(submitted, submitted, exc.code, exc.detail,
                                  feature_order))
            continue
        cands = idx.containing(cmap.ucsc(v.chrom), *v.genomic_span)
        if not cands:
            out[k].append(_reject(v.variant_id, submitted, "NO_MANE_SELECT_CDS",
                                  f"no MANE Select CDS covers {v.variant_id}",
                                  feature_order))
            continue
        staged.append((k, v, cands))

    vep_rows: dict[str, list[VEP.VepRow]] = {}
    if use_vep and staged:
        try:
            vep_rows = VEP.annotate(
                [(f"v{k}", v.chrom, v.pos, v.ref, v.alt) for k, v, _ in staged],
                tag="engine")
        except VEP.VepFailed as exc:
            for k, v, _ in staged:
                out[k].append(_reject(v.variant_id, v.submitted_id, "VEP_FAILED",
                                      str(exc)[:400], feature_order))
            staged = []

    # One pass over the UniProt flat file for the whole batch.
    if staged:
        FA.read_records({acc for acc, _ in
                         (FA.mane_to_uniprot(t.base_id)
                          for _, _, cs in staged for t in cs)
                         if acc})

    for k, v, cands in staged:
        vrows = vep_rows.get(f"v{k}", [])
        if use_vep:
            hits, ambiguous = SCOPE.select_inframe_mane(
                cands, vrows, v.indel_type, VEP.transcript_rows)
            # A disagreeing transcript is refused on its own.
            for tid, detail in ambiguous:
                out[k].append(_reject(v.variant_id, v.submitted_id,
                                      SCOPE.VEP_TRANSCRIPT_AMBIGUOUS,
                                      f"{tid}: {detail}", feature_order))
            if not hits:
                if not ambiguous:
                    only = cands[0]
                    row = VEP.pick_transcript_row(vrows, only.id)
                    if row is None:
                        out[k].append(_reject(
                            v.variant_id, v.submitted_id, "VEP_TRANSCRIPT_ABSENT",
                            ",".join(c.id for c in cands), feature_order))
                    else:
                        want = ("inframe_deletion" if v.indel_type == "deletion"
                                else "inframe_insertion")
                        out[k].append(_reject(
                            v.variant_id, v.submitted_id, "CONSEQUENCE_DISAGREES",
                            f"expected {want} on a MANE Select transcript, VEP "
                            f"says {sorted(row.consequence)}", feature_order))
                continue
        else:
            if len(cands) > 1:
                out[k].append(_reject(
                    v.variant_id, v.submitted_id, "MANE_SELECTION_REQUIRES_VEP",
                    f"{len(cands)} MANE Select CDS regions cover this variant "
                    f"({','.join(c.id for c in cands)}); which of them carry an "
                    f"in-frame consequence cannot be decided with VEP disabled",
                    feature_order))
                continue
            hits = [(cands[0], None)]

        for t, vr in hits:
            row = _one_transcript(v, t, vr, fa, cmap, feature_order, use_vep)
            out[k].append(row)

    # Deterministic row order within one submission.
    for rows_k in out:
        rows_k.sort(key=lambda r: (r.identity["status"] != "OK",
                                   str(r.identity.get("transcript_id") or ""),
                                   str(r.identity.get("status") or ""),
                                   str(r.identity.get("qc") or "")))
    return out


def _one_transcript(v, t, vr, fa, cmap, feature_order, use_vep):
    qc: list[str] = []
    try:
        edit = PE.build(v, t, fa)
    except Rejected as exc:
        return _reject(v.variant_id, v.submitted_id, exc.code, exc.detail,
                       feature_order)

    hgvsp = edit.hgvsp_canonical
    junction = PE.junction_projected(v, t)
    if use_vep:
        tier, detail = SCOPE.crosscheck_protein(
            edit.hgvsp_canonical, edit.wt_protein, edit.mut_protein, vr.hgvsp,
            vr.protein_position, vr.amino_acids,
            lambda h: PE.cross_check_hgvsp(edit, h))
        # A junction projection is kept only on an exact mutant-protein match.
        if junction and tier not in (SCOPE.CROSSCHECK_OK, SCOPE.CROSSCHECK_ALT_OK):
            return _reject(v.variant_id, v.submitted_id,
                           SCOPE.JUNCTION_UNRESOLVED, f"{tier}: {detail}",
                           feature_order)
        if tier == SCOPE.CROSSCHECK_ALT_OK:
            qc.append(SCOPE.CROSSCHECK_ALT_OK)
        elif tier == SCOPE.IDENTITY_NOT_ESTABLISHED:
            return _reject(v.variant_id, v.submitted_id,
                           SCOPE.IDENTITY_NOT_ESTABLISHED, detail, feature_order)
        elif tier != SCOPE.CROSSCHECK_OK:
            return _reject(v.variant_id, v.submitted_id, SCOPE.CROSSCHECK_FAILED,
                           detail, feature_order)
        if junction:
            qc.append(SCOPE.JUNCTION_PROJECTED)
    elif junction:
        return _reject(v.variant_id, v.submitted_id, SCOPE.JUNCTION_UNRESOLVED,
                       "a junction projection cannot be confirmed with VEP "
                       "disabled", feature_order)

    acc, verdict = FA.mane_to_uniprot(t.base_id)
    window = edit.structure_context_residues()

    feats: dict[str, object] = {}
    feats["indel_type"] = v.indel_type
    feats["indel_size_aa"] = v.size_aa
    feats["indel_size_relative"] = FS.indel_size_relative(v, edit)
    feats["dist_to_splice"] = FS.dist_to_splice(v, t)
    feats["cds_length"] = t.cds_length

    loc = FS.local_sequences(edit)
    feats["delta_local_hydropathy"] = FS.delta_local_hydropathy(loc)
    feats["delta_local_entropy"] = FS.delta_local_entropy(loc)
    feats["delta_local_charge"] = FS.delta_local_charge(loc)

    cons = FC.compute(v, edit, t, cmap)
    feats["phylop100_mean"] = cons.phylop_mean
    feats["gerp_mean"] = cons.gerp_mean
    if cons.qc:
        qc.append(f"conservation={cons.qc}")

    con = FK.compute(t.base_id)
    feats["pli"], feats["loeuf"], feats["mis_z"] = con.pli, con.loeuf, con.mis_z
    if not con.joined:
        qc.append("constraint=TRANSCRIPT_NOT_IN_TABLE")

    record = None
    if acc:
        record = FA.read_records({acc}).get(acc)
    ann = FA.compute(record, window)
    feats.update(ann.values)
    if ann.qc:
        qc.append(f"annotation={ann.qc}")

    feats["local_entropy"] = FS.local_entropy(
        edit, record.sequence if record else None)

    st = FT.compute(acc, verdict, window, v.indel_type, v.size_aa)
    feats.update(st.values)
    if st.qc:
        qc.append(f"structure={st.qc}")

    missing = [f for f in feature_order if f not in feats]
    if missing:
        raise RuntimeError(f"engine did not produce {missing}")
    extra = [f for f in feats if f not in feature_order]
    if extra:
        raise RuntimeError(f"engine produced non-contract columns: {extra}")

    return Row({"variant_id": v.variant_id, "submitted_id": v.submitted_id,
                "assembly": "GRCh38", "chrom": v.chrom, "pos": v.pos,
                "ref": v.ref, "alt": v.alt, "gene_symbol": t.gene,
                "transcript_id": t.id, "uniprot_acc": acc,
                "hgvs_p": hgvsp, "status": "OK", "qc": "; ".join(qc)},
               feats)
