#!/usr/bin/env python3
"""Registry of the reference files the engine reads."""
from __future__ import annotations

import os
import pathlib
from importlib import resources


PACKAGED_DATA = pathlib.Path(str(resources.files(__package__) / "_data"))


# External references come from INDELVAR_REFERENCE_DIR and nowhere else.
REFERENCE_ROOT = (pathlib.Path(os.environ["INDELVAR_REFERENCE_DIR"]).expanduser()
                  if os.environ.get("INDELVAR_REFERENCE_DIR") else None)
RAW = REFERENCE_ROOT / "INDELVAR_refs" if REFERENCE_ROOT else None
GRCH38_DIR = REFERENCE_ROOT / "GRCh38" if REFERENCE_ROOT else None


def _ext(base: pathlib.Path | None, *parts: str) -> pathlib.Path | None:
    return base.joinpath(*parts) if base is not None else None


# The packaged tables, by the name the engine knows each one by.
DERIVED = {
    "model": "indelvar_rf.rds",
    "imputation_medians": "imputation_medians.tsv",
    "feature_contract": "feature_contract.tsv",
    "evidence_cutoffs": "evidence_cutoffs.tsv",
    "reference_versions": "reference_versions.tsv",
    "residue_scales": "residue_scales.tsv",
    "tien_maxasa": "tien2013_maxasa.tsv",
    "subcell_map": "subcell_category_map.tsv",
    "contig_map": "grch38_contig_map.tsv",
    "mane_uniprot_map": "mane_uniprot_map.tsv",
    "alphafold_fragment_spans": "alphafold_fragment_spans.tsv",
    "alphafold_member_index": "alphafold_member_index.tsv",
    "runtime_versions": "runtime_versions.tsv",
}


def derived(key: str) -> pathlib.Path:
    try:
        return PACKAGED_DATA / DERIVED[key]
    except KeyError:
        raise SystemExit(f"UNDECLARED_PACKAGED_TABLE:{key}") from None


def _pins() -> dict[str, str]:
    out = {}
    for line in derived("reference_versions").read_text().splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.split("\t")
        if f[0] == "reference":
            continue
        out[f[0]] = f[1]
    return out


PINS = _pins()

# reference key -> (path, pin name)
REFERENCES: dict[str, tuple[pathlib.Path, str]] = {
    "grch38_fasta": (_ext(GRCH38_DIR, "GRCh38.fa"), "grch38_fasta"),
    "grch38_fai": (_ext(GRCH38_DIR, "GRCh38.fa.fai"), "grch38_fasta"),
    "mane_gff": (_ext(RAW, "mane_v1.5", "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"), "mane"),
    "mane_protein": (_ext(RAW, "mane_v1.5", "MANE.GRCh38.v1.5.ensembl_protein.faa.gz"), "mane"),
    "mane_summary": (_ext(RAW, "mane_v1.5", "MANE.GRCh38.v1.5.summary.txt.gz"), "mane"),
    "uniprot_human": (_ext(RAW, "uniprot_2026_02", "uniprot_sprot_human.dat.gz"), "uniprot"),
    "uniprot_subcell": (_ext(RAW, "uniprot_2026_02", "subcell_2026_02.txt"), "uniprot"),
    "gnomad_constraint": (_ext(RAW, "gnomad", "gnomad.v4.0.constraint_metrics.tsv"), "gnomad_constraint"),
    "phylop": (_ext(RAW, "conservation", "hg38.phyloP100way.bw"), "phylop"),
    "gerp": (_ext(RAW, "conservation", "gerp_conservation_scores.homo_sapiens.GRCh38.r111.bw"), "gerp"),
    "alphafold_tar": (_ext(RAW, "alphafold_v6", "UP000005640_9606_HUMAN_v6.tar"), "alphafold"),
    "vep_cache": (_ext(RAW, "vep_116", "cache"), "vep"),
    # derived tables, packaged with the engine
    "residue_scales": (derived("residue_scales"), "n/a"),
    "contig_map": (derived("contig_map"), "n/a"),
    "subcell_map": (derived("subcell_map"), "uniprot"),
    "tien_maxasa": (derived("tien_maxasa"), "n/a"),
    "feature_contract": (derived("feature_contract"), "n/a"),
    "mane_uniprot_map": (derived("mane_uniprot_map"), "mane"),
    "alphafold_member_index": (derived("alphafold_member_index"), "alphafold"),
    "alphafold_fragment_spans": (derived("alphafold_fragment_spans"), "alphafold"),
    "runtime_versions": (derived("runtime_versions"), "n/a"),
    # `_pins()` reads this file directly, before this table exists.
    "reference_versions": (derived("reference_versions"), "n/a"),
}

VEP_IMAGE = ("ensemblorg/ensembl-vep@sha256:"
             "f354dd8d09073e4d943acbbd02f5eb234a9d9e9d444371c1c349910f2123de11")


class ReferenceMissing(RuntimeError):
    pass


def contract_features() -> list[str]:
    """Feature ids in contract order."""
    import csv
    with path("feature_contract").open(newline="") as fh:
        return [r["feature_id"] for r in csv.DictReader(fh, delimiter="\t")]


def _unset(key: str) -> ReferenceMissing:
    return ReferenceMissing(
        f"REFERENCE_ROOT_NOT_SET:{key} needs an external reference set. "
        f"Set INDELVAR_REFERENCE_DIR to the directory holding them.")


def path(key: str) -> pathlib.Path:
    """Resolve a reference file; raises when it is missing."""
    try:
        p, pin = REFERENCES[key]
    except KeyError:
        raise ReferenceMissing(f"UNKNOWN_REFERENCE:{key}") from None
    if p is None:
        raise _unset(key)
    if not p.exists():
        rel = PINS.get(pin, pin)
        where = ("packaged with indelvar" if PACKAGED_DATA in p.parents
                 else "under $INDELVAR_REFERENCE_DIR")
        raise ReferenceMissing(
            f"REFERENCE_MISSING:{key} pin={pin}={rel} file={p.name} "
            f"expected {where}")
    return p
