#!/usr/bin/env python3
"""Annotate variants with the pinned offline Ensembl VEP container."""
from __future__ import annotations

import dataclasses
import os
import pathlib
import subprocess

from . import refs


def _scratch() -> pathlib.Path:
    env = os.environ.get("INDELVAR_SCRATCH")
    if env:
        return pathlib.Path(env).expanduser()
    return pathlib.Path.home() / ".indelvar" / "engine_scratch"


SCRATCH = _scratch()
FIELDS = ("Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,"
          "Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,"
          "Codons,SYMBOL,MANE_SELECT,CANONICAL,HGVSc,HGVSp,STRAND")


@dataclasses.dataclass
class VepRow:
    feature: str
    consequence: set[str]
    hgvsp: str
    protein_position: str
    amino_acids: str


class VepFailed(RuntimeError):
    pass


def _docker_available() -> bool:
    try:
        return subprocess.run(["docker", "info"], capture_output=True,
                              timeout=30).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def annotate(variants: list[tuple[str, int, str, str, str]],
             tag: str = "run") -> dict[str, list[VepRow]]:
    """variants: (id, chrom, pos, ref, alt) in Ensembl naming; returns id -> rows."""
    if not variants:
        return {}
    if not _docker_available():
        raise VepFailed("DOCKER_UNAVAILABLE: the pinned VEP image cannot be run")
    SCRATCH.mkdir(parents=True, exist_ok=True)
    vcf = SCRATCH / f"{tag}.{os.getpid()}.vcf"
    out = SCRATCH / f"{tag}.{os.getpid()}.tsv"
    with vcf.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for vid, chrom, pos, ref, alt in variants:
            fh.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")

    cache = refs.path("vep_cache")
    fasta_dir = refs.path("grch38_fasta").parent
    cmd = [
        "docker", "run", "--rm", "--entrypoint", "vep",
        "-v", f"{cache}:/vepcache:ro",
        "-v", f"{fasta_dir}:/ref:ro",
        "-v", f"{SCRATCH}:/work",
        refs.VEP_IMAGE,
        "--offline", "--cache", "--dir_cache", "/vepcache",
        "--species", "homo_sapiens", "--assembly", "GRCh38",
        "--fasta", "/ref/GRCh38.fa",
        "--input_file", f"/work/{vcf.name}", "--output_file", f"/work/{out.name}",
        "--tab", "--fields", FIELDS,
        "--hgvs", "--mane", "--symbol", "--protein", "--check_ref",
        "--canonical", "--transcript_version", "--force_overwrite", "--no_stats",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if proc.returncode != 0:
        raise VepFailed(f"VEP exited {proc.returncode}: {proc.stderr[-2000:]}")

    rows: dict[str, list[VepRow]] = {}
    header: list[str] = []
    for line in out.read_text().splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line.lstrip("#").split("\t")
            continue
        f = dict(zip(header, line.split("\t")))
        rows.setdefault(f["Uploaded_variation"], []).append(VepRow(
            feature=f.get("Feature", ""),
            consequence=set((f.get("Consequence") or "").split(",")),
            hgvsp=f.get("HGVSp", ""),
            protein_position=f.get("Protein_position", ""),
            amino_acids=f.get("Amino_acids", "")))

    missing = [v[0] for v in variants if v[0] not in rows]
    if missing:
        raise VepFailed(
            "PROBE_COVERAGE_FAILED: VEP returned no row for "
            f"{len(missing)} of {len(variants)} submitted variants "
            f"(first: {missing[0]})")
    return rows


def transcript_rows(rows: list[VepRow], want_enst: str) -> list[VepRow]:
    base = want_enst.split(".")[0]
    return [r for r in rows if r.feature.split(".")[0] == base]


def pick_transcript_row(rows: list[VepRow], want_enst: str) -> VepRow | None:
    got = transcript_rows(rows, want_enst)
    return got[0] if got else None
