#!/usr/bin/env python3
"""The twelve UniProt protein-annotation features."""
from __future__ import annotations

import csv
import dataclasses
import gzip
import re

from . import refs

DOMAIN_TYPES = {"DOMAIN", "REGION", "ZN_FING", "COILED"}
PTM_TYPES = {"MOD_RES", "CARBOHYD", "LIPID", "CROSSLNK"}
ACTIVE_TYPES = {"ACT_SITE", "BINDING"}
WANTED_FT = DOMAIN_TYPES | PTM_TYPES | ACTIVE_TYPES | {"REPEAT", "DISULFID"}

CATEGORIES = ["nucleus", "cytoplasm", "membrane", "mitochondrion", "secreted",
              "other", "unknown"]
SUBCELL_FEATURES = [f"subcell_{c}" for c in CATEGORIES]

ECO = re.compile(r"\{[^}]*\}")
MOLECULE_FORM = re.compile(r"\[[^\]]*\]:\s*")
HYPHEN_WRAP = re.compile(r"-\s+")
FT_LINE = re.compile(r"^FT   (\w+)\s+([<>?\d]+)(?:\.\.([<>?\d]+))?\s*$")

PTM_PROXIMITY = 5


@dataclasses.dataclass
class UniProtRecord:
    accession: str
    sequence: str
    features: list[tuple[str, int, int]]
    subcell_terms: list[str]


@dataclasses.dataclass
class AnnotationResult:
    values: dict[str, object]
    qc: str


def _pos(tok: str) -> int | None:
    """Parse an FT endpoint. An uncertain bound (?, <, >) is not a coordinate."""
    tok = tok.strip()
    if not tok or not tok.isdigit():
        return None
    return int(tok)


_CACHE: dict[str, UniProtRecord] = {}


def read_records(wanted: set[str]) -> dict[str, UniProtRecord]:
    """UniProt records for a batch of accessions, cached per run."""
    missing = {a for a in wanted if a not in _CACHE}
    if not missing:
        return {a: _CACHE[a] for a in wanted if a in _CACHE}
    _CACHE.update(_read_uncached(missing))
    return {a: _CACHE[a] for a in wanted if a in _CACHE}


def _read_uncached(wanted: set[str]) -> dict[str, UniProtRecord]:
    out: dict[str, UniProtRecord] = {}
    acc: str | None = None
    keep = False
    feats: list[tuple[str, int, int]] = []
    seq: list[str] = []
    in_seq = False
    in_sub = False
    sub_block: list[str] = []
    sub_blocks: list[str] = []

    def flush():
        if keep and acc:
            terms: list[str] = []
            for b in sub_blocks:
                terms.extend(structured_terms(b))
            out[acc] = UniProtRecord(acc, "".join(seq), list(feats), terms)

    with gzip.open(refs.path("uniprot_human"), "rt") as fh:
        for line in fh:
            if line.startswith("ID   "):
                flush()
                acc, keep, feats, seq, in_seq = None, False, [], [], False
                in_sub, sub_block, sub_blocks = False, [], []
            elif line.startswith("AC   ") and acc is None:
                acc = line[5:].split(";")[0].strip()
                keep = acc in wanted
            elif not keep:
                continue
            elif line.startswith("CC   -!- SUBCELLULAR LOCATION:"):
                if sub_block:
                    sub_blocks.append(" ".join(sub_block))
                in_sub = True
                sub_block = [line[len("CC   -!- SUBCELLULAR LOCATION:"):].rstrip("\n")]
            elif in_sub and line.startswith("CC       "):
                sub_block.append(line[9:].rstrip("\n"))
            elif line.startswith("CC   ") or line.startswith("FT") or \
                    line.startswith("SQ   ") or line.startswith("//"):
                if in_sub:
                    sub_blocks.append(" ".join(sub_block))
                    in_sub, sub_block = False, []
                if line.startswith("SQ   "):
                    in_seq = True
                elif line.startswith("FT   "):
                    m = FT_LINE.match(line.rstrip("\n"))
                    if m and m.group(1) in WANTED_FT:
                        a = _pos(m.group(2))
                        b = _pos(m.group(3)) if m.group(3) else a
                        if a is not None and b is not None:
                            feats.append((m.group(1), a, b))
            elif in_seq and line.startswith("     "):
                seq.append(line.replace(" ", "").strip())
        flush()
    return out


def structured_terms(block: str) -> list[str]:
    """CV terms of one SUBCELLULAR LOCATION comment; free-text notes are not terms."""
    body = ECO.sub("", block)
    body = body.split("Note=")[0]
    body = MOLECULE_FORM.sub("", body)
    body = " ".join(body.split())
    out = []
    for sentence in body.split("."):
        if not sentence.strip():
            continue
        for part in sentence.split(";"):
            s = part.strip()
            if not s:
                continue
            term = s.split(",")[-1].strip()
            if term:
                out.append(term)
    return out


_CMAP: dict[str, list[str]] | None = None


def category_map() -> dict[str, list[str]]:
    global _CMAP
    if _CMAP is None:
        with refs.path("subcell_map").open(newline="") as fh:
            rows = list(csv.DictReader((l for l in fh if not l.startswith("#")),
                                       delimiter="\t"))
        _CMAP = {r["cv_term"].casefold():
                 (r["category"].split(";") if r["category"] != "other" else ["other"])
                 for r in rows}
    return _CMAP


_MAP: dict[str, tuple[str, str]] | None = None


def mane_to_uniprot(enst_base: str) -> tuple[str | None, str]:
    """(accession, verdict). Only EXACT_MATCH may reach an AlphaFold model."""
    global _MAP
    if _MAP is None:
        m: dict[str, tuple[str, str]] = {}
        with refs.path("mane_uniprot_map").open(newline="") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")),
                                    delimiter="\t"):
                m[r["enst"]] = (r["accession"] or None, r["verdict"])
        _MAP = m
    return _MAP.get(enst_base, (None, "NO_UNIPROT_MAPPING"))


def _intersects(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 <= b1 and b0 <= a1


def compute(record: UniProtRecord | None, window: list[int]) -> AnnotationResult:
    keys = ["n_annotated_regions", "in_repeat", "ptm_count_near",
            "disulfide_nearby", "in_active_site"] + SUBCELL_FEATURES
    if record is None or not window:
        return AnnotationResult({k: None for k in keys},
                                "NO_UNIPROT_RECORD" if record is None else "NO_WINDOW")

    w0, w1 = min(window), max(window)
    p0, p1 = w0 - PTM_PROXIMITY, w1 + PTM_PROXIMITY

    n_regions = 0
    in_repeat = 0
    ptm_count = 0
    disulfide = 0
    active = 0
    for kind, a, b in record.features:
        if kind in DOMAIN_TYPES and _intersects(a, b, w0, w1):
            n_regions += 1
        elif kind == "REPEAT" and _intersects(a, b, w0, w1):
            in_repeat = 1
        elif kind in PTM_TYPES and _intersects(a, b, p0, p1):
            ptm_count += 1
        elif kind == "DISULFID":
            # A disulfide is two residues, not the span between them.
            if p0 <= a <= p1 or p0 <= b <= p1:
                disulfide = 1
        elif kind in ACTIVE_TYPES and _intersects(a, b, w0, w1):
            active = 1

    vals: dict[str, object] = {
        "n_annotated_regions": n_regions,
        "in_repeat": in_repeat,
        "ptm_count_near": ptm_count,
        "disulfide_nearby": disulfide,
        "in_active_site": active,
    }

    cmap = category_map()
    sub = {f: 0 for f in SUBCELL_FEATURES}
    if not record.subcell_terms:
        sub["subcell_unknown"] = 1
    else:
        for t in record.subcell_terms:
            cats = cmap.get(t.casefold()) or cmap.get(HYPHEN_WRAP.sub("-", t).casefold())
            if cats is None:
                return AnnotationResult({k: None for k in keys},
                                        f"SUBCELL_TERM_UNMAPPED: {t!r}")
            for c in cats:
                sub[f"subcell_{c}"] = 1
    vals.update(sub)
    return AnnotationResult(vals, "")
