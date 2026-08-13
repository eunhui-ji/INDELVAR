#!/usr/bin/env python3
"""Per-residue numeric scales and the metrics computed from them."""
from __future__ import annotations

import csv
import math

from . import refs
from .variant import Rejected

HYDROPATHY_SCALE = "kd_hydropathy"


class ChargeDefinitionNotPinned(Rejected):
    def __init__(self, detail: str = ""):
        super().__init__("CHARGE_DEFINITION_NOT_PINNED", detail)


_TABLE: dict[str, dict[str, float]] | None = None


def _load() -> dict[str, dict[str, float]]:
    global _TABLE
    if _TABLE is None:
        out: dict[str, dict[str, float]] = {}
        with refs.path("residue_scales").open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.setdefault(r["scale_id"], {})[r["residue_1"]] = float(r["value"])
        for sid, m in out.items():
            if len(m) != 20:
                raise Rejected("SCALE_INCOMPLETE", f"{sid} has {len(m)} residues, expected 20")
        _TABLE = out
    return _TABLE


def scale(scale_id: str) -> dict[str, float]:
    t = _load()
    try:
        return t[scale_id]
    except KeyError:
        raise Rejected("UNKNOWN_SCALE", f"{scale_id}; known: {sorted(t)}") from None


def pinned_charge_scale_id() -> str:
    for line in refs.path("reference_versions").read_text().splitlines():
        if line.startswith("charge_scale\t"):
            return line.split("\t")[1].strip()
    raise ChargeDefinitionNotPinned(
        "the reference versions table has no charge_scale row")


def charge_scale() -> dict[str, float]:
    return scale(pinned_charge_scale_id())


def mean_over(seq: str, values: dict[str, float]) -> float:
    """Unweighted arithmetic mean of a residue scale over a sequence."""
    if not seq:
        raise Rejected("EMPTY_LOCAL_SEQUENCE", "cannot take a mean over zero residues")
    missing = sorted(set(seq) - set(values))
    if missing:
        raise Rejected("NON_STANDARD_RESIDUE", f"{missing} not in the scale")
    return sum(values[c] for c in seq) / len(seq)


def shannon_bits(seq: str) -> float:
    """Shannon entropy in bits over the 20 standard residues."""
    if not seq:
        raise Rejected("EMPTY_LOCAL_SEQUENCE", "cannot take an entropy over zero residues")
    n = len(seq)
    counts: dict[str, int] = {}
    for c in seq:
        counts[c] = counts.get(c, 0) + 1
    h = 0.0
    for k in counts.values():
        p = k / n
        h -= p * math.log2(p)
    return h
