#!/usr/bin/env python3
"""pli, loeuf, and mis_z, keyed on the MANE transcript."""
from __future__ import annotations

import csv
import dataclasses

from . import refs

PLI_COL = "lof.pLI"
LOEUF_COL = "lof.oe_ci.upper"
MISZ_COL = "mis.z_score"
NULLS = {"", "NA", "NaN", "nan", "None", "."}


@dataclasses.dataclass
class ConstraintResult:
    pli: float | None
    loeuf: float | None
    mis_z: float | None
    joined: bool


_TABLE: dict[str, dict[str, str]] | None = None


def _load() -> dict[str, dict[str, str]]:
    global _TABLE
    if _TABLE is None:
        out: dict[str, dict[str, str]] = {}
        with refs.path("gnomad_constraint").open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("mane_select") != "true":
                    continue
                tx = (r.get("transcript") or "").strip()
                if not tx.startswith("ENST"):
                    continue
                out[tx.split(".")[0]] = r
        _TABLE = out
    return _TABLE


def _num(row: dict[str, str], col: str) -> float | None:
    v = (row.get(col) or "").strip()
    if v in NULLS:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def compute(transcript_base_id: str) -> ConstraintResult:
    row = _load().get(transcript_base_id)
    if row is None:
        return ConstraintResult(None, None, None, joined=False)
    return ConstraintResult(_num(row, PLI_COL), _num(row, LOEUF_COL),
                            _num(row, MISZ_COL), joined=True)
