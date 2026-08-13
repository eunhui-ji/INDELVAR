#!/usr/bin/env python3
"""Calibrated PP3/BP4 cutoffs and the computational-evidence assignment."""
from __future__ import annotations

import csv
import dataclasses
import pathlib
import struct

NO_DISTINCT_BAND = "NO_DISTINCT_BAND_SAME_AS_STRONGER_TIER"

# Abderrazzaq et al.: 1 Supporting, 2 and 3 Moderate, 4 Strong.
STRENGTH = {1: "Supporting", 2: "Moderate", 3: "Moderate", 4: "Strong"}

INDETERMINATE = "Indeterminate"


def _hex_to_double(h: str) -> float:
    return struct.unpack(">d", bytes.fromhex(h))[0]


@dataclasses.dataclass(frozen=True)
class Cutoff:
    indel_type: str
    direction: str
    point: int
    threshold: float
    assignable: bool


class Calibration:
    """The sixteen calibrated tiers, loaded once."""

    def __init__(self, path: pathlib.Path):
        self.path = path
        self.cutoffs: list[Cutoff] = []
        with path.open(newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                dec = float(r["discounted_threshold"])
                exact = _hex_to_double(r["full_precision_discounted_threshold"])
                if dec != exact:
                    raise ValueError(
                        f"calibration table is inconsistent for {r['arm']} "
                        f"{r['point']}: decimal {dec!r} does not equal the "
                        f"recorded full-precision value {exact!r}")
                self.cutoffs.append(Cutoff(
                    # Calibration schema field.
                    indel_type=r["arm"], direction=r["direction"],
                    point=int(r["point"]), threshold=exact,
                    assignable=r["operational_band_status"] != NO_DISTINCT_BAND))
        if len(self.cutoffs) != 16:
            raise ValueError(f"expected 16 calibrated tiers, found {len(self.cutoffs)}")
        self._by = {(c.indel_type, c.direction, c.point): c for c in self.cutoffs}

    def indel_types(self) -> set[str]:
        return {c.indel_type for c in self.cutoffs}

    def _thr(self, indel_type: str, direction: str, point: int) -> float | None:
        c = self._by.get((indel_type, direction, point))
        if c is None or not c.assignable:
            return None
        return c.threshold

    def assign(self, indel_type: str, score: float) -> tuple[str, int]:
        """The strongest tier the score reaches for its indel type, else Indeterminate."""
        if indel_type not in self.indel_types():
            raise ValueError(f"no calibration for indel_type {indel_type!r}")
        for p in (4, 3, 2, 1):
            t = self._thr(indel_type, "pathogenic", p)
            if t is not None and score >= t:
                return f"PP3_{STRENGTH[p]}", p
        for p in (4, 3, 2, 1):
            t = self._thr(indel_type, "benign", -p)
            if t is not None and score <= t:
                return f"BP4_{STRENGTH[p]}", -p
        return INDETERMINATE, 0

    def summary(self) -> list[dict]:
        rows = []
        for c in sorted(self.cutoffs, key=lambda c: (c.indel_type, c.point)):
            label, _ = ((f"PP3_{STRENGTH[abs(c.point)]}", c.point)
                        if c.direction == "pathogenic"
                        else (f"BP4_{STRENGTH[abs(c.point)]}", c.point))
            rows.append({"indel_type": c.indel_type, "evidence": label,
                         "points": c.point, "cutoff": c.threshold,
                         "assignable": c.assignable})
        return rows
