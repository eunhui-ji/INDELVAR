#!/usr/bin/env python3
"""Look up precomputed variant scores."""
from __future__ import annotations

import json
import pathlib
import sqlite3

PROVENANCE_KEYS = ("model_sha256", "feature_contract_sha256",
                   "calibration_sha256")

# The resource stores SQL NULL; a table writes `NA`.
MISSING_TSV = "NA"


def is_missing(v) -> bool:
    return v is None or v == "" or v == MISSING_TSV


class Lookup:
    def __init__(self, db: pathlib.Path):
        self.path = db
        self.conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        self.conn.row_factory = sqlite3.Row
        cols = [r[1] for r in self.conn.execute("PRAGMA table_info(variants)")]
        if not cols:
            raise ValueError(f"{db} has no `variants` table")
        if "variant_id" not in cols or "transcript_id" not in cols:
            raise ValueError(f"{db} `variants` lacks the primary-key columns")
        self.columns = cols

    def provenance(self) -> tuple[dict[str, str] | None, str]:
        """Return the resource provenance and its source."""
        try:
            rows = self.conn.execute(
                "SELECT key, value FROM provenance").fetchall()
            got = {r["key"]: r["value"] for r in rows}
            if all(k in got for k in PROVENANCE_KEYS):
                return {k: got[k] for k in PROVENANCE_KEYS}, "the `provenance` table"
        except sqlite3.Error:
            pass
        schema = self.path.parent / "schema.json"
        if schema.exists():
            try:
                d = json.loads(schema.read_text())
                if all(k in d for k in PROVENANCE_KEYS):
                    return {k: d[k] for k in PROVENANCE_KEYS}, str(schema)
            except (ValueError, OSError):
                pass
        man = self.path.parent / "manifest.tsv"
        if man.exists():
            try:
                d = dict(l.split("\t", 1) for l in man.read_text().splitlines()
                         if "\t" in l)
                if all(k in d for k in PROVENANCE_KEYS):
                    return {k: d[k] for k in PROVENANCE_KEYS}, str(man)
            except OSError:
                pass
        return None, ""

    def meta(self) -> list[tuple]:
        """Every provenance row, in table order, so a duplicate key is visible."""
        try:
            return [(r["key"], r["value"]) for r in
                    self.conn.execute("SELECT key, value FROM provenance")]
        except sqlite3.Error:
            return []

    def row_count(self) -> int:
        return self.conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]

    def get(self, variant_id: str) -> list[dict]:
        rows = self.conn.execute(
            "SELECT * FROM variants WHERE variant_id = ? ORDER BY transcript_id",
            (variant_id,)).fetchall()
        return [dict(r) for r in rows]

    def close(self) -> None:
        self.conn.close()


def open_if_present(db: pathlib.Path | None) -> Lookup | None:
    if db is None or not db.exists():
        return None
    return Lookup(db)
