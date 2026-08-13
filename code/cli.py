#!/usr/bin/env python3
"""Command line interface for scoring GRCh38 in-frame indels."""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import re
import sys

from . import evidence as EV
from . import inputs as IN
from . import lookup as LK
from . import scoring as SC
from importlib.metadata import version as _installed_version

__version__ = _installed_version("indelvar")

DEFAULT_DB = SC.default_db()

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_PREFLIGHT = 3
EXIT_ALL_REJECTED = 4

SCHEMA_VERSION = "1"
GENOME_BUILD = "GRCh38"
DIGEST = re.compile(r"[0-9a-f]{64}")
# What a resource this install will trust states about itself. Each is checked
# for presence, for a single non-empty value, and against what it should say.
# An unknown extra key is left alone.
REQUIRED_PROVENANCE = ("model_sha256", "feature_contract_sha256",
                       "calibration_sha256", "schema_version",
                       "genome_build", "feature_count", "rows")


def _metadata_problems(ctx: SC.Context, meta_rows: list[tuple],
                       n_rows: int) -> list[str]:
    keys = [k for k, _ in meta_rows]
    meta = dict(meta_rows)
    bad: list[str] = []
    unusable: set[str] = set()
    for k in REQUIRED_PROVENANCE:
        if k not in meta:
            bad.append(f"{k} absent")
        elif keys.count(k) > 1:
            bad.append(f"{k} appears {keys.count(k)} times")
        elif meta[k] is None or not str(meta[k]).strip():
            bad.append(f"{k} empty")
        else:
            continue
        unusable.add(k)

    for k in ("model_sha256", "feature_contract_sha256", "calibration_sha256"):
        if k not in unusable and not DIGEST.fullmatch(str(meta[k])):
            bad.append(f"{k} is not a 64-character digest")
    want = {"schema_version": SCHEMA_VERSION, "genome_build": GENOME_BUILD,
            "feature_count": str(len(ctx.feature_order))}
    for k, expected in want.items():
        if k not in unusable and str(meta[k]) != expected:
            bad.append(f"{k} is {meta[k]!r}, not {expected!r}")
    if "rows" not in unusable:
        try:
            stated = int(str(meta["rows"]))
        except ValueError:
            bad.append(f"rows is {meta['rows']!r}, not an integer")
        else:
            if stated != n_rows:
                bad.append(f"rows states {stated}, the table holds {n_rows}")
    return bad


def _write(rows: list[dict], cols: list[str], dest: str) -> None:
    fh = sys.stdout if dest == "-" else open(dest, "w", newline="")
    try:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="raise")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    finally:
        if fh is not sys.stdout:
            fh.close()


def _na_columns(ctx: SC.Context) -> set[str]:
    """Columns whose absent value is written `NA` rather than an empty field."""
    return set(ctx.feature_order) | {"igsr_AC", "igsr_AN", "igsr_AF"}


def _from_lookup(hits: list[dict], sub: IN.Submission, cols: list[str],
                 ctx: SC.Context) -> list[dict]:
    na_cols = _na_columns(ctx)
    out = []
    for h in hits:
        row = {c: "" for c in cols}
        for c in cols:
            if c not in h:
                continue
            v = h[c]
            if c in na_cols and LK.is_missing(v):
                row[c] = LK.MISSING_TSV
            else:
                row[c] = "" if v is None else str(v)
        row["submitted_id"] = sub.submitted_id
        SC.stamp_provenance(row, ctx, "precomputed")
        out.append(row)
    return out


def cmd_score(a: argparse.Namespace) -> int:
    try:
        subs = IN.read_any(a.input, a.variant, a.assembly)
    except IN.InputError as exc:
        print(f"INPUT_ERROR: {exc}", file=sys.stderr)
        return EXIT_USAGE

    try:
        ctx = SC.build_lookup_context()
    except SC.PreflightError as exc:
        print(f"PREFLIGHT_FAILED: {exc}", file=sys.stderr)
        return EXIT_PREFLIGHT
    cols = SC.output_columns(ctx)

    db = LK.open_if_present(None if a.mode == "direct" else a.db)
    if a.mode == "lookup" and db is None:
        print(f"PREFLIGHT_FAILED: MISSING_PRECOMPUTED_DB: {a.db.name} is not "
              "present, and --mode lookup will not fall back to direct scoring. "
              "Set INDELVAR_DB or pass --db.", file=sys.stderr)
        return EXIT_PREFLIGHT

    if db is not None:
        prov, source = db.provenance()
        mine = {"model_sha256": ctx.model_sha,
                "feature_contract_sha256": ctx.contract_sha,
                "calibration_sha256": ctx.calibration_sha}
        if prov is None:
            print(f"PREFLIGHT_FAILED: DB_ARTIFACT_PROVENANCE_UNVERIFIABLE: "
                  f"{a.db.name} states no model, contract or calibration "
                  "digest, and none was found in a schema.json or manifest.tsv "
                  "beside it. Nothing was written.", file=sys.stderr)
            return EXIT_PREFLIGHT
        differ = {k: (prov[k], mine[k]) for k in mine if prov[k] != mine[k]}
        if differ:
            print("PREFLIGHT_FAILED: DB_ARTIFACT_PROVENANCE_MISMATCH: the "
                  "precomputed resource was built with artifacts this install "
                  f"does not have (provenance read from {source}):",
                  file=sys.stderr)
            for k, (theirs, ours) in sorted(differ.items()):
                print(f"  - {k}: resource {theirs[:16]}... vs installed "
                      f"{ours[:16]}...", file=sys.stderr)
            print("  Nothing was written.", file=sys.stderr)
            return EXIT_PREFLIGHT

        # A digest match does not prove the columns are still there.
        absent = [c for c in ctx.feature_order if c not in db.columns]
        if absent:
            print("PREFLIGHT_FAILED: DB_SCHEMA_INCOMPLETE: the precomputed "
                  f"resource is missing {len(absent)} column(s) the feature "
                  f"contract requires: {', '.join(absent)}. Nothing was "
                  "written.", file=sys.stderr)
            return EXIT_PREFLIGHT

    n_dup = IN.count_duplicate_ids(subs)

    norm = SC.Normaliser() if db is not None else None
    if norm is not None and not norm.available:
        print("WARNING: NORMALISATION_UNAVAILABLE: the GRCh38 reference FASTA "
              "could not be opened, so variants are looked up exactly as "
              "submitted, without left-alignment. A non-canonical spelling of a "
              "variant that is present will be reported as absent.",
              file=sys.stderr)
        if norm.reason:
            print(f"         reason: {norm.reason}", file=sys.stderr)
    # (submission index, row), to restore input order.
    ordered: list[tuple[int, dict]] = []
    need_direct: list[tuple[int, IN.Submission]] = []
    for i, s in enumerate(subs):
        if s.reject is not None:
            ordered.append((i, SC.reject_row(s, s.reject[0], s.reject[1], cols)))
            continue
        hits = []
        if db is not None:
            key, bad = norm.key(s)
            if bad is not None:
                ordered.append((i, SC.reject_row(s, bad[0], bad[1], cols)))
                continue
            hits = db.get(key)
        if hits:
            ordered.extend((i, r) for r in _from_lookup(hits, s, cols, ctx))
        elif a.mode == "lookup":
            detail = "no exact normalised match; --mode lookup does not compute"
            if norm is not None and not norm.available:
                detail += "; the submitted spelling was used as the key"
            ordered.append((i, SC.reject_row(s, "NOT_IN_PRECOMPUTED_DB",
                                             detail, cols)))
        else:
            need_direct.append((i, s))

    if need_direct:
        # R is needed only when a variant must be computed.
        problems = SC.preflight_direct()
        if problems:
            print("PREFLIGHT_FAILED: direct scoring needs a runtime that is "
                  "not available:", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            print(f"  {len(need_direct)} of {len(subs)} variant(s) were not in "
                  "the precomputed resource and would have to be computed. "
                  "Nothing was written.", file=sys.stderr)
            return EXIT_PREFLIGHT
        try:
            ctx = SC.build_direct_context()
        except SC.PreflightError as exc:
            print(f"PREFLIGHT_FAILED: {exc}", file=sys.stderr)
            return EXIT_PREFLIGHT
        try:
            ok, rejected = SC.run_engine([s for _, s in need_direct], cols,
                                         use_vep=not a.no_vep)
        except SC.PreflightError as exc:
            print(f"PREFLIGHT_FAILED: {exc}", file=sys.stderr)
            return EXIT_PREFLIGHT
        SC.score_rows([r for _, r in ok], ctx)
        for _, r in ok:
            SC.stamp_provenance(r, ctx, "computed")
        # By index: submitted_id is neither unique nor one-to-one with rows.
        for j, r in ok + rejected:
            ordered.append((need_direct[j][0], r))

    # Stable sort, so the engine's row order is kept.
    ordered.sort(key=lambda p: p[0])
    rows = [r for _, r in ordered]

    # A scored row writes NA for a missing value; a rejected row stays blank.
    na_cols = _na_columns(ctx)
    for r in rows:
        if r["status"] == "OK":
            for c in na_cols & r.keys():
                if LK.is_missing(r[c]):
                    r[c] = LK.MISSING_TSV

    feature_cols = set(ctx.feature_order)
    for r in rows:
        if r["status"] != "OK":
            dirty = [c for c in feature_cols | {"indelvar_score",
                                                "computational_evidence_points"}
                     if r.get(c)]
            if dirty:
                raise RuntimeError(
                    f"rejected row {r['submitted_id']} carries {sorted(dirty)}")

    _write(rows, cols, a.output)
    n_ok = sum(1 for r in rows if r["status"] == "OK")
    pre = sum(1 for r in rows if r.get("scoring_source") == "precomputed")
    comp = sum(1 for r in rows if r.get("scoring_source") == "computed")
    print(f"rows={len(rows)} scored={n_ok} rejected={len(rows) - n_ok} "
          f"(precomputed={pre} computed={comp})", file=sys.stderr)
    if n_dup:
        print(f"note: {n_dup} submitted_id value(s) appear more than once; "
              "every submission still gets its own row", file=sys.stderr)
    if rows and n_ok == 0:
        print("every submitted variant was rejected; see the status column",
              file=sys.stderr)
        return EXIT_ALL_REJECTED
    return EXIT_OK


def cmd_setup(a: argparse.Namespace) -> int:
    try:
        ctx = SC.build_lookup_context()
        lookup_ok, lookup_note = True, ctx.model_sha
    except SC.PreflightError as exc:
        ctx, lookup_ok, lookup_note = None, False, str(exc)

    db = a.db if a.db.exists() else None
    db_prov, db_prov_src, db_prov_state = None, "", "no resource"
    if db is not None and lookup_ok:
        try:
            lk = LK.Lookup(db)
            db_prov, db_prov_src = lk.provenance()
            db_cols, db_meta = list(lk.columns), lk.meta()
            db_rows = lk.row_count()
            lk.close()
        except Exception as exc:                              # noqa: BLE001
            db_prov_state = f"UNREADABLE: {exc}"
        else:
            if db_prov is None:
                db_prov_state = "UNVERIFIABLE (no digest with the resource)"
            elif (db_prov["model_sha256"] == ctx.model_sha
                  and db_prov["feature_contract_sha256"] == ctx.contract_sha
                  and db_prov["calibration_sha256"] == ctx.calibration_sha):
                db_prov_state = f"MATCHES this install (from {db_prov_src})"
            else:
                db_prov_state = f"MISMATCH - built with other artifacts " \
                                f"(from {db_prov_src})"

    # A digest match says which artifacts built the resource, not that the
    # resource is still complete or its metadata still well formed. Both are
    # conditions `score` refuses on, so readiness applies them too.
    if db_prov_state.startswith("MATCHES"):
        absent = [c for c in ctx.feature_order if c not in db_cols]
        bad = _metadata_problems(ctx, db_meta, db_rows)
        if absent:
            db_prov_state = (f"INCOMPLETE - {len(absent)} column(s) the feature "
                             f"contract requires are absent")
        elif bad:
            db_prov_state = f"MALFORMED - {'; '.join(bad)}"

    problems = SC.preflight_direct()
    norm = SC.Normaliser()
    lookup_ready = lookup_ok and (db is None or db_prov_state.startswith("MATCHES"))
    direct_ready = lookup_ok and not problems

    print("INDELVAR setup check")
    print(f"  version                 : {__version__}")
    print(f"  model artifact          : {'OK ' + lookup_note[:16] if lookup_ok else 'MISSING'}")
    print()
    print(f"  LOOKUP readiness        : {'READY' if lookup_ready else 'NOT READY'}"
          "   (no R, no references, no Docker required)")
    print(f"      precomputed resource: "
          f"{db.name if db else 'absent - set INDELVAR_DB or pass --db'}")
    print(f"      resource provenance : {db_prov_state}")
    print(f"      variant normalisation: "
          f"{'full (reference FASTA readable)' if norm.available else 'LIMITED - the submitted spelling is used as the key'}")
    if not norm.available and norm.reason:
        print(f"        - {norm.reason}")
    if lookup_ok:
        print(f"      calibrated tiers    : {len(ctx.calibration.cutoffs)} "
              f"({sum(1 for c in ctx.calibration.cutoffs if c.assignable)} assignable)")
    print()
    print(f"  DIRECT readiness        : {'READY' if direct_ready else 'NOT READY'}"
          f"   (needs R + {' + '.join(p for p, _ in SC.r_requirements())}, "
          "Docker + the pinned VEP image, and the full reference set)")
    for p in problems:
        print(f"      - {p}")

    if a.json:
        print(json.dumps({
            "version": __version__,
            "package_path": str(pathlib.Path(SC.__file__).parent),
            "model_ok": lookup_ok,
            "lookup_ready": lookup_ready,
            "direct_ready": direct_ready,
            "precomputed_db": str(db) if db else None,
            "db_provenance_state": db_prov_state,
            "db_provenance_source": db_prov_src,
            "direct_runtime_problems": problems,
            "normalisation_available": norm.available,
            "normalisation_reason": norm.reason,
        }, indent=2))
    return EXIT_OK if lookup_ok else EXIT_PREFLIGHT


def cmd_cutoffs(a: argparse.Namespace) -> int:
    cal = EV.Calibration(SC.CUTOFFS)
    w = csv.DictWriter(sys.stdout,
                       fieldnames=["indel_type", "evidence", "points",
                                   "cutoff", "assignable"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in cal.summary():
        w.writerow(r)
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="indelvar", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version",
                    version=f"indelvar {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    s = sub.add_parser("score", help="score GRCh38 in-frame indels")
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--input", type=pathlib.Path,
                   help="VCF 4.2+, TSV or CSV carrying CHROM/POS/REF/ALT")
    g.add_argument("--variant", metavar="CHROM:POS:REF:ALT",
                   help="a single variant")
    s.add_argument("--output", default="-", help="output TSV, or - for stdout")
    s.add_argument("--mode", choices=["auto", "lookup", "direct"], default="auto",
                   help="auto: exact precomputed lookup first, then direct. "
                        "lookup: precomputed only. direct: always recompute")
    s.add_argument("--assembly", default="GRCh38",
                   help="assembly of the input when the file does not declare "
                        "one; only GRCh38 is accepted")
    s.add_argument("--db", type=pathlib.Path, default=DEFAULT_DB,
                   help="precomputed SQLite index")
    s.add_argument("--no-vep", action="store_true", help=argparse.SUPPRESS)
    s.set_defaults(func=cmd_score)

    t = sub.add_parser("setup", help="check references and runtimes")
    t.add_argument("--check", action="store_true", default=True,
                   help="report lookup and direct readiness")
    t.add_argument("--json", action="store_true", help="also print JSON")
    t.add_argument("--db", type=pathlib.Path, default=DEFAULT_DB,
                   help="precomputed SQLite index")
    t.set_defaults(func=cmd_setup)

    c = sub.add_parser("cutoffs", help="print the calibrated evidence tiers")
    c.set_defaults(func=cmd_cutoffs)
    return ap


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
