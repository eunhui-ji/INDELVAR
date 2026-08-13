#!/usr/bin/env python3
"""Score variants by lookup or direct computation."""
from __future__ import annotations

import csv
import dataclasses
import hashlib
import os
import pathlib
import shutil
import subprocess
import tempfile

from . import evidence as EV
from ._engine import refs as REFS
# `engine.engine` is NOT imported here: it pulls in the direct-scoring extra.
from ._engine.variant import (Rejected, open_reference,
                              parse as parse_variant, validate_schema)
from .inputs import Submission

PKG = pathlib.Path(__file__).resolve().parent


def _artifact(key: str) -> pathlib.Path:
    return REFS.derived(key)


MODEL = _artifact("model")
MEDIANS = _artifact("imputation_medians")
CONTRACT = _artifact("feature_contract")
CUTOFFS = _artifact("evidence_cutoffs")
RUNTIME_VERSIONS = _artifact("runtime_versions")
R_WRAPPER = PKG / "score_model.R"


def default_db() -> pathlib.Path:
    env = os.environ.get("INDELVAR_DB")
    if env:
        return pathlib.Path(env).expanduser()
    return pathlib.Path.home() / ".indelvar" / "indelvar_precomputed.sqlite"

IDENTITY = ["submitted_id", "variant_id", "assembly", "chrom", "pos", "ref",
            "alt", "gene_symbol", "transcript_id", "uniprot_acc", "hgvs_p",
            "status", "qc"]
TAIL = ["n_features_imputed", "imputed_features", "indelvar_score",
        "computational_evidence", "computational_evidence_points",
        "scoring_source", "model_sha256", "feature_contract_sha256",
        "calibration_sha256"]

# These would read as a completed ACMG classification.
FORBIDDEN_COLUMNS = ("acmg_tier", "acmg_points")


MISSING_R = ("MISSING_R_RUNTIME: Rscript is not on PATH. Computing a score "
             "applies the R randomForest model; lookup does not need it.")
R_CONSUMER = "model scoring"


class PreflightError(RuntimeError):
    pass


def sha256(p: pathlib.Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@dataclasses.dataclass
class Context:
    calibration: EV.Calibration
    feature_order: list[str]
    model_sha: str
    contract_sha: str
    calibration_sha: str
    rscript: str | None = None      # None in a lookup-only context


def build_lookup_context() -> Context:
    for p in (MODEL, CONTRACT, CUTOFFS):
        if not p.exists():
            raise PreflightError(f"MISSING_ARTIFACT: {p.name}")
    return Context(calibration=EV.Calibration(CUTOFFS),
                   feature_order=REFS.contract_features(),
                   model_sha=sha256(MODEL), contract_sha=sha256(CONTRACT),
                   calibration_sha=sha256(CUTOFFS), rscript=None)


def build_direct_context() -> Context:
    ctx = build_lookup_context()
    if not R_WRAPPER.exists():
        raise PreflightError(f"MISSING_ARTIFACT: {R_WRAPPER.name}")
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise PreflightError(MISSING_R)
    ctx.rscript = rscript
    return ctx


DIRECT_PACKAGES = ("numpy", "pyBigWig", "pydssp", "freesasa", "torch", "einops")
DIRECT_EXTRA_HINT = ('pip install "indelvar[direct]"')


def missing_direct_packages() -> list[str]:
    import importlib.util
    return [m for m in DIRECT_PACKAGES if importlib.util.find_spec(m) is None]


# Inputs to the table producers, never opened at run time.
RUNTIME_UNUSED_REFERENCES = {"mane_summary", "uniprot_subcell"}


def r_requirements() -> list[tuple[str, str]]:
    """(package, version) the runtime table requires for scoring."""
    if not RUNTIME_VERSIONS.exists():
        return []
    with RUNTIME_VERSIONS.open(newline="") as fh:
        rows = [r for r in csv.DictReader(
            (l for l in fh if not l.startswith("#")), delimiter="\t")]
    return [(r["package"], r["version"]) for r in rows
            if r["consumed_by"] == R_CONSUMER]


def unsatisfied_r_packages(rscript: str) -> list[str]:
    """Declared R packages R cannot load, or loads at another version."""
    want = r_requirements()
    if not want:
        # An absent table is reported as a missing reference; a present one
        # that declares nothing to check would otherwise pass quietly.
        return ([] if not RUNTIME_VERSIONS.exists()
                else [f"runtime_versions.tsv names no {R_CONSUMER} package"])
    expr = "; ".join(
        f'if (!requireNamespace("{p}", quietly=TRUE) || '
        f'!identical(as.character(packageVersion("{p}")), '
        f'as.character(package_version("{v}")))) cat("{p}={v} ")'
        for p, v in want)
    try:
        r = subprocess.run([rscript, "-e", expr], capture_output=True,
                           text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return [f"{p}={v}" for p, v in want]
    return (r.stdout.split() if r.returncode == 0
            else [f"{p}={v}" for p, v in want])


def preflight_direct() -> list[str]:
    """Problems preventing direct scoring, as a list of messages."""
    problems = []
    absent = missing_direct_packages()
    if absent:
        problems.append(
            f"MISSING_DIRECT_DEPENDENCIES: {', '.join(absent)} are not "
            f"installed. Direct scoring needs them; lookup does not. Install "
            f"with: {DIRECT_EXTRA_HINT}")
    rscript = shutil.which("Rscript")
    if rscript is None:
        problems.append(MISSING_R)
    else:
        absent = unsatisfied_r_packages(rscript)
        if absent:
            problems.append(
                f"R_PACKAGES_UNSATISFIED: {', '.join(absent)}. These are the "
                f"versions runtime_versions.tsv requires, and they are not the "
                f"ones installed in the R that Rscript resolves to.")
    for p in (MODEL, MEDIANS, CONTRACT, CUTOFFS, R_WRAPPER):
        if not p.exists():
            problems.append(f"MISSING_ARTIFACT: {p.name}")

    # Derived from the reference table, so a new reference is checked too.
    for key in sorted(set(REFS.REFERENCES) - RUNTIME_UNUSED_REFERENCES):
        try:
            REFS.path(key)
        except REFS.ReferenceMissing as exc:
            problems.append(str(exc))
    try:
        ok = subprocess.run(["docker", "info"], capture_output=True,
                            timeout=60).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        ok = False
    if not ok:
        problems.append("VEP_RUNTIME_UNAVAILABLE: docker is not usable, so the "
                        f"pinned VEP image {REFS.VEP_IMAGE} cannot be run")
    return problems


def output_columns(ctx: Context) -> list[str]:
    cols = IDENTITY + list(ctx.feature_order) + TAIL
    bad = [c for c in cols if c in FORBIDDEN_COLUMNS]
    if bad:
        raise RuntimeError(f"forbidden output columns present: {bad}")
    return cols


def _blank(cols: list[str]) -> dict[str, str]:
    return {c: "" for c in cols}


class Normaliser:
    """Left-normalises a variant to the resource's key form."""

    def __init__(self) -> None:
        self.available = False
        self.reason = ""
        self._fa = self._cmap = None
        try:
            self._fa, self._cmap = open_reference()
            self.available = True
        except REFS.ReferenceMissing as exc:
            self.reason = str(exc)
        except OSError as exc:
            self.reason = f"{type(exc).__name__}: {exc}"

    def key(self, sub: Submission) -> tuple[str | None, tuple[str, str] | None]:
        """(canonical key, rejection). Both None means 'not normalisable here'."""
        # Schema first, in every mode.
        try:
            validate_schema(sub.chrom, sub.pos, sub.ref, sub.alt)
        except Rejected as exc:
            return None, (exc.code, exc.detail)
        if not self.available:
            return f"{sub.chrom}:{sub.pos}:{sub.ref}:{sub.alt}", None
        try:
            v = parse_variant(sub.chrom, sub.pos, sub.ref, sub.alt,
                              self._fa, self._cmap)
        except Rejected as exc:
            return None, (exc.code, exc.detail)
        except OSError:
            # a reference problem, not an input problem
            return f"{sub.chrom}:{sub.pos}:{sub.ref}:{sub.alt}", None
        return v.variant_id, None


def reject_row(sub: Submission, code: str, detail: str,
               cols: list[str]) -> dict[str, str]:
    """A rejection row: reason only, no features and no score."""
    row = _blank(cols)
    row.update(submitted_id=sub.submitted_id, variant_id="",
               assembly=sub.assembly or "", chrom=sub.chrom, pos=sub.pos,
               ref=sub.ref, alt=sub.alt, status=code, qc=detail,
               scoring_source="")
    return row


def run_engine(subs: list[Submission], cols: list[str],
               use_vep: bool = True) -> tuple[list[tuple[int, dict]],
                                              list[tuple[int, dict]]]:
    """(ok, rejected) as (submission index, row) pairs."""
    try:
        from ._engine import engine as ENG
    except ImportError as exc:
        absent = missing_direct_packages()
        raise PreflightError(
            f"MISSING_DIRECT_DEPENDENCIES: "
            f"{', '.join(absent) if absent else exc.name}. Direct scoring "
            f"needs them; lookup does not. Install with: "
            f"{DIRECT_EXTRA_HINT}") from None
    ok_rows: list[tuple[int, dict]] = []
    rejected: list[tuple[int, dict]] = []
    runnable: list[tuple[int, Submission]] = []
    for i, s in enumerate(subs):
        if s.reject is not None:
            rejected.append((i, reject_row(s, s.reject[0], s.reject[1], cols)))
            continue
        runnable.append((i, s))
    if not runnable:
        return ok_rows, rejected

    try:
        # Position in the outer list is the submission index.
        groups = ENG.run_grouped([(s.chrom, s.pos, s.ref, s.alt)
                                  for _, s in runnable], use_vep=use_vep)
    except REFS.ReferenceMissing as exc:
        raise PreflightError(str(exc)) from None

    if len(groups) != len(runnable):
        raise RuntimeError(
            f"engine returned {len(groups)} groups for {len(runnable)} inputs; "
            "the group-per-input contract was broken")

    for (i, sub), group in zip(runnable, groups):
        for r in group:
            ident = r.identity
            row = _blank(cols)
            row.update({k: ("" if ident.get(k) is None else str(ident.get(k)))
                        for k in IDENTITY if k in ident})
            row["submitted_id"] = sub.submitted_id
            row["assembly"] = "GRCh38"
            if ident["status"] != "OK":
                rejected.append((i, row))
                continue
            for f in r.features:
                v = r.features.get(f)
                row[f] = "" if v is None else str(v)
            ok_rows.append((i, row))
    return ok_rows, rejected


def score_rows(rows: list[dict], ctx: Context) -> None:
    """Adds score, imputation provenance and evidence in place."""
    if not rows:
        return
    with tempfile.TemporaryDirectory(prefix="indelvar_score_") as td:
        tmp = pathlib.Path(td)
        feat_p = tmp / "features.tsv"
        out_p = tmp / "scored.tsv"
        cols = ["variant_id"] + list(ctx.feature_order)
        with feat_p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            for i, r in enumerate(rows):
                rec = {c: r.get(c, "") for c in cols}
                # A duplicate input is legitimate; the join back is positional.
                rec["variant_id"] = f"row{i}"
                w.writerow(rec)
        proc = subprocess.run(
            [ctx.rscript, str(R_WRAPPER), "--features", str(feat_p),
             "--out", str(out_p), "--model", str(MODEL),
             "--contract", str(CONTRACT), "--medians", str(MEDIANS),
             "--runtime-versions", str(RUNTIME_VERSIONS)],
            capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "model scoring failed:\n" + (proc.stderr or proc.stdout))
        with out_p.open(newline="") as fh:
            scored = list(csv.DictReader(fh, delimiter="\t"))
    if len(scored) != len(rows):
        raise RuntimeError(f"scorer returned {len(scored)} rows for {len(rows)}")

    for i, (r, s) in enumerate(zip(rows, scored)):
        if s["variant_id"] != f"row{i}":
            raise RuntimeError(
                f"scorer returned rows out of order at position {i}: "
                f"expected row{i}, got {s['variant_id']}")
        score = float(s["indelvar_score"])
        indel_type = r.get("indel_type", "")
        label, points = ctx.calibration.assign(indel_type, score)
        r["indelvar_score"] = s["indelvar_score"]
        r["n_features_imputed"] = s["n_features_imputed"]
        r["imputed_features"] = s["imputed_features"]
        r["computational_evidence"] = label
        r["computational_evidence_points"] = str(points)
        r["model_sha256"] = ctx.model_sha
        r["feature_contract_sha256"] = ctx.contract_sha
        r["calibration_sha256"] = ctx.calibration_sha


def stamp_provenance(row: dict, ctx: Context, source: str) -> None:
    row["scoring_source"] = source
    row["model_sha256"] = ctx.model_sha
    row["feature_contract_sha256"] = ctx.contract_sha
    row["calibration_sha256"] = ctx.calibration_sha
