#!/usr/bin/env python3
"""The eleven structure features from one AlphaFold fragment."""
from __future__ import annotations

import csv
import dataclasses
import gzip
import importlib.metadata
import importlib.util
import math
import os
import pathlib
import tempfile

import numpy as np

from . import refs

CONTACT_CUTOFF = 8.0
HBOND_ENERGY_CUTOFF = -0.5

# pydssp's __init__ imports torch; the numpy backend is loaded on its own.
C3_ALPHABET = np.array(["-", "H", "E"])
_PYDSSP_NUMPY = None


def _pinned_pydssp_version() -> str:
    """The pydssp version pinned in runtime_versions.tsv."""
    lock = refs.path("runtime_versions")
    for line in lock.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if f[0] == "pydssp":
            return f[1]
    raise RuntimeError(f"pydssp is not pinned in {lock}")


def pydssp_numpy():
    """pydssp's numpy backend, loaded without the torch-importing parent."""
    global _PYDSSP_NUMPY
    if _PYDSSP_NUMPY is not None:
        return _PYDSSP_NUMPY
    want = _pinned_pydssp_version()
    got = importlib.metadata.version("pydssp")
    if got != want:
        raise RuntimeError(f"pydssp {got} is installed; the pinned version is {want}")
    # find_spec on a top-level name locates the package without executing it.
    spec = importlib.util.find_spec("pydssp")
    if spec is None or not spec.origin:
        raise RuntimeError("pydssp is not installed")
    path = pathlib.Path(spec.origin).parent / "pydssp_numpy.py"
    sub = importlib.util.spec_from_file_location("indelvar_pydssp_numpy", path)
    mod = importlib.util.module_from_spec(sub)
    sub.loader.exec_module(mod)
    _PYDSSP_NUMPY = mod
    return mod


def assign_c3(bb: np.ndarray, donor: np.ndarray) -> str:
    onehot = pydssp_numpy().assign(bb, donor_mask=donor)
    return "".join(C3_ALPHABET[np.argmax(onehot, axis=-1)])
BURIED_RSA = 0.20
# Residues with positive Kyte-Doolittle hydropathy.
HYDROPHOBIC = set("ACILMFV")
FLANK = 5

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

STRUCTURE_FEATURES = [
    "af_plddt_mean", "af_helix_frac", "af_sheet_frac", "af_rsa_mean",
    "wcn_mean", "contact_density", "bridging_contacts", "ss_element_break",
    "ss_element_frac_affected",
    "hbond_density", "buried_hydrophobic_fraction",
]


@dataclasses.dataclass
class StructureResult:
    values: dict[str, object]
    qc: str

    @staticmethod
    def absent(qc: str) -> "StructureResult":
        return StructureResult({k: None for k in STRUCTURE_FEATURES}, qc)


_INDEX: dict[tuple[str, int], tuple[int, int, bool]] | None = None
_SPANS: dict[str, list[tuple[int, int, int]]] | None = None


def member_index() -> dict[tuple[str, int], tuple[int, int, bool]]:
    global _INDEX
    if _INDEX is None:
        out = {}
        p = refs.path("alphafold_member_index")
        with p.open(newline="") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")),
                                    delimiter="\t"):
                if r["kind"] != "pdb":
                    continue
                out[(r["accession"], int(r["fragment"]))] = (
                    int(r["offset"]), int(r["size"]), r["gzipped"] == "yes")
        _INDEX = out
    return _INDEX


def fragment_spans(accession: str) -> list[tuple[int, int, int]]:
    """(fragment, begin, end) for an accession, in global UniProt numbering."""
    global _SPANS
    if _SPANS is None:
        by: dict[str, list[tuple[int, int, int]]] = {}
        with refs.path("alphafold_fragment_spans").open(newline="") as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")),
                                    delimiter="\t"):
                by.setdefault(r["accession"], []).append(
                    (int(r["fragment"]), int(r["begin"]), int(r["end"])))
        _SPANS = {a: sorted(v) for a, v in by.items()}
    # The table lists only multi-fragment accessions, so an empty list is F1.
    return _SPANS.get(accession, [])


def select_fragment(w_start: int, w_end: int,
                    frags: list[tuple[int, int, int]]) -> int | None:
    """The fragment containing the whole window, farthest from either edge."""
    cands = [(n, b, e) for n, b, e in frags if b <= w_start and w_end <= e]
    if not cands:
        return None
    return max(cands, key=lambda t: (min(w_start - t[1], t[2] - w_end), -t[0]))[0]


def read_model(accession: str, fragment: int) -> bytes:
    key = (accession, fragment)
    idx = member_index()
    if key not in idx:
        raise KeyError(f"NO_ALPHAFOLD_MODEL:{accession}-F{fragment}")
    offset, size, gz = idx[key]
    with refs.path("alphafold_tar").open("rb") as fh:
        fh.seek(offset)
        blob = fh.read(size)
    return gzip.decompress(blob) if gz else blob


def parse_pdb(text: str):
    """(resnos, resnames1, backbone[N,4,3], plddt[N]); pLDDT is the CA B-factor."""
    order = {"N": 0, "CA": 1, "C": 2, "O": 3}
    res: dict[int, dict] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip()
        if name not in order:
            continue
        num = int(line[22:26])
        r = res.setdefault(num, {"name": line[17:20].strip(),
                                 "xyz": np.full((4, 3), np.nan),
                                 "plddt": np.nan})
        r["xyz"][order[name]] = (float(line[30:38]), float(line[38:46]),
                                 float(line[46:54]))
        if name == "CA":
            r["plddt"] = float(line[60:66])
    nums = sorted(res)
    bb = np.stack([res[n]["xyz"] for n in nums])
    names = [THREE_TO_ONE.get(res[n]["name"], "X") for n in nums]
    plddt = np.array([res[n]["plddt"] for n in nums])
    return np.array(nums), names, bb, plddt


def local_mask(n: int) -> np.ndarray:
    m = ~np.eye(n, dtype=bool)
    m *= ~np.diag(np.ones(n - 1, dtype=bool), k=-1)
    m *= ~np.diag(np.ones(n - 2, dtype=bool), k=-2)
    return m


def maxasa() -> dict[str, float]:
    rows = [l.rstrip("\n") for l in refs.path("tien_maxasa").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    return {r["three_letter"]: float(r["theoretical_allowed"])
            for r in csv.DictReader(rows, delimiter="\t")}


def _elements(c3: str) -> list[tuple[int, int, str]]:
    """Maximal runs of H or E, as 0-based inclusive [a, b] with the state."""
    out, i = [], 0
    while i < len(c3):
        if c3[i] in "HE":
            j = i
            while j + 1 < len(c3) and c3[j + 1] == c3[i]:
                j += 1
            out.append((i, j, c3[i]))
            i = j + 1
        else:
            i += 1
    return out


def compute(accession: str | None, verdict: str, window_global: list[int],
            kind: str, indel_size_aa: int) -> StructureResult:
    if accession is None:
        return StructureResult.absent("NO_UNIPROT_ACCESSION")
    if verdict != "EXACT_MATCH":
        # AlphaFold models the UniProt canonical sequence.
        return StructureResult.absent(f"UNIPROT_SEQUENCE_{verdict}")
    if not window_global:
        return StructureResult.absent("NO_WINDOW")

    w0, w1 = min(window_global), max(window_global)
    frags = fragment_spans(accession)
    if frags:
        fnum = select_fragment(w0, w1, frags)
        if fnum is None:
            return StructureResult.absent("NO_FRAGMENT_CONTAINS_WINDOW")
        begin = next(b for n, b, _ in frags if n == fnum)
    else:
        fnum, begin = 1, 1

    try:
        text = read_model(accession, fnum).decode()
    except KeyError as exc:
        return StructureResult.absent(str(exc))
    nums, names, bb, plddt = parse_pdb(text)
    n = len(nums)
    if n == 0:
        return StructureResult.absent("EMPTY_MODEL")

    # Local index of a global residue. Fragment residues restart at 1.
    def local(g: int) -> int | None:
        i = g - begin
        return i if 0 <= i < n else None

    widx = [local(g) for g in window_global]
    if any(i is None for i in widx):
        return StructureResult.absent("WINDOW_OUTSIDE_FRAGMENT")
    widx = [i for i in widx if i is not None]

    donor = np.ones(n, dtype=float)
    donor[[i for i, a in enumerate(names) if a == "P"]] = 0.0

    if np.isnan(bb).any(axis=(1, 2)).all():
        return StructureResult.absent("NO_COMPLETE_BACKBONE")

    c3 = assign_c3(bb, donor)
    e = pydssp_numpy().get_hbond_map(bb, return_e=True)
    e = e[0] if e.ndim == 3 else e
    B = (e < HBOND_ENERGY_CUTOFF) & local_mask(n) & (donor[:, None] > 0)

    ca = bb[:, 1, :]
    d = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)

    # FreeSASA needs a file.
    import freesasa
    freesasa.setVerbosity(freesasa.silent)
    fd, tmp = tempfile.mkstemp(suffix=".pdb")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        params = freesasa.Parameters({"algorithm": freesasa.LeeRichards,
                                      "probe-radius": 1.4, "n-slices": 20})
        result = freesasa.calc(freesasa.Structure(tmp), params)
        areas = result.residueAreas()
        chain = next(iter(areas))
        areas = areas[chain]
    finally:
        os.unlink(tmp)

    mx = maxasa()
    one_to_three = {v: k for k, v in THREE_TO_ONE.items()}
    rsa = np.full(n, np.nan)
    for i, num in enumerate(nums):
        ra = areas.get(str(num))
        three = one_to_three.get(names[i])
        if ra is not None and three in mx:
            rsa[i] = ra.total / mx[three]

    v: dict[str, object] = {}

    good_plddt = [plddt[i] for i in widx if not math.isnan(plddt[i])]
    v["af_plddt_mean"] = float(np.mean(good_plddt)) if good_plddt else None

    states = [c3[i] for i in widx]
    v["af_helix_frac"] = states.count("H") / len(states)
    v["af_sheet_frac"] = states.count("E") / len(states)

    interior = [i for i in widx if 0 < i < n - 1]
    rs = [rsa[i] for i in interior if not math.isnan(rsa[i])]
    v["af_rsa_mean"] = float(np.mean(rs)) if rs else None

    ok_ca = [i for i in widx if not np.isnan(ca[i]).any()]
    if ok_ca:
        wcn = [float(np.sum(1.0 / np.square(d[i][np.isfinite(d[i])]))) for i in ok_ca]
        v["wcn_mean"] = float(np.mean(wcn))
        v["contact_density"] = float(np.mean([int((d[i] < CONTACT_CUTOFF).sum())
                                              for i in ok_ca]))
    else:
        v["wcn_mean"] = v["contact_density"] = None

    nflank = [local(g) for g in range(w0 - FLANK, w0)]
    cflank = [local(g) for g in range(w1 + 1, w1 + 1 + FLANK)]
    if any(i is None for i in nflank + cflank):
        v["bridging_contacts"] = None
    else:
        v["bridging_contacts"] = int(sum(
            1 for i in widx
            if (d[i][nflank] < CONTACT_CUTOFF).any()
            and (d[i][cflank] < CONTACT_CUTOFF).any()))

    els = _elements(c3)
    lo, hi = min(widx), max(widx)
    containing = [(a, b, s) for a, b, s in els if a < lo and b > hi]
    v["ss_element_break"] = 1 if containing else 0
    touching = [(a, b, s) for a, b, s in els if a <= hi and lo <= b]
    if len(touching) > 1:
        v["ss_element_frac_affected"] = None
    elif len(touching) == 1:
        a, b, _ = touching[0]
        overlap = min(hi, b) - max(lo, a) + 1
        v["ss_element_frac_affected"] = overlap / (b - a + 1)
    else:
        v["ss_element_frac_affected"] = 0.0

    hb = (B.sum(axis=1) + B.sum(axis=0)).astype(float)
    v["hbond_density"] = float(np.mean([hb[i] for i in widx]))

    denom = [i for i in widx if not math.isnan(rsa[i])]
    if denom:
        v["buried_hydrophobic_fraction"] = sum(
            1 for i in denom if rsa[i] < BURIED_RSA and names[i] in HYDROPHOBIC
        ) / len(denom)
    else:
        v["buried_hydrophobic_fraction"] = None

    return StructureResult(v, "")
