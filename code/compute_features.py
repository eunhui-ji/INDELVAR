#!/usr/bin/env python3
"""Compute the 39 INDELVAR features from cached per-protein/per-gene inputs (under
$INDELVAR_ROOT/data); shared by training, the database build, and release scoring.

Usage: compute_features.py <variants.{tsv,parquet}> <out.tsv> [--cache DIR] [--no-conservation]
       compute_features.py --selftest [RAW_COHORT.tsv] [--n N]   # RAW cohort, not train_set.tsv
"""
from __future__ import annotations

import gzip
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Constants
CONTACT_CUTOFF = 8.0
HELIX_TURN = 3.6
HYDROPHOBIC = set("AVILMFWP")
RSA_BURIED_THR = 0.20
SSE_HELIX, SSE_SHEET = "H", "E"
PAD_ENTROPY = 10
PTM_TYPES = {"Modified residue", "Glycosylation", "Lipidation", "Cross-link"}
DOMAIN_TYPES = {"Domain", "Region", "Zinc finger", "Coiled coil"}
THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
AA20 = sorted("ARNDCQEGHILKMFPSTWYV")
AA_IDX = {a: i for i, a in enumerate(AA20)}
SUBCELL_MAP = {
    "Cytoplasm": ["Cytoplasm", "Cytosol", "Cytoskeleton"],
    "Membrane": ["Cell membrane", "Membrane", "Plasma membrane",
                 "Endoplasmic reticulum membrane"],
    "Mitochondrion": ["Mitochondrion"],
    "Nucleus": ["Nucleus", "Nuclear", "Chromosome"],
    "Secreted": ["Secreted", "Extracellular"],
    "Other": ["Endoplasmic reticulum", "Golgi", "Lysosome", "Peroxisome",
              "Vacuole", "Endosome", "Vesicle", "Cilium"],
}
SUBCELL_COLS = [f"subcell_{k}" for k in list(SUBCELL_MAP) + ["Unknown"]]

FEATURE_COLS = [
    "pli", "loeuf", "mis_z", "cds_length",
    "af_plddt_mean", "af_helix_frac", "af_sheet_frac", "af_rsa_mean",
    "wcn_mean", "contact_density", "bridging_contacts",
    "ss_element_break", "ss_element_frac_affected", "helix_register_preserved",
    "hbond_density", "hydrophobic_exposure_risk", "local_ss_class",
    "in_domain", "n_domains", "in_repeat", "ptm_nearby", "ptm_count_near",
    "disulfide_nearby", "in_active_site", "low_plddt_region",
    "indel_size_aa", "indel_type", "indel_size_relative",
    "phylop100_mean", "gerp_mean",
] + SUBCELL_COLS + ["local_entropy", "dist_to_splice"]


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# Per-variant structural feature math
def get_indices(resnos, ps, pe):
    return np.where((resnos >= ps) & (resnos <= pe))[0]


def compute_wcn(ca, idxs):
    if len(idxs) == 0:
        return np.nan
    vals = []
    for i in idxs:
        if np.isnan(ca[i]).any():
            continue
        d2 = np.sum((ca - ca[i]) ** 2, axis=1)
        d2[i] = np.inf
        valid = ~np.isnan(d2) & (d2 > 0)
        if valid.sum() < 1:
            continue
        vals.append(np.sum(1.0 / d2[valid]))
    return float(np.mean(vals)) if vals else np.nan


def compute_contact_density(ca, idxs):
    if len(idxs) == 0:
        return np.nan
    cnts = []
    for i in idxs:
        if np.isnan(ca[i]).any():
            continue
        d = np.linalg.norm(ca - ca[i], axis=1)
        d[i] = np.inf
        cnts.append(int(np.sum(d < CONTACT_CUTOFF)))
    return float(np.mean(cnts)) if cnts else np.nan


def compute_bridging_contacts(ca, idxs, N):
    if len(idxs) == 0:
        return np.nan
    n_flank = np.arange(max(0, idxs[0] - 5), idxs[0])
    c_flank = np.arange(idxs[-1] + 1, min(N, idxs[-1] + 6))
    if len(n_flank) == 0 or len(c_flank) == 0:
        return 0.0
    cnt = 0
    for i in idxs:
        if np.isnan(ca[i]).any():
            continue
        d = np.linalg.norm(ca - ca[i], axis=1)
        if np.any(d[n_flank] < CONTACT_CUTOFF) and np.any(d[c_flank] < CONTACT_CUTOFF):
            cnt += 1
    return float(cnt)


def find_ss_element_span(sse, idx):
    if idx < 0 or idx >= len(sse):
        return None
    s = sse[idx]
    if s not in (SSE_HELIX, SSE_SHEET):
        return None
    a = idx
    while a > 0 and sse[a - 1] == s:
        a -= 1
    b = idx
    while b < len(sse) - 1 and sse[b + 1] == s:
        b += 1
    return (a, b)


def compute_ss_disruption(sse, idxs):
    if len(idxs) == 0:
        return (np.nan, np.nan)
    center = idxs[len(idxs) // 2]
    span = find_ss_element_span(sse, center)
    if span is None:
        return (0, 0.0)
    a, b = span
    elt_len = b - a + 1
    is_break = int((a < idxs[0]) and (b > idxs[-1]))
    intersect = np.sum((idxs >= a) & (idxs <= b))
    frac = float(intersect) / float(elt_len) if elt_len > 0 else np.nan
    return (is_break, frac)


def compute_helix_register(sse, idxs, indel_type, aa_change):
    """|cos(2π·k/3.6)| for helix-internal deletions; NaN for insertions and non-helix windows."""
    if indel_type != "deletion" or len(idxs) == 0:
        return np.nan
    center = idxs[len(idxs) // 2]
    if sse[center] != SSE_HELIX:
        return np.nan
    return float(abs(math.cos(2 * math.pi * aa_change / HELIX_TURN)))


def compute_hbond_density(hbond_map, idxs):
    if len(idxs) == 0 or hbond_map.size == 0:
        return np.nan
    sums = hbond_map[idxs].sum(axis=1) + hbond_map[:, idxs].sum(axis=0)
    return float(np.mean(sums))


def compute_hydrophobic_exposure_risk(resnames, rsa, idxs):
    if len(idxs) == 0:
        return np.nan
    n_risk = 0
    for i in idxs:
        aa1 = THREE2ONE.get(str(resnames[i]).upper(), "X")
        if np.isnan(rsa[i]):
            continue
        if (aa1 in HYDROPHOBIC) and (rsa[i] < RSA_BURIED_THR):
            n_risk += 1
    return float(n_risk) / float(len(idxs))


def compute_local_ss_class(sse, idxs):
    if len(idxs) == 0:
        return np.nan
    chars = sse[idxs]
    n_e = int(np.sum(chars == SSE_SHEET))
    n_h = int(np.sum(chars == SSE_HELIX))
    n_o = len(chars) - n_e - n_h
    if n_e >= n_h and n_e >= n_o:
        return 2.0
    if n_h >= n_o:
        return 1.0
    return 0.0


def _wmean(arr, idxs):
    if len(idxs) == 0:
        return np.nan
    v = np.asarray(arr, dtype=np.float64)[idxs]
    v = v[~np.isnan(v)]
    return float(v.mean()) if v.size else np.nan


# Cached reference data
class Refs:
    """Lazy per-protein/per-gene cache holder."""

    def __init__(self, cache_dir, with_conservation=True):
        c = Path(cache_dir)
        self.npz_dir = c / "features" / "novel_protein_cache"
        self.json_dir = c / "processed" / "uniprot_features"
        self.seq_dir = c / "processed" / "uniprot_seq"
        self.af_dir = c / "raw" / "alphafold_v6_bulk"
        self.gene_lookup = self._load_gene_lookup(c / "precomputed" / "_gene_lookup.csv")
        self.gene_sites = None          # built on first dist_to_splice use
        self.mane_gff = c / "raw" / "mane_gff" / "MANE.GRCh38.v1.5.ensembl_genomic.gff.gz"
        self._prot, self._json, self._seq = {}, {}, {}
        self.with_conservation = with_conservation
        self._bw = None
        self.mane_summary = c / "raw" / "MANE_summary.txt.gz"   # gene -> MANE transcript
        self._txm = None            # {transcript: cds model} (lazy)
        self._g2tx = None           # gene_symbol -> MANE Ensembl transcript id

    @staticmethod
    def _load_gene_lookup(p):
        if not Path(p).exists():
            return {}
        df = pd.read_csv(p)
        return {r.gene_symbol: (r.pli, r.loeuf, r.mis_z, r.cds_length)
                for r in df.itertuples(index=False)}

    def protein(self, uniprot):
        if uniprot in self._prot:
            return self._prot[uniprot]
        npz = self.npz_dir / f"{uniprot}.npz"
        prot = None
        if npz.exists():
            z = np.load(npz, allow_pickle=True)
            n = len(z["resnos"])
            prot = dict(resnos=z["resnos"], resnames=z["resnames"],
                        ca_xyz=z["ca_xyz"], sse=z["sse"], rsa=z["rsa"],
                        hbond_map=z["hbond_map"],
                        plddt=self._plddt(uniprot, n,
                                          z["plddt"] if "plddt" in z.files else None))
        self._prot[uniprot] = prot
        return prot

    def _plddt(self, uniprot, n, npz_plddt):
        """pLDDT from the tier3 parse cache, else the npz array, else NaN."""
        pkl = self.npz_dir.parent / "tier3_protein_cache" / f"{uniprot}.pkl"
        if pkl.exists():
            try:
                import pickle
                pl = np.asarray(pickle.load(open(pkl, "rb")).get("plddt"), dtype=np.float32)
                if pl.size == n:
                    return pl
            except Exception:
                pass
        if npz_plddt is not None and np.isfinite(np.asarray(npz_plddt, dtype=np.float64)).any():
            return npz_plddt
        return np.full(n, np.nan, np.float32)

    def uniprot_json(self, uniprot):
        if uniprot in self._json:
            return self._json[uniprot]
        p = self.json_dir / f"{uniprot}.json"
        obj = None
        if p.exists() and p.stat().st_size >= 100:
            try:
                obj = json.load(open(p))
            except Exception:
                obj = None
        self._json[uniprot] = obj
        return obj

    def sequence(self, acc):
        if acc in self._seq:
            return self._seq[acc]
        s = None
        for ext in (".fasta", ".fa"):
            fp = self.seq_dir / f"{acc}{ext}"
            if fp.exists():
                s = "".join(l.strip() for l in open(fp) if not l.startswith(">")) or None
                if s:
                    break
        if s is None:
            fp = self.af_dir / f"AF-{acc}-F1-model_v6.pdb.gz"
            if fp.exists():
                sq = {}
                with gzip.open(fp, "rt") as fh:
                    for line in fh:
                        if line.startswith("ATOM") and line[12:16].strip() == "CA":
                            try:
                                sq[int(line[22:26])] = THREE2ONE.get(line[17:20].strip(), "X")
                            except Exception:
                                pass
                if sq:
                    s = "".join(sq.get(i, "X") for i in range(1, max(sq) + 1))
        self._seq[acc] = s
        return s

    def splice_sites(self, gene):
        if self.gene_sites is None:
            self.gene_sites = self._build_gene_sites()
        return self.gene_sites.get(gene)

    def _build_gene_sites(self):
        sites = {}
        if not Path(self.mane_gff).exists():
            return sites
        gene_exons = {}
        with gzip.open(self.mane_gff, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "exon":
                    continue
                gname = next((kv[len("gene_name="):] for kv in p[8].split(";")
                              if kv.startswith("gene_name=")), None)
                if gname is None:
                    continue
                gene_exons.setdefault(gname, (p[0], []))[1].append((int(p[3]), int(p[4])))
        for g, (chrom, exons) in gene_exons.items():
            exons = sorted(set(exons))
            gmin = min(e[0] for e in exons)
            gmax = max(e[1] for e in exons)
            ss = []
            for s, e in exons:
                if s != gmin:
                    ss.append(s)
                if e != gmax:
                    ss.append(e)
            sites[g] = (chrom, np.array(sorted(set(ss)), dtype=np.int64))
        return sites

    def txmodel(self, transcript=None, gene=None):
        """CDS genomic-position model for a MANE transcript (insertion conservation).
        Resolves by transcript id or gene; None if unresolved (caller falls back to the
        genomic-span query)."""
        if self._txm is None:
            self._txm, self._g2tx = self._build_txmodels()
        cands = []
        for t in (transcript, self._g2tx.get(gene)):
            if t:
                cands += [t, t.split(".")[0]]
        for t in cands:
            m = self._txm.get(t)
            if m:
                return m
        return None

    def _build_txmodels(self):
        """{transcript: {strand, chrom, cds_pos}} + {gene: MANE transcript}, parsed from
        the MANE genomic GFF CDS features."""
        tx_cds, tx_strand, tx_chrom = {}, {}, {}
        if Path(self.mane_gff).exists():
            with gzip.open(self.mane_gff, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 9 or f[2] != "CDS":
                        continue
                    tid = next((kv[14:] for kv in f[8].split(";")
                                if kv.startswith("transcript_id=")), None)
                    if tid is None:
                        continue
                    tx_cds.setdefault(tid, []).append((int(f[3]), int(f[4])))
                    tx_strand[tid] = f[6]; tx_chrom[tid] = f[0]
        models = {}
        for tid, exons in tx_cds.items():
            strand = tx_strand[tid]; exons = sorted(exons)
            if strand == "-":
                exons = exons[::-1]
            cds = []
            for cs, ce in exons:
                cds.extend(range(cs, ce + 1) if strand == "+" else range(ce, cs - 1, -1))
            models[tid] = {"strand": strand, "chrom": tx_chrom[tid], "cds_pos": cds}
        g2tx = {}
        if Path(self.mane_summary).exists():
            with gzip.open(self.mane_summary, "rt") as fh:
                h = fh.readline().rstrip("\n").split("\t")
                gi, ti = h.index("symbol"), h.index("Ensembl_nuc")
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    g2tx.setdefault(p[gi], p[ti])
        return models, g2tx

    def bigwigs(self):
        """phyloP100 (UCSC) + GERP (Ensembl) bigWig handles; remote by default,
        local copies via env vars. (None, None) if unavailable -> NaN."""
        if not self.with_conservation:
            return None
        if self._bw is None:
            try:
                import pyBigWig
                ph = os.environ.get("INDELVAR_PHYLOP",
                    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw")
                ge = os.environ.get("INDELVAR_GERP",
                    "http://ftp.ensembl.org/pub/release-111/compara/conservation_scores/"
                    "91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw")
                self._bw = (pyBigWig.open(ph), pyBigWig.open(ge))
            except Exception as e:
                log(f"conservation bigWig unavailable ({e}); phylop/gerp -> NaN")
                self._bw = (None, None)
        return self._bw


# Feature-group assemblers
def _structural(prot, ps, pe, indel_type, size):
    if prot is None:
        return {k: np.nan for k in (
            "af_plddt_mean", "af_helix_frac", "af_sheet_frac", "af_rsa_mean",
            "wcn_mean", "contact_density", "bridging_contacts",
            "ss_element_break", "ss_element_frac_affected", "helix_register_preserved",
            "hbond_density", "hydrophobic_exposure_risk", "local_ss_class",
            "low_plddt_region")}
    resnos, sse, rsa, ca = prot["resnos"], prot["sse"], prot["rsa"], prot["ca_xyz"]
    N = len(resnos)
    idxs = get_indices(resnos, ps, pe)
    ssbreak, fraclost = compute_ss_disruption(sse, idxs)
    plddt_mean = _wmean(prot["plddt"], idxs)
    helix_isH = np.asarray([1.0 if s == SSE_HELIX else 0.0 for s in sse])
    sheet_isE = np.asarray([1.0 if s == SSE_SHEET else 0.0 for s in sse])
    # ss_element_frac_affected/local_ss_class: del + ins; helix_register_preserved: del only.
    return dict(
        af_plddt_mean=plddt_mean,
        af_helix_frac=_wmean(helix_isH, idxs),
        af_sheet_frac=_wmean(sheet_isE, idxs),
        af_rsa_mean=_wmean(rsa, idxs),
        wcn_mean=compute_wcn(ca, idxs),
        contact_density=compute_contact_density(ca, idxs),
        bridging_contacts=compute_bridging_contacts(ca, idxs, N),
        ss_element_break=ssbreak,
        ss_element_frac_affected=fraclost,
        helix_register_preserved=compute_helix_register(sse, idxs, indel_type, size),
        hbond_density=compute_hbond_density(prot["hbond_map"], idxs),
        hydrophobic_exposure_risk=compute_hydrophobic_exposure_risk(prot["resnames"], rsa, idxs),
        local_ss_class=compute_local_ss_class(sse, idxs),
        low_plddt_region=int(plddt_mean < 50) if not np.isnan(plddt_mean) else 0,
    )


def _annotation(obj, ps, pe):
    """Tier-4 binaries aggregated over the indel window [ps, pe] (1-based)."""
    out = dict(in_domain=0, n_domains=0, in_repeat=0, ptm_nearby=0,
               ptm_count_near=0, disulfide_nearby=0, in_active_site=0)
    if obj is None:
        return out
    W = 5
    for f in obj.get("features", []):
        t = f.get("type")
        loc = f.get("location", {})
        try:
            s = int(loc["start"]["value"]); e = int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        overlap = (s <= pe and e >= ps)                             # domain/active-site overlap
        near = (ps - W <= s <= pe + W) or (ps - W <= e <= pe + W)   # PTM/disulfide within +/-W
        if t in DOMAIN_TYPES:
            if overlap:
                out["in_domain"] = 1
                out["n_domains"] += 1
        elif t == "Repeat":
            if overlap:
                out["in_repeat"] = 1
        elif t in PTM_TYPES:
            if near:
                out["ptm_nearby"] = 1
                out["ptm_count_near"] += 1
        elif t == "Disulfide bond":
            if near:
                out["disulfide_nearby"] = 1
        elif t in ("Active site", "Binding site"):
            if overlap:
                out["in_active_site"] = 1
    return out


def _classify_subcell(text):
    """Single-label, priority-ordered classification from the first reported location string."""
    if not text:
        return "Unknown"
    t = str(text).lower()
    if "nucle" in t:                       return "Nucleus"
    if "membrane" in t:                    return "Membrane"
    if "secret" in t:                      return "Secreted"
    if "mitochondr" in t:                  return "Mitochondrion"
    if "cytoplasm" in t or "cytosol" in t: return "Cytoplasm"
    if "endoplasm" in t:                   return "Membrane"
    if "golgi" in t:                       return "Membrane"
    if "lysosome" in t:                    return "Membrane"
    if "peroxisome" in t:                  return "Membrane"
    if "extracell" in t:                   return "Secreted"
    return "Other"


def _subcell(obj):
    res = {c: 0 for c in SUBCELL_COLS}
    val = None
    if obj is not None:
        for c in obj.get("comments", []):
            if c.get("commentType") == "SUBCELLULAR LOCATION":
                sl = c.get("subcellularLocations", [])
                if sl:                                  # first location of the first comment
                    val = sl[0].get("location", {}).get("value")
                    break
    res[f"subcell_{_classify_subcell(val)}"] = 1
    return res


def _entropy(seq, ps, pe):
    if seq is None:
        return np.nan
    L = len(seq)
    a = max(1, ps - PAD_ENTROPY)
    b = min(L, pe + PAD_ENTROPY)
    cnt = np.zeros(20, dtype=np.int64)
    for c in seq[a - 1:b]:
        j = AA_IDX.get(c)
        if j is not None:
            cnt[j] += 1
    n = cnt.sum()
    if n < 3:
        return np.nan
    p = cnt[cnt > 0] / n
    return float(round(-(p * np.log2(p)).sum(), 4))


def _dist_to_splice(refs, gene, chrom, pos):
    rec = refs.splice_sites(gene)
    if rec is None or pos is None or (isinstance(pos, float) and np.isnan(pos)):
        return np.nan
    gchrom, sites = rec
    if sites.size == 0 or ("chr" + str(chrom)) != gchrom:
        return np.nan
    i = int(np.clip(np.searchsorted(sites, pos), 0, len(sites) - 1))
    il = max(i - 1, 0)
    return float(min(abs(sites[i] - pos), abs(sites[il] - pos)))


def _conservation(refs, chrom, pos, ref, alt):
    """phyloP100 / GERP mean over the variant's reference nucleotides.
    phyloP uses UCSC 'chrN', GERP uses Ensembl 'N'."""
    bw = refs.bigwigs()
    if bw is None or bw[0] is None:
        return (np.nan, np.nan)
    ph, ge = bw
    pos = int(pos); lr, la = len(str(ref)), len(str(alt))
    if lr > la:
        s, e = pos, pos + lr - 1
    elif lr < la:
        s, e = pos - 1, pos + 1
    else:
        s, e = pos - 1, pos + lr - 1
    if e <= s:
        e = s + 1

    def m(bwh, name):
        try:
            vals = bwh.values(name, s, e)
            arr = [v for v in vals if v is not None and not math.isnan(v)]
            return float(sum(arr) / len(arr)) if arr else np.nan
        except Exception:
            return np.nan
    return (m(ph, chrom if str(chrom).startswith("chr") else "chr" + str(chrom)),
            m(ge, str(chrom).replace("chr", "")))


def _cds_span(model, s, e):
    """Genomic positions (transcript 5'->3') of the CDS nucleotides for residues [s..e].
    Returns (gpos, gmin, gmax, status); status 'ok' only if the span is a single contiguous
    genomic interval (no intron crossing)."""
    cds = model["cds_pos"]; n = len(cds)
    c1 = (s - 1) * 3 + 1; c2 = e * 3
    if c1 < 1 or c2 > n:
        return None, None, None, "cds_oob"
    gpos = cds[c1 - 1:c2]; gmin = min(gpos); gmax = max(gpos)
    if gmax - gmin + 1 != len(gpos):
        return None, None, None, "cross_intron"
    return gpos, gmin, gmax, "ok"


def _conservation_ins_v1(refs, transcript, gene, ps, pe):
    """INSERTION conservation: mean phyloP100 / GERP over the reference codons of the
    junction window [ps..pe]. Same value for either genomic representation. Returns None
    (caller falls back to _conservation) when the CDS model is unavailable."""
    bw = refs.bigwigs()
    if bw is None or bw[0] is None:
        return (np.nan, np.nan)
    model = refs.txmodel(transcript, gene)
    if model is None:
        return None
    ph, ge = bw
    chrom = model["chrom"]
    php = chrom if str(chrom).startswith("chr") else "chr" + str(chrom)
    ghp = str(chrom).replace("chr", "")

    def vals(bwh, name, a, b):
        try:
            return [v for v in bwh.values(name, a - 1, b) if v is not None and not math.isnan(v)]
        except Exception:
            return []
    P, G = [], []
    for res in range(int(ps), int(pe) + 1):
        _gp, gmin, gmax, st = _cds_span(model, res, res)
        if st != "ok":
            continue
        P += vals(ph, php, gmin, gmax)
        G += vals(ge, ghp, gmin, gmax)
    return (float(sum(P) / len(P)) if P else np.nan,
            float(sum(G) / len(G)) if G else np.nan)


# Public API
_POS_DEP_ANNOT = ("in_domain", "n_domains", "in_repeat", "ptm_nearby",
                  "ptm_count_near", "disulfide_nearby", "in_active_site")


def compute_features(variant: dict, refs: Refs) -> dict:
    """Return the 39 model features for one in-frame indel from cached inputs.
    Position-dependent features (structural, annotation, local_entropy,
    indel_size_relative) are NaN when the protein position is unmapped."""
    gene = variant["gene_symbol"]
    up = variant["uniprot"]
    itype = "deletion" if str(variant["indel_type"]).startswith(("del", "deletion")) or \
        variant["indel_type"] == 1 else "insertion"
    size = int(variant["indel_size_aa"])
    try:
        ps = int(variant["protein_pos_start"])
        pe = int(variant["protein_pos_end"])
        pos_ok = True
    except (TypeError, ValueError):
        ps = pe = None
        pos_ok = False

    obj = refs.uniprot_json(up)
    feats = {}

    # tier 2: gene constraint (position-independent)
    g = refs.gene_lookup.get(gene, (np.nan, np.nan, np.nan, np.nan))
    feats.update(pli=g[0], loeuf=g[1], mis_z=g[2], cds_length=g[3])

    # tier 4 subcellular (position-independent)
    feats.update(_subcell(obj))

    # conservation: del = genomic span; ins = junction-window codons (genomic-span fallback)
    if refs.with_conservation and variant.get("pos") is not None:
        cons = None
        if itype == "insertion" and pos_ok:
            cons = _conservation_ins_v1(refs, variant.get("transcript"), gene, ps, pe)
        if cons is None:
            cons = _conservation(refs, variant.get("chrom"), variant.get("pos"),
                                 variant.get("ref"), variant.get("alt"))
        ph, ge = cons
    else:
        ph, ge = (np.nan, np.nan)
    feats.update(phylop100_mean=ph, gerp_mean=ge)

    # dist_to_splice (genomic)
    feats.update(dist_to_splice=_dist_to_splice(refs, gene, variant.get("chrom"),
                                                variant.get("pos")))

    if pos_ok:
        prot = refs.protein(up)
        # tier 3: structural
        feats.update(_structural(prot, ps, pe, itype, size))
        # tier 4: UniProt annotation
        feats.update(_annotation(obj, ps, pe))
        # sequence context
        feats.update(local_entropy=_entropy(refs.sequence(up), ps, pe))
        # indel size relative to protein length
        plen = variant.get("protein_length")
        if plen is None or (isinstance(plen, float) and np.isnan(plen)):
            plen = len(prot["resnos"]) if prot is not None else np.nan
        feats.update(indel_size_relative=(size / plen)
                     if (plen and not (isinstance(plen, float) and np.isnan(plen))) else np.nan)
    else:
        # unmapped position -> position-dependent features are NaN (median-imputed later)
        feats.update(_structural(None, ps, pe, itype, size))   # all-NaN structural dict
        feats.update({k: np.nan for k in _POS_DEP_ANNOT})
        feats.update(local_entropy=np.nan, indel_size_relative=np.nan)

    # indel-intrinsic (always available)
    feats.update(indel_size_aa=size, indel_type=itype)
    return {k: feats.get(k) for k in FEATURE_COLS}


def compute_table(df: pd.DataFrame, refs: Refs) -> pd.DataFrame:
    rows = [compute_features(r._asdict(), refs) for r in df.itertuples(index=False)]
    return pd.DataFrame(rows, columns=FEATURE_COLS)


# CLI
def _read_any(p):
    return pd.read_parquet(p) if str(p).endswith(".parquet") else pd.read_csv(p, sep="\t")


def main(argv):
    if "--selftest" in argv:
        return _selftest(argv)
    args = [a for a in argv if not a.startswith("--")]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    infile, outfile = args[0], args[1]
    cache = next((argv[i + 1] for i, a in enumerate(argv) if a == "--cache"),
                 os.path.join(os.environ.get("INDELVAR_ROOT", "."), "data"))
    refs = Refs(cache, with_conservation="--no-conservation" not in argv)
    df = _read_any(infile)
    log(f"computing 39 features for {len(df):,} variants (cache={cache}) ...")
    feat_tbl = compute_table(df, refs)
    struct = [c for c in ("af_plddt_mean", "af_rsa_mean", "wcn_mean") if c in feat_tbl.columns]
    if len(feat_tbl) and struct and feat_tbl[struct].isna().mean().mean() > 0.9:
        log(f"WARNING: >90% of structural features are NaN: the per-protein caches under "
            f"{cache} look missing/empty (they ship only in the development repo; see README "
            "'Reproduce'). The output will be mostly NaN.")
    out = pd.concat([df.reset_index(drop=True), feat_tbl], axis=1)
    out.to_csv(outfile, sep="\t", index=False)
    log(f"wrote {outfile} ({len(out):,} rows, {out.shape[1]} cols)")


def _selftest(argv):
    """Recompute features for a sample of a RAW cohort and compare to stored values.
    Needs the raw input columns and per-protein caches; pass the raw cohort from
    build_cohort.py (not the computed train_set.tsv)."""
    root = Path(os.environ.get("INDELVAR_ROOT", "."))
    # optional positional raw-cohort path (skip flag values like --n / --cache)
    skip = {i + 1 for i, a in enumerate(argv) if a in ("--n", "--cache")}
    pos = [a for i, a in enumerate(argv)
           if not a.startswith("--") and i not in skip]
    data_path = Path(pos[0]) if pos else (root / "dataset" / "train_data" / "train_set.tsv")
    tf = pd.read_csv(data_path, sep="\t")
    required_raw = ["uniprot", "protein_pos_start", "protein_pos_end"]
    missing = [c for c in required_raw if c not in tf.columns]
    if missing:
        log(f"selftest: {data_path} lacks raw input column(s) {missing}; it is the "
            "computed feature matrix, not a raw cohort. Point --selftest at the raw "
            "cohort from build_cohort.py (with uniprot/protein_pos_start/"
            "protein_pos_end) and ensure the per-protein caches are present.")
        sys.exit(2)
    n = int(next((argv[i + 1] for i, a in enumerate(argv) if a == "--n"), 200))
    cache = root / "data"
    refs = Refs(cache, with_conservation=False)
    sample = tf.sample(min(n, len(tf)), random_state=42)
    log(f"selftest: {len(sample)} variants from {data_path} (conservation skipped)")
    check = [c for c in FEATURE_COLS if c not in
             ("phylop100_mean", "gerp_mean", "pli", "loeuf", "mis_z", "cds_length")]
    nbad = 0
    for r in sample.itertuples(index=False):
        got = compute_features(r._asdict(), refs)
        for c in check:
            a, b = got.get(c), getattr(r, c, None)
            if a is None or b is None:
                continue
            try:
                fa, fb = float(a), float(b)
                if (np.isnan(fa) and np.isnan(fb)):
                    continue
                if abs(fa - fb) > 1e-3:
                    nbad += 1
                    if nbad <= 12:
                        log(f"  MISMATCH {c}: got={fa:.4f} train={fb:.4f}")
            except (TypeError, ValueError):
                if str(a) != str(b):
                    nbad += 1
    log(f"selftest DONE: {nbad} mismatches over {len(sample)}×{len(check)} cells")
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main(sys.argv[1:])
