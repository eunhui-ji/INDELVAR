"""INDELVAR lookup: retrieve the precomputed score, evidence tier, and fold-destabilization
mechanisms (optionally the 39 features) for an in-frame indel, queried by gene (or MANE
transcript) and its HGVS protein-level description (e.g. p.Leu127del or p.Leu127_Ala128insGly).

Usage: indelvar_lookup.py --gene CFTR --hgvsp p.Leu127dup [--transcript ...] [--features]
"""
import argparse, json, os, re, sys
import pandas as pd

AA1 = set("ARNDCQEGHILKMFPSTWYV")
MECH = ["core_packing_loss", "ss_break", "hbond_loss", "hydrophobic_exposure"]

def _count_residues(seq: str) -> int:
    seq = seq.strip()
    if seq.isdigit():
        return int(seq)
    three = re.findall(r"[A-Z][a-z]{2}", seq)
    if three and "".join(three) == seq:
        return len([t for t in three if t != "Ter"])
    if all(c in AA1 for c in seq):
        return len(seq)
    return len(re.findall(r"[A-Z][a-z]{2}", seq)) or len(seq)

def parse_hgvsp(h: str):
    s = re.sub(r"^[^:]*:", "", h.strip())
    s = re.sub(r"^p\.", "", s).strip("()")
    if "delins" in s or "fs" in s or "ext" in s:
        return "oos", {"why": "delins / frameshift / extension not modelled"}
    m = re.fullmatch(r"[A-Za-z]{1,3}(\d+)(?:_[A-Za-z]{1,3}(\d+))?del", s)
    if m:
        a = int(m.group(1)); b = int(m.group(2)) if m.group(2) else a
        return "deletion", {"aa_start": a, "size": b - a + 1}
    m = re.fullmatch(r"[A-Za-z]{1,3}(\d+)(?:_[A-Za-z]{1,3}(\d+))?dup", s)
    if m:
        a = int(m.group(1)); b = int(m.group(2)) if m.group(2) else a
        return "insertion", {"aa_start": b, "size": b - a + 1}
    m = re.fullmatch(r"[A-Za-z]{1,3}(\d+)_[A-Za-z]{1,3}(\d+)ins(.+)", s)
    if m:
        return "insertion", {"aa_start": int(m.group(1)), "size": _count_residues(m.group(3))}
    return "oos", {"why": f"unrecognised / unsupported HGVS p.: {h}"}

def _read(tables, stem, col, val):
    base = val.split(".")[0] if col == "transcript" else None
    pfp = os.path.join(tables, stem + ".parquet")
    if os.path.exists(pfp):
        import pyarrow.dataset as ds, pyarrow.compute as pc
        flt = ds.field(col) == val if base is None else \
              (ds.field(col) == base) | pc.starts_with(ds.field(col), base + ".")
        return ds.dataset(pfp).to_table(filter=flt).to_pandas()
    tfp = os.path.join(tables, stem + ".tsv.gz")
    if os.path.exists(tfp):
        sys.stderr.write(f"[warn] {stem}.parquet not found; scanning {stem}.tsv.gz (slow)\n")
        df = pd.read_csv(tfp, sep="\t")
        return df[df[col].str.split(".").str[0] == base] if base is not None else df[df[col] == val]
    sys.exit(f"ERROR: {stem}.(parquet|tsv.gz) not found in {tables}")

def lookup(gene, transcript, hgvsp, tables, want_features=False):
    kind, key = parse_hgvsp(hgvsp)
    if kind == "oos":
        return {"status": "out_of_scope", "reason": key["why"],
                "note": "INDELVAR covers in-frame deletions and insertions (1-10 aa) on MANE Select only."}
    if not (1 <= key["size"] <= 10):
        return {"status": "out_of_scope", "reason": f"size {key['size']} aa (INDELVAR covers 1-10 aa)"}
    sub = _read(tables, "release_scores", "gene", gene) if gene else \
          _read(tables, "release_scores", "transcript", transcript)
    if sub.empty:
        return {"status": "not_found", "reason": f"gene/transcript not in MANE-based DB: {gene or transcript}"}
    mane = sub.transcript.iloc[0]
    if transcript and transcript.split(".")[0] != mane.split(".")[0]:
        return {"status": "transcript_mismatch", "mane_transcript": mane,
                "reason": f"HGVS p. must be on the MANE Select transcript {mane}; you gave {transcript}. "
                          "Re-run VEP with --mane_select."}
    hit = sub[(sub.indel_type == kind) & (sub.aa_start == key["aa_start"]) &
              (sub.indel_size_aa == key["size"])]
    if hit.empty:
        return {"status": "not_found", "mane_transcript": mane,
                "reason": "no record at that protein position/size (check the HGVS p. is MANE-canonical)"}
    r = hit.iloc[0]
    res = {"status": "ok", "gene": r.gene, "transcript": r.transcript, "type": kind,
           "hgvs_p": r.hgvs_p, "aa_start": int(r.aa_start), "indel_size_aa": int(r.indel_size_aa),
           "indelvar_score": round(float(r.indelvar_score), 3),
           "acmg_tier": r.acmg_tier, "acmg_points": int(r.acmg_points),
           "mechanisms": {m: int(r[m]) for m in MECH}, "n_mechanisms": int(r.n_mechanisms)}
    if kind == "insertion":
        res["scope_note"] = ("score/features depend on position + size only; the inserted "
                             "sequence is not modelled (same score for any in-frame insertion "
                             "of this size at this junction)")
    if want_features:
        res["features"] = _lookup_features(kind, r, tables)
    return res

def _lookup_features(kind, row, tables):
    import pyarrow.dataset as ds
    fp = os.path.join(tables, "release_features.parquet")
    if not os.path.exists(fp):
        return {"error": "release_features.parquet not found (feature table C not deposited)"}
    flt = (ds.field("gene") == row.gene) & (ds.field("indel_type") == kind) & \
          (ds.field("aa_start") == int(row.aa_start)) & (ds.field("indel_size_aa") == int(row.indel_size_aa))
    m = ds.dataset(fp).to_table(filter=flt).to_pandas()
    if m.empty:
        return {"error": "feature row not found (protein-key join)"}
    drop = {"gene", "transcript", "indel_type", "aa_start", "chrom", "pos", "ref", "alt"}
    r0 = m.iloc[0]
    feats = {c: (float(r0[c]) if pd.notna(r0[c]) else None) for c in m.columns if c not in drop}
    feats["indel_type"] = kind
    return feats

def main():
    ap = argparse.ArgumentParser(description="INDELVAR protein-coordinate lookup")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene"); g.add_argument("--transcript")
    ap.add_argument("--hgvsp", required=True, help="MANE-based HGVS p., e.g. p.Leu127dup")
    ap.add_argument("--features", action="store_true", help="also return the 39 model features")
    ap.add_argument("--tables", default=os.environ.get("INDELVAR_TABLES",
                    os.path.join(os.environ.get("INDELVAR_ROOT", "."), "data")),
                    help="directory with the downloaded release_scores / release_features tables")
    a = ap.parse_args()
    print(json.dumps(lookup(a.gene, a.transcript, a.hgvsp, a.tables, a.features), indent=2))

if __name__ == "__main__":
    main()
