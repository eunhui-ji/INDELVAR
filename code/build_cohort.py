#!/usr/bin/env python
"""Build a MANE Select in-frame-indel cohort with Ensembl VEP.

Usage: build_cohort.py {train|test|vus|benign} [--validate-dup-index]
"""
import json, os, re, sys, time, urllib.request, urllib.error
import pandas as pd, numpy as np

ROOT = os.environ.get("INDELVAR_ROOT", os.getcwd())
COH  = os.path.join(ROOT, "data/cohorts")
CACHE_DIR = os.path.join(ROOT, "data/processed/mane_vep_cache")
os.makedirs(CACHE_DIR, exist_ok=True)
VEP_URL = "https://rest.ensembl.org/vep/human/region"

SRC = {  # cohort -> (input postable, output postable, has_label)
    "train": ("postable_train.tsv",   "postable_train_mane.tsv", True),
    "test":  ("postable_test.tsv",   "postable_test_mane.tsv",  True),
    "vus":   ("postable_vus.tsv",    "postable_vus_mane.tsv",   False),
    # gnomAD benign background
    "benign": ("postable_gnomad_benign.tsv", "postable_gnomad_benign_mane.tsv", False),
}

def log(m): print(m, flush=True)

def vcf(v):
    p = v.split("_"); return f"{p[0]} {p[1]} . {p[2]} {p[3]} . . ."

def vep_batch(batch):
    body = json.dumps({"variants": [vcf(v) for v in batch], "mane": 1, "hgvs": 1}).encode()
    req = urllib.request.Request(VEP_URL, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=300))

def vep_all(vids, cache_path):
    """VEP every vid, caching raw records by input-string. Resumable / offline once cached."""
    cache = {}
    if os.path.exists(cache_path):
        cache = {r["input"]: r for r in json.load(open(cache_path))}
        log(f"  cache: {len(cache)} records loaded")
    need = [v for v in vids if vcf(v) not in cache]
    log(f"  need VEP for {len(need)} / {len(vids)} variants")
    for i in range(0, len(need), 200):
        b = need[i:i+200]
        for attempt in range(5):
            try:
                for r in vep_batch(b):
                    cache[r["input"]] = r
                break
            except Exception as e:
                if attempt == 4:
                    log(f"  !! batch {i} failed after retries: {e}")
                else:
                    time.sleep(3 * (attempt + 1))
        if i % 2000 == 0:
            json.dump(list(cache.values()), open(cache_path, "w"))
            log(f"    {i+len(b)}/{len(need)} … cached")
    json.dump(list(cache.values()), open(cache_path, "w"))
    return cache

def mane_tc(rec):
    for tc in rec.get("transcript_consequences", []):
        if tc.get("mane_select"):
            return tc
    return None

# --- HGVSp parsers ---------------------------------------------------------- #
_DUP_RANGE  = re.compile(r"[A-Za-z]{3}(\d+)_[A-Za-z]{3}(\d+)dup")
_DUP_SINGLE = re.compile(r"[A-Za-z]{3}(\d+)dup")
_INS_RANGE  = re.compile(r"[A-Za-z]{3}(\d+)_[A-Za-z]{3}(\d+)ins")
_DEL_RANGE  = re.compile(r"^[A-Za-z]{3}(\d+)_[A-Za-z]{3}(\d+)")
_DEL_SINGLE = re.compile(r"^[A-Za-z]{3}(\d+)")

def _clean_hgvsp(hgvsp):
    if not hgvsp:
        return None
    p = hgvsp.split(":")[-1]
    p = p[2:] if p.startswith("p.") else p
    return p.strip("()")

def parse_del_window(hgvsp):
    """Deleted-residue window [start,end] from the MANE HGVSp."""
    p = _clean_hgvsp(hgvsp)
    if p is None:
        return (None, None)
    m = _DEL_RANGE.match(p)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    m = _DEL_SINGLE.match(p)
    if m:
        return (int(m.group(1)), int(m.group(1)))
    return (None, None)

def junction_L(hgvsp):
    """Junction L (residue 5' of the insertion point, 3'-canonical). Returns (L, kind)."""
    p = _clean_hgvsp(hgvsp)
    if p is None:
        return (None, None)
    m = _DUP_RANGE.search(p)
    if m:
        return (int(m.group(2)), "dup")     # 3'-most duplicated residue
    m = _INS_RANGE.search(p)
    if m:
        return (int(m.group(1)), "ins")     # left flank of the insertion point
    m = _DUP_SINGLE.search(p)
    if m:
        return (int(m.group(1)), "dup")
    return (None, None)


def main(cohort, validate_dup_index=False):
    inp, outp, has_label = SRC[cohort]
    inp_path = os.path.join(COH, inp)
    if not os.path.exists(inp_path):
        sys.exit(f"missing input: {inp_path}\nbuild_cohort.py reads the ClinVar/gnomAD "
                 "postables from the development repo (data/cohorts/), which are not shipped "
                 "in the release; see README 'Reproduce'.")
    base = pd.read_csv(inp_path, sep="\t")
    # ref/alt from variant_id (test/vus postables omit them)
    sp = base.variant_id.str.split("_", expand=True)
    base["chrom"] = sp[0]; base["pos"] = sp[1].astype(int)
    base["ref"] = sp[2];   base["alt"] = sp[3]
    vids = base.variant_id.tolist()
    log(f"[{cohort}] base variants: {len(vids)}")

    cache = vep_all(vids, os.path.join(CACHE_DIR, f"{cohort}.json"))

    # gene/uniprot/protein_length from the MANE transcript + _gene_lookup
    glp = os.path.join(ROOT, "data/precomputed/_gene_lookup.csv")
    gene_lookup = {}
    if os.path.exists(glp):
        gl = pd.read_csv(glp).dropna(subset=["uniprot"]).drop_duplicates("gene_symbol")
        gene_lookup = {row.gene_symbol: (row.uniprot, row.protein_length)
                       for row in gl.itertuples(index=False)}
    else:
        log(f"  WARN: {glp} missing: keeping base-postable gene/uniprot/protein_length")

    rows, kept, dropped_csq, no_mane, no_pos = [], 0, {}, 0, 0
    ins_kind = {"dup": 0, "ins": 0}
    L_clipped = 0                       # window clips at a terminus (L<=1 or L+2>plen)
    gene_reassigned = uniprot_reassigned = 0
    old_ps = base.set_index("variant_id").get("protein_pos_start")
    match_ps = 0; comparable = 0
    for r in base.itertuples(index=False):
        rec = cache.get(vcf(r.variant_id))
        tc = mane_tc(rec) if rec else None
        if tc is None:
            no_mane += 1; continue
        ct = set(tc.get("consequence_terms", []))
        if ct != {"inframe_deletion"} and ct != {"inframe_insertion"}:
            key = ";".join(sorted(ct)); dropped_csq[key] = dropped_csq.get(key, 0) + 1
            continue
        itype = "deletion" if "inframe_deletion" in ct else "insertion"
        # gene/uniprot/protein_length from the MANE transcript
        gene_symbol = tc.get("gene_symbol") or getattr(r, "gene_symbol")
        gl_hit = gene_lookup.get(gene_symbol)
        if gl_hit is not None:
            uniprot, plen = gl_hit[0], gl_hit[1]
        else:
            uniprot, plen = getattr(r, "uniprot"), getattr(r, "protein_length", np.nan)
        if gene_symbol != getattr(r, "gene_symbol"):
            gene_reassigned += 1
        if str(uniprot) != str(getattr(r, "uniprot")):
            uniprot_reassigned += 1
        if itype == "deletion":
            ps, pe = parse_del_window(tc.get("hgvsp"))
        else:
            L, kind = junction_L(tc.get("hgvsp"))
            if L is None:
                no_pos += 1; continue
            ins_kind[kind] += 1
            ps, pe = L - 1, L + 2         # flanking window [L-1 .. L+2]
            if L <= 1 or (pd.notna(plen) and pe > plen):
                L_clipped += 1
        if ps is None:
            no_pos += 1; continue
        row = dict(variant_id=r.variant_id,
                   gene_symbol=gene_symbol,
                   uniprot=uniprot,
                   indel_type=itype,
                   protein_pos_start=ps, protein_pos_end=pe,
                   protein_length=plen,
                   chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alt)
        if has_label:
            row["label"] = getattr(r, "label")
        rows.append(row); kept += 1
        # position sanity vs old (deletions only; insertions intentionally change)
        if itype == "deletion" and old_ps is not None:
            o = old_ps.get(r.variant_id)
            if pd.notna(o):
                comparable += 1; match_ps += int(int(o) == ps)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(COH, outp), sep="\t", index=False)
    log(f"\n[{cohort}] KEPT {kept} / {len(vids)}  -> {outp}")
    log(f"  no MANE tx:        {no_mane}")
    log(f"  no parseable pos/L:{no_pos}")
    log(f"  dropped (non-singleton csq): {sum(dropped_csq.values())}")
    for k, n in sorted(dropped_csq.items(), key=lambda x: -x[1])[:12]:
        log(f"      {n:5d}  {k}")
    log(f"  del start-pos == old: {match_ps}/{comparable} "
        f"({100*match_ps/comparable:.1f}%)" if comparable else "  (no old pos to compare)")
    log(f"  type split: {out.indel_type.value_counts().to_dict()}")
    log(f"  MANE-consistency reassigned: gene {gene_reassigned}, uniprot/plen {uniprot_reassigned}")
    log(f"  insertion junction kind: {ins_kind}  (window clipped at terminus: {L_clipped})")
    if has_label:
        log(f"  label split: {out.label.value_counts().to_dict()}")

    # L-range validity: junction L = ps+1 must be within [1, protein_length]
    ins = out[out.indel_type == "insertion"].copy()
    if len(ins):
        ins["L"] = ins.protein_pos_start.astype(int) + 1
        bad = ins[(ins.L < 1) | (pd.notna(ins.protein_length) & (ins.L > ins.protein_length))]
        log(f"  insertion L in [1,plen]: {len(ins)-len(bad)}/{len(ins)} valid"
            + (f"  (⚠ {len(bad)} out of range)" if len(bad) else ""))

    if validate_dup_index:
        _validate_dup_index(out)


def _validate_dup_index(out):
    """Cross-check HGVSp-derived L against the genomic-left-aligned dup index."""
    idx_fp = os.path.join(ROOT, "data/precomputed/precomputed_insertions_genomic.parquet")
    if not os.path.exists(idx_fp):
        log("  [validate] dup index not found: skipping"); return
    import pyarrow.parquet as pq
    idx = pq.read_table(idx_fp, columns=["chrom", "pos", "ref", "alt",
        "protein_pos_end", "indel_size_aa", "map_status"]).to_pandas()
    idx = idx[idx.map_status == "ok"].copy(); idx["chrom"] = idx.chrom.astype(str)
    idx["size"] = idx.indel_size_aa.astype(int)
    maxpe = idx.groupby(["chrom", "pos", "ref", "alt", "size"])["protein_pos_end"].max()
    ins = out[out.indel_type == "insertion"].copy()
    ins["L"] = ins.protein_pos_start.astype(int) + 1
    ins["size"] = (ins.alt.str.len() - ins.ref.str.len()).abs() // 3
    ins["chrom"] = ins.chrom.apply(lambda c: c[3:] if str(c).startswith("chr") else str(c))
    ins = ins.merge(maxpe.rename("maxpe").reset_index(),
                    on=["chrom", "pos", "ref", "alt", "size"], how="left")
    m = ins[ins.maxpe.notna()]
    if len(m):
        off = m.L.astype(int) - m.maxpe.astype(int)
        log(f"  [validate] dup-index matched {len(m)}/{len(ins)}; "
            f"L==maxpe {100*(off==0).mean():.1f}%; offset>=0 {100*(off>=0).mean():.1f}% "
            f"(3'-shift vs left-align, expected)")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args[0] if args else "train",
         validate_dup_index="--validate-dup-index" in sys.argv)
