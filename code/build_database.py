#!/usr/bin/env python3
"""Assemble the genome-wide Zenodo tables from the precomputed per-type chunks:
release_features.parquet (the 39 features) and release_scores.parquet/.tsv.gz (score, ACMG
tier + points, mechanism flags). Both keyed on (gene, transcript, indel_type, aa_start,
indel_size_aa) so they join. Streams by row group, memory-safe over ~208M rows.

Usage: build_database.py [features|scores|all]   (default all)
"""
import os, gzip, sys
import numpy as np, pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

ROOT = os.environ.get("INDELVAR_ROOT", ".")
PRE = os.path.join(ROOT, "data/precomputed")
GENOMIC = ["chrom", "pos", "ref", "alt"]


def log(m): print(m, flush=True)


def _need(*paths):
    """Clean exit (like build_cohort.py) when a per-type chunk is absent: these ship
    only in the development tree / on Zenodo, not in the release."""
    missing = next((p for p in paths if not os.path.exists(p)), None)
    if missing:
        sys.exit(f"missing input: {missing}\nbuild_database.py assembles the genome-wide "
                 "tables from the precomputed per-type chunks under data/precomputed/, which "
                 "are not shipped in the release (the database is built in the development tree "
                 "and deposited on Zenodo); see README 'Reproduce'.")


# ---- release_features.parquet ----
FMETA = {"gene_symbol", "uniprot", "transcript", "protein_pos_start", "protein_pos_end",
         "ref_seq", "indel_type", "chrom", "pos", "ref", "alt", "map_status"}
FPAIRS = [("features_deletions.parquet", "deletion"), ("features_insertions.parquet", "insertion")]


def build_features():
    _need(*(os.path.join(PRE, src) for src, _ in FPAIRS))
    out = os.path.join(PRE, "release_features.parquet")
    writer = None
    for src, typ in FPAIRS:
        pf = pq.ParquetFile(os.path.join(PRE, src))
        feats = [c for c in pf.schema_arrow.names if c not in FMETA]
        out_cols = ["gene", "transcript", "indel_type", "aa_start", "indel_size_aa"] + \
                   [f for f in feats if f != "indel_size_aa"] + GENOMIC
        log(f"  {typ}: streaming {pf.metadata.num_rows:,} rows ({src}) ...")
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i)
            ps = t.column("protein_pos_start")
            aa_start = ps if typ == "deletion" else pc.add(ps, 1)   # ins: junction L = ps+1
            cols = {"gene": t.column("gene_symbol"), "transcript": t.column("transcript"),
                    "indel_type": pa.array([typ] * t.num_rows, type=pa.string()),
                    "aa_start": pc.cast(aa_start, pa.int64())}
            for f in feats:
                cols[f] = t.column(f)
            for gc in GENOMIC:
                cols[gc] = t.column(gc)
            bt = pa.table({c: cols[c] for c in out_cols})
            if writer is None:
                writer = pq.ParquetWriter(out, bt.schema, compression="zstd")
            writer.write_table(bt)
        log(f"  {typ}: written ({pf.metadata.num_rows:,} rows)")
    writer.close()


# ---- release_scores.parquet + .tsv.gz ----
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
       "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
       "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}
SCOL = ["gene", "transcript", "indel_type", "aa_start", "indel_size_aa", "hgvs_p",
        "indelvar_score", "acmg_tier", "acmg_points",
        "core_packing", "ss_break", "hbond_loss", "hydrophobic_exposure", "n_mechanisms",
        "chrom", "pos", "ref", "alt"]
SPAIRS = [("precomputed_deletions.parquet", "precomputed_deletions_genomic.parquet", "deletion"),
          ("precomputed_insertions.parquet", "precomputed_insertions_genomic.parquet", "insertion")]


def render_hgvsp(df):
    """deletion p.Xaa{ps}(_Yaa{pe})del; insertion (dup segment [ps-k+2..ps+1], junction
    L=ps+1) p.Xaa{L}dup or p.Xaa{ps-k+2}_Yaa{L}dup for k>1."""
    ps = df.protein_pos_start.astype("int64"); pe = df.protein_pos_end.astype("int64")
    k = df.indel_size_aa.astype("int64"); rs = df.ref_seq.fillna("")
    f3 = rs.str[0].map(AA3); l3 = rs.str[-1].map(AA3)
    is_del = df.indel_type.to_numpy() == "deletion"
    out = pd.Series(np.empty(len(df), dtype=object), index=df.index)
    d1 = is_del & (ps.to_numpy() == pe.to_numpy())
    dR = is_del & (ps.to_numpy() != pe.to_numpy())
    out[d1] = "p." + f3[d1] + ps[d1].astype(str) + "del"
    out[dR] = "p." + f3[dR] + ps[dR].astype(str) + "_" + l3[dR] + pe[dR].astype(str) + "del"
    ins = ~is_del
    L = ps + 1; dfirst = ps - k + 2
    i1 = ins & (k.to_numpy() == 1)
    iM = ins & (k.to_numpy() != 1)
    out[i1] = "p." + f3[i1] + L[i1].astype(str) + "dup"
    out[iM] = "p." + f3[iM] + dfirst[iM].astype(str) + "_" + l3[iM] + L[iM].astype(str) + "dup"
    bad = ~rs.str.fullmatch(r"[ARNDCQEGHILKMFPSTWYV]+").fillna(False)
    out[bad.to_numpy()] = ""
    return out.astype(str)


def score_block(prec, gen, typ):
    pp, gp = os.path.join(PRE, prec), os.path.join(PRE, gen)
    n = pq.ParquetFile(pp).metadata.num_rows
    assert pq.ParquetFile(gp).metadata.num_rows == n, f"{typ}: row-count mismatch"
    log(f"  {typ}: assembling {n:,} rows ...")
    p = pq.read_table(pp, columns=["gene_symbol", "transcript", "protein_pos_start",
          "protein_pos_end", "ref_seq", "indel_type", "indel_size_aa", "indelvar_score",
          "indelvar_category", "acmg_points", "mech_compact_core", "mech_ss_break",
          "mech_hbond", "mech_hydrophobic", "n_mechanisms", "uniprot"]).to_pandas()
    g = pq.read_table(gp, columns=["chrom", "pos", "ref", "alt", "uniprot",
          "protein_pos_start", "indel_size_aa"]).to_pandas()
    chk = np.random.RandomState(0).choice(n, size=min(3000, n), replace=False)
    ok = ((p.uniprot.values[chk] == g.uniprot.values[chk]) &
          (p.protein_pos_start.values[chk] == g.protein_pos_start.values[chk]) &
          (p.indel_size_aa.values[chk] == g.indel_size_aa.values[chk])).all()
    assert ok, f"{typ}: precomputed vs genomic NOT row-aligned"
    ps = p.protein_pos_start.astype("int64")
    aa_start = ps if typ == "deletion" else ps + 1     # del: first deleted; ins: junction L
    rel = pd.DataFrame({
        "gene": p.gene_symbol, "transcript": p.transcript, "indel_type": p.indel_type,
        "aa_start": aa_start.astype("int64"), "indel_size_aa": p.indel_size_aa.astype("int64"),
        "hgvs_p": render_hgvsp(p),
        "indelvar_score": p.indelvar_score.astype("float64"),
        "acmg_tier": p.indelvar_category, "acmg_points": p.acmg_points.astype("int64"),
        "core_packing": p.mech_compact_core.astype("int64"), "ss_break": p.mech_ss_break.astype("int64"),
        "hbond_loss": p.mech_hbond.astype("int64"), "hydrophobic_exposure": p.mech_hydrophobic.astype("int64"),
        "n_mechanisms": p.n_mechanisms.astype("int64"),
        "chrom": g.chrom.astype(str), "pos": g.pos.astype("Int64"), "ref": g.ref, "alt": g.alt})[SCOL]
    return rel


def build_scores():
    _need(*(os.path.join(PRE, p) for prec, gen, _ in SPAIRS for p in (prec, gen)))
    out_parquet = os.path.join(PRE, "release_scores.parquet")
    out_tsv = os.path.join(PRE, "release_scores.tsv.gz")
    writer = None; first = True
    for prec, gen, typ in SPAIRS:
        rel = score_block(prec, gen, typ)
        tbl = pa.Table.from_pandas(rel, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(out_parquet, tbl.schema, compression="zstd")
            gz = gzip.open(out_tsv, "wt", compresslevel=6)
        writer.write_table(tbl)
        rel.to_csv(gz, sep="\t", index=False, header=first)
        first = False
        log(f"  {typ}: written ({len(rel):,} rows)")
        del rel, tbl
    writer.close(); gz.close()


def report(name):
    pf = pq.ParquetFile(os.path.join(PRE, name))
    log(f"{name}: {pf.metadata.num_rows:,} rows, {pf.metadata.num_columns} cols")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("features", "all"):
        log("release_features: building unified deletion+insertion feature table ...")
        build_features(); report("release_features.parquet")
    if what in ("scores", "all"):
        log("release_scores: building unified deletion+insertion score table ...")
        build_scores(); report("release_scores.parquet")
    log("build_database done.")
