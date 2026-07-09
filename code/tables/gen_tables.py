"""Table S6: per-bin likelihood ratios on the test set, by indel class and evidence
increment (+0.5 Haldane-Anscombe correction). HGMD-licensed test set not shipped;
place your own test_set_scored.tsv (variant_id, label, indel_type, indelvar) under
dataset/test_data/ to run.
"""
import csv, os, sys

ROOT   = os.environ.get("INDELVAR_ROOT", ".")
SCORED = os.path.join(ROOT, "dataset/test_data/test_set_scored.tsv")
CUTS   = os.path.join(ROOT, "model/evidence_cutoffs.tsv")
CUTS_T = os.path.join(ROOT, "model/cutoffs_local_1sided95.tsv")
OUT    = os.path.join(ROOT, "tables/table_S6_per_bin_LR.csv")
PC     = 0.5   # Haldane-Anscombe continuity correction

if not os.path.exists(SCORED):
    sys.exit("[S6] held-out test set not shipped (HGMD-licensed). Place your own "
             "test_set_scored.tsv under dataset/test_data/ to run this script.")
os.makedirs(os.path.join(ROOT, "tables"), exist_ok=True)

PP3 = ["PP3_Supporting", "PP3_Moderate", "PP3_Moderate3", "PP3_Strong"]
BP4 = ["BP4_Supporting", "BP4_Moderate", "BP4_Moderate3", "BP4_Strong"]

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

# ---- per-type cutoffs + evidence-LR targets ----
cuts, cuts_t = {}, {}
with open(CUTS) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cuts[(r["subset"], r["evidence"])] = r
with open(CUTS_T) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        cuts_t[(r["subset"], r["evidence"])] = r
def rc(sub, ev):  return num(cuts.get((sub, ev), {}).get("cutoff"))
def tg(sub, ev):  return num(cuts_t.get((sub, ev), {}).get("target"))

# ---- scored test set ----
P = {"deletion": [], "insertion": []}
B = {"deletion": [], "insertion": []}
with open(SCORED) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        t = r["indel_type"]
        if t not in P:
            continue
        s = num(r["indelvar"])
        if s is None:
            continue
        path = str(r["label"]) in ("1", "Pathogenic", "pathogenic", "P", "LP")
        (P if path else B)[t].append(s)

def lr_path(nP, nB, NP, NB):  # pathogenic-direction LR (PP3)
    return ((nP + PC) / (NP + PC)) / ((nB + PC) / (NB + PC))
def lr_ben(nP, nB, NP, NB):   # benign-direction LR magnitude (BP4)
    return ((nB + PC) / (NB + PC)) / ((nP + PC) / (NP + PC))

rows = [["Indel type", "Direction", "Increment", "Score interval",
         "n P", "n B", "Per-bin LR", "Target LR"]]
for t, sub in (("deletion", "del"), ("insertion", "ins")):
    NP, NB = len(P[t]), len(B[t])
    # PP3: reachable increments, ascending cutoffs; bin k = [c_k, c_{k+1}), top = [c_K, 1]
    pp = [(i + 1, rc(sub, ev)) for i, ev in enumerate(PP3) if rc(sub, ev) is not None]
    for j, (k, c) in enumerate(pp):
        hi = pp[j + 1][1] if j + 1 < len(pp) else None
        inP = [s for s in P[t] if s >= c and (hi is None or s < hi)]
        inB = [s for s in B[t] if s >= c and (hi is None or s < hi)]
        itv = "[%.3f, %s)" % (c, ("%.3f" % hi) if hi is not None else "1")
        rows.append([t, "PP3 (pathogenic)", "+%d" % k, itv, len(inP), len(inB),
                     round(lr_path(len(inP), len(inB), NP, NB), 2),
                     round(tg(sub, PP3[k - 1]), 2)])
    # BP4: reachable increments, descending cutoffs; bin k = (c_{k+1}, c_k], bottom = (0, c_K]
    bp = [(i + 1, rc(sub, ev)) for i, ev in enumerate(BP4) if rc(sub, ev) is not None]
    for j, (k, c) in enumerate(bp):
        lo = bp[j + 1][1] if j + 1 < len(bp) else 0.0
        inB = [s for s in B[t] if s <= c and s > lo]
        inP = [s for s in P[t] if s <= c and s > lo]
        itv = "(%s, %.3f]" % (("%.3f" % lo) if lo > 0 else "0", c)
        tp = tg(sub, PP3[k - 1])   # symmetric ACMG target magnitude (same strength as +k)
        tgt = round(tp, 2) if tp else ""
        rows.append([t, "BP4 (benign)", "-%d" % k, itv, len(inP), len(inB),
                     round(lr_ben(len(inP), len(inB), NP, NB), 2), tgt])

# order each indel class from the strongest pathogenic (+4) down to the strongest benign (-4)
_trank = {"deletion": 0, "insertion": 1}
rows = [rows[0]] + sorted(rows[1:], key=lambda r: (_trank.get(r[0], 9), -int(r[2])))

with open(OUT, "w", newline="") as f:
    csv.writer(f).writerows(rows)
OUTX = OUT[:-4] + ".xlsx"
try:
    import openpyxl
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "S6_per_bin_LR"
    for r in rows:
        ws.append(r)
    wb.save(OUTX)
    print("wrote %s and %s (%d bins, raw cutoffs, +%.1f continuity)" % (OUT, OUTX, len(rows) - 1, PC))
except ImportError:
    print("wrote %s (%d bins); openpyxl unavailable, .xlsx skipped" % (OUT, len(rows) - 1))
for r in rows[1:]:
    print("  " + " | ".join(str(x) for x in r))
