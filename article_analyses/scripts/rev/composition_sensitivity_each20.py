#!/usr/bin/env python3
"""Composition sensitivity with a per-sample threshold (corrects the pair-total rule of a6_composition_check.R).

A pair enters only if BOTH the flare and the remission sample contain >= 20 nuclei of the population; populations
with < 4 such pairs are not tested. Exact paired Wilcoxon signed-rank test on proportion differences (zeros dropped,
all 2^n sign assignments); Benjamini-Hochberg over the full family of 30 populations (untestable populations
enter as P = 1). The primary all-sample test is recomputed the same way as a check against
composition_tests_recomputed.tsv.

usage: composition_sensitivity_each20.py <results/a6_composition> <out.tsv>
"""
import csv
import itertools
import os
import sys

src, out = sys.argv[1:3]
counts = list(csv.DictReader(open(os.path.join(src, "nuclei_counts_per_subset_per_sample.tsv")), delimiter="\t"))
ref = {r["cell_type"]: r for r in csv.DictReader(open(os.path.join(src, "composition_tests_recomputed.tsv")),
                                                  delimiter="\t") if r["set"] == "all"}
samples = [c for c in counts[0] if c != "cell_type"]
pts = sorted({s.split(".", 1)[1] for s in samples})
tot = {s: sum(int(r[s]) for r in counts) for s in samples}

def exact_p(d):
    d = [x for x in d if x != 0]
    n = len(d)
    if n == 0:
        return 1.0
    order = sorted(range(n), key=lambda i: abs(d[i]))
    rank = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(d[order[j + 1]]) == abs(d[order[i]]):
            j += 1
        for k in range(i, j + 1):
            rank[order[k]] = (i + j) / 2 + 1
        i = j + 1
    obs = sum(rank[i] for i in range(n) if d[i] > 0)
    mid = sum(rank) / 2
    dist = [sum(rank[i] for i in range(n) if s[i]) for s in itertools.product([0, 1], repeat=n)]
    return sum(1 for w in dist if abs(w - mid) >= abs(obs - mid) - 1e-9) / len(dist)

def bh(ps):
    m = len(ps)
    o = sorted(range(m), key=lambda i: ps[i])
    q = [0.0] * m
    prev = 1.0
    for rnk in range(m, 0, -1):
        i = o[rnk - 1]
        prev = min(prev, ps[i] * m / rnk)
        q[i] = prev
    return q

rows = []
for r in counts:
    ct = r["cell_type"]
    d_all = [int(r[f"flare.{p}"]) / tot[f"flare.{p}"] - int(r[f"remission.{p}"]) / tot[f"remission.{p}"] for p in pts]
    keep = [p for p in pts if int(r[f"flare.{p}"]) >= 20 and int(r[f"remission.{p}"]) >= 20]
    d20 = [int(r[f"flare.{p}"]) / tot[f"flare.{p}"] - int(r[f"remission.{p}"]) / tot[f"remission.{p}"] for p in keep]
    rows.append({"cell_type": ct, "p_all": exact_p(d_all), "n_pairs_each20": len(keep),
                 "p_each20": exact_p(d20) if len(keep) >= 4 else 1.0, "testable_each20": len(keep) >= 4,
                 "p_all_reference": float(ref[ct]["p_exact_perm"]) if ct in ref else float("nan")})
for key, col in (("p_all", "fdr_all_30"), ("p_each20", "fdr_each20_30")):
    for r, q in zip(rows, bh([r[key] for r in rows])):
        r[col] = q
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    cols = ["cell_type", "p_all", "p_all_reference", "fdr_all_30", "n_pairs_each20", "testable_each20", "p_each20",
            "fdr_each20_30"]
    w.writerow(cols)
    for r in sorted(rows, key=lambda r: r["p_all"]):
        w.writerow([r[c] for c in cols])
mism = [r["cell_type"] for r in rows if abs(r["p_all"] - r["p_all_reference"]) > 1e-9]
print("families of", len(rows), "| primary P mismatches vs R:", mism or "none")
print("testable with >=20 nuclei in both samples:", sum(r["testable_each20"] for r in rows))
print("FDR<0.05 all:", [r["cell_type"] for r in rows if r["fdr_all_30"] < 0.05])
print("FDR<0.05 each>=20:", [r["cell_type"] for r in rows if r["fdr_each20_30"] < 0.05])
