#!/usr/bin/env python3
"""Global concordance of pseudobulk log2 fold changes, primary (uncorrected) vs SoupX-corrected counts.

Per population, over genes tested in both analyses (primary expression filter in each):
  n_common, Pearson r and Spearman rho of logFC, slope of SoupX logFC on primary logFC (least squares through the origin),
  % same sign (all common genes), and for genes with primary genome-wide FDR<0.05: n, % same sign, % still FDR<0.05,
  median |logFC difference|. Pooled row = all population-gene pairs.

usage: soupx_concordance.py <primary de dir> <soupx de dir> <out dir>
"""
import csv
import math
import os
import sys

prim_dir, soup_dir, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)


def read(f):
    with open(f, encoding="utf-8") as h:
        return {r["gene"]: (float(r["logFC"]), float(r["FDR"])) for r in csv.DictReader(h, delimiter="\t")}, None


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy)


def median(v):
    v = sorted(v)
    n = len(v)
    return float("nan") if n == 0 else (v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2)


def stats(pairs):
    x = [p[0] for p in pairs]
    y = [p[2] for p in pairs]
    sig = [p for p in pairs if p[1] < 0.05]
    return {
        "n_common": len(pairs),
        "pearson_r": pearson(x, y),
        "spearman_rho": pearson(ranks(x), ranks(y)),
        "slope_through_origin": sum(a * b for a, b in zip(x, y)) / sum(a * a for a in x),
        "pct_same_sign_all": 100 * sum((a > 0) == (b > 0) for a, b in zip(x, y)) / len(x),
        "n_primary_fdr05": len(sig),
        "pct_same_sign_primary_fdr05": 100 * sum((p[0] > 0) == (p[2] > 0) for p in sig) / len(sig) if sig else float("nan"),
        "pct_still_fdr05": 100 * sum(p[3] < 0.05 for p in sig) / len(sig) if sig else float("nan"),
        "median_abs_dlogFC_primary_fdr05": median([abs(p[2] - p[0]) for p in sig]),
        "pearson_r_primary_fdr05": pearson([p[0] for p in sig], [p[2] for p in sig]) if len(sig) > 2 else float("nan"),
    }


rows, pooled = [], []
for f in sorted(os.listdir(prim_dir)):
    if not f.endswith("_edgeR_QL.tsv") or not os.path.exists(os.path.join(soup_dir, f)):
        continue
    a, _ = read(os.path.join(prim_dir, f))
    b, _ = read(os.path.join(soup_dir, f))
    pairs = [(a[g][0], a[g][1], b[g][0], b[g][1]) for g in a if g in b]
    if len(pairs) < 10:
        continue
    pooled += pairs
    rows.append({"population": f.replace("_edgeR_QL.tsv", ""), "n_primary": len(a), "n_soupx": len(b), **stats(pairs)})
rows.append({"population": "ALL (pooled)", "n_primary": sum(r["n_primary"] for r in rows), "n_soupx": sum(r["n_soupx"] for r in rows),
             **stats(pooled)})

keys = list(rows[0].keys())
with open(os.path.join(out, "logfc_concordance_primary_vs_soupx.tsv"), "w", encoding="utf-8") as h:
    h.write("\t".join(keys) + "\n")
    for r in rows:
        h.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
for r in rows:
    print(r["population"], r["n_common"], f"r={r['pearson_r']:.3f}", f"rho={r['spearman_rho']:.3f}", f"slope={r['slope_through_origin']:.3f}",
          f"sign={r['pct_same_sign_all']:.1f}%", f"sig={r['n_primary_fdr05']}", f"sigsign={r['pct_same_sign_primary_fdr05']:.1f}%",
          f"still={r['pct_still_fdr05']:.1f}%", f"medD={r['median_abs_dlogFC_primary_fdr05']:.3f}")
