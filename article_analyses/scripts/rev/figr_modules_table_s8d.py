#!/usr/bin/env python3
"""Supplementary Table S8d (transcription factor module tests) as reported in the article.

Input: Table_S8d_FigR_modules.tsv written by figr_supp_tables_v4.py (module tests of figr_module_test.R and
figr_module_test_minlib.R). The article no longer reports the flare-only association analysis, so the three
modules defined on that analysis (module names containing "Table S8") and the consistency-rule column are
removed, and the Benjamini-Hochberg FDR is recomputed over the remaining modules within each FDR scope.

usage: figr_modules_table_s8d.py <Table_S8d_FigR_modules.tsv> <out.tsv>
"""
import csv
import sys

src, out = sys.argv[1:3]
with open(src, encoding="utf-8") as fh:
    m = list(csv.reader(fh, delimiter="\t"))
hdr = m[0]
ix = {c: hdr.index(c) for c in ("Module (motif-supported TF–DORC pairs)", "P, exact within-patient permutation",
                                "FDR", "FDR scope")}
keep = [x for x in m[1:] if "Table S8" not in x[ix["Module (motif-supported TF–DORC pairs)"]]]
assert len(m) - 1 - len(keep) == 3, "expected exactly three modules defined on the flare-only analysis"


def bh(ps):
    n = len(ps)
    o = sorted(range(n), key=lambda i: ps[i])
    q, run = [0.0] * n, 1.0
    for k in range(n - 1, -1, -1):
        run = min(run, ps[o[k]] * n / (k + 1))
        q[o[k]] = run
    return q


for sc in sorted({x[ix["FDR scope"]] for x in keep}):
    rows = [x for x in keep if x[ix["FDR scope"]] == sc]
    for x, q in zip(rows, bh([float(x[ix["P, exact within-patient permutation"]]) for x in rows])):
        x[ix["FDR"]] = repr(q)
        x[ix["FDR scope"]] = f"{len(rows)} modules, " + sc.split(", ", 1)[1]
drop = [i for i, c in enumerate(hdr) if "consistency rule" in c]
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    for x in [hdr] + keep:
        w.writerow([v for i, v in enumerate(x) if i not in drop])
print(f"{len(keep)} modules written to {out}; lowest FDR {min(float(x[ix['FDR']]) for x in keep):.3f}")
