#!/usr/bin/env python3
"""Every SoupX number cited in the article, recomputed from the SoupX results folder (results/ambient_soupx).

usage: soupx_text_numbers.py <soupx results dir> <primary results/pseudobulk dir>
"""
import csv
import os
import re
import statistics as st
import sys

S, P = sys.argv[1:3]


def tsv(f):
    with open(f, encoding="utf-8") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def de(d, pop):
    f = os.path.join(d, "de", pop.replace(" ", "_") + "_edgeR_QL.tsv")
    return {r["gene"]: (float(r["logFC"]), float(r["FDR"])) for r in tsv(f)} if os.path.exists(f) else {}


# contamination
R = tsv(os.path.join(S, "soupx_contamination_per_library_internal.tsv"))
for cond in sorted({r["condition"] for r in R}):
    rr = [r for r in R if r["condition"] == cond]
    rho = [float(r["rho"]) for r in rr]
    ig = [float(r["soup_fraction_immunoglobulin"]) for r in rr]
    rig = [float(r["rho_immunoglobulin"]) for r in rr if r["rho_immunoglobulin"] not in ("NA", "")]
    print(f"{cond}: rho median {st.median(rho):.3f} range {min(rho):.3f}-{max(rho):.3f} | Ig share of soup median {st.median(ig):.4f}"
          f" | rho_IG n={len(rig)}" + (f" median {st.median(rig):.3f} range {min(rig):.3f}-{max(rig):.3f}" if rig else ""))
allrho = sorted(float(r["rho"]) for r in R)
print("highest rho", allrho[-1], "next", allrho[-2], "ratio", round(allrho[-1] / allrho[-2], 2))
pat = {}
for r in R:
    pat.setdefault(r["library"].split("_")[0], {})[r["condition"]] = r
flare_hi = [p for p, v in pat.items() if len(v) == 2 and float(v["MISC"]["soup_fraction_immunoglobulin"]) >
            float([v[c] for c in v if c != "MISC"][0]["soup_fraction_immunoglobulin"])]
print("Ig soup share higher at flare in", len(flare_hi), "of", sum(len(v) == 2 for v in pat.values()), "patients")
for r in sorted(R, key=lambda r: -float(r["rho"]))[:3]:
    print("  top rho", r["condition"], r["rho"], "rho_IG", r["rho_immunoglobulin"], "method", r["method"])

# immunoglobulin genes in T/NK after correction
ig_pat = re.compile(r"^IG[HKL][VDJC]|^IGH[AGMDE][0-9]?$|^IGKC$|^IGLC[0-9]$|^JCHAIN$")
TNK = ["CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM", "CD8 TEM", "Treg", "MAIT", "gdT", "dnT", "NK", "NK CD56bright"]
print("Ig genes passing the filter in T/NK:", sum(1 for p in TNK for g in de(os.path.join(S, "pseudobulk"), p) if ig_pat.search(g)))

# Hallmark camera
cam = tsv(os.path.join(S, "pseudobulk", "camera_hallmark_all.tsv"))
cam0 = tsv(os.path.join(P, "camera_hallmark_all.tsv"))
def sigpops(c, s, d="Up"):
    return sorted(r["cell_type"] for r in c if r["set"] == s and r["Direction"] == d and float(r["FDR"]) < 0.05)
ifng = sigpops(cam, "HALLMARK_INTERFERON_GAMMA_RESPONSE")
print("IFN-g set Up FDR<0.05 in", len(ifng), "of", len({r["cell_type"] for r in cam}), "populations")
for s in ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_E2F_TARGETS"]:
    a, b = sigpops(cam0, s), sigpops(cam, s)
    print(f"{s}: same populations {a == b} (primary {len(a)}, soupx {len(b)})" + ("" if a == b else f" only primary {set(a) - set(b)} only soupx {set(b) - set(a)}"))
down = [(r["cell_type"], r["set"], round(float(r["FDR"]), 3)) for r in cam if r["Direction"] == "Down" and float(r["FDR"]) < 0.05]
print("Hallmark Down FDR<0.05:", down)

# key genes
D = lambda p: de(os.path.join(S, "pseudobulk"), p)
D0 = lambda p: de(P, p)
groups = [(["CD8 TEM", "NK", "gdT"], ["STAT1", "IL6ST", "GBP1", "GZMB", "PRF1"]), (["CD8 TEM", "NK"], ["PRDM1"]),
          (["CD8 TEM"], ["BATF", "HAVCR2", "IFNG"])]
mx, dl = 0, 0
for pops, genes in groups:
    for p in pops:
        for g in genes:
            v = D(p).get(g)
            mx = max(mx, v[1] if v else 9); dl = max(dl, abs(v[0] - D0(p)[g][0]) if v else 9)
print("effector genes max FDR", round(mx, 4))
mx2 = 0
for p in ["CD4 Naive", "CD8 Naive"]:
    for g in ["BATF", "IRF4", "STAT1", "GBP1", "SOCS1", "SOCS3", "IL6ST", "IFNG-AS1"]:
        v = D(p).get(g)
        mx2 = max(mx2, v[1] if v else 9); dl = max(dl, abs(v[0] - D0(p)[g][0]) if v else 9)
print("naive genes max FDR", round(mx2, 4), "| max |dlogFC| over all listed genes", round(dl, 3))
print("STAT4 NK", [round(x, 4) for x in D("NK")["STAT4"]])
for p, g in [("CD8 TEM", "IFNG-AS1"), ("gdT", "NKG7"), ("CD14 Mono", "STAT3")]:
    print("lost?", p, g, [round(x, 4) for x in D(p).get(g, (float("nan"), float("nan")))])
sig = tsv(os.path.join(S, "pseudobulk_followup", "signature_tests_camera.tsv"))
for r in sig:
    if r["cell_type"] == "gdT" and "IL27" in r["set"]:
        print("gdT", r["set"], r["Direction"], round(float(r["FDR"]), 3))
for p, g in [("CD8 TEM", "SOCS1"), ("CD8 TEM", "SOCS3"), ("gdT", "SOCS1"), ("CD8 Naive", "PRDM1"), ("gdT", "BATF"), ("CD8 TEM", "CTLA4"),
             ("CD8 TEM", "LAG3")]:
    print("in primary filter after SoupX?", p, g, g in D(p))
kp = tsv(os.path.join(S, "pseudobulk_followup", "key_gene_panel_pseudobulk.tsv"))
for p, g in [("gdT", "BATF"), ("CD8 TEM", "CTLA4"), ("CD8 TEM", "LAG3")]:
    r = [x for x in kp if x["cell_type"] == p and x["gene"] == g][0]
    print("panel estimate", p, g, round(float(r["logFC"]), 2))
tnf = [r for r in cam if r["cell_type"] == "CD14 Mono" and r["set"] == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"][0]
print("CD14 Mono TNF set", tnf["Direction"], round(float(tnf["FDR"]), 3))
