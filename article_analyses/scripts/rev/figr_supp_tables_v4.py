#!/usr/bin/env python3
"""Supplementary Tables S8d (FigR module tests), S8e (chromVAR families and motifs) and S8f (NK FigR
pairs) with readable headers.

usage: figr_supp_tables_v4.py <results dir> <out dir>
"""
import csv
import os
import sys

res, out = sys.argv[1:3]
CT = {"CD8_TEM": "CD8 TEM", "gdT": "γδ T", "CD8 TEM": "CD8 TEM", "NK": "NK",
      "CD4 TCM": "CD4 TCM", "CD8 TCM": "CD8 TCM", "MAIT": "MAIT"}
MODULE = {"TBX21_all": "TBX21, all pairs", "TBX21_published": "TBX21, pairs in Table S8",
          "published_all": "All pairs in Table S8", "all_pairs": "All pairs"}

def read(*parts):
    with open(os.path.join(res, *parts)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write(name, header, rows):
    with open(os.path.join(out, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)

def frac(a, n):
    return f"{a}/{n}" if n not in ("", "0", "NA") else "not evaluable"

# S8d: module tests. Rule columns from results/consistency_rule; descriptive rho from the module runs.
desc = {(r["cell_type"], r["module"]): r
        for r in read("figr_module", "module_tests.tsv") + read("figr_nk", "module_tests.tsv")}
rows = []
for sub, scope in (("min30", "10 modules, CD8 TEM and γδ T"), ("nk", "7 modules, NK")):
    for r in read("consistency_rule", sub, "modules.tsv"):
        d = desc[(r["cell_type"], r["module"])]
        rows.append([CT[r["cell_type"]], MODULE.get(r["module"], f'{r["module"]}, all pairs'),
                     r["n_pairs"], d["median_rho_flare"], d["median_rho_remission"],
                     f'{d["pairs_flare_higher"]}/{r["n_pairs"]}', r["direction"], r["mean_delta_z"],
                     r["p_perm_module"], r["FDR_module"], scope,
                     frac(r["patients_agree_pooled"], r["patients_evaluable"]),
                     r["median_patient_delta_z"], r["sign_p"], r["p_wilcoxon_patients"],
                     r["consistent"]])
write("Table_S8d_FigR_modules.tsv",
      ["Cell type", "Module (motif-supported TF–DORC pairs)", "Pairs", "Median rho flare",
       "Median rho remission", "Pairs with rho flare > remission", "Direction (pooled)",
       "Mean delta z (pooled)", "P, exact within-patient permutation", "FDR", "FDR scope",
       "Patients agreeing with pooled direction", "Median per-patient delta z", "Sign test P",
       "Wilcoxon P, patients", "Meets consistency rule (exploratory, post hoc)"],
      rows)

# S8e: chromVAR (patient pseudobulk); identical in both consistency_rule runs
mot = {(r["cell_type"], r["feature"]): r for r in read("gdt_chromvar", "family_deviation_flare_vs_remission.tsv")
       + read("gdt_chromvar", "motif_deviation_flare_vs_remission.tsv")}
order = ["CD8 TEM", "gdT", "NK", "CD4 TCM", "CD8 TCM", "MAIT"]
cv = sorted(read("consistency_rule", "min30", "chromvar.tsv"),
            key=lambda r: (order.index(r["cell_type"]), r["level"] != "family", float(r["P"])))
write("Table_S8e_chromVAR.tsv",
      ["Cell type", "Level", "Motif family or motif", "Patient pairs", "Direction", "Delta z (limma)",
       "P", "FDR (within cell type and level)", "Patients with flare > remission",
       "Patients agreeing with direction", "Sign test P", "Meets consistency rule (exploratory, post hoc)"],
      [[CT[r["cell_type"]], r["level"], r["feature"], r["n_pairs"], r["direction"], r["delta_z"],
        r["P"], r["FDR"], mot[(r["cell_type"], r["feature"])]["patients_flare_higher"],
        f'{r["patients_agree"]}/{r["n_pairs"]}', r["sign_p"], r["consistent"]] for r in cv])

# S8f: NK FigR, every motif-supported pair with a defined correlation in both conditions
enr = {(r["Motif"], r["DORC"]): r for r in read("figr_nk", "module_edges.tsv")}
nk = sorted(read("consistency_rule", "nk", "edges_all.tsv"), key=lambda r: float(r["p_perm_pooled"]))
write("Table_S8f_NK_FigR.tsv",
      ["TF", "DORC (target gene)", "Motif enrichment P", "rho flare", "rho remission", "Direction",
       "Delta z (pooled)", "P, exact within-patient permutation", "FDR (2,864 pairs)",
       "Patients with flare > remission", "Patients agreeing with pooled direction", "Sign test P",
       "Meets consistency rule (exploratory, post hoc)"],
      [[r["Motif"], r["DORC"], enr[(r["Motif"], r["DORC"])]["Enrichment.P"], r["rho_flare"],
        r["rho_remission"], r["direction"], r["delta_z"], r["p_perm_pooled"], r["FDR_pooled"],
        frac(enr[(r["Motif"], r["DORC"])]["patients_flare_higher"], r["patients_evaluable"]),
        frac(r["patients_agree_pooled"], r["patients_evaluable"]), r["sign_p"], r["consistent"]]
       for r in nk])
