#!/usr/bin/env python3
"""Apply the reporting rule for flare vs remission differences that do not reach FDR < 0.05.

Rule (post hoc: defined 23 Oct 2026 after the module tests had been run; applied to every row): a difference is reported as "consistent" when
  1. the pooled (nuclei) estimate and the patient-level estimate have the same sign,
  2. the patients agree with that sign beyond chance: exact two-sided sign test p < 0.05
     (8/9 at n = 9, 9/10 at n = 10, 8/8 at n = 8),
  3. modules only: exact Wilcoxon signed-rank p < 0.05 over per-patient mean delta z.
FDR is reported alongside in every row; nothing here is called significant unless FDR < 0.05.

Inputs are the module-test output of figr_module_test.R / figr_module_test_minlib.R (modules and
per-edge patient counts), figr_differential (pooled within-patient permutation p per edge) and
gdt_chromvar (limma per motif on patient pseudobulk; condition 1 = limma direction, condition 3 =
limma P < 0.05).

The expected number of edges passing under no flare-remission difference is n_edges x P(sign test
p < 0.05 | p = 1/2) for each edge's n, i.e. ~3.9% at n = 9.

usage: consistency_rule.py <module_dir> <figr_differential_dir> <chromvar_dir> <out_dir> <label>
"""
import csv
import math
import os
import sys

mod_dir, diff_dir, cv_dir, out_dir, label = sys.argv[1:6]
os.makedirs(out_dir, exist_ok=True)


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def sign_p(k, n):
    """Exact two-sided sign test p for k of n in one direction."""
    if n == 0:
        return float("nan")
    k = max(k, n - k)
    tail = sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def null_rate(n):
    """P(sign test p < 0.05) when each patient is a fair coin."""
    return sum(math.comb(n, k) for k in range(n + 1) if sign_p(k, n) < 0.05) / 2 ** n


def num(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


def write(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# --- modules -------------------------------------------------------------------------------------
mods = read(os.path.join(mod_dir, "module_tests.tsv"))
mods = [r for r in mods if r["patients_flare_higher"] not in ("NA", "")]
for r in mods:
    n = int(r["patients_evaluable"]); up = int(r["patients_flare_higher"])
    pooled = num(r["mean_delta_z"]); pat = num(r["median_patient_delta_z"])
    agree = up if pooled > 0 else n - up
    r["patients_agree_pooled"] = agree
    r["sign_p"] = f"{sign_p(agree, n):.4g}"
    ok = (pooled * pat > 0) and sign_p(agree, n) < 0.05 and agree > n / 2 \
        and num(r["p_wilcoxon_patients"]) < 0.05
    r["direction"] = "flare higher" if pooled > 0 else "flare lower"
    r["consistent"] = "yes" if ok else "no"
mod_cols = ["cell_type", "module", "n_pairs", "direction", "mean_delta_z", "p_perm_module", "FDR_module",
            "patients_evaluable", "patients_agree_pooled", "median_patient_delta_z", "sign_p",
            "p_wilcoxon_patients", "consistent"]
write(os.path.join(out_dir, "modules.tsv"), mods, mod_cols)

# --- edges ---------------------------------------------------------------------------------------
pooled_p = {}
for ct in ("CD8_TEM", "gdT"):
    for r in read(os.path.join(diff_dir, f"{ct}_differential_edges.tsv")):
        pooled_p[(r["cell_type"], r["Motif"], r["DORC"])] = (r["p_perm"], r["FDR"])
edges = read(os.path.join(mod_dir, "module_edges.tsv"))
expected = {}
for r in edges:
    n = int(r["patients_evaluable"]); up = int(r["patients_flare_higher"]); d = num(r["delta_z"])
    agree = up if d > 0 else n - up
    p, fdr = pooled_p.get((r["cell_type"], r["Motif"], r["DORC"]), (r.get("p_perm", "NA"), r.get("FDR", "NA")))
    r.update(p_perm_pooled=p, FDR_pooled=fdr, patients_agree_pooled=agree,
             sign_p=f"{sign_p(agree, n):.4g}", direction="flare higher" if d > 0 else "flare lower")
    r["consistent"] = "yes" if n > 0 and agree > n / 2 and sign_p(agree, n) < 0.05 else "no"
    expected[r["cell_type"]] = expected.get(r["cell_type"], 0) + null_rate(n)
edge_cols = ["cell_type", "Motif", "DORC", "published", "direction", "rho_flare", "rho_remission", "delta_z",
             "p_perm_pooled", "FDR_pooled", "patients_evaluable", "patients_agree_pooled", "sign_p", "consistent"]
write(os.path.join(out_dir, "edges_all.tsv"), edges, edge_cols)
write(os.path.join(out_dir, "edges_consistent.tsv"), [r for r in edges if r["consistent"] == "yes"], edge_cols)

# --- chromVAR ------------------------------------------------------------------------------------
cv = []
for fn, kind in (("family_deviation_flare_vs_remission.tsv", "family"),
                 ("motif_deviation_flare_vs_remission.tsv", "motif")):
    for r in read(os.path.join(cv_dir, fn)):
        n = int(r["n_pairs"]); up = int(r["patients_flare_higher"]); d = num(r["delta_z"])
        agree = up if d > 0 else n - up
        r.update(level=kind, patients_agree=agree, sign_p=f"{sign_p(agree, n):.4g}",
                 direction="flare higher" if d > 0 else "flare lower")
        r["consistent"] = "yes" if agree > n / 2 and sign_p(agree, n) < 0.05 and num(r["P"]) < 0.05 else "no"
        cv.append(r)
cv_cols = ["cell_type", "level", "feature", "direction", "delta_z", "P", "FDR", "n_pairs", "patients_agree",
           "sign_p", "consistent"]
write(os.path.join(out_dir, "chromvar.tsv"), cv, cv_cols)

# --- summary -------------------------------------------------------------------------------------
lines = [f"# Consistency rule applied ({label})", "",
         "Rule: pooled and patient-level estimates have the same sign; patients agree beyond chance",
         "(exact sign test p < 0.05); modules also need Wilcoxon p < 0.05, chromVAR limma P < 0.05.",
         "Exploratory: FDR is shown and none of these rows is called significant unless FDR < 0.05.", "",
         "## Modules", "",
         "| Cell | Module | Direction | Pooled p | FDR | Patients agreeing | Wilcoxon p | Consistent |",
         "|---|---|---|---|---|---|---|---|"]
for r in mods:
    lines.append(f"| {r['cell_type']} | {r['module']} | {r['direction']} | {num(r['p_perm_module']):.3g} | "
                 f"{num(r['FDR_module']):.3g} | {r['patients_agree_pooled']}/{r['patients_evaluable']} | "
                 f"{num(r['p_wilcoxon_patients']):.3g} | {r['consistent']} |")
lines += ["", "## Edges", "",
          "| Cell | Edges tested | Consistent | Expected by chance | Of which in Table S8 | Pooled p < 0.05 too |",
          "|---|---|---|---|---|---|"]
for ct in sorted(expected):
    sub = [r for r in edges if r["cell_type"] == ct]
    c = [r for r in sub if r["consistent"] == "yes"]
    lines.append(f"| {ct} | {len(sub)} | {len(c)} | {expected[ct]:.0f} | "
                 f"{sum(r['published'] == 'TRUE' for r in c)} | {sum(num(r['p_perm_pooled']) < 0.05 for r in c)} |")
lines += ["", "## chromVAR (patient pseudobulk; pairs per row as shown; same table for every FigR run)", "",
          "| Cell | Level | Feature | Direction | Δz | P | FDR | Patients agreeing |", "|---|---|---|---|---|---|---|---|"]
for r in cv:
    if r["consistent"] == "yes":
        lines.append(f"| {r['cell_type']} | {r['level']} | {r['feature']} | {r['direction']} | "
                     f"{num(r['delta_z']):.2f} | {num(r['P']):.3g} | {num(r['FDR']):.3g} | "
                     f"{r['patients_agree']}/{r['n_pairs']} |")
with open(os.path.join(out_dir, "SUMMARY.md"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
