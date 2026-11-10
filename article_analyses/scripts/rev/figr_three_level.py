"""Test the flare-remission difference of each Table S8 edge at three units of analysis, using the
per-condition FigR correlations of the flare and remission runs:

  1. cell level      Fisher z-test, n = nuclei per condition (nuclei treated as independent)
  2. k-NN corrected  Fisher z-test, n = nuclei / k (k = 40, the smoothing neighbourhood)
  3. patient level   exact within-patient permutation from figr_differential_edges.R (joint
                     flare + remission run; motif-supported TF-DORC pairs only)

Also counts TF hubs (edges per TF) among all S8 edges and among motif-supported S8 edges, to check
whether the cell-type architecture described in the Results survives a motif filter.

usage: figr_three_level.py <S8_edge_provenance.tsv> <out_dir> <k>
           <celltype>=<n_flare>,<n_remission>,<published_S8_edge_status.tsv> ...
"""
import csv, sys, math, collections

prov_path, out_dir, k, *specs = sys.argv[1:]
k = int(k)
cts = {}
for s in specs:
    ct, v = s.split("=", 1)
    nf, nr, status = v.split(",", 2)
    cts[ct] = (int(nf), int(nr), status)

def p_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2))

def bh(ps):
    """Benjamini-Hochberg; None entries are skipped."""
    idx = [i for i, p in enumerate(ps) if p is not None]
    o = sorted(idx, key=lambda i: ps[i])
    q, m, n = [None] * len(ps), 1.0, len(idx)
    for r in range(n, 0, -1):
        i = o[r - 1]
        m = min(m, ps[i] * n / r)
        q[i] = m
    return q

def fisher_p(rf, rr, nf, nr):
    z = (math.atanh(rf) - math.atanh(rr)) / math.sqrt(1 / (nf - 3) + 1 / (nr - 3))
    return p_two_sided(z)

prov = collections.defaultdict(dict)
with open(prov_path) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        prov[r["cell_type"]][(r["motif"], r["gene"])] = r

def fmt(x, nd=4):
    return "" if x is None else (f"{x:.{nd}g}" if isinstance(x, float) else str(x))

edge_rows, summary, hubs = [], [], []
for ct, (nf, nr, status_path) in cts.items():
    with open(status_path) as fh:
        status = {(r["motif"], r["gene"]): r for r in csv.DictReader(fh, delimiter="\t")}
    edges = list(prov[ct].values())
    assert len(edges) == len(status), (ct, len(edges), len(status))

    rows = []
    for e in edges:
        key = (e["motif"], e["gene"])
        s = status[key]
        tested = e["remission_status"] == "tested"
        rf = float(e["flare_Corr"])
        rr = float(e["remission_Corr"]) if tested else None
        ep = s["Enrichment.P"]
        motif = None if ep in ("", "NA") else float(ep) < 0.05
        pp = None if s["p_perm"] in ("", "NA") else float(s["p_perm"])
        fdr_pat = None if s["FDR"] in ("", "NA") else float(s["FDR"])
        rows.append(dict(
            cell_type=ct, motif=e["motif"], gene=e["gene"], gene_set=s["gene_set"],
            rho_flare=rf, rho_remission=rr, remission_tested=tested, motif_support=motif,
            p_cell=fisher_p(rf, rr, nf, nr) if tested else None,
            p_knn=fisher_p(rf, rr, nf / k, nr / k) if tested else None,
            p_patient=pp, FDR_patient=fdr_pat, patient_status=s["status"]))
    for lvl in ("cell", "knn"):
        q = bh([r[f"p_{lvl}"] for r in rows])
        for r, x in zip(rows, q):
            r[f"FDR_{lvl}"] = x
    edge_rows += rows

    def count(sub, col):
        return sum(1 for r in sub if r[col] is not None and r[col] < 0.05)
    tested = [r for r in rows if r["remission_tested"]]
    tested_motif = [r for r in tested if r["motif_support"]]
    summary.append(dict(
        cell_type=ct, nuclei_flare=nf, nuclei_remission=nr, S8_edges=len(rows),
        tested_in_remission=len(tested), not_tested_in_remission=len(rows) - len(tested),
        sig_cell=count(tested, "FDR_cell"), sig_knn=count(tested, "FDR_knn"),
        motif_supported=sum(1 for r in rows if r["motif_support"]),
        motif_supported_tested=len(tested_motif),
        motif_sig_cell=count(tested_motif, "FDR_cell"), motif_sig_knn=count(tested_motif, "FDR_knn"),
        patient_tested=sum(1 for r in rows if r["FDR_patient"] is not None),
        motif_sig_patient=count(rows, "FDR_patient")))

    all_c = collections.Counter(r["motif"] for r in rows)
    mot_c = collections.Counter(r["motif"] for r in rows if r["motif_support"])
    rank_all = {tf: i + 1 for i, (tf, _) in enumerate(all_c.most_common())}
    rank_mot = {tf: i + 1 for i, (tf, _) in enumerate(mot_c.most_common())}
    for tf, n in all_c.most_common():
        hubs.append(dict(cell_type=ct, TF=tf, edges_S8=n, rank_S8=rank_all[tf],
                         edges_motif_supported=mot_c.get(tf, 0), rank_motif_supported=rank_mot.get(tf)))

cols = ["cell_type", "motif", "gene", "gene_set", "rho_flare", "rho_remission", "remission_tested",
        "motif_support", "p_cell", "FDR_cell", "p_knn", "FDR_knn", "p_patient", "FDR_patient",
        "patient_status"]
with open(f"{out_dir}/S8_edge_three_level_tests.tsv", "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(cols)
    for r in edge_rows:
        w.writerow([fmt(r[c]) for c in cols])

with open(f"{out_dir}/SUMMARY.tsv", "w") as out:
    w = csv.DictWriter(out, fieldnames=list(summary[0]), delimiter="\t")
    w.writeheader(); w.writerows(summary)

with open(f"{out_dir}/TF_hubs.tsv", "w") as out:
    w = csv.DictWriter(out, fieldnames=list(hubs[0]), delimiter="\t")
    w.writeheader()
    for h in hubs:
        w.writerow({c: fmt(v) for c, v in h.items()})

for s in summary:
    print("\t".join(f"{a}={b}" for a, b in s.items()))
