"""Check that Table S8 edges are reproduced by the per-cell-type FigR outputs, and record
which edges were testable on the remission side.

usage: figr_provenance.py <table_S8.tsv> <out_dir> <celltype>=<flare_figR_d>,<remission_figR_d> ...
"""
import csv, sys, collections

s8_path, out_dir, *pairs = sys.argv[1:]
runs = {}
for p in pairs:
    ct, files = p.split("=", 1)
    runs[ct] = files.split(",")

s8 = collections.defaultdict(list)
with open(s8_path) as fh:
    next(fh)
    for r in csv.DictReader(fh, delimiter="\t"):
        s8[r["cell_type"]].append(r)

def load(path, keys):
    """Return {(Motif, DORC): row} for the requested keys plus the TF and DORC universes."""
    hits, tfs, dorcs = {}, set(), set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tfs.add(r["Motif"]); dorcs.add(r["DORC"])
            k = (r["Motif"], r["DORC"])
            if k in keys:
                hits[k] = r
    return hits, tfs, dorcs

summary = []
with open(f"{out_dir}/S8_edge_provenance.tsv", "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["cell_type", "motif", "gene", "S8_corr", "S8_pval", "flare_Corr", "flare_Corr.P",
                "flare_match", "remission_status", "remission_Corr", "remission_Corr.P"])
    for ct, (fl, rm) in runs.items():
        edges = s8[ct]
        keys = {(e["motif"], e["gene"]) for e in edges}
        fh_, ftf, fdorc = load(fl, keys)
        rh_, rtf, rdorc = load(rm, keys)
        c = collections.Counter()
        for e in edges:
            k = (e["motif"], e["gene"])
            f = fh_.get(k)
            match = bool(f) and abs(float(f["Corr"]) - float(e["corr"])) < 1e-3 \
                and abs(float(f["Corr.P"]) - float(e["pval"])) <= max(1e-6, 1e-3 * float(e["pval"]))
            if k in rh_:
                rs = "tested"
            elif e["motif"] not in rtf and e["gene"] not in rdorc:
                rs = "TF_and_DORC_absent"
            elif e["motif"] not in rtf:
                rs = "TF_absent"
            elif e["gene"] not in rdorc:
                rs = "DORC_absent"
            else:
                rs = "pair_absent"
            r = rh_.get(k, {})
            w.writerow([ct, e["motif"], e["gene"], e["corr"], e["pval"], f["Corr"] if f else "",
                        f["Corr.P"] if f else "", match, rs, r.get("Corr", ""), r.get("Corr.P", "")])
            c["match" if match else ("missing" if not f else "value_differs")] += 1
            c[rs] += 1
            if rs == "tested":
                c["rem_P<0.05" if float(r["Corr.P"]) < 0.05 else "rem_P>=0.05"] += 1
        summary.append((ct, len(edges), len(ftf), len(rtf), len(ftf - rtf), len(rtf - ftf),
                        len(fdorc), len(rdorc), dict(c)))

with open(f"{out_dir}/SUMMARY.tsv", "w") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["cell_type", "S8_edges", "flare_TFs", "remission_TFs", "TFs_flare_only",
                "TFs_remission_only", "flare_DORCs", "remission_DORCs", "counts"])
    for s in summary:
        w.writerow(list(s[:-1]) + [";".join(f"{k}={v}" for k, v in sorted(s[-1].items()))])
print(open(f"{out_dir}/SUMMARY.tsv").read())
