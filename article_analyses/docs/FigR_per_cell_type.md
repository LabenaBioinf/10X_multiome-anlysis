# Per-cell-type FigR runs (input of Table S8)

`FigR_networks_MISC.R` (flare) and `FigR_networks_Control.R` (remission) with their SLURM wrappers.
Both split the harmonised object by `celltype.condition`, run cisTopic (1–50 topics), k-NN smoothing,
`runGenePeakcorr` (hg38, 100 background peaks, pvalZ ≤ 0.05), `getDORCScores` and `runFigRGRN`, and
write `{cell type}_{condition}_figR_d.txt`. Differences: the flare script sets k = min(40, max(5,
0.1 × nuclei)) — 40 for every population used here — and the remission script sets NA DORC scores to 0
(no NAs occurred). Both keep only TFs whose genes are also DORCs before `runFigRGRN`, so TF sets
differ between conditions. CD8 TEM: runs of October 2024; γδ T: runs of February 2026.
`scripts/rev/figr_provenance.py` checks that every Table S8 edge is reproduced from these outputs.
