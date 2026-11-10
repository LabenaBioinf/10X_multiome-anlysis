# Pre-declaration: chromVAR extension to naive and regulatory populations

Written 6 November 2026, before the run. Declared post hoc extension of the 23 October 2026
chromVAR analysis (results/gdt_chromvar/MOTIFS.md).

Question: does accessibility at STAT-family motifs (STAT1, STAT4 cisBP motifs) rise at flare in
naive lymphocytes, as it does in CD8 TEM, gdT and CD4 TCM cells?

Added populations: CD4 Naive, CD8 Naive, B naive, Treg. CD4 TEM is not testable (2 pairs with
>= 20 nuclei in both samples) and is not added.

Unchanged from the original analysis: input (ATAC pseudobulk, >= 20 nuclei per sample), peak
filtering, cisBP motif set, families, chromVAR background (200 iterations, seed 20261016), limma
~ patient + condition on deviation z, BH within cell type (families and single motifs separately),
minimum four complete pairs. The original six populations are re-run in the same job as a
reproducibility check; their results must match results/gdt_chromvar.

Primary readout for this extension: STAT1 and STAT4 single motifs in the four added populations
(FDR < 0.05, BH within cell type over all single motifs, as before). Families are reported as
context. Results are reported whatever their direction or significance.

Script: scripts/rev/chromvar_naive_extension.R; job: scripts/rev/chromvar_naive_extension.sb.
