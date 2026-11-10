################################################################################
##
##  FigR_filter_flare_only_networks.R
##
##  PURPOSE:
##    Retain only significant TF→DORC regulatory edges from FigR (runFigRGRN)
##    that are exclusively present during MIS-C flare and absent in remission samples of the same children ("Control" in
##    file names). Edges are stratified into three biologically motivated
##    gene sets: IFN-γ pathway, TGF-β pathway, and Checkpoint/Exhaustion.
##    Analysis is performed independently for CD8+ TEM and γδ T cells.
##
##  INPUT:
##    FigR runFigRGRN output tables (one per cell type × condition):
##      - {cell_type}_MISC_figR_d.txt
##      - {cell_type}_Control_figR_d.txt
##
##  OUTPUT (per cell type):
##    - {ct}_all_merged_{gene_set}.xlsx        All edges from both conditions
##    - {ct}_only_MISC_{gene_set}.xlsx         MIS-C–exclusive edges
##    - {ct}_only_Control_{gene_set}.xlsx      Control-exclusive edges
##    - {ct}_shared_{gene_set}.xlsx            Edges present in both conditions
##    - {ct}_network_MISC_{gene_set}.html/rds  FigR network visualization (MIS-C only)
##    - {ct}_network_Control_{gene_set}.html/rds  FigR network (Control only)
##    - Supplementary_Table_FigR_{ct}.xlsx     Combined table for publication
##    - {ct}_summary_statistics.txt            Edge counts and QC metrics
##
##  CHANGELOG vs original script:
##    - FIXED: TGF-β repressor subsetting used Control instead of MISC (line 182)
##    - FIXED: Checkpoint repressor subsetting used Control instead of MISC (line 280)
##    - FIXED: TGF-β network visualization used IFN-γ variable (lines 217-222)
##    - REVISED: Gene lists aligned with v3 30-panel transcriptomic report
##    - ADDED:  Checkpoint gene set (replaces old Exhaustion markers)
##    - ADDED:  γδ T cell support with cell-type-specific gene lists
##    - ADDED:  Summary statistics and QC reporting
##    - ADDED:  Reproducible function-based workflow (DRY principle)
##
################################################################################


# ══════════════════════════════════════════════════════════════════════════════
# 0. SETUP
# ══════════════════════════════════════════════════════════════════════════════

# Libraries
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(dplyr)
library(Cairo)
library(FigR)
library(grid)
library(tidyr)
library(openxlsx)

# ── USER CONFIGURATION ──────────────────────────────────────────────────────
# Modify these paths to match your directory structure

base_dir       <- "/path/to/FigR_results"
figR_input_dir <- file.path(base_dir, "cell_types_FigR_CD8_TEM_gdT")
output_dir     <- file.path(base_dir, "CORRECTED_FigR_results")

# Cell types to process
cell_types <- c("CD8 TEM", "gdT")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Set seed for reproducibility
set.seed(42)


# ══════════════════════════════════════════════════════════════════════════════
# 1. GENE SET DEFINITIONS
# ══════════════════════════════════════════════════════════════════════════════
#
# Gene sets define the DORC target genes whose regulatory TF networks we want
# to extract. Each gene was selected based on:
#   (a) membership in the corresponding biological pathway
#   (b) relevance to CD8+ TEM or γδ T cell biology
#   (c) differential expression or DORC status in MIS-C vs controls
#   (d) centrality in the v3 30-panel transcriptomic report findings
#
# For cross-cell-type comparison, we define SHARED CORE gene sets used in both
# cell types, with cell-type-specific additions where biology diverges.
#
# ══════════════════════════════════════════════════════════════════════════════


# ── 1a. IFN-γ PATHWAY ───────────────────────────────────────────────────────
# Rationale: The IFN-γ pathway is the dominant transcriptional feature of
# MIS-C, with IFNG upregulated across all 3 producer populations (CD8+ TEM,
# γδ T, NK) and IFN-γ response signatures activated in 26/26 cell types.
# This gene set captures the production axis (IFNG, IFNG-AS1, effector TFs),
# the signaling cascade (STAT1, JAK2, SOCS1), and key downstream ISGs.
#
# IFNG-AS1 is included because our analysis revealed a paradoxical pattern:
# IFNG is upregulated while IFNG-AS1 (its cis-antisense positive regulator)
# is downregulated in 15/26 cell types — resolved at the chromatin level
# by TBX21 dual regulation (opens IFNG, closes IFNG-AS1 chromatin).

IFNG_genes_CORE <- c(
  # IFN-γ production axis (v3: Panel 1 - IFNG Production)
  "IFNG",          # Effector cytokine; CD8 TEM FC=+1.12***, γδT FC=+0.98
  "IFNG-AS1",      # Antisense lncRNA — paradox target; CD8 FC=-1.49***, γδT FC=-1.61***
  "BATF",          # AP-1 effector TF; CD8 FC=+1.33***, γδT FC=+2.40*** (strongest)
  "PRDM1",         # Blimp-1, terminal differentiation; CD8 FC=+1.78***, γδT FC=+1.05***
  
  # IFN-γ signaling cascade (v3: Panel 2 - IFNG Signaling)
  "STAT1",         # Core IFN-γ signal transducer; CD8 FC=+0.79***, γδT FC=+1.14***
  "STAT4",         # IL-12/IL-27 signal — DOWN = IL-12 bypass; CD8 FC=-0.44***, γδT FC=-0.69***
  "JAK2",          # IFN-γ receptor kinase; CD8 FC=+0.63***, γδT FC=+0.45
  "IRF1",          # IFN-γ primary response TF; CD8 FC=-0.15
  "IRF9",          # ISGF3 component; CD8 FC=+0.77***, γδT FC=+0.70**
  "SOCS1",         # Direct IFN-γ/STAT1 negative feedback; CD8 FC=+1.52***, γδT FC=+1.74***
  "SOCS3",         # Cytokine signaling feedback; CD8 FC=+1.29***, γδT FC=+1.93***
  
  # IFN-γ response genes / ISGs (v3: Panels 3,4 - IFNG Response, ISGs)
  "GZMB",          # Granzyme B, cytotoxic effector; CD8 FC=+1.54***, γδT FC=+1.46***
  "GBP1",          # Top ISG; CD8 FC=+1.66***, γδT FC=+1.76***
  "GBP4",          # ISG; CD8 FC=+1.17***, γδT FC=+1.22***
  "PSMB8",         # Immunoproteasome subunit; CD8 FC=+1.46***, γδT FC=+1.40***
  "ISG15",         # Pan-ISG; CD8 FC=+2.33***, γδT FC=+1.32***
  
  # Effector TFs (v3: Panel 5 - Effector TFs)
  "ZEB2",          # Terminal effector TF; CD8 FC=+0.14
  "ID2",           # Effector vs memory balance; CD8 FC=+0.65***, γδT FC=+0.51***
  "IRF4",          # TCR-driven effector fate; CD8 FC=+0.93*, γδT FC=+0.89
  "IRF8",          # IFN-γ amplifier; CD8 FC=+0.90, γδT FC=+1.39
  
  # IFN-γ–induced checkpoint ligand (links IFN-γ to exhaustion)
  "CD274"          # PD-L1; CD8 FC=+0.77, γδT FC=+0.76
)

# CD8+ TEM-specific additions
# Rationale: TBX21, BHLHE40, EOMES, RUNX3 are adaptive CTL effector TFs
# that are upregulated or functionally relevant in CD8+ TEM but not in γδ T.
# TBX21 in particular is UP in CD8 TEM (+0.36) but DOWN in γδT (−0.39),
# making it cell-type-discriminating.
IFNG_genes_CD8_TEM <- c(
  IFNG_genes_CORE,
  "TBX21",         # T-bet; CD8 FC=+0.36 vs γδT FC=-0.39 (OPPOSITE)
  "BHLHE40",       # DEC1, IFN-γ co-regulator; CD8 FC=+0.46* vs γδT FC=-0.14
  "EOMES",         # CD8 effector/memory TF; CD8 FC=+0.43, γδT: not detected
  "RUNX3"          # CTL lineage TF; CD8: DORC, γδT FC=-0.32
)

# γδ T cell-specific additions
# Rationale: γδ T cells use a partially overlapping but distinct effector
# program. They lack TBX21/EOMES dependence and instead rely on innate-like
# cytotoxicity (NKG7, PRF1, GZMA) and stronger IL-27 receptor signaling
# (IL6ST/gp130, STAT3). CX3CR1 marks their terminal effector differentiation.
IFNG_genes_gdT <- c(
  IFNG_genes_CORE,
  "NKG7",          # Cytotoxic granule; γδT FC=+1.67*** (strongest)
  "PRF1",          # Perforin; γδT FC=+1.28***
  "GZMA",          # Granzyme A; γδT FC=+1.30***
  "KLRD1",         # CD94/NKG2 NK receptor; γδT FC=+0.87***
  "IL6ST",         # gp130 / IL-27 receptor; γδT FC=+1.10*** (signature)
  "STAT3",         # IL-27 → STAT3 arm; γδT FC=+0.58***
  "CX3CR1"         # Terminal effector homing receptor; γδT FC=+1.19**
)


# ── 1b. TGF-β PATHWAY ───────────────────────────────────────────────────────
# Rationale: TGF-β signaling is profoundly suppressed in MIS-C across
# lymphoid populations. TGFB1, SMAD7, SKIL, SKI, and PMEPA1 are all
# downregulated, indicating collapse of both ligand production and
# intracellular SMAD signaling. This gene set captures the entire
# canonical pathway from ligand/receptors through SMADs to target TFs.
#
# Notably, γδ T cells show the DEEPEST TGF-β suppression (SMAD7 −1.44***,
# SKIL −1.26***), making cross-cell-type comparison of TF regulators
# particularly informative.
#
# IDENTICAL gene list for both cell types — the TGF-β pathway is the same
# signaling cascade regardless of cell type. Differences emerge in which
# TFs regulate these DORCs, not which DORCs to examine.

TGFb_genes <- c(
  # Ligand and receptors (v3: Panel 6 - TGFb Ligands/Receptors)
  "TGFB1",         # Core ligand; CD8 FC=-0.56***, γδT FC=-0.76***
  "TGFBR1",        # Type I receptor (ALK5); CD8 FC=+0.17, γδT FC=-0.14
  "TGFBR2",        # Type II receptor; CD8 FC=+0.22, γδT FC=+0.52***
  "TGFBR3",        # Type III receptor (betaglycan); CD8 FC=+0.30***
  "BAMBI",         # Pseudoreceptor / TGF-β decoy inhibitor
  
  # SMAD signaling (v3: Panel 7 - TGFb SMAD Cascade)
  "SMAD2",         # Receptor-SMAD; DORC in CD8 TEM (poised chromatin)
  "SMAD3",         # Receptor-SMAD; DORC in CD8 TEM (poised chromatin)
  "SMAD4",         # Common-mediator SMAD; essential signal integrator
  "SMAD7",         # Inhibitory SMAD; CD8 FC=-0.58, γδT FC=-1.44*** (deepest)
  "SKI",           # SMAD corepressor; CD8 FC=-0.52***
  "SKIL",          # SMAD corepressor (SnoN); CD8 FC=-0.79***, γδT FC=-1.26***
  "PMEPA1",        # TGF-β–induced feedback; CD8 FC=-1.37, γδT FC=-1.45
  
  # TGF-β target TFs (v3: Panel 8 - TGFb Target TFs)
  "ZEB1",          # EMT-TF / TGF-β effector; CD8 FC=-0.43***
  "ID3",           # TGF-β–induced, inhibitor of differentiation
  "FOXP3",         # Treg TF — not transcribed in CD8 TEM but DORC present
                   # (accessible chromatin = poised regulatory state)
  
  # TGF-β–associated regulatory genes
  "LTBP2",         # ECM-bound TGF-β regulator; DORC in CD8 TEM
  "RAB31"          # Vesicle trafficking / TGF-β presentation; CD8 FC=-1.75
)


# ── 1c. CHECKPOINT / EXHAUSTION ──────────────────────────────────────────────
# Rationale: A central finding of our analysis is that MIS-C CD8+ TEM cells
# show checkpoint receptor upregulation (CTLA4, HAVCR2, LAG3) WITHOUT
# canonical exhaustion TF activation (TOX flat, NR4A1/2/3 downregulated).
# This "activated but not exhausted" phenotype is critical for understanding
# why the IFN-γ response is not self-limited.
#
# γδ T cells present an even more striking pattern: checkpoint receptors are
# barely detectable (PDCD1 ND, HAVCR2 +0.15 ns) despite parallel effector
# activation, suggesting fundamentally different exhaustion trajectories.
#
# This gene set includes inhibitory receptors, exhaustion-driving TFs, and
# counter-regulators (TCF7, BCL6, MYC) that oppose the exhaustion program.

Checkpoint_genes_CORE <- c(
  # Inhibitory checkpoint receptors
  "PDCD1",         # PD-1; CD8 FC=+1.28, γδT: not detected
  "CTLA4",         # CTLA-4; CD8 FC=+2.40***, γδT FC=+1.95
  "HAVCR2",        # TIM-3; CD8 FC=+1.34***, γδT FC=+0.15 (absent!)
  "TIGIT",         # TIGIT; CD8 FC=+0.12, γδT FC=+0.36
  "LAG3",          # LAG-3; CD8 FC=+1.82***, γδT FC=+1.52***
  "ENTPD1",        # CD39, ectonucleotidase; CD8 FC=+0.15, γδT FC=-0.11
  
  # Exhaustion-driving TFs
  "TOX",           # Master exhaustion TF; CD8 FC=-0.11, γδT FC=-0.12 (NOT engaged)
  "NR4A1",         # Nur77, NFAT-driven exhaustion; CD8 FC=-0.32
  "NR4A2",         # Nurr1; CD8 FC=-0.70***, γδT FC=-0.96*** (DOWN in both)
  "NR4A3",         # NOR-1; CD8 FC=-0.67**, γδT FC=-0.67**
  
  # NFAT axis (exhaustion vs activation switch)
  "NFATC1",        # NFAT1; CD8 FC=-0.22, γδT FC=-0.51
  "NFATC2",        # NFAT2, partners with TOX; CD8 FC=+0.15, γδT FC=-0.23
  
  # Anti-exhaustion / counter-regulators
  "TCF7",          # TCF-1, stem-like maintenance; CD8 FC=-0.73***, γδT FC=-0.29
  "MYC",           # Metabolic counter-regulator; CD8 FC=+1.30*, γδT FC=+1.17**
  "BCL6"           # Represses exhaustion program; CD8 FC=+0.37, γδT FC=+0.46
)

# CD8+ TEM-specific additions
# Rationale: CD8 TEM shows classical exhaustion features (CD38 high,
# TIM-3/LAG-3 high) and terminal exhaustion markers (LAYN, CD101) as DORCs.
Checkpoint_genes_CD8_TEM <- c(
  Checkpoint_genes_CORE,
  "CD38",          # Activation/exhaustion marker; CD8 FC=+2.17***
  "TNFRSF9",       # 4-1BB, co-stimulatory activation; CD8 FC=+1.06
  "LGALS9",        # Galectin-9, TIM-3 ligand; CD8 FC=+1.29**
  "LAYN",          # Layilin, terminal exhaustion marker; DORC in CD8 TEM
  "CD101"          # CD101/IGSF2, terminal exhaustion; DORC in CD8 TEM
)

# γδ T cell-specific additions
# Rationale: γδ T cells escape classical exhaustion. ICOS (co-stimulatory)
# is strongly upregulated, and KLRC2 (NKG2C inhibitory receptor) represents
# an NK-like inhibitory axis distinct from PD-1/TIM-3.
Checkpoint_genes_gdT <- c(
  Checkpoint_genes_CORE,
  "CD38",          # Activation marker; γδT FC=+1.14***
  "KLRC2",         # NKG2C, inhibitory NK receptor; γδT FC=-0.42
  "ICOS"           # Co-stimulatory receptor; γδT FC=+1.16***
)


# ══════════════════════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

#' Extract condition-exclusive TF→DORC regulatory edges for a gene set
#'
#' @param MISC_full     Full FigR output (figR.d) for MIS-C condition
#' @param Control_full  Full FigR output (figR.d) for Control condition
#' @param gene_list     Character vector of DORC target genes to filter
#' @param gene_set_name Label for the gene set (e.g., "IFNG", "TGFb")
#' @param p_threshold   P-value threshold for FigR correlation (default: 0.05)
#' @return List with: MISC_only, Control_only, shared, merged, summary stats
extract_flare_only <- function(MISC_full, Control_full, gene_list, gene_set_name,
                                p_threshold = 0.05) {
  
  # ── Step 1: Filter to gene set DORCs with significant TF–DORC correlation
  MISC_figR   <- subset(MISC_full,    DORC %in% gene_list & Corr.P < p_threshold)
  Control_figR <- subset(Control_full, DORC %in% gene_list & Corr.P < p_threshold)
  
  # ── Step 2: Split by regulatory direction
  # CRITICAL FIX: Both activators AND repressors must come from the CORRECT condition
  MISC_act    <- subset(MISC_figR,    Corr >= 0)
  MISC_rep    <- subset(MISC_figR,    Corr < 0)     # ← FIXED: was Control in original
  Control_act <- subset(Control_figR, Corr >= 0)
  Control_rep <- subset(Control_figR, Corr < 0)
  
  # ── Step 3: Identify condition-exclusive edges via anti_join on TF–DORC pairs
  # MIS-C–only: present in MIS-C, absent in Control
  act_only_MISC <- anti_join(MISC_act, Control_act, by = c("DORC", "Motif"))
  rep_only_MISC <- anti_join(MISC_rep, Control_rep, by = c("DORC", "Motif"))
  
  # Control-only: present in Control, absent in MIS-C
  act_only_Control <- anti_join(Control_act, MISC_act, by = c("DORC", "Motif"))
  rep_only_Control <- anti_join(Control_rep, MISC_rep, by = c("DORC", "Motif"))
  
  # Shared: present in both conditions
  act_shared <- inner_join(MISC_act, Control_act, by = c("DORC", "Motif"),
                           suffix = c("_MISC", "_Control"))
  rep_shared <- inner_join(MISC_rep, Control_rep, by = c("DORC", "Motif"),
                           suffix = c("_MISC", "_Control"))
  
  # ── Step 4: Combine activators + repressors
  MISC_only    <- bind_rows(act_only_MISC, rep_only_MISC)
  Control_only <- bind_rows(act_only_Control, rep_only_Control)
  shared       <- bind_rows(act_shared, rep_shared)
  
  # Add gene set label
  MISC_only$gene_set    <- gene_set_name
  Control_only$gene_set <- gene_set_name
  
  # ── Step 5: Create merged table (all edges from both conditions)
  merged <- bind_rows(
    MISC_figR %>% mutate(Condition = "MISC"),
    Control_figR %>% mutate(Condition = "Control")
  ) %>% mutate(gene_set = gene_set_name)
  
  # ── Step 6: Summary statistics
  genes_queried   <- length(gene_list)
  genes_as_DORCs  <- length(unique(c(MISC_figR$DORC, Control_figR$DORC)))
  genes_not_DORCs <- setdiff(gene_list, unique(c(MISC_figR$DORC, Control_figR$DORC)))
  
  summary <- data.frame(
    gene_set         = gene_set_name,
    genes_queried    = genes_queried,
    genes_detected   = genes_as_DORCs,
    genes_not_DORC   = paste(genes_not_DORCs, collapse = ", "),
    MISC_total       = nrow(MISC_figR),
    MISC_activators  = nrow(MISC_act),
    MISC_repressors  = nrow(MISC_rep),
    Control_total    = nrow(Control_figR),
    Control_act      = nrow(Control_act),
    Control_rep      = nrow(Control_rep),
    MISC_only_edges  = nrow(MISC_only),
    MISC_only_act    = nrow(act_only_MISC),
    MISC_only_rep    = nrow(rep_only_MISC),
    Control_only_edges = nrow(Control_only),
    shared_edges     = nrow(shared),
    stringsAsFactors = FALSE
  )
  
  # ── Print summary to console
  cat(sprintf("\n=== %s ===\n", gene_set_name))
  cat(sprintf("  Genes queried: %d | Detected as DORCs: %d\n", genes_queried, genes_as_DORCs))
  if (length(genes_not_DORCs) > 0) {
    cat(sprintf("  Not detected as DORCs: %s\n", paste(genes_not_DORCs, collapse = ", ")))
  }
  cat(sprintf("  MIS-C edges:    %d total (%d activators, %d repressors)\n",
              nrow(MISC_figR), nrow(MISC_act), nrow(MISC_rep)))
  cat(sprintf("  Control edges:  %d total (%d activators, %d repressors)\n",
              nrow(Control_figR), nrow(Control_act), nrow(Control_rep)))
  cat(sprintf("  MIS-C–only:     %d (%d act, %d rep)\n",
              nrow(MISC_only), nrow(act_only_MISC), nrow(rep_only_MISC)))
  cat(sprintf("  Control-only:   %d (%d act, %d rep)\n",
              nrow(Control_only), nrow(act_only_Control), nrow(rep_only_Control)))
  cat(sprintf("  Shared:         %d edges\n", nrow(shared)))
  
  return(list(
    MISC_only    = MISC_only,
    Control_only = Control_only,
    shared       = shared,
    merged       = merged,
    summary      = summary
  ))
}


#' Plot and save FigR network for a set of regulatory edges
#'
#' @param edges     Data frame of TF→DORC edges (from extract_flare_only)
#' @param file_html Path for HTML output
#' @param file_rds  Path for RDS output
#' @param score_cut Score cutoff for plotfigRNetwork (default: 0)
plot_and_save_network <- function(edges, file_html, file_rds, score_cut = 0) {
  if (nrow(edges) == 0) {
    cat("  [SKIP] No edges to plot\n")
    return(NULL)
  }
  
  net <- FigR::plotfigRNetwork(
    edges,
    score.cut  = score_cut,
    TFs        = unique(edges$Motif),
    weight.edges = TRUE
  )
  r2d3::save_d3_html(net, file_html)
  saveRDS(net, file = file_rds)
  cat(sprintf("  Saved network: %s\n", basename(file_html)))
  return(net)
}


#' Format the Supplementary Table for publication
#'
#' @param misc_only_list List of MISC_only data frames (one per gene set)
#' @param cell_type_label Label for the cell type
#' @return Combined data frame formatted for Supplementary Table
format_supp_table <- function(misc_only_list, cell_type_label) {
  
  combined <- bind_rows(misc_only_list)
  
  # Select and rename columns for publication
  supp <- combined %>%
    transmute(
      motif       = Motif,
      corr        = round(Corr, 4),
      pval        = round(Corr.P, 6),
      gene        = DORC,
      gene_set    = gene_set,
      logp_corr   = round(-log10(Corr.P), 4)
    ) %>%
    arrange(gene_set, gene, desc(abs(corr)))
  
  supp$cell_type <- cell_type_label
  
  return(supp)
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. MAIN ANALYSIS LOOP
# ══════════════════════════════════════════════════════════════════════════════

# Container for combined supplementary table
all_supp_tables <- list()

for (ct in cell_types) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("  PROCESSING: %s\n", ct))
  cat(strrep("=", 80), "\n")
  
  # ── 3a. Set cell-type-specific gene lists ─────────────────────────────────
  if (ct == "CD8 TEM") {
    IFNG_genes       <- IFNG_genes_CD8_TEM
    Checkpoint_genes <- Checkpoint_genes_CD8_TEM
    ct_label <- "CD8_TEM"
  } else if (ct == "gdT") {
    IFNG_genes       <- IFNG_genes_gdT
    Checkpoint_genes <- Checkpoint_genes_gdT
    ct_label <- "gdT"
  } else {
    # For any other cell type, use core lists
    IFNG_genes       <- IFNG_genes_CORE
    Checkpoint_genes <- Checkpoint_genes_CORE
    ct_label <- gsub(" ", "_", ct)
  }
  # TGF-β is identical for all cell types
  TGFb_genes_ct <- TGFb_genes
  
  # ── 3b. Create cell-type output directory ─────────────────────────────────
  ct_out <- file.path(output_dir, ct_label)
  dir.create(ct_out, recursive = TRUE, showWarnings = FALSE)
  
  # ── 3c. Load FigR results ─────────────────────────────────────────────────
  misc_file    <- file.path(figR_input_dir, paste0(ct, "_MISC_figR_d.txt"))
  control_file <- file.path(figR_input_dir, paste0(ct, "_Control_figR_d.txt"))
  
  if (!file.exists(misc_file) | !file.exists(control_file)) {
    cat(sprintf("  [WARNING] FigR files not found for %s, skipping.\n", ct))
    cat(sprintf("    Expected: %s\n", misc_file))
    cat(sprintf("    Expected: %s\n", control_file))
    next
  }
  
  MISC_figR.d_full    <- read.delim(misc_file)
  Control_figR.d_full <- read.delim(control_file)
  
  cat(sprintf("  Loaded MISC:    %d total TF→DORC edges\n", nrow(MISC_figR.d_full)))
  cat(sprintf("  Loaded Control: %d total TF→DORC edges\n", nrow(Control_figR.d_full)))
  
  # ── 3d. Process each gene set ────────────────────────────────────────────
  
  # IFN-γ pathway
  res_IFNG <- extract_flare_only(
    MISC_figR.d_full, Control_figR.d_full,
    IFNG_genes, "IFNG"
  )
  
  # TGF-β pathway
  res_TGF <- extract_flare_only(
    MISC_figR.d_full, Control_figR.d_full,
    TGFb_genes_ct, "TGFb"
  )
  
  # Checkpoint / Exhaustion
  res_CKP <- extract_flare_only(
    MISC_figR.d_full, Control_figR.d_full,
    Checkpoint_genes, "Checkpoint"
  )
  
  # ── 3e. Save per-gene-set Excel files ─────────────────────────────────────
  cat("\n  Saving Excel files...\n")
  
  for (gs_name in c("IFNG", "TGFb", "Checkpoint")) {
    res <- switch(gs_name, IFNG = res_IFNG, TGFb = res_TGF, Checkpoint = res_CKP)
    
    # All edges (merged MIS-C + Control)
    write.xlsx(res$merged,
               file.path(ct_out, sprintf("%s_all_merged_%s.xlsx", ct_label, gs_name)))
    
    # MIS-C–only edges
    write.xlsx(res$MISC_only,
               file.path(ct_out, sprintf("%s_only_MISC_%s.xlsx", ct_label, gs_name)))
    
    # Control-only edges
    if (nrow(res$Control_only) > 0) {
      write.xlsx(res$Control_only,
                 file.path(ct_out, sprintf("%s_only_Control_%s.xlsx", ct_label, gs_name)))
    }
    
    cat(sprintf("    %s: %d MIS-C–only, %d Control-only, %d shared\n",
                gs_name, nrow(res$MISC_only), nrow(res$Control_only), nrow(res$shared)))
  }
  
  # ── 3f. Save network visualizations ───────────────────────────────────────
  cat("\n  Generating networks...\n")
  
  for (gs_name in c("IFNG", "TGFb", "Checkpoint")) {
    res <- switch(gs_name, IFNG = res_IFNG, TGFb = res_TGF, Checkpoint = res_CKP)
    
    # MIS-C–only network (FIXED: each gene set uses its OWN data)
    plot_and_save_network(
      res$MISC_only,
      file.path(ct_out, sprintf("%s_network_MISC_%s.html", ct_label, gs_name)),
      file.path(ct_out, sprintf("%s_network_MISC_%s.rds", ct_label, gs_name))
    )
    
    # Control-only network
    plot_and_save_network(
      res$Control_only,
      file.path(ct_out, sprintf("%s_network_Control_%s.html", ct_label, gs_name)),
      file.path(ct_out, sprintf("%s_network_Control_%s.rds", ct_label, gs_name))
    )
  }
  
  # ── 3g. Build Supplementary Table for this cell type ──────────────────────
  cat("\n  Building Supplementary Table...\n")
  
  supp <- format_supp_table(
    list(res_IFNG$MISC_only, res_TGF$MISC_only, res_CKP$MISC_only),
    ct_label
  )
  
  write.xlsx(supp, file.path(ct_out, sprintf("Supplementary_Table_FigR_%s.xlsx", ct_label)))
  cat(sprintf("  Supplementary Table: %d total edges\n", nrow(supp)))
  
  # Store for combined output
  all_supp_tables[[ct_label]] <- supp
  
  # ── 3h. Save summary statistics ───────────────────────────────────────────
  summary_all <- bind_rows(res_IFNG$summary, res_TGF$summary, res_CKP$summary)
  write.table(summary_all,
              file.path(ct_out, sprintf("%s_summary_statistics.txt", ct_label)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(sprintf("\n  ✓ %s complete. All outputs in: %s\n", ct, ct_out))
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. COMBINED SUPPLEMENTARY TABLE (ALL CELL TYPES)
# ══════════════════════════════════════════════════════════════════════════════

if (length(all_supp_tables) > 0) {
  combined_supp <- bind_rows(all_supp_tables)
  write.xlsx(combined_supp,
             file.path(output_dir, "Supplementary_Table_FigR_ALL.xlsx"))
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("  COMBINED SUPPLEMENTARY TABLE\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("  Total edges: %d\n", nrow(combined_supp)))
  for (ct_name in names(all_supp_tables)) {
    ct_sub <- combined_supp[combined_supp$cell_type == ct_name, ]
    cat(sprintf("\n  %s:\n", ct_name))
    for (gs in unique(ct_sub$gene_set)) {
      gs_sub <- ct_sub[ct_sub$gene_set == gs, ]
      cat(sprintf("    %s: %d edges (%d act, %d rep)\n",
                  gs, nrow(gs_sub), sum(gs_sub$corr >= 0), sum(gs_sub$corr < 0)))
    }
  }
  cat(sprintf("\n  Saved: %s\n", file.path(output_dir, "Supplementary_Table_FigR_ALL.xlsx")))
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. SESSION INFO
# ══════════════════════════════════════════════════════════════════════════════

cat("\n\n")
cat(strrep("=", 80), "\n")
cat("  SESSION INFO\n")
cat(strrep("=", 80), "\n")
cat(sprintf("  Date: %s\n", Sys.time()))
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  FigR version: %s\n", packageVersion("FigR")))
cat(sprintf("  Output directory: %s\n", output_dir))

sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n  ✓ Analysis complete.\n\n")

