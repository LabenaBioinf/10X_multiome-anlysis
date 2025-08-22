# Retaine only significant gene-peak correlations exclusively present 
# during MIS-C flare, removing any shared correlations

# Libries
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(dplyr)
library(Cairo)
library(FigR)
library(grid)
library(tidyr)
library(openxlsx)



# IFN-γ Pathway: Effector cytokine, transcriptional regulators, downstream mediators
IFNG_genes <- c(
  "IFNG",          # Effector cytokine
  "IFNG-AS1",      # Enhancer lncRNA of IFNG transcription
  "STAT1", "STAT4", # Key signal transducers in IFNG pathway
  "NFATC1",        # Calcium-activated TF involved in early activation and exhaustion
  "EOMES", "TBX21", # Master regulators of CD8+ T effector fate
  "RUNX3", "BATF", "FOSB", # Cooperating transcription factors in Th1 and cytotoxic programming
  "PRDM1", "TCF7L2",       # Terminal differentiation and memory/naïve regulation
  "CXCL9", "CXCL10", "GZMB", # Effector molecules induced by IFN-γ
  "IRF1", "IRF4", "IRF8", # Interferon-response TFs, regulatory modulators
  "BHLHE40", "GATA3", "ZEB2", # TFs influencing differentiation and exhaustion tone
  "ID2", "ID3",             # Negative regulators of terminal differentiation
  "SOCS2", "SOCS3",          # Feedback inhibition of cytokine signaling
  "GBP1", "IFITM1", "PARP9"
)

# TGF-β Pathway: Canonical and non-canonical regulators of immune suppression and tissue homeostasis
TGFb_genes <- c(
  "TGFB1", "TGFBR1", "TGFBR2", # Ligand and receptor complex
  "SMAD2", "SMAD3", "SMAD7",   # Canonical intracellular effectors (SMAD7 is inhibitory)
  "SKI", "SKIL",               # Repressors of SMAD activity
  "PMEPA1", "LTBP2", "THBS1",  # Feedback and ECM-associated regulators
  "RAB31",                     # Trafficking-related; modulates surface presentation
  "ZEB1",                      # Epigenetic and EMT-linked modulator of TGF-β output
  "IL10", "IL1RN",             # Additional immunoregulatory cytokines
  "FOXP3"                      # Master TF of regulatory T cells (Treg)
)

# T cell Exhaustion Markers: Checkpoint molecules and exhaustion-associated TFs
Exhaustion_markers <- c(
  "PDCD1",  # PD-1 checkpoint
  "CTLA4", "HAVCR2", "TIGIT", "LAG3",  # Other co-inhibitory receptors
  "TOX", "NR4A1", "NR4A2", "NR4A3",    # Key transcriptional drivers of exhaustion
  "EOMES", "TBX21", "BATF", "IRF4",    # Transcriptional regulators involved in exhaustion or effector fate
  "CD101", "LAYN", "ENTPD1",          # Markers of terminal exhaustion
  "NFATC1", "NFATC2",                 # NFATs regulate exhaustion programs (with TOX)
  "BCL6",                             # Represses exhaustion program (inverse regulator)
  "MYC"                               # Metabolic counter-regulator of exhaustion
)


############ CD8 TEM
MISC_figR.d_full <- read.delim("~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/cell_types_figR_d/CD8 TEM_MISC_figR_d.txt")
Control_figR.d_full <- read.delim("~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/cell_types_figR_d/CD8 TEM_Control_figR_d.txt")


#FOR SELECTED GENES
merged_data_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_select_IFNG_merged_data_060125.xlsx"
merged_data_only_MISC_Control_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_select_IFNG_only_MISC_060125.xlsx"
network_MISC_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_select_IFNG_network_060125.html"
network_MISC_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_select_IFNG_network_MISC_060125.rds"
network_Control_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_select_IFNG_Control_network_060125.html"
network_Control_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_select_IFNG_network_Control_060125.rds"


# MISC_figR.d <- subset(MISC_figR.d_6, DORC == c('IFNG', 'IFNG-AS1)' & Corr.P < 0.05)
MISC_figR.d <- subset(MISC_figR.d_full, DORC %in% IFNG_genes & Corr.P < 0.05)
Control_figR.d <- subset(Control_figR.d_full, DORC %in% IFNG_genes & Corr.P < 0.05)

# Add a Condition column to each subset
MISC_figR.d_out <- MISC_figR.d %>% mutate(Condition = "MISC")
Control_figR.d_out <- Control_figR.d %>% mutate(Condition = "Control")

# Merge the two data frames
merged_data <- bind_rows(MISC_figR.d_out, Control_figR.d_out)

# Write the merged data to an Excel file
write.xlsx(merged_data, merged_data_file)

## MISC
MISC_sig_activator <- subset(MISC_figR.d, Corr >= 0)
MISC_sig_repressor <- subset(MISC_figR.d, Corr <= 0)

## Control
#### inspect data
Control_sig_activator <- subset(Control_figR.d, Corr >= 0)
Control_sig_repressor <- subset(Control_figR.d, Corr <= 0)


#### Merged Control and MISC

## Sig Activators table
Activators_sig_merged_table <- merge(MISC_sig_activator, Control_sig_activator, by = c("DORC", "Motif"), all = TRUE)
Activators_sig_only_MISC <- anti_join(MISC_sig_activator, Control_sig_activator, by = c("DORC", "Motif"))
Activators_sig_only_Control <- anti_join(Control_sig_activator, MISC_sig_activator, by = c("DORC", "Motif"))
Activators_sig_both <- inner_join(Control_sig_activator, MISC_sig_activator, by = c("DORC", "Motif"))

## Sig Represors table
Repressors_sig_merged_table <- merge(MISC_sig_repressor, Control_sig_repressor, by = c("DORC", "Motif"), all = TRUE)
Repressors_sig_only_MISC <- anti_join(MISC_sig_repressor, Control_sig_repressor, by = c("DORC", "Motif"))
Repressors_sig_only_Control <- anti_join(Control_sig_repressor, MISC_sig_repressor, by = c("DORC", "Motif"))
Repressors_sig_both <- inner_join(Control_sig_repressor, MISC_sig_repressor, by = c("DORC", "Motif"))

Merged_signific_acti_repres_only_MISC <- rbind(Activators_sig_only_MISC, Repressors_sig_only_MISC)
Merged_signific_acti_repres_only_Control <- rbind(Activators_sig_only_Control, Repressors_sig_only_Control)

# Add a Condition column to each subset
Merged_signific_acti_repres_only_MISC_out <- Merged_signific_acti_repres_only_MISC %>% mutate(Condition = "MISC")
Merged_signific_acti_repres_only_Control_out <- Merged_signific_acti_repres_only_Control %>% mutate(Condition = "Control")

# Merge the two data frames
# merged_data_only_MISC_Control <- bind_rows(Merged_signific_acti_repres_only_MISC_out, Merged_signific_acti_repres_only_Control_out)
write.xlsx(Merged_signific_acti_repres_only_MISC_out, merged_data_only_MISC_Control_file)


# plot network of associations
network_MISC <- FigR::plotfigRNetwork(Merged_signific_acti_repres_only_MISC,
                                      score.cut = 0,
                                      TFs = unique(Merged_signific_acti_repres_only_MISC$Motif),
                                      weight.edges = TRUE)
r2d3::save_d3_html(network_MISC, network_MISC_file)
saveRDS(network_MISC, file = network_MISC_object)

# plot network of associations
network_Control <- FigR::plotfigRNetwork(Merged_signific_acti_repres_only_Control,
                                         score.cut = 0,
                                         TFs = unique(Merged_signific_acti_repres_only_Control$Motif),
                                         weight.edges = TRUE)
# r2d3::save_d3_html(network_Control, network_Control_file)
saveRDS(network_Control, file = network_Control_object)



#################################################
#           TGF Beta pathway
################################################
# genes_from_TGF_beta_CD8_TEM_Cells <- c("RAB31",
#                                        "PMEPA1",
#                                        "SKI",
#                                        "THBS1",
#                                        "TGFB1",
#                                        "SKIL",
#                                        "LTBP2",
#                                        "SMAD2",
#                                        "SMAD3",
#                                        "SMAD7",
#                                        "TGFBR1")


merged_data_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_select_TGF_merged_data_TGF_052025.xlsx"
merged_data_only_MISC_file_TGF = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_TGF_only_MISC_TGF_052025.xlsx"
network_MISC_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_TGF_network_052025.html"
network_MISC_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_TGF_network_MISC_052025.rds"
network_Control_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_TGF_Control_network_052025.html"
network_Control_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_TGF_network_Control_052025.rds"



# MISC_figR.d <- subset(MISC_figR.d_6, DORC == c('IFNG', 'IFNG-AS1)' & Corr.P < 0.05)
MISC_figR.d_TGF <- subset(MISC_figR.d_full, DORC %in% TGFb_genes & Corr.P < 0.05)
Control_figR.d_TGF <- subset(Control_figR.d_full, DORC %in% TGFb_genes & Corr.P < 0.05)

# Add a Condition column to each subset
MISC_figR.d_TGF_out <- MISC_figR.d_TGF %>% mutate(Condition = "MISC")
Control_figR.d_TGF_out <- Control_figR.d_TGF %>% mutate(Condition = "Control")

# Merge the two data frames
merged_data_TGF <- bind_rows(MISC_figR.d_TGF_out, Control_figR.d_TGF_out)

# Write the merged data to an Excel file
write.xlsx(merged_data_TGF, merged_data_file)

## MISC
MISC_sig_activator_TGF <- subset(MISC_figR.d_TGF, Corr >= 0)
MISC_sig_repressor_TGF <- subset(Control_figR.d_TGF, Corr <= 0)

## Control
#### inspect data
Control_sig_activator_TGF <- subset(Control_figR.d_TGF, Corr >= 0)
Control_sig_repressor_TGF <- subset(Control_figR.d_TGF, Corr <= 0)


#### Merged Control and MISC

## Sig Activators table
Activators_sig_merged_table_TGF <- merge(MISC_sig_activator_TGF, Control_sig_activator_TGF, by = c("DORC", "Motif"), all = TRUE)
Activators_sig_only_MISC_TGF <- anti_join(MISC_sig_activator_TGF, Control_sig_activator_TGF, by = c("DORC", "Motif"))
Activators_sig_only_Control_TGF <- anti_join(Control_sig_activator_TGF, MISC_sig_activator_TGF, by = c("DORC", "Motif"))
Activators_sig_both_TGF <- inner_join(Control_sig_activator_TGF, MISC_sig_activator_TGF, by = c("DORC", "Motif"))

## Sig Represors table
Repressors_sig_merged_table_TGF <- merge(MISC_sig_repressor_TGF, Control_sig_repressor_TGF, by = c("DORC", "Motif"), all = TRUE)
Repressors_sig_only_MISC_TGF <- anti_join(MISC_sig_repressor_TGF, Control_sig_repressor_TGF, by = c("DORC", "Motif"))
Repressors_sig_only_Control_TGF <- anti_join(Control_sig_repressor_TGF, MISC_sig_repressor_TGF, by = c("DORC", "Motif"))
Repressors_sig_both_TGF <- inner_join(Control_sig_repressor_TGF, MISC_sig_repressor_TGF, by = c("DORC", "Motif"))

Merged_signific_acti_repres_only_MISC_TGF <- rbind(Activators_sig_only_MISC_TGF, Repressors_sig_only_MISC_TGF)
Merged_signific_acti_repres_only_Control_TGF <- rbind(Activators_sig_only_Control_TGF, Repressors_sig_only_Control_TGF)

# Add a Condition column to each subset
# Merged_signific_acti_repres_only_MISC_out <- Merged_signific_acti_repres_only_MISC %>% mutate(Condition = "MISC")
# Merged_signific_acti_repres_only_Control_out <- Merged_signific_acti_repres_only_Control %>% mutate(Condition = "Control")

# Merge the two data frames
#merged_data_only_MISC_TGFl <- bind_rows(Merged_signific_acti_repres_only_MISC_out, Merged_signific_acti_repres_only_Control_out)
write.xlsx(Merged_signific_acti_repres_only_MISC_TGF, merged_data_only_MISC_file_TGF)


# plot network of associations
network_MISC <- FigR::plotfigRNetwork(Merged_signific_acti_repres_only_MISC,
                                      score.cut = 0,
                                      TFs = unique(Merged_signific_acti_repres_only_MISC$Motif),
                                      weight.edges = TRUE)
r2d3::save_d3_html(network_MISC, network_MISC_file)
saveRDS(network_MISC, file = network_MISC_object)

# plot network of associations
network_Control <- FigR::plotfigRNetwork(Merged_signific_acti_repres_only_Control,
                                         score.cut = 0,
                                         TFs = unique(Merged_signific_acti_repres_only_Control$Motif),
                                         weight.edges = TRUE)
r2d3::save_d3_html(network_Control, network_Control_file)
saveRDS(network_Control, file = network_Control_object)



##############################################################################
#           Exostion markers
##############################################################################


# exostion_markers <- c("PDCD1",
#                       "HAVCR2",
#                       "LAG3",
#                       "TIGIT",
#                       "CTLA4",
#                       "TOX",
#                       "NR4A1",
#                       "NR4A2",
#                       "NR4A3",
#                       "EOMES",
#                       "TBX21",
#                       "BATF",
#                       "IRF4",
#                       "TOX")


merged_data_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_select_EXOSTION_M_merged_data_EXOSTION_M_060125.xlsx"
merged_data_only_MISC_file_EXOSTION_M = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_EXOSTION_M_only_MISC_EXOSTION_M_060125.xlsx"
network_MISC_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_EXOSTION_M_network_060125.html"
network_MISC_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_EXOSTION_M_network_MISC_060125.rds"
network_Control_file = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_EXOSTION_M_Control_network_060125.html"
network_Control_object = "~/Documents/MISC/FIGR_RESULTS_1_CELL_TYPE/CD8_TEM_EXOSTION_M_network_Control_060125.rds"



# MISC_figR.d <- subset(MISC_figR.d_6, DORC == c('IFNG', 'IFNG-AS1)' & Corr.P < 0.05)
MISC_figR.d_EXOSTION_M <- subset(MISC_figR.d_full, DORC %in% Exhaustion_markers & Corr.P < 0.05)
Control_figR.d_EXOSTION_M <- subset(Control_figR.d_full, DORC %in% Exhaustion_markers & Corr.P < 0.05)

# Add a Condition column to each subset
MISC_figR.d_EXOSTION_M_out <- MISC_figR.d_EXOSTION_M %>% mutate(Condition = "MISC")
Control_figR.d_EXOSTION_M_out <- Control_figR.d_EXOSTION_M %>% mutate(Condition = "Control")

# Merge the two data frames
merged_data_EXOSTION_M <- bind_rows(MISC_figR.d_EXOSTION_M_out, Control_figR.d_EXOSTION_M_out)

# Write the merged data to an Excel file
write.xlsx(merged_data_EXOSTION_M, merged_data_file)

## MISC
MISC_sig_activator_EXOSTION_M <- subset(MISC_figR.d_EXOSTION_M, Corr >= 0)
MISC_sig_repressor_EXOSTION_M <- subset(Control_figR.d_EXOSTION_M, Corr <= 0)

## Control
#### inspect data
Control_sig_activator_EXOSTION_M <- subset(Control_figR.d_EXOSTION_M, Corr >= 0)
Control_sig_repressor_EXOSTION_M <- subset(Control_figR.d_EXOSTION_M, Corr <= 0)


#### Merged Control and MISC

## Sig Activators table
Activators_sig_merged_table_EXOSTION_M <- merge(MISC_sig_activator_EXOSTION_M, Control_sig_activator_EXOSTION_M, by = c("DORC", "Motif"), all = TRUE)
Activators_sig_only_MISC_EXOSTION_M <- anti_join(MISC_sig_activator_EXOSTION_M, Control_sig_activator_EXOSTION_M, by = c("DORC", "Motif"))
Activators_sig_only_Control_EXOSTION_M <- anti_join(Control_sig_activator_EXOSTION_M, MISC_sig_activator_EXOSTION_M, by = c("DORC", "Motif"))
Activators_sig_both_EXOSTION_M <- inner_join(Control_sig_activator_EXOSTION_M, MISC_sig_activator_EXOSTION_M, by = c("DORC", "Motif"))

## Sig Represors table
Repressors_sig_merged_table_EXOSTION_M <- merge(MISC_sig_repressor_EXOSTION_M, Control_sig_repressor_EXOSTION_M, by = c("DORC", "Motif"), all = TRUE)
Repressors_sig_only_MISC_EXOSTION_M <- anti_join(MISC_sig_repressor_EXOSTION_M, Control_sig_repressor_EXOSTION_M, by = c("DORC", "Motif"))
Repressors_sig_only_Control_EXOSTION_M <- anti_join(Control_sig_repressor_EXOSTION_M, MISC_sig_repressor_EXOSTION_M, by = c("DORC", "Motif"))
Repressors_sig_both_EXOSTION_M <- inner_join(Control_sig_repressor_EXOSTION_M, MISC_sig_repressor_EXOSTION_M, by = c("DORC", "Motif"))

Merged_signific_acti_repres_only_MISC_EXOSTION_M <- rbind(Activators_sig_only_MISC_EXOSTION_M, Repressors_sig_only_MISC_EXOSTION_M)
Merged_signific_acti_repres_only_Control_EXOSTION_M <- rbind(Activators_sig_only_Control_EXOSTION_M, Repressors_sig_only_Control_EXOSTION_M)

write.xlsx(Merged_signific_acti_repres_only_MISC_EXOSTION_M, merged_data_only_MISC_file_EXOSTION_M)


