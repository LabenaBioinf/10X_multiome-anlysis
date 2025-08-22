# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(writexl)
library(stringr)
library(grid)
library(cowplot)
library(png)


###### Functions

# Calculate the proportion of each cell type
calculate_proportions <- function(data) {
  total_MISC <- sum(data$MISC)  # Total number of cells in MISC group
  total_Control <- sum(data$Control)  # Total number of cells in Control group
  
  data %>%
    group_by(predicted.id) %>%
    summarise(MISC_proportion = sum(MISC) / total_MISC,
              Control_proportion = sum(Control) / total_Control)
}


# Function to sum CD4/CD8 subtypes and calculate their proportions
sum_and_calculate_proportions <- function(data) {
  data %>%
    # Filter for CD4 and CD8 types and categorize them
    filter(predicted.id %in% c(cd4_types, cd8_types)) %>%
    mutate(cell_type_group = ifelse(predicted.id %in% cd4_types, "CD4", "CD8")) %>%
    # Sum the counts for MISC and Control separately
    group_by(cell_type_group) %>%
    summarise(
      MISC_total = sum(MISC, na.rm = TRUE),
      Control_total = sum(Control, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Calculate the total cells in each group
    mutate(
      total_MISC = sum(MISC_total),
      total_Control = sum(Control_total)
    ) %>%
    # Calculate the proportion of CD4 and CD8 cells
    mutate(
      MISC_proportion = MISC_total / total_MISC,
      Control_proportion = Control_total / total_Control
    ) %>%
    select(cell_type_group, MISC_proportion, Control_proportion)
}

# Function to find common genes in the same pathway across clusters
get_common_genes <- function(data) {
  data %>%
    group_by(pathway) %>%
    # Split the leadingEdge genes and collapse them into a single list per pathway
    summarise(all_genes = list(unlist(str_split(leadingEdge, ",")))) %>%
    # Remove duplicates from the gene lists for each pathway
    mutate(common_genes = sapply(all_genes, function(genes) paste(unique(genes), collapse = ", "))) %>%
    select(pathway, common_genes) %>%  # Select relevant columns
    ungroup()
}


######### Code
# Load the datasets
cell_type_matrix_P2 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P2.csv")
cell_type_matrix_P3 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P3.csv")
cell_type_matrix_P4 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P4.csv")
cell_type_matrix_P5 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P5.csv")
cell_type_matrix_P7 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P7.csv")
cell_type_matrix_P8 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P8.csv")
cell_type_matrix_P10 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P10.csv")
cell_type_matrix_P11 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P11.csv")
cell_type_matrix_P12 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P12.csv")
cell_type_matrix_P13 <- read.csv("~/Documents/MISC/cell_counts/cell_type_matrix_P13.csv")


# Apply the function to both datasets
proportions_p2 <- calculate_proportions(cell_type_matrix_P2)
proportions_p3 <- calculate_proportions(cell_type_matrix_P3)
proportions_p4 <- calculate_proportions(cell_type_matrix_P4)
proportions_p5 <- calculate_proportions(cell_type_matrix_P5)
proportions_p7 <- calculate_proportions(cell_type_matrix_P7)
proportions_p8 <- calculate_proportions(cell_type_matrix_P8)
proportions_p10 <- calculate_proportions(cell_type_matrix_P10)
proportions_p11 <- calculate_proportions(cell_type_matrix_P11)
proportions_p12 <- calculate_proportions(cell_type_matrix_P12)
proportions_p13 <- calculate_proportions(cell_type_matrix_P13)


# Merge the datasets on relevant columns (assuming common columns are `sample_id` and `predicted.id`)
combined_proportions <- bind_rows(proportions_p2,
                                  proportions_p3,
                                  proportions_p4,
                                  proportions_p5,
                                  proportions_p7,
                                  proportions_p8,
                                  proportions_p10,
                                  proportions_p11,
                                  proportions_p12,
                                  proportions_p13) %>%
  group_by(predicted.id)


# Perform Wilcoxon test across all datasets
test_results <- combined_proportions %>%
  summarise(
    median_MISC = median(MISC_proportion),
    IQR_MISC = IQR(MISC_proportion),
    median_Control = median(Control_proportion),
    IQR_Control = IQR(Control_proportion),
    p_value = wilcox.test(MISC_proportion, Control_proportion, paired = TRUE)$p.value
  )

# Adjust p-values using Benjamini-Hochberg (BH) correction
test_results <- test_results %>%
  mutate(adj_p_value = p.adjust(p_value, method = "BH"))

write_xlsx(test_results, path = "~/Documents/MISC/cell_counts/count_proportions_test_results.xlsx")


# Filter for significant changes after adjustment
adj_significance_threshold <- 0.05
significant_adj_results <- test_results %>%
  filter(adj_p_value < adj_significance_threshold)

write_xlsx(test_results, path = "~/Documents/MISC/cell_counts/count_proportions_significant_test_results.xlsx")



### Calculate ratio CD4/CD8
cd4_types <- c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM")
cd8_types <- c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM")

# Apply the function to both datasets
cd4_cd8_proportions_p2 <- sum_and_calculate_proportions(cell_type_matrix_P2)
cd4_cd8_proportions_p3 <- sum_and_calculate_proportions(cell_type_matrix_P3)
cd4_cd8_proportions_p4 <- sum_and_calculate_proportions(cell_type_matrix_P4)
cd4_cd8_proportions_p5 <- sum_and_calculate_proportions(cell_type_matrix_P5)
cd4_cd8_proportions_p7 <- sum_and_calculate_proportions(cell_type_matrix_P7)
cd4_cd8_proportions_p8 <- sum_and_calculate_proportions(cell_type_matrix_P8)
cd4_cd8_proportions_p10 <- sum_and_calculate_proportions(cell_type_matrix_P10)
cd4_cd8_proportions_p11 <- sum_and_calculate_proportions(cell_type_matrix_P11)
cd4_cd8_proportions_p12 <- sum_and_calculate_proportions(cell_type_matrix_P12)
cd4_cd8_proportions_p13 <- sum_and_calculate_proportions(cell_type_matrix_P13)


cd4_cd8_combined_proportions <- bind_rows(cd4_cd8_proportions_p2, 
                                          cd4_cd8_proportions_p3,
                                          cd4_cd8_proportions_p4,
                                          cd4_cd8_proportions_p5,
                                          cd4_cd8_proportions_p7,
                                          cd4_cd8_proportions_p8,
                                          cd4_cd8_proportions_p10,
                                          cd4_cd8_proportions_p11,
                                          cd4_cd8_proportions_p12,
                                          cd4_cd8_proportions_p13) 

# Add a row identifier to make merging easier later
cd4_cd8_combined_proportions <- cd4_cd8_combined_proportions %>%
  mutate(row_id = rep(1:(n()/2), each = 2))

# Pivot the data to separate CD4 and CD8 proportions for both MISC and Control
cd4_cd8_wide <- cd4_cd8_combined_proportions %>%
  pivot_wider(names_from = cell_type_group, values_from = c(MISC_proportion, Control_proportion))

# Calculate the CD4/CD8 ratios for MISC and Control
cd4_cd8_wide <- cd4_cd8_wide %>%
  mutate(
    MISC_CD4_CD8_ratio = MISC_proportion_CD4 / MISC_proportion_CD8,
    Control_CD4_CD8_ratio = Control_proportion_CD4 / Control_proportion_CD8
  )

# Combine the ratios into one column for comparison
cd4_cd8_ratios <- cd4_cd8_wide %>%
  pivot_longer(cols = c(MISC_CD4_CD8_ratio, Control_CD4_CD8_ratio), 
               names_to = "Group", values_to = "CD4_CD8_ratio") %>%
  mutate(Group = ifelse(Group == "MISC_CD4_CD8_ratio", "MISC", "Control"))

# Calculate median and IQR for each group
summary_stats <- cd4_cd8_ratios %>%
  group_by(Group) %>%
  summarise(
    Median = median(CD4_CD8_ratio),
    IQR = IQR(CD4_CD8_ratio),
    .groups = 'drop'
  )

wilcox_test_result <- wilcox.test(CD4_CD8_ratio ~ Group, data = cd4_cd8_ratios)
p_value <- wilcox_test_result$p.value
summary_stats <- summary_stats %>% mutate(p_value = p_value)
write_xlsx(summary_stats, "~/Documents/MISC/cell_counts/cd4_cd8_porpotions.xlsx")