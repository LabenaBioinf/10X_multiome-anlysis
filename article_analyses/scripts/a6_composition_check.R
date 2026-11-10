# A6 (reviewer checklist): recompute paired composition tests from the per-patient
# cell-count CSVs used for Table S4 / S4a, and check exact vs normal-approximation
# p-values, zero/tie handling, BH family size, and absolute counts per subset.
# Input: per-patient count tables cell_type_matrix_P*.csv
#   columns Batch, predicted.id, MISC (= flare), Control (= remission)
# Output: directory given as the second argument

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(in_dir, pattern = "^cell_type_matrix_P[0-9]+\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
names(d) <- c("patient", "cell_type", "flare", "remission")

types    <- sort(unique(d$cell_type))
patients <- unique(d$patient)
patients <- patients[order(as.integer(sub("P", "", patients)))]

# Full patient x cell type grid; a missing row means 0 nuclei in both samples.
grid <- expand.grid(patient = patients, cell_type = types, stringsAsFactors = FALSE)
g <- merge(grid, d, all.x = TRUE)
g$flare[is.na(g$flare)] <- 0
g$remission[is.na(g$remission)] <- 0
g$tot_f <- ave(g$flare, g$patient, FUN = sum)
g$tot_r <- ave(g$remission, g$patient, FUN = sum)
g$p_f <- g$flare / g$tot_f
g$p_r <- g$remission / g$tot_r

# Absolute counts table (checklist: new supplementary table).
counts_wide <- reshape(g[, c("patient", "cell_type", "flare", "remission")],
                       idvar = "cell_type", timevar = "patient", direction = "wide")
write.table(counts_wide, file.path(out_dir, "nuclei_counts_per_subset_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
totals <- unique(g[, c("patient", "tot_f", "tot_r")])
write.table(totals, file.path(out_dir, "nuclei_totals_per_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

test_set <- function(pts, label, min_pair = 0) {
  res <- lapply(types, function(ct) {
    x <- g[g$cell_type == ct & g$patient %in% pts, ]
    x <- x[match(pts, x$patient), ]
    if (min_pair > 0) x <- x[x$flare + x$remission >= min_pair, ]
    if (nrow(x) < 5) return(NULL)
    dif <- x$p_f - x$p_r
    n_zero <- sum(dif == 0)
    n_ties <- sum(duplicated(abs(dif[dif != 0])))
    w_def  <- suppressWarnings(wilcox.test(x$p_f, x$p_r, paired = TRUE))
    w_appr <- suppressWarnings(wilcox.test(x$p_f, x$p_r, paired = TRUE, exact = FALSE, correct = TRUE))
    # Exact permutation (sign-flip) p on the Wilcoxon statistic, valid with zeros/ties
    nz <- dif[dif != 0]
    if (length(nz) == 0) p_perm <- 1 else {
      r <- rank(abs(nz)); obs <- sum(r[nz > 0])
      s <- as.matrix(expand.grid(rep(list(c(0, 1)), length(nz))))
      null <- s %*% r
      mu <- sum(r) / 2
      p_perm <- mean(abs(null - mu) >= abs(obs - mu) - 1e-9)
    }
    min_cells <- min(c(x$flare, x$remission))
    data.frame(set = label, cell_type = ct, n_pairs = nrow(x),
               n_zero_diff = n_zero, n_tied_abs = n_ties,
               min_nuclei_in_any_sample = min_cells,
               n_samples_lt10 = sum(c(x$flare, x$remission) < 10),
               median_flare = median(x$p_f), median_remission = median(x$p_r),
               p_default = w_def$p.value, p_approx = w_appr$p.value, p_exact_perm = p_perm)
  })
  res <- do.call(rbind, res)
  res$fdr_default_30 <- p.adjust(res$p_default, "BH")
  res$fdr_exact_30   <- p.adjust(res$p_exact_perm, "BH")
  res[order(res$p_exact_perm), ]
}

all10 <- test_set(patients, "all")
# The two pairs with the most unbalanced flare nuclei counts (702 and 11,692).
excl  <- test_set(setdiff(patients, c("P8", "P11")), "excl_msP6_msP8")
# Minimum-nuclei sensitivity: drop a pair when flare + remission < 20 nuclei.
thr20 <- test_set(patients, "all_pair_ge20", min_pair = 20)
out <- rbind(all10, excl, thr20)
write.table(out, file.path(out_dir, "composition_tests_recomputed.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("patients:", patients, "\n")
cat("cell types:", length(types), "\n")
cat("min exact two-sided p, n=10:", 2 / 2^10, " n=8:", 2 / 2^8, "\n")
print(totals)
options(width = 200)
print(all10[, c("cell_type", "n_zero_diff", "n_tied_abs", "min_nuclei_in_any_sample",
                "n_samples_lt10", "p_default", "p_approx", "p_exact_perm",
                "fdr_default_30", "fdr_exact_30")], row.names = FALSE, digits = 4)
cat("\n== pair >= 20 nuclei sensitivity\n")
print(thr20[, c("cell_type", "n_pairs", "p_exact_perm", "fdr_exact_30")], row.names = FALSE, digits = 4)
cat("\n== excluding the two most unbalanced pairs\n")
print(excl[, c("cell_type", "n_pairs", "p_default", "p_exact_perm", "fdr_exact_30")], row.names = FALSE, digits = 4)
