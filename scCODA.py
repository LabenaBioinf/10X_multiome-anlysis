#!/usr/bin/env python3
"""Runs scCODA

This script allows the user to run scCODA compositional data analysis on scRNA-seq cell count data.

Input files are `cell_type_count.txt` files, generated using `SC_GEX-ATAC.R` script, and saved in the `/results/per_patient/[patient_id]` directory.
The input files are expected to be in tab separated matrix format, with cells as rows, and samples as columns, e.g:
                Control     MISC
B memory        113         12
B naive         150         113
CD14 Mono       13          40

This script requires that `os`, `pandas`, `sccoda`, and `matplotlib` are installed within the Python
environment you are running this script in.
"""


# needed modules
import argparse
import os
import pandas as pd
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
import matplotlib.pyplot as plt


# import command line arguments
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--input_directory', help='Path to the directory you ran SC_GEX-ATAC.R script in.', required=True)

args = parser.parse_args()

###############################################################################

# get sample IDs
samples = []
for entry in os.scandir('./results/per_patient'):
    if entry.is_dir():
        samples.append(entry.name)

# get input files
input_data = []
for sample in samples:
    file_path = args.input_directory + '/results/per_patient/' + sample + '/cell_type_count.txt'
    file = pd.read_csv(file_path, sep = '\t')
    file = file.transpose()
    file.index.name = 'Condition'
    file.reset_index(inplace = True)
    file.insert(loc = 0, column = 'Sample', value = sample)
    file.insert(loc = 2, column = 'Full_ID', value = sample + '_' + file.Condition)
    input_data.append(file)

# merge input files
results = pd.concat(input_data)
results.reset_index(inplace = True, drop = True)
results.fillna(0, inplace = True)
results[results.columns[3:]] = results[results.columns[3:]].astype('int64')

###############################################################################

# convert data to anndata object
sccoda_data = dat.from_pandas(results, covariate_columns = ['Sample', 'Condition', 'Full_ID'])
print(sccoda_data)

# plot boxplot
os.mkdir('./results/per_patient/cell_composition')
viz.boxplots(sccoda_data, feature_name = 'Condition', figsize=(12, 5))
plt.savefig('./results/per_patient/cell_composition/cell_composition_boxplot.png')
#plt.show()

# plot stacked barplot
viz.stacked_barplot(sccoda_data, feature_name = 'Condition', figsize=(8, 9))
plt.savefig('./results/per_patient/cell_composition/cell_composition_barplot.png', bbox_inches = 'tight')
#plt.show()

# finding a potential reference cell type
viz.rel_abundance_dispersion_plot(data = sccoda_data, abundant_threshold = 0.9)
plt.savefig('./results/per_patient/cell_composition/cell_composition_find_reference.png')
#plt.show()

###############################################################################

# set up the model
sccoda_model = mod.CompositionalAnalysis(sccoda_data, formula = 'Condition', reference_cell_type = 'automatic')

# run MCMC
sccoda_results = sccoda_model.sample_hmc()
sccoda_results = sccoda_model.sample_nuts()
# acceptance rate should be between 0.4 and 0.9

# check results
sccoda_results.summary()
sccoda_results.summary_extended()
sccoda_results.credible_effects()

