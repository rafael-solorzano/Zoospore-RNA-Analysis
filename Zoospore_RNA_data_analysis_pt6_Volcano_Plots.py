import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from os.path import join as pjoin

input_folder = r'input' 
temp_folder = r'temp'
output_folder = r'output'

# Helpful website: https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html

# Color Notes:
# # color for lightblue3: #9AC0CD
# # color for gray40: #666666
# # color for indianred3: #CD5555

# Use font Roboto

"""
Functions
"""
def generate_dge_set(df, dir, reg_colname_up, reg_colname_down):
    """
    Generate a set of differentially expressed genes based on the following criteria.

    Inputs:
    df: pandas dataframe, dge summary data
    dir: up, down, neither
    reg_colname_up: str, column name to describe if a gene is differentially regulated in zoosp (values are 0 or 1)
    reg_colname_down: str, column name to describe if a gene is differentially regulated in mat (values are 0 or 1)

    Outputs:
    dge_set: pandas dataframe, filtered dge summary data for only differentially expressed genes
    """
    # copy the dataframe and filter for only rows with a value in the column reg_colname
    dge_set = df.copy()
    if dir == "up":
        dge_set = dge_set[dge_set[reg_colname_up] == 1]
    elif dir == "down":
        dge_set = dge_set[dge_set[reg_colname_down] == 1]
    elif dir == "neither":
        dge_set = dge_set[(dge_set[reg_colname_up] == 0) & (dge_set[reg_colname_down] == 0)]
    else:
        print("Error: invalid direction. Use 'up' or 'down'.")
        return None
    return dge_set

def generate_volcano_plot(up_tfs, down_tfs, ns_tfs, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, plot_name, output_folder, dp_size=100, dp_alpha=0.5, ax_label_size=30, ax_tick_size=30, ax_tick_len=10, ax_tick_width=2, dash_lens=[10,10], leg_size=25, legend=True):
    plt.scatter(up_tfs[lfc_colname], -np.log10(up_tfs[pval_colname]), color='#CD5555', label='upregulated in zoospores', s=dp_size, alpha=dp_alpha)
    plt.scatter(down_tfs[lfc_colname], -np.log10(down_tfs[pval_colname]), color='#9AC0CD', label='downregulated in zoospores', s=dp_size, alpha=dp_alpha)
    plt.scatter(ns_tfs[lfc_colname], -np.log10(ns_tfs[pval_colname]), color='#666666', label='not significant', s=dp_size, alpha=dp_alpha)

    plt.axhline(-np.log10(pval_cutoff), color='grey', linestyle='--', dashes=dash_lens)
    plt.axvline(lfc_cutoff, color='grey', linestyle='--', dashes=dash_lens)
    plt.axvline(-lfc_cutoff, color='grey', linestyle='--', dashes=dash_lens)

    plt.xlabel('log2 fold-change', fontsize=ax_label_size)
    plt.ylabel('-log10(q)', fontsize=ax_label_size)

    plt.xticks(fontsize=ax_tick_size)
    plt.yticks(fontsize=ax_tick_size)
    plt.tick_params(axis='both', which='major', length = ax_tick_len, width = ax_tick_width)
    
    if legend:
        plt.legend(fontsize=leg_size)

    # Save plot in output folder
    plt.savefig(pjoin(*[output_folder, plot_name]))
    plt.close()
    return

"""
Filenames for import
"""
input_dge_filename = "Supplementary Dataset DGE_summary draft 2.xlsx"
# Notes: 

"""
Column names for volcano plot inputs
"""
# # logFC column = 'log2FC'
lfc_colname = 'log2FC'
# # p-value column = 'padj'
pval_colname = 'padj'
# # Mat TPM avg column = 'mat_tpm_avg'
# mat_tpm_colname = 'mat_tpm_avg'
# # Zoospore TPM avg column = 'zoosp_tpm_avg'
# zoosp_tpm_colname = 'zoosp_tpm_avg'
mat_upreg_colname = "mat_upreg"
zoosp_upreg_colname = "zoosp_upreg"

# """
# Significance cutoffs
# """
pval_cutoff = 0.05
lfc_cutoff = 1
# tpm_cutoff = 1


"""
Import data
"""
dge_summary = pd.read_excel(pjoin(*[input_folder, input_dge_filename]), sheet_name="DGE summary")
dge_scaffoldins = pd.read_excel(pjoin(*[input_folder, input_dge_filename]), sheet_name="Scaffoldins")
dge_tfs = pd.read_excel(pjoin(*[input_folder, input_dge_filename]), sheet_name="TFs")
dge_sms = pd.read_excel(pjoin(*[input_folder, input_dge_filename]), sheet_name="Core Secondary Metabolite Genes")

dge_cazymes = pd.read_excel(pjoin(*[input_folder, input_dge_filename]), sheet_name="DGE summary")
# filter dge_cazymes for only rows with a value in column "CAZyme_defline"
dge_cazymes = dge_cazymes[dge_cazymes["CAZyme_defline"].notnull()]


"""
Make volcano plot for transcription factors
"""
# Get today's date
today = pd.to_datetime('today').strftime('%Y-%m-%d')
plot_name = "Volcano_plot_TFs_" + today + ".png"
up_tfs = generate_dge_set(dge_tfs, "up", zoosp_upreg_colname, mat_upreg_colname)
down_tfs = generate_dge_set(dge_tfs, "down", zoosp_upreg_colname, mat_upreg_colname)
# ns = not significant
ns_tfs = generate_dge_set(dge_tfs, "neither", zoosp_upreg_colname, mat_upreg_colname)

# plot
plt.figure(figsize=(10,10))
# set data point size: 10
dp_size = 100
# set data point transparency: 0.5
dp_alpha = 0.5
# set axis label size: 20
ax_label_size = 30
# set tick label size: 20
ax_tick_size = 30
# set tick length and width
ax_tick_len = 10
ax_tick_width = 2
# increase dash lenght and spacing
dash_lens = [10,10]
# set dashed line for p-value cutoff: 0.05
# set dashed lines for log2FC cutoffs: 1 and -1
# set axis y title: '-log10(q)'
# set axis x title: 'log2 fold-change'
# create legend
# legend font size: 20
leg_size = 25

generate_volcano_plot(up_tfs, down_tfs, ns_tfs, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, plot_name, output_folder, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend=False)



"""
Make volcano plot for transcription factors
"""








