import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from os.path import join as pjoin
from adjustText import adjust_text

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

def generate_volcano_plot(up_df, down_df, ns_df, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size):
    plt.scatter(up_df[lfc_colname], -np.log10(up_df[pval_colname]), color='#CD5555', label='upregulated in zoospores', s=dp_size, alpha=dp_alpha)
    plt.scatter(down_df[lfc_colname], -np.log10(down_df[pval_colname]), color='#9AC0CD', label='downregulated in zoospores', s=dp_size, alpha=dp_alpha)
    plt.scatter(ns_df[lfc_colname], -np.log10(ns_df[pval_colname]), color='#666666', label='not significant', s=dp_size, alpha=dp_alpha)

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

    # Add a title
    plt.title(title, fontsize=title_size)

    return

def label_points_volcano_plot_proteinIDs(reg_df, pval_colname, num_top_points, lfc_colname):
    # Filter up_df for the top 10 most significant genes (by p-value)
    reg_df = reg_df.sort_values(by=pval_colname).head(num_top_points)

    texts = []
    for i,r in reg_df.iterrows():
        texts.append(plt.text(r[lfc_colname], -np.log10(r[pval_colname]), 'proteinID \n' + str(r['proteinID']), fontsize=20, ha='center', va='center'))

    adjust_text(texts,arrowprops=dict(arrowstyle='-', color='black',alpha=0.5), force_text=0.1, force_points=0.1, expand_points=(1, 1), expand_text=(1, 1), lim=100)
    return

def organize_volcano_plot_inputs(dge_df, zoosp_upreg_colname='zoosp_upreg', mat_upreg_colname='mat_upreg',lfc_colname='log2FC', pval_colname='padj', pval_cutoff=0.05, lfc_cutoff=1, dp_size=100, dp_alpha=0.5, ax_label_size=30, ax_tick_size=30, ax_tick_len=10, ax_tick_width=2, dash_lens=[10,10], leg_size=25, legend=True, title="Volcano Plot", title_size=30):
    up_df = generate_dge_set(dge_df, "up", zoosp_upreg_colname, mat_upreg_colname)
    down_df = generate_dge_set(dge_df, "down", zoosp_upreg_colname, mat_upreg_colname)
    # ns = not significant
    ns_df = generate_dge_set(dge_df, "neither", zoosp_upreg_colname, mat_upreg_colname)

    generate_volcano_plot(up_df, down_df, ns_df, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)

    # Fit the border to the plot
    plt.tight_layout()
    return up_df, down_df

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

# Get today's date
today = pd.to_datetime('today').strftime('%Y-%m-%d')


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
Common Plot Settings
"""
# Default plot settings
# Considered putting a break in the y-axis, but decided against it
# bax = brokenaxes(ylims=((-5, 50), (80, 100)), hspace=.05)
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
title_size = 30


"""
Make volcano plot for transcription factors
"""
plot_name = "Volcano_plot_TFs_" + today + ".png"
# Determine the number of top points to label in either direction
num_top_points_up = 4
num_top_points_down = 0
title = "Transcription Factors"
legend = True

plt.figure(figsize=(10,10))
up_tfs, down_tfs = organize_volcano_plot_inputs(dge_tfs, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# Label points
label_points_volcano_plot_proteinIDs(up_tfs, pval_colname, num_top_points_up, lfc_colname)
label_points_volcano_plot_proteinIDs(down_tfs, pval_colname, num_top_points_down, lfc_colname)

# Save plot in output folder, use dpi=300
plt.savefig(pjoin(*[output_folder, plot_name]), dpi=300)
plt.show()
plt.close()


# """
# Make volcano plot for core secondary metabolite genes
# """
# plot_name = "Volcano_plot_SMS_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 0
# num_top_points_down = 0
# title = "Core Secondary Metabolite Genes"


# # plot









