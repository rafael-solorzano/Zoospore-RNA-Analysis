import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from os.path import join as pjoin
from adjustText import adjust_text
import upsetplot

input_folder = r'input' 
temp_folder = r'temp'
output_folder = r'output'

# Helpful website: https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html

# Color Notes:
# # color for lightblue3: #9AC0CD
# # color for gray40: #666666
# # color for indianred3: #CD5555

# Use font Roboto
plt.rcParams['font.sans-serif'] = 'Roboto'
plt.rcParams['font.family'] = 'sans-serif'


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
    plt.scatter(up_df[lfc_colname], -np.log10(up_df[pval_colname]), color='#CD5555', label='Upregulated in Zoospores', s=dp_size, alpha=dp_alpha)
    plt.scatter(down_df[lfc_colname], -np.log10(down_df[pval_colname]), color='#9AC0CD', label='Upregulated in Mats', s=dp_size, alpha=dp_alpha)
    plt.scatter(ns_df[lfc_colname], -np.log10(ns_df[pval_colname]), color='#666666', label='Not Significant', s=dp_size, alpha=dp_alpha)

    plt.axhline(-np.log10(pval_cutoff), color='grey', linestyle='--', dashes=dash_lens)
    plt.axvline(lfc_cutoff, color='grey', linestyle='--', dashes=dash_lens)
    plt.axvline(-lfc_cutoff, color='grey', linestyle='--', dashes=dash_lens)

    plt.xlabel('Log2 Fold-Change', fontsize=ax_label_size)
    plt.ylabel('-Log10(q)', fontsize=ax_label_size)

    plt.xticks(fontsize=ax_tick_size)
    plt.yticks(fontsize=ax_tick_size)
    plt.tick_params(axis='both', which='major', length = ax_tick_len, width = ax_tick_width)
    
    if legend:
        plt.legend(fontsize=leg_size, loc="best")

    # Add a title
    plt.title(title, fontsize=title_size)

    return

def label_points_volcano_plot_proteinIDs(reg_df, pval_colname, num_top_points, lfc_colname):
    # Filter up_df for the top n most significant genes (by p-value)
    reg_df = reg_df.sort_values(by=pval_colname).head(num_top_points)

    texts = []
    for i,r in reg_df.iterrows():
        texts.append(plt.text(r[lfc_colname], -np.log10(r[pval_colname]), 'proteinID \n' + str(r['proteinID']), fontsize=20, ha='center', va='center'))

    adjust_text(texts,arrowprops=dict(arrowstyle='-', color='black',alpha=0.5, shrinkA=5), force_text=0.1, force_points=0.1, expand_points=(1, 1), expand_text=(1, 1), lim=100)
    return

def label_points_volcano_plot_proteinIDs_plus_note(reg_df_up, reg_df_down, pval_colname, num_top_points_up, num_top_points_down, lfc_colname, note_colname):
    # Filter up_df for the top n most significant genes (by p-value)
    reg_df_up = reg_df_up.sort_values(by=pval_colname).head(num_top_points_up)
    reg_df_down = reg_df_down.sort_values(by=pval_colname).head(num_top_points_down)

    texts = []
    for i,r in reg_df_up.iterrows():
        texts.append(plt.text(r[lfc_colname], -np.log10(r[pval_colname]), str(r[note_colname] + ' proteinID\n' + str(r['proteinID'])+'\n'), fontsize=20, ha='center', va='center'))

    for i,r in reg_df_down.iterrows():
        texts.append(plt.text(r[lfc_colname], -np.log10(r[pval_colname]), str(r[note_colname] + ' proteinID\n' + str(r['proteinID'])+'\n'), fontsize=20, ha='center', va='center'))

    # use adjust_text to prevent overlapping text
    adjust_text(texts,arrowprops=dict(arrowstyle='-', color='black',alpha=0.5, shrinkA=5), force_text=0.1, force_points=0.1, expand_points=(1, 1), expand_text=(1, 1), lim=100)
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

# GH48
dge_gh48 = dge_cazymes.copy()
# filter for rows containing the string "Glycoside Hydrolase Family 48" in the column "CAZyme_defline"
dge_gh48 = dge_gh48[dge_gh48["CAZyme_defline"].str.contains("Glycoside Hydrolase Family 48")]

# GH1
dge_gh1 = dge_cazymes.copy()
# filter for rows containing the string "Glycoside Hydrolase Family 1" in the column "CAZyme_defline"
dge_gh1 = dge_gh1[dge_gh1["CAZyme_defline"].str.contains("Glycoside Hydrolase Family 1")]

# GH3
dge_gh3 = dge_cazymes.copy()
# filter for rows containing the string "Glycoside Hydrolase Family 3" in the column "CAZyme_defline"
dge_gh3 = dge_gh3[dge_gh3["CAZyme_defline"].str.contains("Glycoside Hydrolase Family 3")]

# CBM18
dge_cbm18 = dge_cazymes.copy()
# filter for rows containing the string "Carbohydrate-binding module family 18" in the column "CAZyme_defline"
dge_cbm18 = dge_cbm18[dge_cbm18["CAZyme_defline"].str.contains("Carbohydrate-Binding Module Family 18 ")]

# GT8
# filter for rows containing the string "Glycosyltransferase Family 8" in the column "CAZyme_defline"
dge_gt8 = dge_cazymes.copy()
dge_gt8 = dge_gt8[dge_gt8["CAZyme_defline"].str.contains("Glycosyltransferase Family 8")]

# oxidative stress genes, filter the column "KEGG_pathway" for strings containing "Glutathione metabolism"
dge_oxid = dge_summary.copy()
# remove rows with nan values in the column "KEGG_pathway"
dge_oxid = dge_oxid[dge_oxid["KEGG_pathway"].notnull()]
dge_oxid = dge_oxid[dge_oxid["KEGG_pathway"].str.contains("Glutathione metabolism")]

# Proteins with dockerin domain
# Filter dge_summary for rows where column "CAZyme_description" contains either "DOC2" or "CBM10", as this is how I have defined proteins with dockerin domains
dge_dockerin = dge_summary.copy()
# filter for no nan values in the column "CAZyme_description"
dge_dockerin = dge_dockerin[dge_dockerin["CAZyme_description"].notnull()]
dge_dockerin = dge_dockerin[dge_dockerin["CAZyme_description"].str.contains("DOC2|CBM10")]

# Proteins with any glycoside hydrolase, carbohydrate esterase, or polysaccharide lyase annotation (CAZyme_defline contains one of the following strings: "Polysaccharide Lyase", "Glycoside Hydrolase", "Carbohydrate Esterase")
dge_cat_cazymes = dge_summary.copy()
# remove nan values in the column "CAZyme_defline"
dge_cat_cazymes = dge_cat_cazymes[dge_cat_cazymes["CAZyme_defline"].notnull()]
dge_cat_cazymes = dge_cat_cazymes[dge_cat_cazymes["CAZyme_defline"].str.contains("Polysaccharide Lyase|Glycoside Hydrolase|Carbohydrate Esterase")]

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
dpi_val = 1200

# """
# Transcription Factors Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# # """
# # Make volcano plot for transcription factors
# # """
# plot_name = "Volcano_plot_TFs_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 4
# num_top_points_down = 0
# title = "Transcription Factors"
# legend = False

# plt.figure(figsize=(10,10))
# up_tfs, down_tfs = organize_volcano_plot_inputs(dge_tfs, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_tfs, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_tfs, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """


# """
# Secondary Metabolite Core Genes Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_SMs_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 3
# num_top_points_down = 5
# title = "Core Secondary Metabolite Genes"
# legend = False

# plt.figure(figsize=(10,10))
# up_sms, down_sms = organize_volcano_plot_inputs(dge_sms, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs_plus_note(up_sms, down_sms, pval_colname, num_top_points_up, num_top_points_down, lfc_colname, note_colname='SM_cluster_type')


# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# Oxidative Stress Genes Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_oxidative_stress_genes_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 0
# num_top_points_down = 0
# title = "Oxidative Stress Genes"
# legend = False

# plt.figure(figsize=(10,10))
# up_oxid, down_oxid = organize_volcano_plot_inputs(dge_oxid, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_oxid, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_oxid, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# Catabolic CAZymes Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_catabolic_cazymes_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 9
# num_top_points_down = 4
# title = "Catabolic CAZymes"
# legend = False

# plt.figure(figsize=(10,10))
# up_cat_cazymes, down_cat_cazymes = organize_volcano_plot_inputs(dge_cat_cazymes, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_cat_cazymes, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_cat_cazymes, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()

# """
# Scaffoldins Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_scaffoldins_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 0
# num_top_points_down = 5
# # title = "Proteins with Scaffoldin Domains"
# title = "Scaffoldin-containing Proteins"
# legend = False

# plt.figure(figsize=(10,10))
# up_scaffoldins, down_scaffoldins = organize_volcano_plot_inputs(dge_scaffoldins, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_scaffoldins, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_scaffoldins, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """


# """
# GH48 Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_GH48_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 5
# num_top_points_down = 0
# title = "Glycoside Hydrolase Family 48"
# # title = "GH48"
# legend = False

# plt.figure(figsize=(10,10))
# up_gh48, down_gh48 = organize_volcano_plot_inputs(dge_gh48, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_gh48, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_gh48, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# GH1 Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_GH1_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 5
# num_top_points_down = 2
# title = "Glycoside Hydrolase Family 1"
# # title = "GH1"
# legend = False

# plt.figure(figsize=(10,10))
# up_gh1, down_gh1 = organize_volcano_plot_inputs(dge_gh1, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_gh1, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_gh1, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# GH3 Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_GH3_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 5
# num_top_points_down = 2
# title = "Glycoside Hydrolase Family 3"
# # title = "GH3"
# legend = False

# plt.figure(figsize=(10,10))
# up_gh3, down_gh3 = organize_volcano_plot_inputs(dge_gh3, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_gh3, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_gh3, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# CBM18 volcano plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_CBM18_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 8
# num_top_points_down = 4

# title = "Carbohydrate-binding module family 18"
# # title = "CBM18"
# legend = False

# plt.figure(figsize=(10,10))
# up_cbm18, down_cbm18 = organize_volcano_plot_inputs(dge_cbm18, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_cbm18, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_cbm18, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """


# """
# GT8 Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_GT8_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 0
# num_top_points_down = 6

# title = "Glycosyltransferase Family 8"
# # title = "GT8"
# legend = False

# plt.figure(figsize=(10,10))
# up_gt8, down_gt8 = organize_volcano_plot_inputs(dge_gt8, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_gt8, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_gt8, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# Proteins with Dockerin Domain Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_dockerin_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 8
# num_top_points_down = 5
# # title = "Proteins with Dockerin Domains"
# title = "Dockerin-containing Proteins"
# legend = False

# plt.figure(figsize=(10,10))
# up_dockerin, down_dockerin = organize_volcano_plot_inputs(dge_dockerin, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_dockerin, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_dockerin, pval_colname, num_top_points_down, lfc_colname)

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """

# """
# All Genes Volcano Plot
# """
# # *** Remove comments to create volcano Plot ***
# plot_name = "Volcano_plot_all_genes_" + today + ".png"
# # Determine the number of top points to label in either direction
# num_top_points_up = 0
# num_top_points_down = 0
# title = ""
# legend = True

# plt.figure(figsize=(10,10))
# up_all, down_all = organize_volcano_plot_inputs(dge_summary, zoosp_upreg_colname, mat_upreg_colname, lfc_colname, pval_colname, pval_cutoff, lfc_cutoff, dp_size, dp_alpha, ax_label_size, ax_tick_size, ax_tick_len, ax_tick_width, dash_lens, leg_size, legend, title, title_size)
# # Label points
# label_points_volcano_plot_proteinIDs(up_all, pval_colname, num_top_points_up, lfc_colname)
# label_points_volcano_plot_proteinIDs(down_all, pval_colname, num_top_points_down, lfc_colname)

# # For axis number labels, have only the value for 0, middle point, and max value
# plt.xticks(ticks=[-10, 0, 10, 20])
# plt.yticks(ticks=[0, 50, 100])

# # Change axis labels to be larger
# label_size_all_genes = 40
# plt.xlabel('Log2 Fold-Change', fontsize=label_size_all_genes)
# plt.ylabel('-Log10(q)', fontsize=label_size_all_genes)
# plt.xticks(fontsize=label_size_all_genes)
# plt.yticks(fontsize=label_size_all_genes)

# # Make legend a little smaller
# plt.legend(fontsize=20, loc="best")
# # re-fit borders
# plt.tight_layout()

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, plot_name]), dpi=1200)
# plt.show()
# plt.close()
# """
# """


# """
# Generate just the legend as output
# """
# plt.gca().set_axis_off()
# plt.scatter([], [], color='#CD5555', label='Upregulated in Zoospores', s=dp_size, alpha=dp_alpha)
# plt.scatter([], [], color='#9AC0CD', label='Upregulated in Mats', s=dp_size, alpha=dp_alpha)
# plt.scatter([], [], color='#666666', label='Not Significant', s=dp_size, alpha=dp_alpha)
# # Add legend and center it
# plt.legend(fontsize=23, loc="center")

# # Save plot in output folder
# plt.savefig(pjoin(*[output_folder, "Volcano_plot_legend_" + today + ".png"]), dpi=dpi_val)
# plt.show()
# plt.close()
# """
# """


"""
Upset Plot for Carbohydrate Transporter and L-Arabinose Isomerase Annotations
"""
# Create figure
fig = plt.figure(figsize=(10,10))

# Create a venn diagram based on the number of proteinIDs in each list
# l-arabinose isomerase, KEGG, 89 in list
list_l_arabinose_isomerase = [674304, 637955, 697352, 706567, 707593, 423452, 669724, 669727, 669728, 455201, 669729, 669731, 669730, 669734, 669732, 673828, 673830, 395261, 675891, 675896, 419902, 662078, 458855, 666731, 461955, 461961, 644752, 699024, 461981, 464541, 665248, 461988, 500901, 389810, 456374, 461498, 456894, 388295, 385233, 388307, 706771, 388311, 434397, 675044, 675051, 675053, 703727, 705775, 669432, 703224, 703742, 465211, 460610, 411466, 513355, 460621, 448343, 450398, 455007, 450401, 514405, 701804, 500078, 673646, 649073, 411508, 662408, 465803, 455566, 455570, 378777, 458146, 455592, 381867, 395181, 708013, 409010, 701877, 393142, 393146, 267210, 700375, 709089, 679404, 674296, 674297, 674298, 637949, 674303]
# Predicted transporter (major facilitator superfamily), KOG, 133 in list
list_mfs = [674304, 637955, 172552, 697352, 706567, 707593, 423452, 669724, 669727, 669728, 455201, 669729, 669730, 669731, 669732, 669734, 673828, 673830, 708645, 675891, 675896, 419902, 506942, 518207, 662078, 451678, 403043, 458855, 666731, 704109, 461955, 703108, 461961, 28810, 28813, 664205, 644752, 699024, 702612, 461981, 464541, 465056, 665248, 28836, 461988, 500901, 392368, 389810, 708792, 461498, 456894, 451270, 388295, 28874, 337105, 388307, 706771, 388311, 434397, 334558, 122593, 675044, 376042, 675051, 675053, 703727, 705775, 669432, 703224, 703742, 698624, 641812, 679224, 465211, 460610, 411466, 513355, 460621, 456533, 448343, 664919, 324443, 450398, 455007, 450401, 514405, 701804, 500078, 673646, 649073, 411508, 514936, 670074, 662408, 465803, 455566, 455570, 506774, 378777, 669083, 670107, 458146, 410533, 698279, 455592, 381867, 708013, 409010, 701877, 393142, 393146, 700865, 700866, 700867, 619460, 267210, 452043, 637949, 709586, 709587, 700375, 412635, 630237, 701919, 267232, 261604, 397797, 679404, 674296, 674297, 674298, 395261, 674303]
# Sugar transporter, IPR, 88 in list
list_sugar_transporter = [674304, 637955, 697352, 706567, 707593, 423452, 669724, 669727, 669728, 455201, 669729, 669731, 669732, 669734, 669730, 673828, 673830, 395261, 675891, 675896, 419902, 662078, 458855, 666731, 461955, 461961, 644752, 699024, 461981, 464541, 665248, 461988, 389810, 461498, 456894, 388295, 388307, 706771, 388311, 434397, 122593, 675044, 675051, 675053, 703727, 705775, 669432, 703224, 703742, 465211, 460610, 411466, 513355, 460621, 448343, 450398, 455007, 450401, 514405, 701804, 500078, 673646, 649073, 411508, 514936, 662408, 465803, 455566, 455570, 378777, 458146, 410533, 455592, 381867, 708013, 409010, 701877, 393142, 393146, 267210, 700375, 630237, 679404, 674296, 674297, 674298, 637949, 674303]

list_carbohydrate_transport = [674304, 461955, 662408, 461961, 465803, 706567, 707593, 423452, 461981, 464541, 461988, 673828, 673830, 674296, 708013, 389810, 516915, 675891, 393142, 675896, 461498, 393146, 465211, 662078, 460610, 388295, 513355, 516940, 460621, 388307, 706771, 388311, 434397, 122593, 675044, 514405, 675051, 679404, 675053, 673646, 705775, 649073, 514936, 674297, 674298, 674303]
data_annot = upsetplot.from_contents({'L-arabinose isomerase (KEGG)': list_l_arabinose_isomerase, 'Predicted transporter MFS (KOG)': list_mfs, 'Sugar transporter (InterPro)': list_sugar_transporter, 'Carbohydrate transporter (GO)': list_carbohydrate_transport})
ax_dict2 = upsetplot.UpSet(data_annot, subset_size='count').plot()

plt.tight_layout()

# Save the ax_dict2 figure as a png file
plt.savefig(pjoin(*[output_folder, 'upsetplot_carbohydrate_transport_arabinose_isomerase.png']) , dpi=600)
plt.show()