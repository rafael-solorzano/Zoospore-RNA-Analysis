"""
Zoospore RNA Manuscript: RNA Data Analysis for CAZyme and Cellulosome Annotations

@author: Lazarina Butkovich
"""
import os
from os.path import join as pjoin
import pandas as pd
import numpy as np
import scipy.stats

"""
This script is designed to further format and analyze differential regulation of CAZyme-annotated and cellulosome-annotated genes in the comparison of N. californiae zoospore and fungal mat samples. 
"""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Functions
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def check_ID_column(pd_df, id_name="proteinID"):
    """
    The function check_ID_column checks if the ID column id_name is present in the pandas dataframe, pd_df.

    inputs:
    pd_df: pandas dataframe
    id_name: name of id column, default set to 'proteinID' (string)

    output:
    boolean: True if id_name column in pd_df, exit if not
    """
    if id_name in pd_df.columns:
        return
    else:
        print("Error: no expected " + id_name + " column in pd_df")
        exit()

def make_proteinID_annot_dict(pd_annot, annot_col, id_name="proteinID"):
    """
    The function make_proteinID_annot_dict makes a dictionary with proteinIDs as keys and annotations as values from the annot_col column of the pandas dataframe pd_annot. The function removes duplicate annotations for a given ID. Multiple annotations for a single ID are all kept and separated by a comma in a single value string. Each proteinID key will therefore only have 1 value string. The purpose of this formatting is to prepare the data for addition to the summary dataframe, with the single string as a value in a column.

    inputs:
    pd_annot: pandas dataframe with annotations
    - must have id_name (string) column
    annot_col: specific column name (string) with annotations
    id_name: name of id column, default set to 'proteinID' (string)

    output:
    proteinID_annot_dict: dictionary with IDs as keys and annotations listed (sep=', ') as the values. Remove duplicate annotations for a given ID.
    """
    proteinID_annot_dict = {}
    # if no proteinID (or id_name) column in pd_annot, print error message
    check_ID_column(pd_annot, id_name)
    for i in range(len(pd_annot)):
        proteinID = pd_annot[id_name][i]
        annot = pd_annot[annot_col][i]
        if proteinID in proteinID_annot_dict:
            # if annot is different from proteinID_annot_dict[proteinID] ...
            # (I don't want duplicates of same annotation values)
            if annot not in proteinID_annot_dict[proteinID]:
                # ... then append annot to proteinID_annot_dict[proteinID]
                proteinID_annot_dict[proteinID].append(annot)
        else:
            proteinID_annot_dict[proteinID] = [annot]
    for key in proteinID_annot_dict:
        if len(proteinID_annot_dict[key]) > 1:
            vals = list(map(str, proteinID_annot_dict[key]))
            proteinID_annot_dict[key] = ",".join(vals)
        else:
            proteinID_annot_dict[key] = proteinID_annot_dict[key][0]
    return proteinID_annot_dict

def add_to_df(df, X_annot_pd, annot_cols_list, shared_col='proteinID'):
    """
    The function add_to_df adds annotation columns from the pandas dataframe X_annot_pd to the pandas dataframe df. The function returns df with the new columns. For example, df is the DGE_summary dataframe and X_annot_pd could be the KOG annotations dataframe.

    inputs: 
    df: pandas dataframe with proteinIDs
    X_annot_pd: pandas dataframe with proteinIDs and annotations
    annot_cols_list: list of column names in X_annot_pd to add to DGE_summary (list of strings)
    shared_col: column name in df and X_annot_pd that is shared (string) default set to proteinID

    output: df with new columns (pandas dataframe)
    """
    X_dicts_list = []
    for col in annot_cols_list:
        X_dicts_list.append(make_proteinID_annot_dict(X_annot_pd, col, shared_col))

    # Add annotations to df
    annot_col_num = 0
    for annot_dict in X_dicts_list:
        df[annot_cols_list[annot_col_num]] = df[shared_col].map(annot_dict)
        annot_col_num += 1

    return df

def fisher_exact_verify(fisher_a, fisher_b, fisher_c, fisher_d, p_val_original):
    """
    The function fisher_exact_verify is called in the functions that run the Fisher Exact Tests. It tests that if an additional gene were observed to be upregulated, that the Fisher exact test would return an even lower p-value (ie: you're testing if the result is on the correct side of the two-sided test).

    inputs: a, b, c, d from original Fisher Exact test, p_value from original Fisher Exact test
    outputs: Boolean (true or false) of whether the p-value decreases or increases
    """
    # Suppose that an additional gene in the specific annotation group was upregulated instead of not upregulated. Would the p-value decrease or increase?
    # Descriptions of a, b, c, d:
    # a = # genes upregulated in that annotation group
    # b = # genes upregulated not in that annotation group
    # c = # genes in that annotation group that are not upregulated
    # d = # genes not upregulated and not in that annotation group
    # We get that a increases by 1, b stays the same, c decreases by 1, d stays the same
    if fisher_c == 0:
        # if c = 0, then there are no genes in that annotation group that are not upregulated. So, we are on the correct side of the two-sided test. We don't want to try to run the fisher exact test with a negative number.
        return True
    oddsratio, p_val = scipy.stats.fisher_exact([[fisher_a + 1, fisher_b], [fisher_c - 1, fisher_d]], alternative='two-sided')
    if p_val < p_val_original:
        return True
    else:
        return False

def fisher_exact_test(fisher_a, fisher_b, fisher_c, fisher_d, p_sig=0.05):
    """
    The function fisher_exact_test performs the Statistical Fisher's Exact Test for 2x2 contingency tables. It returns the probability of the observed distribution, assuming that the null hypothesis is true.

    inputs: a, b, c, d for Fisher's Exact test.

    outputs: 
    p_val: p-value of whether the result indicates significant upregulation of the gene annotation or not
    sig: boolean of whether the result is significant or not
    """
    oddsratio, p_val = scipy.stats.fisher_exact([[fisher_a, fisher_b], [fisher_c, fisher_d]], alternative='two-sided')
    sig = fisher_exact_verify(fisher_a, fisher_b, fisher_c, fisher_d, p_val)
    return p_val, sig

def fisher_exact_test_run_for_group(annot_X_pd, num_genes_upreg, num_genes_total, upreg_col_name, pval_cutoff=0.05, annot_genes_min=10, gene_count_col_name='proteinIDs count'):
    """
    The function fisher_exact_test_run_for_group acquires a dataframe of the results of the Fisher Exact test for each group of genes in the input annotation dataframe.

    inputs:
    annot_X_pd: annotation dataframe (ie: annot_KOG_defline_pd)
    num_genes_upreg: number of genes upregulated overall (calculate this number from DGE_summary before calling this function)
    num_genes_total: number of genes total identified for the reference genome(calculate this number from DGE_summary before calling this function)
    upreg_col_name: name of the column in the annotation dataframe that contains the number of genes upregulated in the group (ie: 'zoosp upreg count' or 'mat upreg count'
    pval_cutoff: p-value cutoff for significance (default = 0.05)
    annot_genes_min: minimum number of genes in the annotation group for the test to be conclusive (default = 10)
    gene_count_col_name: name of the column in the annotation dataframe that contains the number of genes in the group (ie: default 'proteinIDs count')

    output: 
    df: pandas dataframe of the results of the Fisher Exact test for each group of genes in the input annotation dataframe. (columns = original columns of annot_X_pd + ['Fisher_a','Fisher_b','Fisher_c','Fisher_d','Fisher_p-value','Fisher_sig_verify','Fisher_above_min_num_genes_cutoff'])
    """
    # Make a new pandas dataframe with the same first columns and values as annot_X_pd
    df = pd.DataFrame(columns=annot_X_pd.columns)
    for col in annot_X_pd.columns:
        df[col] = annot_X_pd[col]
    Fisher_a = []
    Fisher_b = []
    Fisher_b0 = num_genes_upreg
    Fisher_c = []
    Fisher_d = []
    Fisher_d0 = num_genes_total
    Fisher_p_val = []
    Fisher_sig_verify = []
    Fisher_above_min_num_genes_cutoff = []
    for idx in range(len(annot_X_pd)):
        Fisher_a.append(annot_X_pd.iloc[idx][upreg_col_name])
        Fisher_b.append(Fisher_b0 - Fisher_a[idx])
        Fisher_c0 = annot_X_pd.iloc[idx][gene_count_col_name]
        Fisher_c.append(Fisher_c0 - Fisher_a[idx])
        Fisher_d.append(Fisher_d0 - Fisher_a[idx] - Fisher_b[idx] - Fisher_c[idx])
        p_val, sig = fisher_exact_test(Fisher_a[idx], Fisher_b[idx], Fisher_c[idx], Fisher_d[idx])
        if p_val >= pval_cutoff:
            sig = False # Fisher_sig_verify indicates if p_val<0.05 AND if the result is on the correct side of the two-sided test
        Fisher_p_val.append(p_val)
        Fisher_sig_verify.append(sig)
        if Fisher_c0 > annot_genes_min:
            Fisher_above_min_num_genes_cutoff.append(True)
        else:
            Fisher_above_min_num_genes_cutoff.append(False)

    df['Fisher_a'] = Fisher_a
    df['Fisher_b'] = Fisher_b
    df['Fisher_c'] = Fisher_c
    df['Fisher_d'] = Fisher_d
    df['Fisher_p-value'] = Fisher_p_val
    df['Fisher_sig_verify'] = Fisher_sig_verify
    df['Fisher_above_min_num_genes_cutoff'] = Fisher_above_min_num_genes_cutoff
    return df

def count_non_NaN_in_df(df,col_name):
    """
    Count the number of values in the col_name column in df that do not have NaN value
    """
    return df[col_name].count()

def fisher_filter_sort(df):
    """
    The function fisher_filter_sort filters the input df for Fisher_sig_verify=True and Fisher_above_min_num_genes_cutoff=True
    Sort by Fisher_p-value
    input: df from the fisher_start function process
    output: filtered and sorted df
    """
    # Check that df has columns Fisher_sig_verify, Fisher_above_min_num_genes_cutoff, and Fisher_p-value
    if 'Fisher_sig_verify' not in df.columns:
        print('ERROR: df does not have column Fisher_sig_verify')
        return
    if 'Fisher_above_min_num_genes_cutoff' not in df.columns:
        print('ERROR: df does not have column Fisher_above_min_num_genes_cutoff')
        return
    if 'Fisher_p-value' not in df.columns:
        print('ERROR: df does not have column Fisher_p-value')
        return
    # Filter and sort
    df = df[(df['Fisher_sig_verify'] == True) & (df['Fisher_above_min_num_genes_cutoff'] == True)]
    df = df.sort_values(by=['Fisher_p-value'])
    return df

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Changeable values    **************************************
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Import data to append orthologs to
input_folder = r'input' 
temp_folder = r'temp'
output_folder = r'output'

"""
Inputs
"""
# Manually Added Inputs:
# 1) DGE Output Summary:
DGE_summary_filename = "DGE_Summary_Fisher_main_annotations.xlsx"

# Inputs from Previous Scripts (deposited in temp folder)
# 1) dbCAN2 CAZyme predictions:
# dbCAN2 (https://doi.org/10.1093/nar/gky418) is a tool for annotating CAZYme annotations. I predominately used CAZYme annotations from JGI's Mycocosm, but you can also look at dbCAN2 CAZyme predicitons. This would be more helpful for annotating de novo transcripts from this dataset.
# List dbCAN2 data filenames:
dbCAN2_filenames = ['G1_part1_dbCAN_output','G1_part2_dbCAN_output','G1_part3_dbCAN_output','G1_part4_dbCAN_output','G1_part5_dbCAN_output','G1_part6_dbCAN_output']
dbCAN2_file_ext = ".xlsx"

# 2) Cellulosome proteomics data and annotations
# In order to match proteomics data to the genes listed in DGE Summary, we need to compare sequences by BLASTp. This can be accomplished by running a local BLASTp of the proteomics data amino acid sequences against the amino acid sequences for the genes listed in DGE Summary, which can be downloaded from JGI MycoCosm. The BLASTp match cutoffs are as follows: e-value < 1E-5, percent identity > 90%, and “query coverage per HSP” > 75%. Ensure that the BLASTp results (imported here as 'cellulosomes_BLASTp_results.csv') have only one hit per gene and that there are no repeats of the same gene in the BLASTp results (remove duplicates).
cellulosomes_BLASTp_results_filename = "cellulosomes_BLASTp_results.csv"
cellulosome_annot_filename = "G1_cellulosomes_proteomics.csv"


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Consolidate dbCAN2 results
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Import CAZyme annotation data and combine dbCAN2 data into one dataframe
dbCAN2_df = pd.DataFrame()
for filename in dbCAN2_filenames:
    # import dbCAN2 data
    dbCAN2_file = pd.read_excel(pjoin(*[input_folder, filename + dbCAN2_file_ext]), sheet_name=filename)
    # create a column in dbCAN2_file as the first column
    proteinIDs = []
    for i in range(len(dbCAN2_file)):
        # proteinID is the value in the Gene ID column between the 2nd and 3rd "|"
        proteinIDs.append(int(dbCAN2_file['Gene ID'][i].split("|")[2]))
    dbCAN2_file.insert(0, "proteinID", proteinIDs)
    # append dbCAN2_file to dbCAN2_df using pd.concat
    dbCAN2_df = pd.concat([dbCAN2_df, dbCAN2_file], ignore_index=True)

# remove duplicate proteinID rows
dbCAN2_df = dbCAN2_df.drop_duplicates(subset='proteinID', keep='first')

# Sort dbCAN2_df by proteinID
dbCAN2_df = dbCAN2_df.sort_values(by=['proteinID'])

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Append dbCAN2 results to DGE_summary and create output excel
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Except for proteinID and GeneID, rename dbCAN2 columns to start with 'dbCAN2_'
dbCAN2_df = dbCAN2_df.rename(columns={'proteinID':'proteinID','Gene ID':'dbCAN2_GeneID','EC#':'dbCAN2_EC#','HMMER':'dbCAN2_HMMER','dbCAN_sub':'dbCAN2_dbCAN_sub','DIAMOND':'dbCAN2_DIAMOND','Signalp':'dbCAN2_Signalp','#ofTools':'dbCAN2_#ofTools'})

# Import DGE_summary
DGE_summary = pd.read_excel(pjoin(*[temp_folder, DGE_summary_filename]), sheet_name='DGE_summary')

# Add dbCAN2 results to DGE_summary 
new_cols_dbCAN2 = ['dbCAN2_EC#','dbCAN2_HMMER','dbCAN2_dbCAN_sub','dbCAN2_DIAMOND','dbCAN2_Signalp','dbCAN2_#ofTools']
DGE_summary = add_to_df(DGE_summary, dbCAN2_df, new_cols_dbCAN2)
# Move new_cols_dbCAN2 to after the 12th column in DGE_summary
cols = DGE_summary.columns.tolist()
cols = cols[:12] + new_cols_dbCAN2 + cols[12:len(cols)-len(new_cols_dbCAN2)]
DGE_summary = DGE_summary[cols]

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Sheet with all proteinIDs from cellulosome_df
name_out = 'DGE_summary_CAZymes_dbCAN2_cellulosomes'
file_path_out = pjoin(*[output_folder, name_out + '.xlsx'])
# https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')

# Output DGE_summary to excel file without index column
DGE_summary.to_excel(writer, sheet_name='DGE_summary',index=False)

# Output dbCAN2_df to excel file without index column
dbCAN2_df.to_excel(writer, sheet_name='dbCAN2_all_results',index=False)

# Filter DGE_summary for rows with dbCAN2_#ofTools values of 2 or higher
DGE_summary = DGE_summary[DGE_summary['dbCAN2_#ofTools'] >= 2]

# Filter DGE_summary for rows with any value in KOG_defline column
DGE_summary = DGE_summary[DGE_summary['KOG_defline'].notnull()]

# Make a table of the values in KOG_defline with counts of how often they appear
dbCAN2_KOG_defline_counts = DGE_summary['KOG_defline'].value_counts()


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Add Transcriptomic Data to Cellulosome-proteomics table and export
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Import cellulosomes_BLASTp_results.csv from output_folder using pjoin
cellulosomes_BLASTp_results = pd.read_csv(pjoin(input_folder,cellulosomes_BLASTp_results_filename))
# Cellulosome-associated:
cellulosome_df = pd.read_csv(pjoin(*[input_folder, cellulosome_annot_filename]))
# Remove rows with value "#VALUE!" in proteinID column
cellulosome_df = cellulosome_df[cellulosome_df['proteinID'] != '#VALUE!']
# Make proteinID column values integers
cellulosome_df['proteinID'] = cellulosome_df['proteinID'].astype(int)
# Reset index
cellulosome_df = cellulosome_df.reset_index(drop=True)

# Reset DGE_Summmary index
DGE_summary = DGE_summary.reset_index(drop=True)

# Add DGE_summary columns 1-12 and 39-42 to cellulosome_df
# Make list of column names to add to cellulosome_df
new_cols_DGE = DGE_summary.columns.tolist()
new_cols_DGE = new_cols_DGE[1:12] + new_cols_DGE[39:42]
cellulosome_df = add_to_df(cellulosome_df, DGE_summary, new_cols_DGE)

# Rearrange columns. Make the column order be 1, 5, 13-27, 2-4, 6-12
cols = cellulosome_df.columns.tolist()
cols = cols[0:1] + cols[4:5] + cols[12:26] + cols[1:4] + cols[5:11]
cellulosome_df = cellulosome_df[cols]

# Rename proteinID header to proteinID_original
cellulosome_df = cellulosome_df.rename(columns={'proteinID':'proteinID_original'})

cellulosome_df.to_excel(writer, sheet_name='All cellulosome',index=False)

# Add a sheet in the excel for the cellulosomes_BLASTp_results
cellulosomes_BLASTp_results.to_excel(writer, sheet_name='cellulosomes_BLASTp_results',index=False)

writer.close()