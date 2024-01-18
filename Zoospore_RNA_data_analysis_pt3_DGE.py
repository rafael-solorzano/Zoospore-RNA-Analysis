"""
Zoospore RNA Manuscript: RNA Data Analysis DESeq2

@author: Lazarina Butkovich
"""
import os
from os.path import join as pjoin
import pandas as pd
import numpy as np
import scipy.stats
import time
from bioinfokit import analys, visuz
start = time.time()

"""
Script description

This script is designed to analyze the RNA-seq data for N. californiae zoospore vs fungal mat samples. 

"""

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Functions
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def check_proteinID_column(pd_df):
    """
    input:
    pd_df: pandas dataframe

    output:
    boolean: True if proteinID column in pd_df, False if not
    """
    if "proteinID" in pd_df.columns:
        return
    else:
        print("Error: no proteinID column in pandas dataframe")
        exit()

def make_proteinID_annot_dict(pd_annot, annot_col):
    """
    inputs:
    pd_annot: pandas dataframe with annotations
    - must have proteinID column
    annot_col: specific column name (string) with annotations

    output:
    proteinID_annot_dict: dictionary with proteinIDs as keys and annotations as values. Remove duplicate annotations for a given proteinID.
    """
    proteinID_annot_dict = {}
    # if no proteinID column in pd_annot, print error message
    check_proteinID_column(pd_annot)
    for i in range(len(pd_annot)):
        proteinID = pd_annot['proteinID'][i]
        annot = pd_annot[annot_col][i]
        if proteinID in proteinID_annot_dict:
            # if annot is different from proteinID_annot_dict[proteinID] ...
            # I don't want duplicates of same annotation values
            if annot not in proteinID_annot_dict[proteinID]:
                # append annot to proteinID_annot_dict[proteinID]
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

def check_proteinIDs_with_multiple_annotations(X_annot_pd, annot_col_main):
    """
    inputs: 
    X_annot_pd: pandas dataframe with proteinIDs and annotations
    annot_col_main: column name in X_annot_pd to search for multiple annotations for same proteinID

    output: X_dup_proteinIDs_unique_annot list of proteinIDs with multiple annotations
    """
    # Objective: look for proteinIDs that have multiple X annotations, and add these proteinIDs to lists
    # for rows in X_annot_pd with the same proteinID, add proteinID to X_dup_proteinIDs list
    X_dup_proteinIDs = []
    for proteinID in X_annot_pd['proteinID']:
        if len(X_annot_pd[X_annot_pd['proteinID'] == proteinID]) > 1:
            X_dup_proteinIDs.append(proteinID)
    # for rows in X_annot with the same proteinID and different annotation, add proteinID to X_dup_proteinIDs_unique_annot list
    X_dup_proteinIDs_unique_annot = []
    for proteinID in X_dup_proteinIDs:
        if len(X_annot_pd[(X_annot_pd['proteinID'] == proteinID) & (X_annot_pd[annot_col_main].duplicated() == False)]) > 1:
            X_dup_proteinIDs_unique_annot.append(proteinID)
    return X_dup_proteinIDs_unique_annot

def add_to_df(df, X_annot_pd, annot_cols_list, shared_col='proteinID'):
    """
    inputs: 
    df: pandas dataframe with proteinIDs
    X_annot_pd: pandas dataframe with proteinIDs and annotations
    annot_cols_list: list of column names in X_annot_pd to add to DGE_summary (list of strings)

    output: df with new columns (pandas dataframe)
    """
    X_dicts_list = []
    for col in annot_cols_list:
        X_dicts_list.append(make_proteinID_annot_dict(X_annot_pd, col))

    # Add annotations to df
    annot_col_num = 0
    for annot_dict in X_dicts_list:
        df[annot_cols_list[annot_col_num]] = df[shared_col].map(annot_dict)
        annot_col_num += 1

    return df

def make_annot_proteinID_pd(pd_annot, DGE_summary, annot_col):
    """
    inputs:
    pd_annot: pandas dataframe with annotations (ie: KOG_annot)
    - must have proteinID column
    DGE_summary: input the DGE_summary dataframe that you develop in this script; therefore, this function can only be run after DGE_summary is fully created
    annot_col: specific column name (string) to use as keys in output annot_dict

    output:
    annot_dict: dictionary with annotations as keys and proteinIDs as values. Remove duplicate proteinIDs for a given annotation. 
    pd_annot_out: output pandas dataframe with annotation column, proteinIDs column, count of those proteinIDs, counts of "upreg in zoospores" column (info from DGE_summary), counts of "upreg in mats" column (info from DGE_summary)
    # The counts are for performing Fisher Exact test
    """
    annot_dict = {}
    # if no proteinID column in pd_annot, print error message
    check_proteinID_column(pd_annot)
    # iterate through pd_annot, fetching the values in annot_col column and using these values as keys for annot_dict. So, make a list of the consolidated values in annot_col to use as keys for annot_dict.
    annot_keys_list = []
    for i in range(len(pd_annot)):
        annot = pd_annot[annot_col][i]
        if annot not in annot_keys_list:
            annot_keys_list.append(annot)
    # iterate through annot_keys_list, making each value a key in annot_dict
    for annot in annot_keys_list:
        annot_dict[annot] = []
    # iterate through pd_annot, fetching the proteinIDs in proteinID column and using these as values  for annot_dict. For each proteinID, see the annot_col value for that proteinID, and append the proteinID to the list of values for that key in annot_dict. 
    for i in range(len(pd_annot)):
        proteinID = pd_annot['proteinID'][i]
        annot = pd_annot[annot_col][i]
        annot_dict[annot].append(proteinID)
    # remove duplicate proteinIDs for a given annotation
    for annot in annot_dict:
        annot_dict[annot] = list(set(annot_dict[annot]))
        # note: By interating through pd_annot rather than DGE_summary, each row should have 1 proteinID and 1 annotation, whereas DGE_summary can have multiple annotations listed per proteinID (row).
    
    # make a pandas dataframe with the following columns: annotation, proteinIDs (list), count of those proteinIDs, zoosp upreg count column (info from DGE_summary), mat upreg count column (info from DGE_summary)
    header = [annot_col, "proteinIDs", "proteinIDs count", "zoosp upreg count", "mat upreg count"]
    rows = []
    for annot in annot_dict:
        row = [annot]
        proteinIDs = annot_dict[annot]
        proteinIDs_count = 0
        zoosp_upreg_count = 0
        mat_upreg_count = 0
        proteinIDs_in_DGE_summary = []
        for proteinID in proteinIDs:
            # if proteinID is in DGE_summary, then add to zoosp_upreg_count or mat_upreg_count
            if proteinID in DGE_summary['proteinID'].tolist():
                # get the row index for the proteinID in DGE_summary
                row_index = DGE_summary.index[DGE_summary['proteinID'] == proteinID].tolist()[0]
                # get the value for the "zoosp_upreg" column and add to zoosp_upreg_count; value should be 0 or 1
                zoosp_upreg_count += DGE_summary['zoosp_upreg'][row_index]
                # get the value for the "mat_upreg" column and add to mat_upreg_count; value should be 0 or 1
                mat_upreg_count += DGE_summary['mat_upreg'][row_index]
                # Only add proteinIDs to "proteinIDs count" that exist in DGE_summary
                proteinIDs_count += 1
                proteinIDs_in_DGE_summary.append(proteinID)
        proteinIDs = ", ".join(sorted(list(map(str, proteinIDs))))
        row.append(proteinIDs_in_DGE_summary)
        row.append(proteinIDs_count)
        row.append(zoosp_upreg_count)
        row.append(mat_upreg_count)
        rows.append(row)

    pd_annot_out = pd.DataFrame(rows, columns = header)
    pd_annot_out.set_index(annot_col, inplace=True)

    # For list in proteinIDs column, convert to string, with a comma separating each proteinID
    pd_annot_out['proteinIDs'] = pd_annot_out['proteinIDs'].apply(lambda x: ', '.join(map(str, x)))

    # Sort pd_annot_out by annot_col
    pd_annot_out.sort_index(inplace=True)
    return pd_annot_out

def make_annot_proteinID_pd_for_keywords_list(pd_annot, DGE_summary, annot_col,annot_keys_list):
    """
    inputs:
    pd_annot: pandas dataframe with annotations (ie: KOG_annot)
    - must have proteinID column
    DGE_summary: input the DGE_summary dataframe that you develop in this script; therefore, this function can only be run after DGE_summary is fully created
    annot_col: specific column name (string) to use as keys in output annot_dict
    annot_keys_list: list of keywords to look for in annotations

    output:
    annot_dict: dictionary with annotations as keys and proteinIDs as values. Remove duplicate proteinIDs for a given annotation. 
    pd_annot_out: output pandas dataframe with annotation column, proteinIDs column, count of those proteinIDs, counts of "upreg in zoospores" column (info from DGE_summary), counts of "upreg in mats" column (info from DGE_summary)
    # The counts are for performing Fisher Exact test
    """
    annot_dict = {}
    # if no proteinID column in pd_annot, print error message
    check_proteinID_column(pd_annot)
    # iterate through pd_annot, fetching the values in annot_col column and using these values as keys for annot_dict. So, make a list of the consolidated values in annot_col to use as keys for annot_dict.
    # iterate through annot_keys_list, making each value a key in annot_dict
    for annot in annot_keys_list:
        annot_dict[annot] = []
    # iterate through pd_annot, fetching the proteinIDs in proteinID column and using these as values  for annot_dict. For each proteinID, see the annot_col value for that proteinID, and append the proteinID to the list of values for that key in annot_dict. 
    for i in range(len(pd_annot)):
        proteinID = pd_annot['proteinID'][i]
        annot = pd_annot[annot_col][i]
        for key in annot_keys_list:
            if key in annot:
                annot_dict[key].append(proteinID)
        # *** note, this part above is different from function  make_annot_proteinID_pd *** 
    # remove duplicate proteinIDs for a given annotation
    for annot in annot_dict:
        annot_dict[annot] = list(set(annot_dict[annot]))
        # note: By interating through pd_annot rather than DGE_summary, each row should have 1 proteinID and 1 annotation, whereas DGE_summary can have multiple annotations listed per proteinID (row).
    
    # make a pandas dataframe with the following columns: annotation, proteinIDs (list), count of those proteinIDs, zoosp upreg count column (info from DGE_summary), mat upreg count column (info from DGE_summary)
    header = [annot_col, "proteinIDs", "proteinIDs count", "zoosp upreg count", "mat upreg count"]
    rows = []
    for annot in annot_dict:
        row = [annot]
        proteinIDs = annot_dict[annot]
        proteinIDs_count = 0
        zoosp_upreg_count = 0
        mat_upreg_count = 0
        proteinIDs_in_DGE_summary = []
        for proteinID in proteinIDs:
            # if proteinID is in DGE_summary, then add to zoosp_upreg_count or mat_upreg_count
            if proteinID in DGE_summary['proteinID'].tolist():
                # get the row index for the proteinID in DGE_summary
                row_index = DGE_summary.index[DGE_summary['proteinID'] == proteinID].tolist()[0]
                # get the value for the "zoosp_upreg" column and add to zoosp_upreg_count; value should be 0 or 1
                zoosp_upreg_count += DGE_summary['zoosp_upreg'][row_index]
                # get the value for the "mat_upreg" column and add to mat_upreg_count; value should be 0 or 1
                mat_upreg_count += DGE_summary['mat_upreg'][row_index]
                # Only add proteinIDs to "proteinIDs count" that exist in DGE_summary
                proteinIDs_count += 1
                proteinIDs_in_DGE_summary.append(proteinID)
        proteinIDs = ", ".join(sorted(list(map(str, proteinIDs))))
        row.append(proteinIDs_in_DGE_summary)
        row.append(proteinIDs_count)
        row.append(zoosp_upreg_count)
        row.append(mat_upreg_count)
        rows.append(row)

    pd_annot_out = pd.DataFrame(rows, columns = header)
    pd_annot_out.set_index(annot_col, inplace=True)

    # For list in proteinIDs column, convert to string, with a comma separating each proteinID
    pd_annot_out['proteinIDs'] = pd_annot_out['proteinIDs'].apply(lambda x: ', '.join(map(str, x)))

    # Sort pd_annot_out by annot_col
    pd_annot_out.sort_index(inplace=True)
    return pd_annot_out

def pd_annot_expand(pd_annot, col_to_expand_upon, separator):
    """
    inputs:
    pd_annot: pandas dataframe with annotations (ie: KOG_annot)
    col_to_expand_upon: specific column name (string) to expand upon (ie: contains multiple different annotations for a single proteinID)
    separator: string that separates the values in the column to expand upon
    """
    # iterate over pd_annot column col_to_expand_upon, splitting the values in that column by separator. Make new rows for each value in the split list.
    # example: pd_annot = pd.DataFrame({'proteinID': ['P1', 'P2', 'P3'], 'KOG': ['K1;K2;K3', 'K4;K5', 'K6']})
    # pd_annot_expand(pd_annot, 'KOG', ';')
    # output: pd_annot = pd.DataFrame({'proteinID': ['P1', 'P1', 'P1', 'P2', 'P2', 'P3'], 'KOG': ['K1', 'K2', 'K3', 'K4', 'K5', 'K6']})
    # columns of pd_annot_expanded are "proteinID" and col_to_expand_upon
    pd_annot_expanded = pd.DataFrame(columns = ['proteinID', col_to_expand_upon])
    for i in range(len(pd_annot)):
        proteinID = pd_annot['proteinID'][i]
        annot = pd_annot[col_to_expand_upon][i]
        annot_list = annot.split(separator)
        for annot in annot_list:
            # concatenate dataframes, but remove identical rows
            pd_annot_expanded = pd.concat([pd_annot_expanded, pd.DataFrame({'proteinID': [proteinID], col_to_expand_upon: [annot]})], ignore_index = True).drop_duplicates()
    #Sort pd_annot_expanded by col_to_expand_upon
    pd_annot_expanded = pd_annot_expanded.sort_values(by = [col_to_expand_upon])
    return pd_annot_expanded

def fisher_exact_verify(fisher_a, fisher_b, fisher_c, fisher_d, p_val_original):
    """
    Test that if an additional gene were observed to be upregulated, that the Fisher exact test would return an even lower p-value (ie: you're testing if the result is on the correct side of the two-sided test)

    inputs: a, b, c, d from original Fisher Exact test, p_value from original Fisher Exact test
    outputs: true or false of whether the p-value decreases or increases
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
    """Statistical Fisher's exact test for 2x2 contingency tables.
    Returns the probability of the observed distribution, assuming that the null
    hypothesis is true.

    inputs: a, b, c, d for  Fisher Exact test.
    outputs: p-value and true or false of whether the result indicates significant upregulation of the gene annotation or not
    """
    oddsratio, p_val = scipy.stats.fisher_exact([[fisher_a, fisher_b], [fisher_c, fisher_d]], alternative='two-sided')
    sig = fisher_exact_verify(fisher_a, fisher_b, fisher_c, fisher_d, p_val)
    return p_val, sig

def fisher_exact_test_run_for_group(annot_X_pd, num_genes_upreg, num_genes_total, upreg_col_name, pval_cutoff=0.05, annot_genes_min=10, gene_count_col_name='proteinIDs count'):
    """
    Use this function to acquire a dataframe of the results of the Fisher Exact test for each group of genes in the input annotation dataframe.

    inputs:
    annot_X_pd: annotation dataframe (ie: annot_KOG_defline_pd)
    num_genes_upreg: number of genes upregulated overall (calculate this number from DGE_summary before calling this function)
    num_genes_total: number of genes total identified for the reference genome(calculate this number from DGE_summary before calling this function)
    upreg_col_name: name of the column in the annotation dataframe that contains the number of genes upregulated in the group (ie: 'zoosp upreg count' or 'mat upreg count'
    pval_cutoff: p-value cutoff for significance (default = 0.05)
    annot_genes_min: minimum number of genes in the annotation group for the test to be conclusive (default = 10)
    gene_count_col_name: name of the column in the annotation dataframe that contains the number of genes in the group (ie: default 'proteinIDs count')


    output: pandas dataframe of the results of the Fisher Exact test for each group of genes in the input annotation dataframe. (columns = original columns of annot_X_pd + ['Fisher_a','Fisher_b','Fisher_c','Fisher_d','Fisher_p-value','Fisher_sig_verify','Fisher_above_min_num_genes_cutoff'])
    
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

def fisher_start(df, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total, zoosp_col = 'zoosp upreg count', mat_col = 'mat upreg count'):
    # Run fisher_exact_test_run_for_group for zoosp and then for mat. Return both dfs
    df_zoosp = fisher_exact_test_run_for_group(df, num_genes_upreg_zoosp, num_genes_total, zoosp_col)
    df_mat = fisher_exact_test_run_for_group(df, num_genes_upreg_mat, num_genes_total, mat_col)
    return df_zoosp, df_mat

def count_non_NaN_in_df(df,col_name):
    """
    Count the number of values in the col_name column in df that do not have NaN value
    """
    return df[col_name].count()

def fisher_filter_sort(df):
    """
    Filter df for Fisher_sig_verify=True and Fisher_above_min_num_genes_cutoff=True
    Sort by Fisher_p-value
    input: df like Fisher_KOG_mat_upreg
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

def extract_proteinIDs_from_orthologs_column(ortho_list, gf_name):
    # extract protein IDs from orthologs column
    # ortho_list: list of orthologs
    # gf_name: name of the GF
    # return: list of protein IDs, where multiple IDs in a column value are separated by ", " and list of counts of those proteinIDs

    # proteinIDs are written as "jgi|" + gf_name + "|" + proteinID + "|"
    # extract proteinIDs from ortho_list
    final_proteinIDs = []
    counts_list = []
    for proteinIDs_in_string in ortho_list:
        if proteinIDs_in_string == "none":
            final_proteinIDs.append("none")
            counts_list.append(0)
        else:
            final_proteinID_str = ""
            # ProteinIDs in Ortho format are separated by ","
            split_str = proteinIDs_in_string.split(",")
            num = 0
            for substr in split_str:
                if "jgi|" + gf_name + "|" in substr:
                    proteinID = substr.split("|")[2]
                    if num < len(split_str) - 1:
                        final_proteinID_str += proteinID + ", "
                    else:
                        final_proteinID_str += proteinID
                    num += 1
            final_proteinIDs.append(final_proteinID_str)
            counts_list.append(len(split_str))

    return final_proteinIDs, counts_list

def print_stats(df, col, keyword, zoosp_upreg_col="zoosp_upreg",mat_upreg_col="mat_upreg"):
    # Print the number of rows in df that contain keyword in col
    print("Number of proteinIDs with " + keyword + " in " + col + ": " + str(len(df[df[col].str.contains(keyword)])))
    # print sum of values in zoosp_upreg_col for rows that contain keyword in col
    print("Sum of values in zoosp_upreg_col for rows that contain " + keyword + " in " + col + ": " + str(df[df[col].str.contains(keyword)][zoosp_upreg_col].sum()))
    # print sum of values in mat_upreg_col for rows that contain keyword in col
    print("Sum of values in mat_upreg_col for rows that contain " + keyword + " in " + col + ": " + str(df[df[col].str.contains(keyword)][mat_upreg_col].sum()))
    return

def return_class_stats(df, col, keyword, zoosp_upreg_col="zoosp_upreg",mat_upreg_col="mat_upreg"):
    # Print the number of rows in df that contain keyword in col
    CAZyme_total = len(df[df[col].str.contains(keyword)])
    # print sum of values in zoosp_upreg_col for rows that contain keyword in col
    CAZyme_zoosp_upreg = df[df[col].str.contains(keyword)][zoosp_upreg_col].sum()
    # print sum of values in mat_upreg_col for rows that contain keyword in col
    CAZyme_mat_upreg = df[df[col].str.contains(keyword)][mat_upreg_col].sum()
    return CAZyme_total, CAZyme_zoosp_upreg, CAZyme_mat_upreg

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Changeable values    **************************************
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
input_folder = r'input'
temp_folder = r'temp' 
output_folder = r'output'
# RNA sample groupings
mat_samples = ["HHCCW", "HHCCX", "HHCCY", "HHCGA", "HHCGB", "HHCGG", "HHCGH", "HHCGN", "HHCGO", "HHCGP", "HHCGT", "HHCGU", "HHCGW", "HHCGX", "HHCGY"]
zoosp_samples = ["HHCGC", "HHCGS"]
# log2FC cutoff (absolute value)
log2FC_cutoff = 1
# p-value cutoff
pval_cutoff = 0.05
# TPM cutoff
tpm_cutoff = 1
# For orthofinder alignment
GF_NAME = "Neosp1"

    
"""
File names for import
"""
# Manually Added Inputs: 
# 1) Downloaded MycoCosm Annotations: 
    # Multiple annotation files:
KOG_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_KOG.tsv"
GO_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_GO.tsv"
IPR_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_IPR.tsv"
KEGG_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_KEGG.tsv"
# Additional MycoCosm Annotations:
# 2) Secondary metabolites:
SM_annot_filename = "Neosp1_SMs_orthologs.xlsx"
# 3) OrthoFinder data:
Orthologs_filename = "Orthogroups_GF.tsv"
# 4) CAZymes with selected dbCAN2 predictions:
CAZyme_annot_filename = "G1_cazymes_with_dbCAN2.csv"
# 5) Cellulosome-associated proteinIDs:
cellulosome_annot_filename = "G1_cellulosomes_proteomics.csv"

# Inputs from Previous Scripts (deposited in temp folder)
# 1) Zoospore vs fungal mat DESeq2 data:
deseq2_filename = "deseq2_output.csv"
    # From my DESeq2 analysis in RStudio
    # Columns: proteinID, baseMean, mat_vs_zoosp_log2FC, lfcSE, stat, pvalue, mat_vs_zoosp_Padj
    # R v4.2.2
    # Bioconductor v3.16
    # DESeq2 v1.36.0
# 2) DESeq2 Normalized counts:
deseq2_norm_counts_filename = "deseq2_normalized_counts_labeled.csv"
    # Columns: proteinID, ~all sample names~,zoosp_avg_DESeq2_norm_cts, mat_avg_DESeq2_norm_cts, zoosp_var_DESeq2_norm_cts, mat_var_DESeq2_norm_cts, log2FC_check
# 3) Counts:
counts_filename = "counts_RNAseq_updated.csv"
    # From "Zoospore_data_cleanup_pipeline.py"
# 4) TPM Counts:
tpm_filename = "tpm_counts_RNAseq_updated.csv" 
    # From "Zoospore_data_cleanup_pipeline.py"


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Import files
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
deseq2 = pd.read_csv(pjoin(*[temp_folder, deseq2_filename]), sep=',')
deseq2_cts = pd.read_csv(pjoin(*[temp_folder, deseq2_norm_counts_filename]), sep=',')
counts = pd.read_csv(pjoin(*[temp_folder, counts_filename]), sep=',')
tpm = pd.read_csv(pjoin(*[temp_folder, tpm_filename]), sep=',')
KOG_annot = pd.read_csv(pjoin(*[input_folder, KOG_annot_filename]), sep='\t')
GO_annot = pd.read_csv(pjoin(*[input_folder, GO_annot_filename]), sep='\t')
IPR_annot = pd.read_csv(pjoin(*[input_folder, IPR_annot_filename]), sep='\t')
KEGG_annot = pd.read_csv(pjoin(*[input_folder, KEGG_annot_filename]), sep='\t')
# SigP_annot = pd.read_csv(pjoin(*[input_folder, SigP_annot_filename]), sep='\t')
SM_annot = pd.read_excel(pjoin(*[input_folder, SM_annot_filename]), sheet_name='SMs')
df_OrthoFinder = pd.read_csv(pjoin(*[input_folder, Orthologs_filename]), sep="\t", header=0)
CAZyme_annot = pd.read_csv(pjoin(*[input_folder, CAZyme_annot_filename]), sep=',')

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Align dataframes for DGE_Summary
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
Initialize DGE_summary dataframe
"""
# Make proteinID values int or float64 (for decimals, use ~double)
deseq2['proteinID'] = deseq2['proteinID'].astype(int)
deseq2 = deseq2.sort_values(by=['proteinID'])
# Some cond_1v2_log2FC values are "NA" (not a number), so I need to convert them to NaN
deseq2['mat_vs_zoosp_log2FC'] = pd.to_numeric(deseq2['mat_vs_zoosp_log2FC'], errors='coerce')
deseq2['mat_vs_zoosp_log2FC'] = deseq2['mat_vs_zoosp_log2FC'].astype(np.float64)
deseq2['mat_vs_zoosp_Padj'] = pd.to_numeric(deseq2['mat_vs_zoosp_Padj'], errors='coerce')
deseq2['mat_vs_zoosp_Padj'] = deseq2['mat_vs_zoosp_Padj'].astype(np.float64)

# Create DGE_summary dataframe and add DESeq2 data
DGE_summary = deseq2[['proteinID']].copy()
DGE_summary['log2FC'] = deseq2[['mat_vs_zoosp_log2FC']]
DGE_summary['padj'] = deseq2[['mat_vs_zoosp_Padj']]
# Create a column sig that is TRUE if padj<pval_cutoff and FALSE otherwise
DGE_summary['sig'] = DGE_summary['padj'] < pval_cutoff

# Sort DGE_summary by proteinID number
DGE_summary = DGE_summary.sort_values(by=['proteinID'])

# reset index
DGE_summary = DGE_summary.reset_index(drop=True)


"""
Counts
"""
# Make proteinID column by taking number values after 2nd "|" in GeneID
counts['proteinID'] = counts['GeneID'].str.split('|').str[2]
# Make proteinID column into int
counts['proteinID'] = counts['proteinID'].astype(int)
# Sort counts by proteinID number
counts = counts.sort_values(by=['proteinID'])
# Calculate counts averages for sample groups, rounding to 3 decimal places
mat_counts_avg = counts[mat_samples].mean(axis=1).round(3)
zoosp_counts_avg = counts[zoosp_samples].mean(axis=1).round(3)
counts['mat_counts_avg'] = mat_counts_avg
counts['zoosp_counts_avg'] = zoosp_counts_avg
# Align counts dataframe with DGE_summary dataframe based on proteinID
DGE_summary = DGE_summary.merge(counts[['proteinID', 'mat_counts_avg', 'zoosp_counts_avg']], on='proteinID', how='left')


"""
TPM Counts
"""
# Make proteinID column by taking number values after 2nd "|" in GeneID
tpm['proteinID'] = tpm['GeneID'].str.split('|').str[2]
# Make proteinID column into int
tpm['proteinID'] = tpm['proteinID'].astype(int)
# Sort tpm by proteinID number
tpm = tpm.sort_values(by=['proteinID'])
# Calculate tpm averages for sample groups, rounding to 3 decimal places
mat_tpm_avg = tpm[mat_samples].mean(axis=1).round(3)
zoosp_tpm_avg = tpm[zoosp_samples].mean(axis=1).round(3)
tpm['mat_tpm_avg'] = mat_tpm_avg
tpm['zoosp_tpm_avg'] = zoosp_tpm_avg
# Align tpm dataframe with DGE_summary dataframe based on proteinID
DGE_summary = DGE_summary.merge(tpm[['proteinID', 'mat_tpm_avg', 'zoosp_tpm_avg']], on='proteinID', how='left')


"""
DESeq2-normalized counts
"""
# Make deseq2_cts proteinID column values int
deseq2_cts['proteinID'] = deseq2_cts['proteinID'].astype(int)
# Sort deseq2_cts by proteinID number
deseq2_cts = deseq2_cts.sort_values(by=['proteinID'])
# Reset index
deseq2_cts = deseq2_cts.reset_index(drop=True)
# From deseq2_cts, add columns to DGE_summary
DGE_summary['mat_avg_DESeq2_normalized_cts'] = deseq2_cts[['mat_avg_DESeq2_norm_cts']]
DGE_summary['zoosp_avg_DESeq2_normalized_cts'] = deseq2_cts[['zoosp_avg_DESeq2_norm_cts']]


"""
DGE info
"""
# Use cutoffs to determine if proteinID is significantly upregulated in mat or in zoospores
# Note that log2FC < 0 means that the proteinID is upregulated in zoospores
DGE_summary['mat_upreg'] = np.where((DGE_summary['log2FC'] > log2FC_cutoff) & (DGE_summary['padj'] < pval_cutoff) & (DGE_summary['mat_tpm_avg'] > tpm_cutoff), 1, 0)
DGE_summary['zoosp_upreg'] = np.where((DGE_summary['log2FC'] < -log2FC_cutoff) & (DGE_summary['padj'] < pval_cutoff) & (DGE_summary['zoosp_tpm_avg'] > tpm_cutoff), 1, 0)


"""
Export this basic DGE summary
"""
# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

filename_out = "DGE_summary_output_basic.xlsx"
file_path_out = pjoin(*[output_folder, filename_out])

writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')

# Write each dataframe to a different sheet (with no index column)
DGE_summary.to_excel(writer, sheet_name='DGE_summary',index=False)

# Close the Pandas Excel writer and output the Excel file.
# 10/12/23 note: .save() may give an error; use close() instead, since .save() is deprecated
writer.close()


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Align dataframes for DGE_Summary cont'd: Annotations
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
KOG annotations
"""
# Objective: Align KOG annotations by proteinID in DGE_summary and proteinId in KOG_annot
# Important: Some proteinIDs have multiple KOG annotations
# columns of interest: kogid, kogdefline, kogClass, and kogGroup

# rename proteinId to proteinID
KOG_annot = KOG_annot.rename(columns={'proteinId': 'proteinID'})

# if rows in KOG_annot have all identical values, remove duplicates
KOG_annot = KOG_annot.drop_duplicates()

# rename KOG_annot columns
KOG_annot = KOG_annot.rename(columns={'kogid': 'KOG_id', 'kogdefline': 'KOG_defline', 'kogClass': 'KOG_class', 'kogGroup': 'KOG_group'})

# KOG_dup_proteinIDs list of proteinIDs with multiple kogdefline annotations
KOG_dup_proteinIDs = check_proteinIDs_with_multiple_annotations(KOG_annot, 'KOG_defline')

DGE_summary = add_to_df(DGE_summary, KOG_annot, ['KOG_id', 'KOG_defline', 'KOG_class', 'KOG_group'])


"""
GO terms
"""
# Objective: Align GO annotations by proteinID in DGE_summary and #proteinId in GO_annot
# Important: Some proteinIDs have multiple GO annotations
# columns of interest: gotermId, goName, gotermType, goAcc

# rename #proteinId to proteinID
GO_annot = GO_annot.rename(columns={'#proteinId': 'proteinID'})

# if rows in GO_annot have identical values, remove duplicates
GO_annot = GO_annot.drop_duplicates()

# rename GO_annot columns
GO_annot = GO_annot.rename(columns={'gotermId': 'GO_id', 'goName': 'GO_name', 'gotermType': 'GO_type', 'goAcc': 'GO_acc'})

# GO_dup_proteinIDs list of proteinIDs with multiple GO annotations
GO_dup_proteinIDs = check_proteinIDs_with_multiple_annotations(GO_annot, 'GO_name')

DGE_summary = add_to_df(DGE_summary, GO_annot, ['GO_id', 'GO_name', 'GO_type', 'GO_acc'])


"""
IPR (InterPro) annotations
"""
# Objective: Align InterPro annotations by proteinID in DGE_summary and #proteinId in IPR_annot
# Important: Some proteinIDs have multiple InterPro annotations
# columns of interest: iprId, iprDesc, domainDb, domainId, domainDesc, numHits, score

# rename #proteinId to proteinID
IPR_annot = IPR_annot.rename(columns={'#proteinId': 'proteinID'})

# if rows in IPR_annot have identical values, remove duplicates
IPR_annot = IPR_annot.drop_duplicates()

# rename IPR_annot columns
IPR_annot = IPR_annot.rename(columns={'iprId': 'IPR_id', 'iprDesc': 'IPR_desc', 'domainDb': 'domain_db', 'domainId': 'domain_id', 'domainDesc': 'domain_desc', 'numHits': 'num_hits', 'score': 'IPR_score'})

DGE_summary = add_to_df(DGE_summary, IPR_annot, ['IPR_id', 'IPR_desc', 'domain_db', 'domain_id', 'domain_desc', 'num_hits', 'IPR_score'])


"""
KEGG annotations
"""
# Objective: Align KEGG annotations by proteinID in DGE_summary and #proteinId in KEGG_annot
# Important: Some proteinIDs have multiple KEGG annotations
# columns of interest: ecNum, definition, catalyticActivity, pathway, pathway_class, pathway_type

# rename #proteinId to proteinID
KEGG_annot = KEGG_annot.rename(columns={'#proteinId': 'proteinID'})

# if rows in KEGG_annot have identical values, remove duplicates
KEGG_annot = KEGG_annot.drop_duplicates()

# rename KEGG_annot columns
KEGG_annot = KEGG_annot.rename(columns={'ecNum': 'KEGG_ecNum', 'definition': 'KEGG_definition', 'catalyticActivity': 'KEGG_catalyticActivity', 'pathway': 'KEGG_pathway', 'pathway_class': 'KEGG_pathway_class', 'pathway_type': 'KEGG_pathway_type'})

DGE_summary = add_to_df(DGE_summary, KEGG_annot, ['KEGG_ecNum', 'KEGG_definition', 'KEGG_catalyticActivity', 'KEGG_pathway', 'KEGG_pathway_class', 'KEGG_pathway_type'])


"""
CAZymes
"""
# Objective: Align KEGG annotations by proteinID in DGE_summary and proteinID in CAZyme_annot
# Important: Some proteinIDs could have multiple CAZyme annotations
# columns of interest: description, modelnotes, defline

# if rows in CAZyme_annot have identical values, remove duplicates
CAZyme_annot = CAZyme_annot.drop_duplicates()

# rename CAZyme_annot columns
CAZyme_annot = CAZyme_annot.rename(columns={'description': 'CAZyme_description', 'modelnotes': 'CAZyme_modelnotes', 'defline': 'CAZyme_defline'})

DGE_summary = add_to_df(DGE_summary, CAZyme_annot, ['CAZyme_description', 'CAZyme_modelnotes', 'CAZyme_defline'])


"""
Secondary Metabolites (SMs)
"""
# Objective: Align SM annotations by proteinID in DGE_summary and proteinID in SM_annot
# columns of interest: Cluster Id, Cluster Type, Scaffold, Core

# if rows in SM_annot have identical values, remove duplicates
SM_annot = SM_annot.drop_duplicates()

# rename SM_annot columns
SM_annot = SM_annot.rename(columns={'Cluster Id': 'SM_cluster_id', 'Cluster Type': 'SM_cluster_type', 'Scaffold': 'SM_scaffold', 'Core': 'SM_core'})

# reset index for SM_annot
SM_annot = SM_annot.reset_index(drop=True)

# remove rows with nan proteinIDs
SM_annot = SM_annot[SM_annot['proteinID'].notna()]

# SM_annot SM_core values to TRUE or FALSE
# to-do

DGE_summary = add_to_df(DGE_summary, SM_annot, ['SM_cluster_id', 'SM_cluster_type', 'SM_scaffold', 'SM_core'])

print('Finished initializing DGE Summary, took %.2f seconds' % (time.time() - start))
start = time.time()

"""
OrthoFinder data
"""
# Extract the header name of the GF ortholog columns
anasp1_Ortho_header = df_OrthoFinder.columns[1]
caecom1_Ortho_header = df_OrthoFinder.columns[2]
neosp1_Ortho_header = df_OrthoFinder.columns[3]
pirfi3_Ortho_header = df_OrthoFinder.columns[4]

Ortho_headers = [anasp1_Ortho_header, caecom1_Ortho_header, neosp1_Ortho_header, pirfi3_Ortho_header]

if GF_NAME == "Anasp1":
    target_GF_Ortho_header = anasp1_Ortho_header
elif GF_NAME == "Caecom1":
    target_GF_Ortho_header = caecom1_Ortho_header
elif GF_NAME == "Neosp1":
    target_GF_Ortho_header = neosp1_Ortho_header
elif GF_NAME == "Pirfi3":
    target_GF_Ortho_header = pirfi3_Ortho_header
else:
    print("Error: GF name not recognized")

# Some values in df_OrthoFinder have multiple proteinIDs.
# The gene name is formatted like "jgi|GF_NAME|PROTEINID#|..."
# where GF_name is "Anasp1", "Caecom1", "Neosp1", or "Pirfi3"
# These names can be gotten also from the header, as text before the first "_"
# PROTEINID# is the proteinID number
# Gene names are separated by ", "

# Create a new column in DGE_summary for ortholog name to search. The format will be "jgi|" + GF_NAME + "|" + PROTEINID + "|"
DGE_summary["proteinID_Ortho_format"] = "jgi|" + GF_NAME + "|" + DGE_summary["proteinID"].astype(str) + "|"
# proteinIDs_Ortho_format = source column from source annotations data. These are proteinIDs formatted to a substring of how the gene appears in the OrthoFinder data:
proteinIDs_Ortho_format = DGE_summary["proteinID_Ortho_format"].tolist()
# For each row in DGE_summary, search the target_GF_Ortho_header column in df_OrthoFinder for the substring in the proteinID_Ortho_format column

# Make a list GF_ortholog_searchfor the df_orhologs[target_GF_Ortho_header] column.
# GF_ortholog_search= GF column-to-search from OrthoFinder data:
GF_ortholog_search = df_OrthoFinder[target_GF_Ortho_header].tolist()

# These are the lists we will use to fill the new columns in DGE_summary
# Check that the list for this GF matches the expected ()
anasp1_all_orthologs = df_OrthoFinder[anasp1_Ortho_header].tolist()
caecom1_all_orthologs = df_OrthoFinder[caecom1_Ortho_header].tolist()
neosp1_all_orthologs = df_OrthoFinder[neosp1_Ortho_header].tolist()
pirfi3_all_orthologs = df_OrthoFinder[pirfi3_Ortho_header].tolist()

# Make a list of all orthogroups
all_orthogroups = df_OrthoFinder["Orthogroup"].tolist()

anasp1_source_orthologs = []
caecom1_source_orthologs = []
neosp1_source_orthologs = []
pirfi3_source_orthologs = []

index_Ortho = []
source_orthogroups = []
for source_gene in proteinIDs_Ortho_format:
    for ortho_description in GF_ortholog_search:
        index = 0
        filled = 0
        # if ortho_description is not NaN:
        if isinstance(ortho_description, str):
            if source_gene in ortho_description:
                filled = 1
                # append index position of ortho_description to index_Ortho
                index = GF_ortholog_search.index(ortho_description)
                index_Ortho.append(index)
                # Grab all the ortholog info at that index; will later add to DGE_summary as new columns
                # anasp1 column was previously searched for nans and replaced with "none", but now I have to do this with the other GF
                if isinstance(anasp1_all_orthologs[index], str):
                    anasp1_source_orthologs.append(anasp1_all_orthologs[index])
                else:
                    anasp1_source_orthologs.append("none")
                # anasp1_source_orthologs.append(anasp1_all_orthologs[index])
                # if value is type string
                if isinstance(caecom1_all_orthologs[index], str):
                    caecom1_source_orthologs.append(caecom1_all_orthologs[index])
                else:
                    caecom1_source_orthologs.append("none")
                if isinstance(neosp1_all_orthologs[index], str):
                    neosp1_source_orthologs.append(neosp1_all_orthologs[index])
                else:
                    neosp1_source_orthologs.append("none")
                if isinstance(pirfi3_all_orthologs[index], str):
                    pirfi3_source_orthologs.append(pirfi3_all_orthologs[index])
                else:
                    pirfi3_source_orthologs.append("none")
                # Also grab orthogroup name
                source_orthogroups.append(all_orthogroups[index])
                break
    if filled == 0:
        # If no match was found, append "none" too all relevant lists!
        index_Ortho.append("none")
        anasp1_source_orthologs.append("none")
        caecom1_source_orthologs.append("none")
        neosp1_source_orthologs.append("none")
        pirfi3_source_orthologs.append("none")
        source_orthogroups.append("none")

# Add new columns to DGE_summary
DGE_summary["Orthogroup"] = source_orthogroups
DGE_summary["anasp1_ortholog"] = anasp1_source_orthologs
DGE_summary["caecom1_ortholog"] = caecom1_source_orthologs
DGE_summary["neosp1_ortholog"] = neosp1_source_orthologs
DGE_summary["pirfi3_ortholog"] = pirfi3_source_orthologs

# Make more simplified versions of the ortholog columns, where proteinIDs are extracted an separated by ", "
DGE_summary["S4 orthologs"], DGE_summary["S4 ortholog counts"] = extract_proteinIDs_from_orthologs_column(anasp1_source_orthologs, "Anasp1")
DGE_summary["CC orthologs"], DGE_summary["CC ortholog counts"] = extract_proteinIDs_from_orthologs_column(caecom1_source_orthologs, "Caecom1")
DGE_summary["G1 orthologs"], DGE_summary["G1 ortholog counts"] = extract_proteinIDs_from_orthologs_column(neosp1_source_orthologs, "Neosp1")
DGE_summary["PF orthologs"], DGE_summary["PF ortholog counts"] = extract_proteinIDs_from_orthologs_column(pirfi3_source_orthologs, "Pirfi3")

print('Finished adding Orthologs to DGE Summary, took %.2f seconds' % (time.time() - start))
start = time.time()

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Annotation Dictionaries
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
KOG Class
"""
annot_KOG_class_pd = make_annot_proteinID_pd(KOG_annot, DGE_summary, 'KOG_class')
# sort by proteinIDs count
annot_KOG_class_pd = annot_KOG_class_pd.sort_values(by=['proteinIDs count'], ascending=False)

"""
KOG annotations
"""
annot_KOG_defline_pd = make_annot_proteinID_pd(KOG_annot, DGE_summary, 'KOG_defline')
# sort by proteinIDs count
annot_KOG_defline_pd = annot_KOG_defline_pd.sort_values(by=['proteinIDs count'], ascending=False)

"""
GO terms
"""
annot_GO_name_pd = make_annot_proteinID_pd(GO_annot, DGE_summary, 'GO_name')

"""
IPR (InterPro) annotations
"""
annot_IPR_desc_pd = make_annot_proteinID_pd(IPR_annot, DGE_summary, 'IPR_desc')

"""
KEGG annotations
"""
annot_KEGG_desc_pd = make_annot_proteinID_pd(KEGG_annot, DGE_summary, 'KEGG_definition')
annot_KEGG_pathway_pd = make_annot_proteinID_pd(KEGG_annot, DGE_summary, 'KEGG_pathway')
annot_KEGG_pathway_class_pd = make_annot_proteinID_pd(KEGG_annot, DGE_summary, 'KEGG_pathway_class')

"""
CAZymes
"""
# Note for CAZyme deflines; some deflines contain multiple annotations, separated by " / ". Use pd_annot_expand function to expand these deflines into separate rows
CAZyme_defline_pd_expanded = pd_annot_expand(CAZyme_annot, 'CAZyme_defline', ' / ')
# needs more cleaning up; some deflines have " protein" and others don't, but the annotation meaning is the same
# iterate over rows in pd_CAZyme_defline_expanded; if defline contains " protein", remove " protein" from defline.
for index, row in CAZyme_defline_pd_expanded.iterrows():
    if ' protein' in row['CAZyme_defline']:
        CAZyme_defline_pd_expanded.at[index, 'CAZyme_defline'] = row['CAZyme_defline'].replace(' protein', ' ')
# if two rows deflines differ only by a ' ', combine the rows

annot_CAZyme_defline_pd = make_annot_proteinID_pd(CAZyme_defline_pd_expanded, DGE_summary, 'CAZyme_defline')

# Add column for ratio of 'zoosp upreg count' to 'proteinIDs count'
annot_CAZyme_defline_pd['zoosp upreg ratio'] = annot_CAZyme_defline_pd['zoosp upreg count'] / annot_CAZyme_defline_pd['proteinIDs count']
# Round ratio values to 2 decimal places
annot_CAZyme_defline_pd['zoosp upreg ratio'] = annot_CAZyme_defline_pd['zoosp upreg ratio'].round(2)

# Add column for ratio of 'mat upreg count' count' to 'proteinIDs count'
annot_CAZyme_defline_pd['mat upreg ratio'] = annot_CAZyme_defline_pd['mat upreg count'] / annot_CAZyme_defline_pd['proteinIDs count']
# Round ratio values to 2 decimal places
annot_CAZyme_defline_pd['mat upreg ratio'] = annot_CAZyme_defline_pd['mat upreg ratio'].round(2)


"""
Key CAZyme Classes
"""
# Make annot_pd for key CAZyme classes
CAZyme_class_list = ["Carbohydrate Esterase", "Carbohydrate-Binding Module", "Glycoside Hydrolase", "Glycosyltransferase", "Polysaccharide Lyase"]
annot_key_CAZyme_classes_pd = make_annot_proteinID_pd_for_keywords_list(CAZyme_defline_pd_expanded, DGE_summary, 'CAZyme_defline',CAZyme_class_list)

"""
Secondary Metabolites (SMs)
"""
# Note, it is not necessary to run Fisher exact tests for SMs. I will not perform Fisher Exact test for any SM annotations, since the annotation groups are too small.


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Summary statistics
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Number of proteinIDs in genome
num_genes_total = len(DGE_summary.index) #(=19968)

# Number of proteinIDs with KOG annotations
num_genes_total_KOG = count_non_NaN_in_df(DGE_summary, 'KOG_defline') #(=12643)
# Number of proteinIDs with GO annotations
num_genes_total_GO = count_non_NaN_in_df(DGE_summary, 'GO_name') #(=3618)
# Number of proteinIDs with IPR annotations
num_genes_total_IPR = count_non_NaN_in_df(DGE_summary, 'IPR_desc') #(=12375)
# Number of proteinIDs with KEGG annotations
num_genes_total_KEGG = count_non_NaN_in_df(DGE_summary, 'KEGG_definition') #(=4727)
# Number of proteinIDs with KEGG pathway annotations
num_genes_total_KEGG_pathway = count_non_NaN_in_df(DGE_summary, 'KEGG_pathway') #(=4727)
# Number of proteinIDs with KEGG pathway class annotations
num_genes_total_KEGG_pathway_class = count_non_NaN_in_df(DGE_summary, 'KEGG_pathway_class') #(=4727)
# Number of proteinIDs with CAZyme annotations
num_genes_total_CAZyme = count_non_NaN_in_df(DGE_summary, 'CAZyme_defline') #(=1450)

# Number of proteinIDs upregulated in zoospores
num_genes_upreg_zoosp = len(DGE_summary[DGE_summary['zoosp_upreg'] == True].index) #(=3138)
# Number of proteinIDs upregulated in mats
num_genes_upreg_mat = len(DGE_summary[DGE_summary['mat_upreg'] == True].index) #(=2515)

# Number of proteinIDs with log2FC<0 and padj<0.05
num_genes_upreg_zoosp_p_val_sig_only = len(DGE_summary[(DGE_summary['log2FC'] < 0) & (DGE_summary['padj'] < pval_cutoff)].index) #(=4060 checked in excel DGE_summary output)

# Number of proteinIDs with log2FC>0 and padj<0.05
num_genes_upreg_mat_p_val_sig_only = len(DGE_summary[(DGE_summary['log2FC'] > 0) & (DGE_summary['padj'] < pval_cutoff)].index) #(=3922)

# Number of proteinIDs detected in any zoospore
num_genes_detected_zoosp = len(DGE_summary[DGE_summary['zoosp_tpm_avg'] > tpm_cutoff].index) #(=10305)
# Number of proteinIDs detected in any mat
num_genes_detected_mat = len(DGE_summary[DGE_summary['mat_tpm_avg'] > tpm_cutoff].index) #(=12170)
# Number of proteinIDs detected in any of both zoospores and mats
num_genes_detected_zoosp_and_mat = len(DGE_summary[(DGE_summary['zoosp_tpm_avg'] > tpm_cutoff) & (DGE_summary['mat_tpm_avg'] > tpm_cutoff)].index) #(=9657)
# Number of proteinIDs detected in either zoospores and mats
num_genes_detected_zoosp_or_mat = len(DGE_summary[(DGE_summary['zoosp_tpm_avg'] > tpm_cutoff) | (DGE_summary['mat_tpm_avg'] > tpm_cutoff)].index) #(=12818)
# Number of proteinIDs detected in neither zoospores nor mats
num_genes_detected_neither = len(DGE_summary[(DGE_summary['zoosp_tpm_avg'] <= tpm_cutoff) & (DGE_summary['mat_tpm_avg'] <= tpm_cutoff)].index) #(=7150)

# Number of genes upregulated in Zoospores but zoosp_tpm_avg<mat_tpm_avg (likely artefact of sample sizes)
num_genes_upreg_zoosp_TPM_zoosp_avg_less_TPM_mat_avg = len(DGE_summary[(DGE_summary['zoosp_upreg'] == True) & (DGE_summary['zoosp_tpm_avg'] < DGE_summary['mat_tpm_avg'])].index) #(=159)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fisher Exact Test 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# My previous experience with Fisher Exact test involved manually inputting values into a 2x2 table using PRISM (two-sided Fisher Exact test)
# Objective: use Python package scipy.stats.fisher_exact to perform automated Fisher Exact tests on DGE_summary dataframe
# Null hypothesis: observations were sampled at random from the 2 populations

# Create lists for columns in Fisher_summary dataframes : Fisher_mat_upreg_p_val, Fisher_zoosp_upreg_p_val

# Fisher_summary dataframes are a collection of dataframes, where each has a column of unique annotations, a column with lists of proteinIDs with that annotation, counts of how many of those proteinIDs are sig upreg in mats and in zoospores

# Create 2x2 contingency table
# 4 squares of the table: a, b, c, and d, all non-negative integers
# a = number of upreg with specific annotation
# b = number of total genes upregulated - a
# c = number of total genes with specific annotation - a
# d = number of total genes - a - b - c

fisher_p_value_cutoff = 0.05
min_annot_size = 10 # minimum number of proteinIDs with specific annotation, to avoid positive hits largely due to small sample size. min_annot_size = 10 is arbitrary, so you can adjust this number. Others have used small sample sizes with Fisher's Exact Test, but some sources mark these results as questionable (https://influentialpoints.com/Training/Fishers_exact_test_use_and_misuse.htm)
num_genes_total = len(DGE_summary.index) # note, num_genes_total includes unannotated genes (=19968)

# The following functions are used to run Fisher Exact tests for each annotation type: fisher_start, fisher_exact_test_run_for_group, fisher_exact_test, fisher_exact_verify

"""
KOG Class Fisher Exact tests
"""
Fisher_KOG_class_zoosp_upreg_unfiltered, Fisher_KOG_class_mat_upreg_unfiltered = fisher_start(annot_KOG_class_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)
# sort pandas dataframe Fisher_KOG_class_zoosp_upreg by increasing Fisher_p-value
Fisher_KOG_class_zoosp_upreg_unfiltered = Fisher_KOG_class_zoosp_upreg_unfiltered.sort_values(by=['Fisher_p-value'], ascending=True)
# sort pandas dataframe Fisher_KOG_class_mat_upreg by increasing Fisher_p-value
Fisher_KOG_class_mat_upreg_unfiltered = Fisher_KOG_class_mat_upreg_unfiltered.sort_values(by=['Fisher_p-value'], ascending=True)

"""
KOG Fisher Exact tests
"""
Fisher_KOG_zoosp_upreg_unfiltered, Fisher_KOG_mat_upreg_unfiltered = fisher_start(annot_KOG_defline_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for KOG annotations, took %.2f seconds' % ((time.time() - start)))
start = time.time()

"""
GO Fisher Exact tests
"""
Fisher_GO_zoosp_upreg_unfiltered, Fisher_GO_mat_upreg_unfiltered = fisher_start(annot_GO_name_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for GO annotations, took %.2f seconds' % ((time.time() - start)))
start = time.time()

"""
IPR Fisher Exact tests
"""
Fisher_IPR_zoosp_upreg_unfiltered, Fisher_IPR_mat_upreg_unfiltered = fisher_start(annot_IPR_desc_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for IPR annotations, took %.2f seconds' % ((time.time() - start)))
start = time.time()
"""
KEGG Fisher Exact tests
"""
Fisher_KEGG_zoosp_upreg_unfiltered, Fisher_KEGG_mat_upreg_unfiltered = fisher_start(annot_KEGG_desc_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for KEGG annotations, took %.2f seconds' % ((time.time() - start)))
start = time.time()

"""
KEGG Pathway Fisher Exact tests
"""
Fisher_KEGG_pathway_zoosp_upreg_unfiltered, Fisher_KEGG_pathway_mat_upreg_unfiltered = fisher_start(annot_KEGG_pathway_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for KEGG pathway annotations, took %.2f seconds' % ((time.time() - start)))

"""
KEGG Pathway Class Fisher Exact tests
"""
Fisher_KEGG_pathway_class_zoosp_upreg_unfiltered, Fisher_KEGG_pathway_class_mat_upreg_unfiltered = fisher_start(annot_KEGG_pathway_class_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

print('Finished running Fisher\'s Exact Tests for KEGG pathway class annotations, took %.2f seconds' % ((time.time() - start)))

"""
CAZyme Fisher Exact tests
"""
Fisher_CAZyme_zoosp_upreg_unfiltered, Fisher_CAZyme_mat_upreg_unfiltered = fisher_start(annot_CAZyme_defline_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)

"""
CAZyme Key Classes
"""
# Run Fisher stats for key CAZyme classes
Fisher_key_CAZyme_classes_zoosp_upreg, Fisher_key_CAZyme_classes_mat_upreg = fisher_start(annot_key_CAZyme_classes_pd, num_genes_upreg_zoosp, num_genes_upreg_mat, num_genes_total)
# Note, I will not filter these results
Fisher_key_CAZyme_classes_zoosp_upreg.sort_values(by=['Fisher_p-value'])
Fisher_key_CAZyme_classes_mat_upreg.sort_values(by=['Fisher_p-value'])

print('Finished running Fisher\'s Exact Tests for CAZyme annotations, took %.2f seconds' % ((time.time() - start)))
start = time.time()


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Filter Fisher Exact Data for Significant Results
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Filter for Fisher_sig_verify=True and Fisher_above_min_num_genes_cutoff=True
# Sort by increasing p-value
# KOG Class
Fisher_KOG_class_mat_upreg = fisher_filter_sort(Fisher_KOG_class_mat_upreg_unfiltered)
Fisher_KOG_class_zoosp_upreg = fisher_filter_sort(Fisher_KOG_class_zoosp_upreg_unfiltered)
# KOG
Fisher_KOG_mat_upreg = fisher_filter_sort(Fisher_KOG_mat_upreg_unfiltered)
Fisher_KOG_zoosp_upreg = fisher_filter_sort(Fisher_KOG_zoosp_upreg_unfiltered)
# GO
Fisher_GO_mat_upreg = fisher_filter_sort(Fisher_GO_mat_upreg_unfiltered)
Fisher_GO_zoosp_upreg = fisher_filter_sort(Fisher_GO_zoosp_upreg_unfiltered)
# IPR
Fisher_IPR_mat_upreg = fisher_filter_sort(Fisher_IPR_mat_upreg_unfiltered)
Fisher_IPR_zoosp_upreg = fisher_filter_sort(Fisher_IPR_zoosp_upreg_unfiltered)
# KEGG
Fisher_KEGG_mat_upreg = fisher_filter_sort(Fisher_KEGG_mat_upreg_unfiltered)
Fisher_KEGG_zoosp_upreg = fisher_filter_sort(Fisher_KEGG_zoosp_upreg_unfiltered)
# KEGG Pathway
Fisher_KEGG_pathway_mat_upreg = fisher_filter_sort(Fisher_KEGG_pathway_mat_upreg_unfiltered)
Fisher_KEGG_pathway_zoosp_upreg = fisher_filter_sort(Fisher_KEGG_pathway_zoosp_upreg_unfiltered)
# KEGG Pathway Class
Fisher_KEGG_pathway_class_mat_upreg = fisher_filter_sort(Fisher_KEGG_pathway_class_mat_upreg_unfiltered)
Fisher_KEGG_pathway_class_zoosp_upreg = fisher_filter_sort(Fisher_KEGG_pathway_class_zoosp_upreg_unfiltered)
# CAZyme
Fisher_CAZyme_mat_upreg = fisher_filter_sort(Fisher_CAZyme_mat_upreg_unfiltered)
Fisher_CAZyme_zoosp_upreg = fisher_filter_sort(Fisher_CAZyme_zoosp_upreg_unfiltered)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
General Statistics and Tables to Export
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
General Stats
"""
stat_description = []
stat_value = []

stat_description+=['Number of proteinIDs in genome']
stat_value+=[num_genes_total]

stat_description+=['Number of proteinIDs with KOG annotations']
stat_value+=[num_genes_total_KOG]

stat_description+=['Number of proteinIDs with GO annotations']
stat_value+=[num_genes_total_GO]

stat_description+=['Number of proteinIDs with IPR annotations']
stat_value+=[num_genes_total_IPR]

stat_description+=['Number of proteinIDs with KEGG annotations']
stat_value+=[num_genes_total_KEGG]

stat_description+=['Number of proteinIDs with MycoCosm CAZyme annotations']
stat_value+=[num_genes_total_CAZyme]

stat_description+=['Total number of proteinIDs upregulated in zoospores']
stat_value+=[num_genes_upreg_zoosp]

stat_description+=['Total number of proteinIDs upregulated in mats']
stat_value+=[num_genes_upreg_mat]

stat_description+=['Total number of proteinIDs upregulated in zoospores* (q<0.05)']
stat_value+=[num_genes_upreg_zoosp_p_val_sig_only]

stat_description+=['Total number of proteinIDs upregulated in mats* (q<0.05)']
stat_value+=[num_genes_upreg_mat_p_val_sig_only]

stat_description+=['Number of proteinIDs detected in any zoospore']
stat_value+=[num_genes_detected_zoosp]

stat_description+=['Number of proteinIDs detected in any mat']
stat_value+=[num_genes_detected_mat]

stat_description+=['Number of proteinIDs detected in both zoospores and mats']
stat_value+=[num_genes_detected_zoosp_and_mat]

stat_description+=['Number of proteinIDs detected in either zoospores or mats']
stat_value+=[num_genes_detected_zoosp_or_mat]

stat_description+=['Number of proteinIDs detected in neither zoospores nor mats']
stat_value+=[num_genes_detected_neither]

stat_description+=['Number of genes upregulated in Zoospores but zoosp_tpm_avg<mat_tpm_avg']
stat_value+=[num_genes_upreg_zoosp_TPM_zoosp_avg_less_TPM_mat_avg]

# Make pd dataframe out of lists stat_description and stat_value
percent_of_total = [100*stat_value[i]/num_genes_total for i in range(len(stat_value))]
# format values of percent_of_total as percent to 1 decimal place
percent_of_total = ['%.1f%%' % percent_of_total[i] for i in range(len(percent_of_total))]
stats_pd = pd.DataFrame({'Statistic': stat_description, 'Number of Genes': stat_value, 'Percent of Total Genes from Reference Genome': percent_of_total})

# remove index from stats_pd
stats_pd = stats_pd.reset_index(drop=True)


"""
Volcano Plot
"""
# Values from DGE_summary
# Volcano plot y-axis = -log10(padj)
# Volcano plot x-axis = log2FC
# each gene (row) is a point on the plot, so there should be num_genes_total number of points
# Color points red that have: log2FC<1, padj<0.05, and (1) if mat_upreg=1, then mat_tpm_avg>1 or (2) if zoosp_upreg=1, then zoosp_tpm_avg>1
# Color points blue that have: log2FC>1, padj<0.05, and (1) if mat_upreg=1, then mat_tpm_avg>1 or (2) if zoosp_upreg=1, then zoosp_tpm_avg>1

# Make volcano plot
# installed bioinfokit: https://github.com/reneshbedre/bioinfokit
# about volcano plots: https://www.reneshbedre.com/blog/volcano.html

# Copy DGE_summary to volcano_df
volcano_df = DGE_summary.copy()

# Remove values from volcano_df that have blanks/NA in the log2FC column
volcano_df = volcano_df[volcano_df['log2FC'].notna()]
# Remove values from volcano_df that have blanks/NA in the padj column
volcano_df = volcano_df[volcano_df['padj'].notna()]

# Reset index
volcano_df = volcano_df.reset_index(drop=True)

# valpha = transparency of points
# q = p_adjusted
# When show=False, the plot gets saved as a .png image in the same folder as the script 
# color for lightblue3: #9AC0CD
# color for gray40: #666666
# color for indianred3: #CD5555
current_folder = os.getcwd()
os.chdir(output_folder)
visuz.GeneExpression.volcano(df=volcano_df, lfc='log2FC', pv='padj', show=False, plotlegend=True, legendpos='upper right', axtickfontname="Roboto", axlabelfontsize=12, axtickfontsize=12, legendanchor=(1.01, 1.01), color=["#9AC0CD","#666666","#CD5555"], valpha=0.5, axxlabel='log2-fold change', axylabel='-log10(q)', sign_line=True,legendlabels=['upregulated in mats*','no significant regulation','upregulated in zoospores*'],lfc_thr=(log2FC_cutoff,log2FC_cutoff),pv_thr=(pval_cutoff,pval_cutoff))
os.chdir(current_folder)


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Export files
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""
Output Statistics and Data Tables
"""
filename_out = "DGE_Statistics_and_Data_Tables.xlsx"
file_path_out = pjoin(*[output_folder, filename_out])

writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')
stats_pd.to_excel(writer, sheet_name='General',index=True)

Fisher_KOG_class_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KOG_class_zoosp_upreg',index=True)
Fisher_KOG_class_mat_upreg.to_excel(writer, sheet_name='Fisher_KOG_class_mat_upreg',index=True)

# Close the Pandas Excel writer and output the Excel file.
writer.close()

"""
Write each dataframe to a different sheet in output excel
"""
filename_out = "DGE_Summary_output_main.xlsx"
file_path_out = pjoin(*[output_folder, filename_out])

# https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')

# Write each dataframe to a different sheet (with no index column)
DGE_summary.to_excel(writer, sheet_name='DGE_summary',index=False)
Fisher_KOG_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KOG_zoosp_upreg',index=True)
Fisher_KOG_mat_upreg.to_excel(writer, sheet_name='Fisher_KOG_mat_upreg',index=True)
Fisher_GO_zoosp_upreg.to_excel(writer, sheet_name='Fisher_GO_zoosp_upreg',index=True)
Fisher_GO_mat_upreg.to_excel(writer, sheet_name='Fisher_GO_mat_upreg',index=True)
Fisher_IPR_zoosp_upreg.to_excel(writer, sheet_name='Fisher_IPR_zoosp_upreg',index=True)
Fisher_IPR_mat_upreg.to_excel(writer, sheet_name='Fisher_IPR_mat_upreg',index=True)
Fisher_KEGG_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_zoosp_upreg',index=True)
Fisher_KEGG_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_mat_upreg',index=True)
Fisher_KEGG_pathway_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathway_zoosp_upreg',index=True)
Fisher_KEGG_pathway_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathway_mat_upreg',index=True)
Fisher_KEGG_pathway_class_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathclass_z_upreg',index=True)
Fisher_KEGG_pathway_class_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathclass_m_upreg',index=True)
Fisher_CAZyme_zoosp_upreg.to_excel(writer, sheet_name='Fisher_CAZyme_zoosp_upreg',index=True)
Fisher_CAZyme_mat_upreg.to_excel(writer, sheet_name='Fisher_CAZyme_mat_upreg',index=True)

# Close the Pandas Excel writer and output the Excel file.
writer.close()

# Also write DGE_summary_output to the temp folder, for access to downstream scripts (copy-pasted above section and changed output_folder to temp_folder)
file_path_out = pjoin(*[temp_folder, filename_out])

# https://xlsxwriter.readthedocs.io/example_pandas_multiple.html
writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')

# Write each dataframe to a different sheet (with no index column)
DGE_summary.to_excel(writer, sheet_name='DGE_summary',index=False)
Fisher_KOG_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KOG_zoosp_upreg',index=True)
Fisher_KOG_mat_upreg.to_excel(writer, sheet_name='Fisher_KOG_mat_upreg',index=True)
Fisher_GO_zoosp_upreg.to_excel(writer, sheet_name='Fisher_GO_zoosp_upreg',index=True)
Fisher_GO_mat_upreg.to_excel(writer, sheet_name='Fisher_GO_mat_upreg',index=True)
Fisher_IPR_zoosp_upreg.to_excel(writer, sheet_name='Fisher_IPR_zoosp_upreg',index=True)
Fisher_IPR_mat_upreg.to_excel(writer, sheet_name='Fisher_IPR_mat_upreg',index=True)
Fisher_KEGG_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_zoosp_upreg',index=True)
Fisher_KEGG_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_mat_upreg',index=True)
Fisher_KEGG_pathway_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathway_zoosp_upreg',index=True)
Fisher_KEGG_pathway_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathway_mat_upreg',index=True)
Fisher_KEGG_pathway_class_zoosp_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathclass_z_upreg',index=True)
Fisher_KEGG_pathway_class_mat_upreg.to_excel(writer, sheet_name='Fisher_KEGG_pathclass_m_upreg',index=True)
Fisher_CAZyme_zoosp_upreg.to_excel(writer, sheet_name='Fisher_CAZyme_zoosp_upreg',index=True)
Fisher_CAZyme_mat_upreg.to_excel(writer, sheet_name='Fisher_CAZyme_mat_upreg',index=True)

# Close the Pandas Excel writer and output the Excel file.
writer.close()


"""
Output CAZymes defline table
"""
filename_out = "DGE_Summary_CAZymes_output.xlsx"
file_path_out = pjoin(*[output_folder, filename_out])

# First, output annot_CAZyme_defline_pd sorted by proteinIDs count descending
annot_CAZyme_defline_pd = annot_CAZyme_defline_pd.sort_values(by=['proteinIDs count'], ascending=False)

writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')
annot_CAZyme_defline_pd.to_excel(writer, sheet_name='CAZymes_defline_count_sorted',index=True)

# Second output annot_CAZyme_defline_pd sorted by zoosp upreg ratio descending
# kind='mergesort': ties from the sort are broken based on previous order (in this case, descending total count)
annot_CAZyme_defline_pd = annot_CAZyme_defline_pd.sort_values(by=['zoosp upreg ratio'], kind='mergesort', ascending=False)

annot_CAZyme_defline_pd.to_excel(writer, sheet_name='CAZymes_defline_zoosp_upreg',index=True)

# reset sort by count to allow mergesort in sort for 'mat upreg ratio'
annot_CAZyme_defline_pd = annot_CAZyme_defline_pd.sort_values(by=['proteinIDs count'], ascending=False)

# Third, output annot_CAZyme_defline_pd sorted by mat upreg ratio descending
annot_CAZyme_defline_pd = annot_CAZyme_defline_pd.sort_values(by=['mat upreg ratio'], kind='mergesort', ascending=False)

annot_CAZyme_defline_pd.to_excel(writer, sheet_name='CAZymes_defline_mat_upreg',index=True)

# Fourth, output in alphabetical order of CAZyme_defline
annot_CAZyme_defline_pd = annot_CAZyme_defline_pd.sort_values(by=['CAZyme_defline'], ascending=True)
annot_CAZyme_defline_pd.to_excel(writer, sheet_name='CAZymes_defline_name_ordered',index=True)

# Include Fisher results for CAZyme key classes
Fisher_key_CAZyme_classes_zoosp_upreg.to_excel(writer, sheet_name='Fisher_key_CAZyme_zoosp',index=True)

Fisher_key_CAZyme_classes_mat_upreg.to_excel(writer, sheet_name='Fisher_key_CAZyme_mat',index=True)

# Close the Pandas Excel writer and output the Excel file.
writer.close()

"""
Output All Annot Pd dataframes
"""
filename_out = "All_annotation_dataframes.xlsx"
file_path_out = pjoin(*[temp_folder, filename_out])

writer = pd.ExcelWriter(file_path_out, engine='xlsxwriter')

annot_KOG_defline_pd.to_excel(writer, sheet_name='KOG_defline',index=True)
annot_GO_name_pd.to_excel(writer, sheet_name='GO_name',index=True)
annot_IPR_desc_pd.to_excel(writer, sheet_name='IPR_desc',index=True)
annot_KEGG_desc_pd.to_excel(writer, sheet_name='KEGG_desc',index=True)
annot_KEGG_pathway_pd.to_excel(writer, sheet_name='KEGG_pathway',index=True)
annot_KEGG_pathway_class_pd.to_excel(writer, sheet_name='KEGG_pathway_class',index=True)
annot_CAZyme_defline_pd.to_excel(writer, sheet_name='CAZymes_defline',index=True)

# Close the Pandas Excel writer and output the Excel file.
writer.close()