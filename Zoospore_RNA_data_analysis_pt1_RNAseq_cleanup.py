"""
Zoospore RNA Manuscript: RNA Data Analysis RNAseq Cleanup

@author: Lazarina Butkovich

Features include:
    Modified proteinIDs FASTA with removal of genes with duplicate amino acid sequence
    Consolidation of counts for dupes
    Consolidation of TPM counts for dupes
"""
import pandas as pd
import os
from os.path import join as pjoin

"""
Functions
"""
def parse_fasta_file(f):
    """
    Parses a fasta file into a dictionary of 'proteinID':'amino acid sequence'

    input: 
    f: a file handle for an input fasta file

    output:
    parsed_seqs: a dictionary of 'proteinID':'amino acid sequence'
    """
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if curr_seq_id is not None:
                parsed_seqs[curr_seq_id] = ''.join(curr_seq)

            curr_seq_id = line[1:]
            curr_seq = []
            continue
        curr_seq.append(line)
    #Add the final sequence to the dict
    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
    return parsed_seqs

def remove_proteinIDs_with_duplicate_aa_seqs(counts, dupes, prefix, index_id):
    """
    Removes proteinIDs with duplicate amino acid sequences from counts dataframe and consolidates counts for proteinIDs with identical amino acid sequences

    inputs: 
    counts: original counts dataframe. ProteinIDs are in the index_id column, written as prefix+str(proteinID)
    dupes: dictionary of proteinIDs with duplicate amino acid sequences
    prefix: a string prefix to the proteinID value in the counts dataframe

    output: 
    counts: a dataframe of counts with proteinIDs with duplicate amino acid sequences removed and counts consolidated for proteinIDs with identical amino acid sequences
    """
    counts = counts.set_index(index_id)
    # Get header names 
    samples = list(counts.columns.values)
    # Iterate over duped_values and make changes to counts_old
    for k,v in dupes.items():
        protID_keep = prefix + str(k)
        for val in v: #note that v is a set, so its values don't have an index
            #skip key value
            if str(val) != str(k):
                protID_remove = prefix + str(val)
                for col in samples:
                    # set col,protID_keep count to sum of keep and remove counts
                    counts[col][protID_keep] = counts.at[protID_keep,col] + counts.at[protID_remove,col]
                #remove protID_remove row
                counts = counts.drop(protID_remove)
    return counts


"""
Values
""" 
input_folder = r'input' 
temp_folder = r'temp'

# Input filenames
deflines_in_filename = 'Neosp1_FilteredModels5_deflines_pre.fasta'
counts_in_filename = 'counts_RNAseq_original.txt'
tpm_counts_in_filename = 'tpm_counts_RNAseq_original.txt'

# Prefix for proteinIDs in counts_in_filename file
fasta_prefix = 'jgi.p|Neosp1|'
counts_id_col_name = 'GeneID'

# Output filenames
deflines_out_filename = 'Neosp1_FilteredModels5_deflines_post.fasta'
counts_out_filename = 'counts_RNAseq_updated.csv'
tpm_counts_out_filename = 'tpm_counts_RNAseq_updated.csv'
duped_values_out_filename = 'Neosp1_FilteredModels5_deflines_duped_proteinIDs_sorted.fasta'


"""
Input files
"""
# Parse protein sequences into a dictionary of 'proteinID':'amino acid sequence'
f = open(pjoin(*[input_folder, deflines_in_filename]))
parsed_seqs = parse_fasta_file(f)
f.close()

# Input RNAseq raw counts data file as a pandas dataframe. These counts values will be consolidated for proteinIDs with identical amino acid sequences
counts_original = pd.read_csv(pjoin(*[input_folder, counts_in_filename]),sep='\t')

# Input RNAseq TPM counts data file as a pandas dataframe. These counts values will be consolidated for proteinIDs with identical amino acid sequences
tpm_counts_original = pd.read_csv(pjoin(*[input_folder,tpm_counts_in_filename]), sep="\t")

"""
Determine proteinIDs with duplicate amino acid sequences
"""
components = []
any_dups = []
dup_dict = {}

# Reverse the dictionary key:value from the parsed_seqs dictionary. This should consolidate the amino acid sequences with multiple proteinID names and the values are the proteinIDs that match the sequence 
rev_parsed_seqs = {}
for key, value in parsed_seqs.items():
    rev_parsed_seqs.setdefault(value,set()).add(key) #https://stackoverflow.com/questions/20672238/find-dictionary-keys-with-duplicate-values
# duped_values dictionary describes (i) which proteinIDs have identical amino acid sequences, (ii)) which proteinIDs were removed and (iii) to what proteinID their data was consolidated to
duped_values = {} 
# new_parsed_seqs dictionary has proteinID key and amino acid sequence values for the updated FASTA file contents
new_parsed_seqs={}
# Convert string proteinIDs to number values for sorting
for key, value in rev_parsed_seqs.items():
    k = sorted(list(map(int, value)))
    new_parsed_seqs[k[0]]=key
    if len(k)>1: # if sequence with duplicated proteinIDs:
        duped_values[k[0]]=k # duped_values values include the key proteinID

duped_values_sorted = {}
for key,value in duped_values.items():
    vals = sorted(value)
    k = vals[0]
    duped_values_sorted[k]=vals


"""
Consolidate counts and TPM counts data for proteinIDs with identical amino acid sequences
"""
counts = remove_proteinIDs_with_duplicate_aa_seqs(counts_original, duped_values_sorted, fasta_prefix, counts_id_col_name)            
tpm_counts = remove_proteinIDs_with_duplicate_aa_seqs(tpm_counts_original, duped_values_sorted, fasta_prefix, counts_id_col_name)            


"""
Output files
"""
# Manual analysis of the replicates using Pearson Correlation Coefficients shows that zoospore replicate HHCCZ is an outlier. For the purposes of this analysis, we will remove this replicate from the counts and TPM counts dataframes.
counts = counts.drop(['HHCCZ'], axis=1)
tpm_counts = tpm_counts.drop(['HHCCZ'], axis=1)

# Create temp folder if it doesn't exist
if not os.path.exists(temp_folder):
    os.makedirs(temp_folder)

# Write FASTA file with consolidated proteinIDs
fastA_file = open(pjoin(*[temp_folder, deflines_out_filename]), 'w') 
for k,v in new_parsed_seqs.items():
    fastA_file.write(f">{k}\n{v}\n")
fastA_file.close()

# Write duped_values_sorted dictionary to file
f2 = open(pjoin(*[temp_folder, duped_values_out_filename]),'w')
for k,v in duped_values_sorted.items():
    f2.write(f"{k}: {v}\n")
f2.close()

# Write counts csv file with consolidated proteinIDs
counts.to_csv(path_or_buf=pjoin(*[temp_folder, counts_out_filename]), sep=',')
# Write TPM csv counts file with consolidated proteinIDs
tpm_counts.to_csv(path_or_buf=pjoin(*[temp_folder, tpm_counts_out_filename]), sep=',')