import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from os.path import join as pjoin

input_folder = r'input' 
temp_folder = r'temp'
output_folder = r'output'

"""
Filenames for import
"""
input_dge_filename = "Supplementary Dataset DGE_summary draft 2.xlsx"


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


