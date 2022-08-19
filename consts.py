import os
import pandas as pd
import json


muted_columns = ['Tissue', 'fraq', 'Group', 'Dataset.Type', 'Freq']

base_dir = os.path.dirname(os.path.abspath(__file__))
genes_data = json.load(open('/huge/bulk/ENSEMBLE/genes.json', 'r'))
motifs_data = pd.read_csv(f'{base_dir}/../data/new_splicing_factors_symbols.tsv', sep='\t', index_col=0)
tcga_sfs = pd.read_csv('/huge/bulk/TCGA/TCGA-COMBINED/combined/sfs_FPKM.tsv', sep='\t', index_col=0).T
