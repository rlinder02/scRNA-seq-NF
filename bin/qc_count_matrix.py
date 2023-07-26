#!/usr/bin/env python 
# 4-space indented, v1.0.0
# File name: qc_count_matrix.py
# Description: This script can be used to calculate various QC metrics from scRNAseq data.
# Author: Robert Linder
# Date: 2023-07-20

import os
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'
import argparse
import re
from functools import wraps
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric
from sklearn.decomposition import PCA
import leidenalg as leiden

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Input 10x Genomics-style count matrix")
	parser.add_argument("count_matrix_folder", type=str, help="Folder that contains the 10x Genomics count matrix")
	args = parser.parse_args()
	return args

def pre_process(adata):
    """Filters and normalizes the data"""
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)  
    sc.tl.pca(adata, svd_solver='arpack')
    return adata

def clustering(adata):
    """Clusters the data"""
    #computing the neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata)
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata)
    return adata

def dea(adata):
    """Finds differentially expressed genes"""
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='gene_ranks.pdf')
    #get a table with scores and groups
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    score_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']})
    #compare to a single cluster
    sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
    adata.write('data.h5ad', compression='gzip')
    score_df.to_csv('data.csv')
    print(os.getcwd()) 

def main():
    inputs = parse_args()
    anndata = sc.read_10x_mtx(inputs.count_matrix_folder, var_names = 'gene_symbols',
                        cache=True)
    print(f"The dimensions of this dataset are {anndata.shape}")
    preprocessed_data = pre_process(anndata)
    clustered_data = clustering(preprocessed_data)
    diff_data = dea(clustered_data)

if __name__ == "__main__":
	main()