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
import scvi
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric
from sklearn.decomposition import PCA
from scipy.sparse import csr_matrix
import leidenalg as leiden

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Input 10x Genomics-style count matrix")
	parser.add_argument("count_matrix_folder", type=str, help="Folder that contains the 10x Genomics count matrix")
	args = parser.parse_args()
	return args

def pre_process(adata):
    """Filters count matrices from individual samples, starting with doublet removal"""
    adata_raw = adata
    sc.pp.filter_genes(adata_raw, min_cells=3)
    sc.pp.highly_variable_genes(adata_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
    scvi.model.SCVI.setup_anndata(adata_raw)
    vae = scvi.model.SCVI(adata_raw)
    vae.train()
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    df.index = df.index.map(lambda x: x[:-2])
    df['dif'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    sc.pp.filter_cells(adata, min_genes=200)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    return adata

def combine_samples(sample_ids):
     """Concatenates count matrices from multiple samples into a single count matrix, then filter, normalize, transform, and continue to process the data"""
     sc.pp.normalize_total(adata, target_sum=1e4)
     sc.pp.log1p(adata)
     sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
     adata = adata[:, adata.var.highly_variable]
     sc.pp.scale(adata, max_value=10)  
     sc.tl.pca(adata, svd_solver='arpack')
     sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50, save='PC_ranking.pdf')

def clustering(adata):
    """Clusters the data"""
    ## computing the neighborhood graph using 20 PCs as default
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.pl.umap(adata, save='UMAP.pdf')
    sc.tl.leiden(adata) 
    sc.pl.umap(adata, color=['leiden'], save='UMAP_Leiden_overlay.pdf') 
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata)
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata)
    return adata

def dea(adata):
    """Finds differentially expressed genes"""
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='gene_ranks.pdf')
    ## get a table with scores and groups
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    score_df = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']})
    ## compare to a single cluster
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