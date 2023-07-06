#!/bin/env python

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scbi
from sklearn.linear_model import LogisticRegressionCV
from sklearn import metrics


def score_geneset(adata, geneset, colname, method='pca'):
    if method == 'pca':
        testset = adata[:, geneset]
        sc.pp.scale(testset)
        adata.obs[colname] = sc.tl.pca(testset, copy=True, n_comps=1).obsm['X_pca'][:, 0]
        return testset
    if method == 'score_genes':
        testset = adata.copy()
        sc.pp.scale(testset)
        sc.tl.score_genes(testset, geneset, score_name=colname)
        adata.obs[colname] = testset.obs[colname]


def validate_score(adata, cell_type, ct_col, score, scoring='accruracy'):
    adata.obs['is_type'] = False
    adata.obs.loc[adata.obs['ct_col'] == cell_type, 'is_type'] = True
    X = np.reshape(adata.obs[score].array, (-1, 1))
    y = adata.obs.is_type
    clf = LogisticRegressionCV(cv=5, random_state=0, class_weight='balanced', scoring=scoring).fit(X, adata.obs.is_type)
    return clf.score(X, y)


#def plot_violin_overview(adata, cell_types, colname_ct, colname_fibro, colname_macro, colname_epi):
#    plot_adata = adata[adata.obs['colname_ct'].isin(cell_types)]
#    f, axes = plt.subplots(3, 1)
#    sns.violinplot(x='colname_ct', y='colname_fibro', data=plot_adata.obs, scale='width', ax=axes[0], order=cell_types)
#    axes[0].set_xticks([])
#    axes[0].set_xlabel("")
#    sns.violinplot(x='colname_ct', y='colname_macro', data=plot_adata.obs, scale='width', ax=axes[1], order=cell_types)
#    axes[1].set_xticks([])
#    axes[1].set_xlabel("")
#    sns.violinplot(x='colname_ct', y='colname_epi', data=plot_adata.obs, scale='width', ax=axes[2], order=cell_types)
#    axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, horizontalalignment='right')


#import diffxpy.api as de


# def run_diff_unc(adata):
#     sc.pp.filter_genes(adata, min_cells=50)
#     test_tmp = de.test.wald(
#         data=adata.layers['counts'],
#         formula_loc="~ 1 + cond_test",
#         # as_numeric=['size_factors'],
#         coef_to_test=["cond_test[T.control]"],
#         sample_description=adata.obs,
#         # constraints_loc={'identifier':'cond'},
#         noise_model='nb',
#         gene_names=adata.var_names,
#         size_factors='size_factors'
#     )
#     return test_tmp


from matplotlib import colors

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7, 0.8, 1))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)


def plot_markers(adata):
    print('Macrophages')
    sc.pl.umap(adata, color=['Apoe', 'Mrc1', 'Marco', 'Mertk'], cmap=mymap, size=20)
    print('Monocytes')
    sc.pl.umap(adata, color=['Cd14', 'Vcan'], cmap=mymap, size=20)
    print('T-Cells')
    sc.pl.umap(adata, color=['Cd3e', 'Cd4', 'Cd8a', 'Cd3d'], cmap=mymap, size=40)
    print('B-Cells')
    sc.pl.umap(adata, color=['Cd19', 'Cd79a'], cmap=mymap, size=40)
    print('AT2')
    sc.pl.umap(adata, color=['Muc1', 'Sftpc', 'Sftpd', 'Lcn2'], cmap=mymap, size=30)
    print('AT1')
    sc.pl.umap(adata, color=['Vegfa'], cmap=mymap, size=30)
    print('Krt8')
    sc.pl.umap(adata, color='Krt8', cmap=mymap, size=30)
    print('Endothelial')
    sc.pl.umap(adata, color='Pecam1', cmap=mymap, size=30)
    print('Fibroblasts')
    sc.pl.umap(adata, color='Col1a2', cmap=mymap, size=30)
    print('Myofibroblasts')
    sc.pl.umap(adata, color=['Col3a1', 'Cthrc1', 'Postn', 'Spp1', 'Tnc', 'S100a6', 'Ccl2'], cmap=mymap, size=30)