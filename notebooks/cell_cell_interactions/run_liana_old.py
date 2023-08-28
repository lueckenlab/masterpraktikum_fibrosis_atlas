import sys
from os import walk

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import liana as li
from liana.method import rank_aggregate
import decoupler as dc

import anndata2ri
import rpy2
from rpy2 import robjects
from rpy2.robjects import r
import random

anndata2ri.activate()

robjects.r('''
    suppressPackageStartupMessages({
        library(reticulate)
        library(ggplot2)
        library(tidyr)
        library(dplyr)
        library(purrr)
        library(tibble)
    })
    library("nichenetr", lib="/home/d/danilina/mambaforge/envs/scanpy_r/lib/R/library")
''')


# Helper function to obtain sufficiently expressed genes
from functools import reduce


def get_expressed_genes(adata, cell_type, expr_prop):
    # calculate proportions
    temp = adata[adata.obs["manual_celltype_annotation"] == cell_type, :]
    a = temp.X.getnnz(axis=0) / temp.X.shape[0]
    stats = (
        pd.DataFrame({"genes": temp.var_names, "props": a})
        .assign(cell_type=cell_type)
        .sort_values("genes")
    )

    # obtain expressed genes
    stats = stats[stats["props"] >= expr_prop]
    expressed_genes = stats["genes"].values

    return expressed_genes


robjects.r('''
    # Increase timeout threshold
    options(timeout=600)

    # Load PK
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
''')


directory = sys.argv[1]
filenames = next(walk(sys.argv[1]), (None, None, []))[2]
#print(directory, filenames)

for file in filenames:
    adata = sc.read(str(directory)+"/"+file)

    condition = [x for x in ['bleomycin', 'bleo', 'Bleo', 'asbestos'] if x in list(adata.obs["condition"].cat.categories)][0]
    control = [x for x in ['saline', 'healthy', 'UT', 'control'] if x in list(adata.obs["condition"].cat.categories)][0]

    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    # log1p normalize the data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # make sure dc won't throw errors
    adata.obs["batch"] = adata.obs["batch"].astype("str")
    adata.obs["manual_celltype_annotation"] = adata.obs["manual_celltype_annotation"].astype("str")

    # prepare pseudobulk for later
    # Get pseudo-bulk profile
    pdata = dc.get_pseudobulk(
        adata,
        sample_col="batch",
        groups_col="manual_celltype_annotation",
        min_prop=0.1,
        min_smpls=3,
        layer="counts",
    )
    # Storing the raw counts
    pdata.layers["counts"] = pdata.X.copy()
    # Does PC1 captures a meaningful biological or technical fact?
    pdata.obs["lib_size"] = pdata.X.sum(1)
    # Normalize
    sc.pp.normalize_total(pdata, target_sum=None)
    sc.pp.log1p(pdata)

    logFCs, pvals = dc.get_contrast(
        pdata,
        group_col="manual_celltype_annotation",
        condition_col="condition",
        condition=condition,
        reference=control,
        method="t-test",
    )
    # format results
    deg = dc.format_contrast_results(logFCs, pvals)
    
    # change back for liana later
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    adata.obs["manual_celltype_annotation"] = adata.obs["manual_celltype_annotation"].astype("category")

    for cond in [condition, control]:
        adata_cond = adata[adata.obs["condition"]==cond].copy()

        # run liana consensus
        print("Running rank_aggregate on "+file[:-5]+", "+cond)
        rank_aggregate(
            adata_cond, groupby="manual_celltype_annotation", resource_name = 'mouseconsensus',
            return_all_lrs=True, use_raw=False, verbose=True)
        liana_res = adata_cond.uns["liana_res"].drop_duplicates(["ligand_complex", "receptor_complex"]).sort_values(["magnitude_rank", "specificity_rank"],)
        liana_res.to_csv("./results/"+file[:-5]+"_"+cond+"_liana.csv")
        adata_cond.write("../../../data/liana_anndatas/"+file[:-5]+"_"+cond+"_liana.h5ad", compression='gzip')
        fig = li.pl.dotplot(
            adata=adata_cond,
            colour="magnitude_rank",
            size="specificity_rank",
            inverse_colour=True,  # we inverse sign since we want small p-values to have large sizes
            inverse_size=True,
            # since the rank_aggregate can also be interpreted as a probability distribution
            # we can again filter them according to their specificity significance
            # yet here the interactions are filtered according to
            # how consistently highly-ranked is their specificity across the methods
            filterby="specificity_rank",
            filter_lambda=lambda x: x <= 0.05,
            # again, we can also further order according to magnitude
            orderby="magnitude_rank",
            orderby_ascending=True,  # prioritize those with lowest values
            top_n=20,  # and we want to keep only the top 20 interactions
            figure_size=(46, 26),
            #size_range=(1, 6),
            return_fig=True
        )   
        fig.save("./results/"+file[:-5]+"_"+cond+"_liana.png", dpi=500, limitsize=False)


        # take top 20 interacting pairs of celltypes from liana results
        top_res = liana_res.drop_duplicates(["source", "target"])[:20]
        sender_celltypes = list(top_res.drop_duplicates("source").source)
        receiver_celltypes = list(top_res.drop_duplicates("target").target)

        print("preparing to run NicheNet on interacting celltypes")

        # create baseline and DEG with decoupler
        sender_expressed = reduce(
                np.union1d,
                [
                    get_expressed_genes(adata_cond, cell_type=cell_type, expr_prop=0.1)
                    for cell_type in sender_celltypes
                ],
            )
        receiver_expressed = reduce(
                np.union1d,
                [
                    get_expressed_genes(adata_cond, cell_type=cell_type, expr_prop=0.1)
                    for cell_type in receiver_celltypes
                ],
            )

        print(f"running NicheNet on {sender_celltypes} - {receiver_celltypes} pair")
        robjects.globalenv['sender_expressed'] = sender_expressed
        robjects.globalenv['receiver_expressed'] = receiver_expressed
        robjects.r('''
            # get ligands and receptors in the resource
            ligands <- lr_network %>% pull(from) %>% unique()
            receptors <- lr_network %>% pull(to) %>% unique()

            # only keep the intersect between the resource and the data
            expressed_ligands <- intersect(ligands, sender_expressed)
            expressed_receptors <- intersect(receptors, receiver_expressed)

            # filter the network to only include ligands for which both the ligand and receptor are expressed
            potential_ligands <- lr_network %>% 
            filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
            pull(from) %>% unique()
        ''')

        # only keep the receiver cell type(s)
        deg = deg[np.isin(deg["contrast"], receiver_celltypes)]
        # define background of sufficiently expressed genes
        background_genes = deg["name"].values

        # only keep significant and positive DE genes
        deg = deg[(deg["pvals"] <= 0.05) & (deg["logFCs"] > 1)]
        # get geneset of interest
        geneset_oi = deg["name"].values

        print("getting ligand activities")
        robjects.globalenv['geneset_oi'] = geneset_oi
        robjects.globalenv['background_genes'] = background_genes
        robjects.r('''
            ligand_activities <- predict_ligand_activities(geneset = geneset_oi, 
                                                        background_expressed_genes = background_genes,
                                                        ligand_target_matrix = ligand_target_matrix,
                                                        potential_ligands = potential_ligands)

            ligand_activities <- ligand_activities %>% 
            arrange(-aupr) %>% 
            mutate(rank = rank(desc(aupr)))
                    ''')
        ligand_activities = robjects.r["ligand_activities"]

        robjects.r('''
            top_ligands <- ligand_activities %>%
            top_n(15, aupr) %>% 
            arrange(-aupr) %>%
            pull(test_ligand) %>%
            unique()

            # get regulatory potentials
            ligand_target_potential <- map(top_ligands,
                                        ~get_weighted_ligand_target_links(.x,
                                                                            geneset = geneset_oi,
                                                                            ligand_target_matrix = ligand_target_matrix,
                                                                            n = 500)
                                        ) %>%
                bind_rows() %>% 
                drop_na()
                
            # prep for visualization
            active_ligand_target_links <- 
            prepare_ligand_target_visualization(ligand_target_df = ligand_target_potential, 
                                                ligand_target_matrix = ligand_target_matrix)

            # order ligands & targets
            order_ligands <- intersect(top_ligands,
                                    colnames(active_ligand_target_links)) %>% rev() %>% make.names()
            order_targets <- ligand_target_potential$target %>%
            unique() %>% 
            intersect(rownames(active_ligand_target_links)) %>%
            make.names()
            rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>%
            make.names() # make.names() for heatmap visualization of genes like H2-T23
            colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>%
            make.names() # make.names() for heatmap visualization of genes like H2-T23

            vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>%
            t()
                
            # convert to dataframe, and then it's returned to py
            vis_ligand_target <- vis_ligand_target %>%
                as.data.frame() %>%
                rownames_to_column("ligand") %>%
                as_tibble()
        ''')
        vis_ligand_target = robjects.r["vis_ligand_target"]

        # convert dot to underscore and set ligand as index
        vis_ligand_target["ligand"] = vis_ligand_target["ligand"].replace("\.", "_", regex=True)
        vis_ligand_target.set_index("ligand", inplace=True)
        # keep only columns where at least one gene has a regulatory potential >= 0.05
        vis_ligand_target = vis_ligand_target.loc[
            :, vis_ligand_target[vis_ligand_target >= 0.05].any()
        ]

        fig, ax = plt.subplots(1, 1, figsize=(15, 5))
        sns.heatmap(vis_ligand_target, xticklabels=True, ax=ax)
        fig.save("./results/"+file[:-5]+"_"+cond+"_"+sender_celltypes+"_"+receiver_celltypes+"_heatmap.png", dpi=500, limitsize=False)

        ligand_oi = ligand_activities.head(5)["test_ligand"].values

        fig = li.pl.dotplot(
            adata=adata_cond,
            colour="lr_means",
            size="cellphone_pvals",
            inverse_size=True,  # we inverse sign since we want small p-values to have large sizes
            # We choose only the cell types which we wish to plot
            source_labels=sender_celltypes,
            target_labels=receiver_celltypes,
            # keep only those ligands
            filterby="ligand_complex",
            filter_lambda=lambda x: np.isin(x, ligand_oi),
            # as this type of methods tends to result in large numbers
            # of predictions, we can also further order according to
            # expression magnitude
            orderby="magnitude_rank",
            orderby_ascending=False,  # we want to prioritize those with highest expression
            top_n=25,  # and we want to keep only the top 25 interactions
            figure_size=(23, 13),
            #size_range=(1, 6),
            return_fig=True
        )
        fig.save("./results/"+file[:-5]+"_"+cond+"_"+sender_celltypes+"_"+receiver_celltypes+"_ligands.png", dpi=500, limitsize=False)




        break
    #decide sender and reciever per cond or overall ? where deg genes?
    #break
