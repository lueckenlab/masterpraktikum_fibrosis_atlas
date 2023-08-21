import sys
from os import walk
from rpy2 import robjects

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import liana as li
from liana.method import rank_aggregate


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

    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    # log1p normalize the data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # prepare pseudobulk for later
    # Get pseudo-bulk profile
    pdata = dc.get_pseudobulk(
        adata,
        sample_col="batch",
        groups_col="manual_celltype_annotation",
        min_prop=0.1,
        min_smpls=3,
        layer="raw_counts",
    )
    

    for cond in list(adata.obs["condition"].cat.categories):
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
        sources = list(top_res.source)
        targets = list(top_res.target)

        print("preparing to run NicheNet on interacting celltypes")
        for i in range (20):
            sender_celltypes = sources[i]
            receiver_celltypes = targets[i]

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
            robjects.r('''
                -i sender_expressed -i receiver_expressed
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
        break
    #decide sender and reciever per cond or overall ? where deg genes?
    #break
