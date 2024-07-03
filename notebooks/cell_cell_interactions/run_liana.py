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


directory = sys.argv[1]
filenames = next(walk(sys.argv[1]), (None, None, []))[2]
print(directory, filenames)

#filenames = filenames[-1:] + filenames[:2]
print(filenames)

for file in filenames:
    adata = sc.read(str(directory)+"/"+file)

    condition = [x for x in ['bleomycin', 'bleo', 'Bleo', 'asbestos'] if x in list(adata.obs["condition"].cat.categories)][0]
    control = [x for x in ['saline', 'healthy', 'UT', 'control'] if x in list(adata.obs["condition"].cat.categories)][0]

    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    # log1p normalize the data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # make sure the format is correct
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
    break

