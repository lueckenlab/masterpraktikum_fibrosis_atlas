import warnings


warnings.filterwarnings("ignore")


import matplotlib.pyplot as plt

import seaborn as sns

import scanpy as sc

import pandas as pd

import numpy as np

import random

import sc_toolbox

import anndata


sc.settings.verbosity = 0


def aggregate_and_filter(

    adata,

    cell_identity,

    donor_key="sample",

    condition_key="condition_2",

    cell_identity_key="final_annotation",

    obs_to_keep=[],  # which additional metadata to keep, e.g. gender, age, etc.

    replicates_per_patient=1,

):

    # subset adata to the given cell identity

    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()

    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])


    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")

    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):

        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")            

        adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]

        # create replicates for each donor

        indices = list(adata_donor.obs_names)

        random.shuffle(indices)

        indices = np.array_split(np.array(indices), replicates_per_patient)

        for i, rep_idx in enumerate(indices):

            adata_replicate = adata_donor[rep_idx]

            # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information

            agg_dict = {gene: "sum" for gene in adata_replicate.var_names}

            for obs in obs_to_keep:

                agg_dict[obs] = "first"

            # create a df with all genes, donor and condition info

            df_donor = pd.DataFrame(adata_replicate.X.A)

            df_donor.index = adata_replicate.obs_names

            df_donor.columns = adata_replicate.var_names

            df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])

            # aggregate

            df_donor = df_donor.groupby(donor_key).agg(agg_dict)

            df_donor[donor_key] = donor

            df.loc[f"donor_{donor}_{i}"] = df_donor.loc[donor]

    print("\n")

    # create AnnData object from the df

    adata_cell_pop = sc.AnnData(

        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)

    )

    return adata_cell_pop


adata = sc.read("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/all_datasets_merged/complete_merged_public_and_galapagos_integrated_scvi_neighbors_umap_withfinalanno_zwischenstand_small.h5ad")


#small_adata = sc.pp.subsample(adata,fraction=0.01, random_state=42, copy=True)


#small_adata.write("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/all_datasets_merged/complete_merged_public_and_galapagos_integrated_scvi_neighbors_umap_withfinalanno_zwischenstand_small.h5ad")


#print("new adata was made")

#adata = small_adata.copy()
#print("adata was copied")
print(adata)


obs_to_keep = ["dataset", "batch", "sample", "condition", "condition_2", "fibrotic/control",

               "coarse_harmonized_anno", "harmonized_anno", "final_annotation"]   # obs columns that we want to keep

print("adata was read")


# process first cell type separately...

cell_type = adata.obs["final_annotation"].cat.categories[0]

print(

    f'Processing {cell_type} (1 out of {len(adata.obs["final_annotation"].cat.categories)})...'

)

adata_pb = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)

for i, cell_type in enumerate(adata.obs["final_annotation"].cat.categories[1:]):

    print(

        f'Processing {cell_type} ({i+2} out of {len(adata.obs["final_annotation"].cat.categories)})...'

    )

    adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)

    adata_pb = anndata.concat([adata_pb, adata_cell_type])   # --> adata_pb.concatenate(adata_cell_type) does not work!


print("processing ended, now doing normalisation")

adata_pb.layers['counts'] = adata_pb.X.copy()


sc.pp.normalize_total(adata_pb, target_sum=1e6)

sc.pp.log1p(adata_pb)

sc.pp.pca(adata_pb)


print("normalization ended, now doing something with lib size")


adata_pb.obs["lib_size"] = np.sum(adata_pb.layers["counts"], axis=1)

adata_pb.obs["lib_size"] = adata_pb.obs["lib_size"].astype(float)   # --> not in example but important

adata_pb.obs["log_lib_size"] = np.log(adata_pb.obs["lib_size"])


print ("plotting pca and saving it (hopefully)")


#pcafig = sc.pl.pca(adata_pb, color=["dataset", "batch", "sample", "condition", "condition_2", "fibrotic/control",

                           "coarse_harmonized_anno"], show=False)

#pcafig.savefig("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/masterpraktikum_fibrosis_atlas/notebooks/differential_gene_expression_analysis/outputfigs/pcafigures.png")


adata_pb.X = adata_pb.layers['counts'].copy()

print(np.max(adata_pb.X))   # check if .X contains raw counts


print("plot sum of raw counts")


column_sums = np.sum(adata_pb.X, axis=1)

plt.hist(column_sums, bins=90)   # plot the sums of the raw counts as a sanity check

plt.gca().set(title='Pseudobulk - raw counts histogram', ylabel='Frequency')

plt.savefig("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/masterpraktikum_fibrosis_atlas/notebooks/differential_gene_expression_analysis/outputfigs/raw_counts_histogram_pseudobulk.png", dpi = 600)


print("adata X conversion and saving adata")

# very important conversions as adata_pb.write() does not work otherwise 

adata_pb.X = np.vstack(adata_pb.X[:, :]).astype(float)

adata_pb.layers['counts'] = np.vstack(adata_pb.layers['counts'][:, :]).astype(float)


adata_pb.write("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/all_datasets_merged/pseudobulk_merged_data_for_diffEx_edgeR_condition-2_withoutSchiller_finalanno.h5ad")

print("everything worked :)")
