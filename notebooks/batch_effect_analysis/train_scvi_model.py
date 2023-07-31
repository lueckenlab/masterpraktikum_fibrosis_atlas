# import packages
import scanpy as sc
import scvi
from rich import print
from scib_metrics.benchmark import Benchmarker
from scvi.model.utils import mde
from scvi_colab import install
import tqdm as notebook_tqdm
import pandas as pd

# set scvi seed to reproduce results
scvi.settings.seed = 0

# load adata object
adata = sc.read("/home/h/hollenberg/MaPra/merged_data.h5ad")

# setup scvi
scvi.model.SCVI.setup_anndata(adata, layer="raw_counts", batch_key="batch")

# init model 
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

# train model
vae.train()

# save trained model
# the folder must not exist
vae.save("/home/h/hollenberg/MaPra/scvi_model/")

# safe embedding to obs
adata.obsm["X_scVI"] = vae.get_latent_representation()

# calculate Clusters
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata)

adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

# Add the category 'unknown' to the existing categories
adata.obs['manual_celltype_annotation'] = adata.obs['manual_celltype_annotation'].cat.add_categories('unknown')

# Replace NaN values with 'unknown'
adata.obs['manual_celltype_annotation'].fillna('unknown', inplace=True)

# setup scvi
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="manual_celltype_annotation",
    unlabeled_category="unknown",
)

# train model
lvae.train(max_epochs=20, n_samples_per_label=100)

# safe model
# the folder must not exist
lvae.save("/home/h/hollenberg/MaPra/scanvi_model/")

# safe embedding to obs
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

adata.obsm["X_mde_scanvi"] = mde(adata.obsm["X_scANVI"])

# Compute integration metrics
bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key="manual_celltype_annotation",
    embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI"],
    n_jobs=-1,
)
bm.benchmark()


# check if outdir exist, create dir if it does not
outdir = "/home/h/hollenberg/MaPra/benchmark_results/"
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# plot results
bm.plot_results_table(min_max_scale=False, save_dir = outdir)

# store results in a dataframe 
df = bm.get_results(min_max_scale=False)
df.to_csv('/home/h/hollenberg/MaPra/benchmark_results/results.csv', index=True)



