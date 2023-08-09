import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import pertpy as pt
import mudata as mu
import yaml
import sys
import os

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod
from sccoda.util import data_visualization as viz

config_file = sys.argv[1]

with open(config_file, 'r') as configfile:
    config = yaml.safe_load(configfile)

# path of partly preprocessed HLCA file used in this notebook:
adata_path = config['adata_path']
add_annotation = config['add_annotation']
cell_type_identifier=config['cell_type_identifier']
sample_identifier=config['sample_identifier']
covariate_obs=config['covariate_obs']
group1 = config['group1']
group2 = config['group2']
reference_cell_type=config['reference_cell_type']
fdr = config['fdr']
lfc_thresh = config['lfc_thresh']
dir_out = config['dir_out']

adata = sc.read(adata_path)
# load annotation if not in adata
if add_annotation is not None:
    adata.obs = pd.read_csv(add_annotation)
    

meta = adata.obs[[covariate_obs, sample_identifier]].drop_duplicates().set_index(sample_identifier)
scc_dat = dat.from_scanpy(
    adata,
    cell_type_identifier=cell_type_identifier,
    sample_identifier=sample_identifier,
    covariate_df=meta
)
scc_dat
sccoda_mod = mod.CompositionalAnalysis(
    scc_dat,
    formula=covariate_obs,
    reference_cell_type=reference_cell_type,
)
sccoda_res = sccoda_mod.sample_hmc()
#scc_dat.write(output_adata)

# build sccoda model
sccoda_model = pt.tl.Sccoda()
sccoda_data = sccoda_model.load(adata, type="cell_level", generate_sample_level=True, 
                                cell_type_identifier=cell_type_identifier,
                                sample_identifier=sample_identifier, covariate_obs=[covariate_obs])

# new combi 
group1_group2 = group1 + "_" + group2
sccoda_data.mod[group1_group2] = sccoda_data["coda"][sccoda_data["coda"].obs[covariate_obs].isin([group1, group2])].copy()


sc.settings.set_figure_params(
    dpi=100,
    color_map='plasma',
    dpi_save=200,
    vector_friendly=True,
    frameon=False,
    fontsize=10,
    figsize=(8,6),
    format='png',
)
# set thresholds
sccoda_res.set_fdr(est_fdr=fdr)

# write model summary
with open(os.path.join(dir_out, "sccoda_summary_" + adata_path.split("/")[-1].split(".")[0] + ".txt"), 'w') as f:
    print(sccoda_res.summary_extended(), file=f)
    cred_effects = sccoda_res.credible_effects()
    print(cred_effects, file=f)

    # plot significant results
    with PdfPages(os.path.join(dir_out, "sccoda_" + adata_path.split("/")[-1].split(".")[0] + ".pdf")) as pdf:

        plt.rcParams['figure.figsize'] = (16, 10)
        pt.pl.coda.boxplots(sccoda_data, modality_key=group1_group2, feature_name=covariate_obs, add_dots=True)
        pdf.savefig(bbox_inches='tight')
        plt.close()

        covariates = cred_effects.index.get_level_values('Covariate').unique()
        for covariate in covariates:
            plt.rcParams['figure.figsize'] = (8,6)
            effects_df = sccoda_res.effect_df.loc[covariate]
            effects_df = effects_df[np.abs(effects_df["log2-fold change"]) >= lfc_thresh]
            effects_df = effects_df.loc[cred_effects[covariate]].reset_index()
            effects_df = effects_df.sort_values("log2-fold change", ascending=False)

            print(effects_df, file=f)

            err_min = np.abs(effects_df['HDI 3%'] - effects_df['Final Parameter'])
            err_max = np.abs(effects_df['HDI 97%'] - effects_df['Final Parameter'])
    
            if effects_df.shape[0] == 0:
                print(f'skip {covariate}')
                continue
    
            # Final parameter
            # plt.bar(
            #     data=effects_df,
            #     x="Cell Type",
            #     height="Final Parameter",
            # )
            plt.errorbar(
                data=effects_df,
                x="Cell Type",
                y="Final Parameter",
                yerr=[err_min, err_max],
                fmt='o',
            )
            plt.axhline(y=0, color='black', linestyle='-')
            plt.title(f'Final Parameter {covariate} FDR={fdr}')
            plt.xticks(rotation=90)
            plt.ylabel ('Final Parameter')
            pdf.savefig(bbox_inches='tight')
            plt.close()

            # significance
            sns.barplot(
                data=effects_df,
                x="Cell Type",
                y="log2-fold change",
                order=effects_df['Cell Type'],
            ).set(title=f'{covariate} FDR={fdr}')
            plt.xticks(rotation=90)
            plt.ylabel ('Log2 Fold Change')
            pdf.savefig(bbox_inches='tight')
            plt.close()
    
    
        # add metadata
        d = pdf.infodict()
        d['Title'] = f'scCODA on {adata_path.split("/")[-1]}'