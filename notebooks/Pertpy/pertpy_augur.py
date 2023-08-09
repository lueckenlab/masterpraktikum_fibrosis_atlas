#!/usr/bin/python3

import warnings

warnings.filterwarnings("ignore")
import scanpy as sc
import pertpy as pt

import json
import pickle


adata = sc.read("/lustre/groups/ml01/workspace/daniel.michaela.masterpraktikum23/all_datasets_merged/complete_merged_public_and_galapagos_harmonized_doublet.h5ad")


datasets = [
    'galapagos_bleo',
    'galapagos_rad', 'misharin', 'peyser', 'schiller',
       'tsukui', 'xie']
for x in datasets: 
    print(x)
    subset = adata[adata.obs["dataset"].isin([x])]
    subset.obs["cell_type"] = subset.obs["harmonized_anno"] 
    subset.obs["label"] = subset.obs["fibrotic/control"]

    ag_rfc = pt.tl.Augur("random_forest_classifier")
    loaded_data = ag_rfc.load(subset)

    h_adata, h_results = ag_rfc.predict(loaded_data, subsample_size=20, n_threads=4)
    #save summary metrics 
    h_results["summary_metrics"].to_pickle("/home/icb/leonie.pohl/masterpraktikum_fibrosis_atlas/notebooks/Pertpy/Output_augur/summary_metrics_"+x+".pkl")  


    
    with open("/home/icb/leonie.pohl/masterpraktikum_fibrosis_atlas/notebooks/Pertpy/Output_augur/all_metrics"+x+".json", 'wb') as fp:
        pickle.dump(h_results, fp)
    
    del h_adata.uns["augurpy_results"]
    # save hdata 
    h_adata.write("/home/icb/leonie.pohl/data/augur_"+x+".h5ad")

    # save lolipop plots
    #lollipop = pt.pl.ag.lollipop(h_results)

    #important_features = pt.pl.ag.important_features(h_results)
    
