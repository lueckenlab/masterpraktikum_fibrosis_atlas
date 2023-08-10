library(anndata)
library(edgeR)
library(stringr)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(cowplot)

source("DiffEx_functions.R")   # import helper functions

# 1. Define paths:
# ----------------
# TODO: get the paths via arguments?
# path to pseudobulk object
path_adata <- "/Users/emmaschonner/Desktop/tmp/merged_data_for_diffEx_edgeR_condition.h5ad"
# path to where the plots should be saved
path_out <- "/Users/emmaschonner/Desktop/tmp/DiffEx_merged_public"

# 2. Prepare the data:
# --------------------
# read in the pseudobulk object
adata_pb <- read_h5ad(path_adata)
print(paste0("Dimensions of the annData object: ", dim(adata_pb)[1], " ", dim(adata_pb)[2]))

# convert adata_pb to a SCE (fit_model() does not work otherwise)
sce_pb <- AnnData2SCE(adata_pb, X_name = NULL, layers = TRUE, uns = TRUE, 
                      var = TRUE, obs = TRUE, varm = TRUE, obsm = TRUE,
                      varp = TRUE, obsp = TRUE, raw = FALSE, skip_assays = FALSE,
                      hdf5_backed = TRUE, verbose = NULL)

# 3. Perform DiffEx for each cell type:
# -------------------------------------
# iterate through the cell types of the coarse harmonized annotation and
# calculate differentially expressed genes for each
all_cell_types <- unique(adata_pb$obs$coarse_harmonized_anno)   # all types
augur_cell_types <- list("Fibroblasts", "Interstitial_macrophages", "Epithelial")
milo_cell_types <- list()

# choose the cell type list that you want to use
for (ct in augur_cell_types){
  ct_print <- str_replace(ct, "_", " ")   # will be used for plot titles
  print(ct_print)
  
  # 3.1. subset sce_pb to one cell type
  adata_ct <- sce_pb[, colData(sce_pb)$coarse_harmonized_anno == ct]
  
  # 3.2. fit the model
  # condition-2 ... condition where 'saline' is renamed to 'untreated'
  outs <- fit_model(adata_ = adata_ct, fit_on = "condition-2")
  fit <- outs$fit
  y <- outs$y
  
  # 3.3. plot overviews (MDS and BVC plot)
  # TODO: does not save plots --> needed?
  plot_MDS_and_BVC_plots(cell_type=ct_print, y=y, path_out=path_out)
  
  # 3.4. create gene table for each condition combination
  # TODO: use only asbestos/untreated, bleomycin/untreated, ...?
  #       -> at the moment we use every possible condition combination
  contrasts_tmp <- list()
  for (colname in colnames(y$design)){
    if (startsWith(colname, "condition")){
      contrasts_tmp <- append(contrasts_tmp, colname)
    }
  }
  # create tuple out of contrasts_tmp
  contrast_tuples <- t(combn(contrasts_tmp,2))
  for (i in 1:(length(contrast_tuples)/2)){
    # make a String out of the tuple
    cont_tmp <- paste(contrast_tuples[i,1], contrast_tuples[i,2], sep = "-")
    myContrast <- makeContrasts(toString(cont_tmp), levels = y$design)
    qlf <- glmQLFTest(fit, contrast = myContrast)
    # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
    tt_tmp <- topTags(qlf, n = Inf)
    tt_tmp <- tt_tmp$table
    tt_tmp$PValue_adj <- p.adjust(tt_tmp$PValue, "BH")   # calculate padj
    
    # 3.5. plot p-value distributions
    png(paste(path_out, '/merged_public_p-vals_', ct,'_', cont_tmp, '.png', sep = ""),
        width = 9, height = 5, units = "in", res = 1200)
    par(mfrow=c(1,2))
    hist(tt_tmp$PValue, xlab = 'p-value', main = "")
    hist(tt_tmp$PValue_adj, xlab = 'adjusted p-value', main = "")
    mtext(paste("P-value distribution for", ct_print, "(",cont_tmp, ")"),
          side = 3, line = - 2, outer = TRUE)
    dev.off()
    
    # 3.6. Vulcano plot for each condition combination per cell type
    print(EnhancedVolcano(tt_tmp, lab = rownames(tt_tmp), x = 'logFC', y = 'PValue',
                          title = paste("Vulcano plot for", ct_print, "(",cont_tmp, ")")))
    ggsave(
      paste('merged_public_Vulcano_', ct,'_', cont_tmp, '.png', sep = ""),
      plot = last_plot(), device = NULL, path = path_out, scale = 1, width = 30,
      height = 30, units = "cm", dpi = 300, limitsize = TRUE, bg = NULL
    )
  }
}
