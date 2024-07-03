library(anndata)
library(edgeR)
library(stringr)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(cowplot)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
library(ggupset)
library(gridExtra)

# TODO: update path !
source("/Users/emmaschonner/RStudio - tmp/DiffEx_functions.R")   # import helper functions

# 1. Define paths:
# ----------------
# path to pseudobulk object
path_adata <- "/Users/emmaschonner/Desktop/tmp/pseudobulk_merged_data_for_diffEx_edgeR_condition-2_finalanno.h5ad"
# path to where the plots and data frames should be saved
path_out <- "/Users/emmaschonner/Desktop/tmp/_new"

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
# iterate through the cell types of the final harmonized annotation and
# calculate differentially expressed genes for each
all_cell_types <- unique(adata_pb$obs$final_annotation)   # all types
# TODO: update cell types if needed
leonie_cell_types <- list("Alveolar_macrophages", "Endothelial",
                          "Aerocytes", "Transitioning_epithelial")

# choose the cell type list that you want to use
for (ct in leonie_cell_types){
  ct_print <- str_replace(ct, "_", " ")   # will be used for plot titles
  print("")
  print("")
  print(ct_print)
  print("")
  print("")
  
  # 3.1. subset sce_pb to one cell type
  adata_ct <- sce_pb[, colData(sce_pb)$final_annotation == ct]
  
  # 3.2. fit the model
  # condition-2 ... condition where 'saline' is renamed to 'untreated'
  outs <- fit_model(adata_ = adata_ct, fit_on = "condition-2")
  fit <- outs$fit
  y <- outs$y
  
  # 3.3. plot overviews (MDS and BVC plot)
  # TODO: does not save plots --> needed?
  #get_MDS_and_BVC_plots(cell_type=ct_print, y=y, path_out=path_out)
  
  # 3.4. create gene table for each condition combination
  #      (only asbestos/untreated, bleomycin/untreated, ...)
  condition_list <- grep("condition", colnames(y$design), value = TRUE)
  condition_list <- grep("conditionuntreated", condition_list, value = TRUE, invert = TRUE)
  for (cond in condition_list){
    cont_tmp <- paste(cond, "conditionuntreated", sep = "-")
    cont_tmp_print <- gsub("condition", "condition ", cont_tmp)
    myContrast <- makeContrasts(toString(cont_tmp), levels = y$design)
    qlf <- glmQLFTest(fit, contrast = myContrast)
    # get all of the DE genes
    tt_tmp <- topTags(qlf, n = Inf)
    tt_tmp <- tt_tmp$table
    
    write.csv(tt_tmp, paste0(path_out, "/", ct,'_', cont_tmp,".csv"), row.names=TRUE)
    
    # 3.5. plot p-value distribution
    hist(tt_tmp$PValue, xlab = 'p-value', main = "")
    
    # 3.6. Vulcano plot for each condition combination per cell type
    print(EnhancedVolcano(tt_tmp, lab = rownames(tt_tmp), x = 'logFC', y = 'PValue',
                          title = paste("Vulcano plot for ", ct_print, " (", cont_tmp_print, ")", sep = "")))
    ggsave(
      paste('Vulcano_', ct, '_', cont_tmp, '.png', sep = ""),
      plot = last_plot(), device = NULL, path = path_out, scale = 1, width = 30,
      height = 30, units = "cm", dpi = 300, limitsize = TRUE, bg = NULL
    )
    
    # 3.7. GO
    get_GO_analysis_plots(ct, ct_print, cont_tmp, tt_tmp)
  }
}
