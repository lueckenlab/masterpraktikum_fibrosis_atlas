fit_model <- function(adata_, fit_on){
  if (toupper(fit_on) == toupper("condition")){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group = colData(adata_)$condition)
    condition <- colData(adata_)$condition   # saline, asbestos, untreated, ...
  }
  else if (toupper(fit_on) == toupper("condition-2")){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group = colData(adata_)$condition_2)
    condition <- colData(adata_)$condition_2   # saline changed to untreated
  }
  else if (toupper(fit_on) == toupper("fibrotic/control")){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "X"), group = colData(adata_)$fibrotic.control)
    condition <- colData(adata_)$fibrotic.control   # fibrosis and untreated
  }
  else {
    # end program in this case
    stop(paste0(fit_on, "is not possible!")) 
  }
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  keep <- filterByExpr(y)   # keeps genes that have at least min.count (default: 10) reads in a worthwhile number samples
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  # normalize
  y <- calcNormFactors(y)
  # create the design matrix
  project <- colData(adata_)$dataset       # schiller, misharin, peyser, ...
  design <- model.matrix(~ 0 + condition + project)
  print(colnames(design))   # printing columns from the design matrix
  # estimate dispersion
  y <- estimateDisp(y, design = design)
  # fit the model
  fit <- glmQLFit(y, design)
  return(list("fit"=fit, "design"=design, "y"=y))
}


plot_MDS_and_BVC_plots <- function(cell_type, y, path_out){
  # MDS plot:
  plotMDS(y, labels = NULL, col = as.numeric(y$samples$group),
          title=paste("MDS plot of the public merged data for", cell_type))
  plotMDS(y, labels = NULL, col = as.numeric(y$samples$group), 
          pch = as.numeric(y$samples$cell), cex = 2,
          title=paste("MDS plot of the public merged data for", cell_type))
  # Biological Coefficient of Variation (BCV) plot:
  plotBCV(y, title=paste("BVC plot of the public merged data for", cell_type))
  ggsave(
    paste('merged_public_BVC_', cell_type, '.png', sep = ""),
    plot = last_plot(), device = NULL, path = path_out, scale = 1, width = 30,
    height = 30, units = "cm", dpi = 300, limitsize = TRUE, bg = NULL
  )
}
