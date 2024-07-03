library(VennDiagram)
library(UpSetR)
library(grid)

augur_genes <- list("Hbb-bs", "Scgb1a1", "Crip1", "Rag1", "S100a4", "Napsa",
                    "H2-Q7", "Tmsb10", "Ctss", "Ccl6", "Maf", "Il7r", "Ccnd2",
                    "Ccr2", "Rgcc", "Thy1", "Wfdc17", "Fyb", "Ccl9", "Lgals3")


# Analyzed cell types:
# Alveolar_macrophages, Endothelial, Aerocytes, Transitioning_epithelial

cell_type <- "Endothelial"

df_1 <- read.csv('/Users/emmaschonner/Desktop/tmp/DiffEx_merged_public/Transitioning_epithelial_conditionasbestos-conditionuntreated.csv')
df_2 <- read.csv('/Users/emmaschonner/Desktop/tmp/DiffEx_merged_public/Transitioning_epithelial_conditionbleomycin-conditionuntreated.csv')
df_3 <- read.csv('/Users/emmaschonner/Desktop/tmp/DiffEx_merged_public/Transitioning_epithelial_conditionradiation-conditionuntreated.csv')

df_1_signif <- df_1[df_1$logFC >= 2 | df_1$logFC <= -2, ]   # filter genes based on FC
df_1_signif <- arrange(df_1_signif, desc(abs(logFC)))
genes_1_signif <- df_1_signif$X

df_2_signif <- df_2[df_2$logFC >= 2 | df_2$logFC <= -2, ]   # filter genes based on FC
df_2_signif <- arrange(df_2_signif, desc(abs(logFC)))
genes_2_signif <- df_2_signif$X

df_3_signif <- df_3[df_3$logFC >= 2 | df_3$logFC <= -2, ]   # filter genes based on FC
df_3_signif <- arrange(df_3_signif, desc(abs(logFC)))
genes_3_signif <- df_3_signif$X

print(intersect(augur_genes, genes_1_signif))
print(intersect(augur_genes, genes_2_signif))
print(intersect(augur_genes, genes_3_signif))

# Venn
#venn.diagram(
#  x = list(genes_1_signif, genes_2_signif, genes_3_signif),
#  category.names = c("asbestos" , "bleomycin", "radiation"),
#  filename = '/Users/emmaschonner/Desktop/tmp/14_venn_diagramm.png',
#  output=TRUE
#)
 
# Upside --> better visualization than Venn
listInput <- list(augur = augur_genes, asbestos = c(df_1_signif$X), bleomycin = c(df_2_signif$X),
                  radiation = c(df_3_signif$X))
upset(fromList(listInput), order.by = "freq", )
grid.text(paste("Upset plot for ", cell_type, sep = ""),
          x = 0.65, y=0.95, gp=gpar(fontsize=15))