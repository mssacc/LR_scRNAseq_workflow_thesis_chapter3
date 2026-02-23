###Single-cell level DE analysis

#Load required libraries
library(Seurat)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(patchwork)

#Load gene RDS
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Rejoin gene datasets after integration
all_samples_integrated <- JoinLayers(all_samples_integrated)

#Define identities by fertility and cell cluster
all_samples_integrated$seurat_clusters.fertility <- paste(all_samples_integrated$seurat_clusters, all_samples_integrated$fertility, sep = "_")
Idents(all_samples_integrated) <- "seurat_clusters.fertility"

#Perform single-cell differential expression
pre_unc_mark <- FindMarkers(all_samples_integrated, ident.1 = "0_Infertile", ident.2 = "0_Fertile", verbose = FALSE)
unc_mark <- FindMarkers(all_samples_integrated, ident.1 = "1_Infertile", ident.2 = "1_Fertile", verbose = FALSE)
cil_mark <- FindMarkers(all_samples_integrated, ident.1 = "2_Infertile", ident.2 = "2_Fertile", verbose = FALSE)
sec_mark <- FindMarkers(all_samples_integrated, ident.1 = "3_Infertile", ident.2 = "3_Fertile", verbose = FALSE)
pre_cil_mark <- FindMarkers(all_samples_integrated, ident.1 = "4_Infertile", ident.2 = "4_Fertile", verbose = FALSE)
prolif_mark <- FindMarkers(all_samples_integrated, ident.1 = "5_Infertile", ident.2 = "5_Fertile", verbose = FALSE)

#Generate volcano plots
pu <- EnhancedVolcano(pre_unc_mark, lab = rownames(pre_unc_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Pre-Unciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#F8766D"),
                      caption = NULL, legendPosition = 'bottom')
u <- EnhancedVolcano(unc_mark, lab = rownames(unc_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Unciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#ABA300"),
                      caption = NULL, legendPosition = 'bottom')
c <- EnhancedVolcano(cil_mark, lab = rownames(cil_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Ciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#0CB702"),
                      caption = NULL, legendPosition = 'bottom')
s <- EnhancedVolcano(sec_mark, lab = rownames(sec_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Pseudobulk Secretory DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#00BFC4"),
                      caption = NULL, legendPosition = 'bottom')
pc <- EnhancedVolcano(pre_cil_mark, lab = rownames(pre_cil_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14, 
                      title = "Pre-Ciliated DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#849AFF"),
                      caption = NULL, legendPosition = 'bottom')
p <- EnhancedVolcano(prolif_mark, lab = rownames(prolif_mark),
                      x = 'avg_log2FC', y = 'p_val_adj',
                      pCutoff = 0.05, FCcutoff = 0.585,
                      pointSize = 2, labSize = 4, axisLabSize = 14,  
                      title = "Proliferative DEGs - Infertile vs Fertile", titleLabSize = 18, subtitle = NULL,
                      col = c("#717171","#717171","#717171","#FF61CC"),
                      caption = NULL, legendPosition = 'bottom')

(wrap_plots(pu, u, c, s, pc, p, ncol = 3)) +
  plot_annotation(theme = theme(plot.title = element_text(size = 18, hjust = 0.5)))
