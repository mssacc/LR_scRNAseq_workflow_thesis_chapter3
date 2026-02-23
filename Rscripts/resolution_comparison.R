###Assess cell clustering and assign cell subtypes

#Load required libraries
library(ggraph)
library(clustree)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)

#Load seurat object RDS
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")


##Clustree analysis to help visualize cell splitting at different resolutions
    #Gene Clustree
    all_samples_clustree <- FindClusters(all_samples_integrated, reduction = "umap", resolution = seq(0.05, 0.6, 0.05), dims = 1:10)
    clustree(all_samples_clustree, prefix="RNA_snn_res.")



##Assess endo epi marker gene expression globally
    #Marker gene sets
    marker_lists <- list(
      Proliferative = c("MKI67","CDK1","CKS1B","H2AFZ","PCNA","HMGB2","PCLAF",
                        "PTTG1","STMN1","TOP2A","CENPF","TYMS","AURKB"),
      Ciliated = c("FOXJ1","PIFO","TPPP3","RSPH1","DNAAF1","C1orf194","C9orf24",
                   "C20orf85","CAPSL","C5orf49","ZMYND10","SNTN"),
      Pre_Ciliated = c("CCNO","CDC20B","MUC12"),
      Glandular = c("FOXA2","MUC1","CXCL14","DPP4","ABCG1","PAEP","SCGB2A2",
                    "IER3","MUC16","GPX3","SLPI"),
      Luminal = c("PTGS1","CLDN22","IL6","LEFTY1","LGR5","VTCN1","CRISP3","PAX2","SULT1E1"),
      Pre_Unciliated = c("SULT1E1","PTGS1","ABCG1"),
      Secretory = c("PAEP","SPP1","CXCL14","DPP4","C2CD4A","GPX3","MUC1",
                    "GADD45A","IGFBP2","PAX8","S100P"),
      Unciliated = c("AKAP12","ADGRF1","CCL20","CCND2","CD9","CLEC4E","COL4A1","CRYAB",
                     "CXCL2","CYP3A5","DUSP5","GCNT3","GDF15","IDO1","KLF5",
                     "KRT7","KRT8","KRT18","KRT19","KRT23","LAMB3","MAP1B","MUC16",
                     "PERP","PLA2G16","PMEPA1","S100A6","SDC4","SLCO4A1",
                     "SLPI","TACSTD2","TIMP2","TINAGL1","TSPAN1","UCA1"),
      Stem = c("FUT4","SOX9","AXIN2","CDH2"))
  
    #Ensure genes exist
    marker_lists <- lapply(marker_lists, function(genes) {
      genes[genes %in% rownames(all_samples_integrated)]})
    
    #Calculate module scores
    for (name in names(marker_lists)) {
      all_samples_integrated <- AddModuleScore(
        all_samples_integrated,
        features = list(marker_lists[[name]]),
        name = name)}
    
    #Extract UMAP + scores
    umap_data <- as.data.frame(Embeddings(all_samples_integrated, reduction = "umap"))
    
    for (name in names(marker_lists)) {
      umap_data[[name]] <- all_samples_integrated@meta.data[[paste0(name, "1")]]}
    
    #Create plots (independent scales)
    plot_list <- list()
    
    for (name in names(marker_lists)) {
      
      p <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = .data[[name]])) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_viridis(option = "plasma") +
        labs(color = "Score",
             x = "UMAP 1", y = "UMAP 2",
             title = paste(name, "Marker Gene Expression")) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
      
      plot_list[[name]] <- p}
    
    #Combine plots
    wrap_plots(plot_list)



###Compare grouping across resolutions
    #3 groupings
    DimPlot(object = all_samples_integrated, reduction = "umap",
            cols = c("#F8766D","#F8766D","#0CB702","#F8766D","#0CB702","#FF61CC"))
    
    #4 groupings
    DimPlot(object = all_samples_integrated, reduction = "umap",
            cols = c("#F8766D","#F8766D","#0CB702","#F8766D","#849AFF","#FF61CC"))



##Compare cell clustering across resolutions
    DefaultAssay(all_samples_integrated) <- "RNA"
    DimPlot(all_samples_integrated, reduction="umap", group.by = "seurat_clusters")
    
    #Change resolution
    all_samples_integrated <- FindClusters(all_samples_integrated, resolution = 0.15)
    p1 <- DimPlot(all_samples_integrated, reduction="umap", group.by = "seurat_clusters")
    
    all_samples_integrated <- FindClusters(all_samples_integrated, resolution = 0.2)
    p2 <- DimPlot(all_samples_integrated, reduction="umap", group.by = "seurat_clusters")
    
    all_samples_integrated <- FindClusters(all_samples_integrated, resolution = 0.25)
    p3 <- DimPlot(all_samples_integrated, reduction="umap", group.by = "seurat_clusters")
    
    all_samples_integrated <- FindClusters(all_samples_integrated, resolution = 0.3)
    p4 <- DimPlot(all_samples_integrated, reduction="umap", group.by = "seurat_clusters")
    
    (p1 + p2 + p3 + p4) + plot_layout(ncol = 4)



###Heatmap of top 10 genes per cluster across resolutions
    #Define the resolutions you want to test
    resolutions <- c(0.2, 0.25, 0.3)
    
    #Store heatmaps in a list
    heatmap_list <- list()
    
    for (res in resolutions) {
      
      #Re-cluster at this resolution
      all_samples_integrated <- FindClusters(all_samples_integrated, resolution = res)
      
      #Find all markers for this clustering
      markers <- FindAllMarkers(all_samples_integrated, assay = "RNA", 
                                logfc.threshold = 0.585, min.pct = 0.2, only.pos = FALSE) %>% 
                                filter(p_val_adj < 0.05)
      
      #Select top 10 genes per cluster (positive logFC > 1)
      top_genes <- markers %>%
        group_by(cluster) %>%
        filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup()
      
      #Generate heatmap
      hm <- DoHeatmap(all_samples_integrated, features = top_genes$gene, assay = "RNA") +
        ggtitle(paste0("Resolution ", res)) +
        theme(text = element_text(size = 14))
      
      #Store in list
      heatmap_list[[paste0("res_", res)]] <- hm}
    
    #Combine plots (optional, adjust layout as needed)
    wrap_plots(heatmap_list, ncol = 3)
