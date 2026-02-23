###Integration method comparison

#Load required libraries
library(SeuratObject)
library(Seurat)
library(harmony)
library(future)

#Load seurat object RDS
all_samples <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples.rds")

###Harmony integration
    all_samples_integrated <- IntegrateLayers(object=all_samples,
                                              method=HarmonyIntegration,
                                              orig.reduction="pca", new.reduction="integrated",
                                              verbose=FALSE)
    all_samples_integrated <- FindNeighbors(all_samples_integrated, dims=1:10, reduction="integrated")
    all_samples_integrated <- FindClusters(all_samples_integrated, resolution=0.25, cluster.name="integrated_cluster")
    all_samples_integrated <- RunUMAP(all_samples_integrated, reduction="integrated", dims=1:10, seed.use=123)
    pca <- RunPCA(all_samples_integrated)
    #Visualize data
    DimPlot(all_samples_integrated, reduction="umap")
    DimPlot(all_samples_integrated, reduction="umap", group.by=c("integrated_cluster", "orig.ident", "fertility"))


###CCA integration
    all_samples_integrated_2 <- IntegrateLayers(object=all_samples,
                                              method=CCAIntegration,
                                              orig.reduction="pca", new.reduction="integrated",
                                              verbose=FALSE)
    all_samples_integrated_2 <- FindNeighbors(all_samples_integrated_2, dims=1:10, reduction="integrated")
    all_samples_integrated_2 <- FindClusters(all_samples_integrated_2, resolution=0.25, cluster.name="integrated_cluster")
    all_samples_integrated_2 <- RunUMAP(all_samples_integrated_2, reduction="integrated", dims=1:10, seed.use=123)
    pca <- RunPCA(all_samples_integrated_2)
    #Visualize data
    DimPlot(all_samples_integrated_2, reduction="umap")
    DimPlot(all_samples_integrated_2, reduction="umap", group.by=c("integrated_cluster", "orig.ident", "fertility"))

    
#RPCA integration
    #Increase limit to 30 GiB
    options(future.globals.maxSize = 30 * 1024^3)
    
    all_samples_integrated_3 <- IntegrateLayers(object=all_samples,
                                                method=RPCAIntegration,
                                                orig.reduction="pca", new.reduction="integrated",
                                                verbose=FALSE)
    all_samples_integrated_3 <- FindNeighbors(all_samples_integrated_3, dims=1:10, reduction="integrated")
    all_samples_integrated_3 <- FindClusters(all_samples_integrated_3, resolution=0.25, cluster.name="integrated_cluster")
    all_samples_integrated_3 <- RunUMAP(all_samples_integrated_3, reduction="integrated", dims=1:10, seed.use=123)
    pca <- RunPCA(all_samples_integrated_3)
    #Visualize data
    DimPlot(all_samples_integrated_3, reduction="umap")
    DimPlot(all_samples_integrated_3, reduction="umap", group.by=c("integrated_cluster", "orig.ident", "fertility"))
