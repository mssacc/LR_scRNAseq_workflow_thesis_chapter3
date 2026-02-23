###Single-cell level DE analysis

#Load required libraries
library(DESeq2)
library(Seurat)
library(EnhancedVolcano)
library(dittoSeq)
library(dplyr)
library(ggplot2)
library(patchwork)

#Load gene RDS
all_samples_integrated <- readRDS("/data/gpfs/projects/punim1901/flames_v2/seurat_workspace/all_samples_integrated.rds")

#Rejoin gene datasets after integration
all_samples_integrated <- JoinLayers(all_samples_integrated)


#DEG analysis between fertility status - bulk
    #Add column denoting fertility status and cell subtype
    pseudo_all_genes <- AggregateExpression(all_samples_integrated, assays="RNA", 
                                            return.seurat=T, 
                                            group.by = c("fertility", "endo.ID", "cell_type"))
    #Run if doing by fertility
    Idents(pseudo_all_genes) <- "fertility"
    
    #Run analysis
    bulk <- FindMarkers(pseudo_all_genes, ident.1="Infertile", ident.2="Fertile", test.use = "DESeq2")
    
    #List of significant DEGs
    upregulated <- subset(bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    downregulated <- subset(bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
        
    #Define colours
        keyvals <- ifelse(
        bulk$avg_log2FC < -0.585, 'red',
        ifelse(bulk$avg_log2FC > 0.585, 'green4', 'black'))
        keyvals[is.na(keyvals)] <- 'black'
        names(keyvals)[keyvals == 'green4'] <- 'Upregulated'
        names(keyvals)[keyvals == 'red'] <- 'Downregulated'
        names(keyvals)[keyvals == 'black'] <- 'No Significance'
    
    #Volcano plot with all DEGs combined
    EnhancedVolcano(bulk, lab=rownames(bulk),
                    x='avg_log2FC', y='p_val_adj',
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, labSize = 0,
                    colCustom = keyvals,
                    title = "Differentially Expressed Genes - Fertile vs Infertile", subtitle = NULL,
                    legendPosition = 'bottom')


#DEG analysis within the same cell subtype between conditions
    #Add column denoting fertility status and cell subtype
    pseudo_all_genes <- AggregateExpression(all_samples_integrated, assays="RNA", 
                                            return.seurat=T, 
                                            group.by = c("fertility", "endo.ID", "cell_type"))
    #Run if doing by individual cluster
    pseudo_all_genes$celltype.fertility <- paste(pseudo_all_genes$cell_type, pseudo_all_genes$fertility, sep="_")
    Idents(pseudo_all_genes) <- "celltype.fertility"
    
    #Run analysis
    c0_bulk <- FindMarkers(pseudo_all_genes, ident.1="0_Infertile", ident.2="0_Fertile", test.use = "DESeq2")
    c1_bulk <- FindMarkers(pseudo_all_genes, ident.1="1_Infertile", ident.2="1_Fertile", test.use = "DESeq2")
    c2_bulk <- FindMarkers(pseudo_all_genes, ident.1="2_Infertile", ident.2="2_Fertile", test.use = "DESeq2")
    c3_bulk <- FindMarkers(pseudo_all_genes, ident.1="3_Infertile", ident.2="3_Fertile", test.use = "DESeq2")
    c4_bulk <- FindMarkers(pseudo_all_genes, ident.1="4_Infertile", ident.2="4_Fertile", test.use = "DESeq2")
    c5_bulk <- FindMarkers(pseudo_all_genes, ident.1="5_Infertile", ident.2="5_Fertile", test.use = "DESeq2")
    
    #List of significant DEGs
    c0_up <- subset(c0_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c0_down <- subset(c0_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    c1_up <- subset(c1_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c1_down <- subset(c1_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    c2_up <- subset(c2_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c2_down <- subset(c2_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    c3_up <- subset(c3_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c3_down <- subset(c3_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    c4_up <- subset(c4_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c4_down <- subset(c4_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    c5_up <- subset(c5_bulk, subset= avg_log2FC>0.585 & p_val_adj<0.05)
    c5_down <- subset(c5_bulk, subset= avg_log2FC<(-0.585) & p_val_adj<0.05)
    
    #Merging up and downregulated genes into 1 dataset
    c0_DEGs <- rbind(c0_up,c0_down)
    c1_DEGs <- rbind(c1_up,c1_down)
    c2_DEGs <- rbind(c2_up,c2_down)
    c3_DEGs <- rbind(c3_up,c3_down)
    c4_DEGs <- rbind(c4_up,c4_down)
    c5_DEGs <- rbind(c5_up,c5_down)
    
    #Merging all datasets
    #Add column with cell subtype
    c0_bulk['cell_type']='pre_unciliated'
    c1_bulk['cell_type']='unciliated'
    c2_bulk['cell_type']='ciliated'
    c3_bulk['cell_type']='secretory'
    c4_bulk['cell_type']='pre_ciliated'
    c5_bulk['cell_type']='proliferative'
    
    #Convert gene names into a column
    c0_bulk <- tibble::rownames_to_column(c0_bulk, "gene")
    c1_bulk <- tibble::rownames_to_column(c1_bulk, "gene")
    c2_bulk <- tibble::rownames_to_column(c2_bulk, "gene")
    c3_bulk <- tibble::rownames_to_column(c3_bulk, "gene")
    c4_bulk <- tibble::rownames_to_column(c4_bulk, "gene")
    c5_bulk <- tibble::rownames_to_column(c5_bulk, "gene")
    
    #Volcano plots divided by cluster
    #Pre-Unciliated cells
    p1 <- EnhancedVolcano(c0_bulk, lab=c0_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "DEGs - Pre-Unciliated Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#F8766D'),
                    legendPosition = 'bottom',
                    xlim=c(-2,6), ylim=c(0,38),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    titleLabSize = 24)
    p1
    #Unciliated Cells
    p2 <- EnhancedVolcano(c1_bulk, lab=c1_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "Differentially Expressed Genes - Unciliated Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#ABA300'),
                    legendPosition = 'bottom',
                    xlim=c(-2,6), ylim=c(0,30),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    max.overlaps = 20,
                    titleLabSize = 24)
    #Ciliated Cells
    p3 <- EnhancedVolcano(c2_bulk, lab=c2_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "Differentially Expressed Genes - Ciliated Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#0CB702'),
                    legendPosition = 'bottom',
                    xlim=c(-3,6), ylim=c(0,30),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    titleLabSize = 24)
    #Secretory Cells
    p4 <- EnhancedVolcano(c3_bulk, lab=c3_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "Differentially Expressed Genes - Secretory Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#00BFC4'),
                    legendPosition = 'bottom',
                    xlim=c(-2,5), ylim=c(0,25),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    titleLabSize = 24)
    #Pre-Ciliated Cells
    p5 <- EnhancedVolcano(c4_bulk, lab=c4_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "Differentially Expressed Genes - Pre-Ciliated Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#849AFF'),
                    legendPosition = 'bottom',
                    xlim=c(-3,6), ylim=c(0,20),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    titleLabSize = 24)
    #Proliferative Cells
    p6 <- EnhancedVolcano(c5_bulk, lab=c5_bulk$gene, 
                    x='avg_log2FC', y='p_val_adj',
                    title = "Differentially Expressed Genes - Proliferative Cells", subtitle = NULL,
                    pCutoff=0.05, FCcutoff=0.26303,
                    pointSize = 3, 
                    col=c('black','black','black','#FF61CC'),
                    legendPosition = 'bottom',
                    xlim=c(-3,5), ylim=c(0,15),
                    drawConnectors=TRUE, widthConnectors = 0.5,
                    arrowheads = FALSE,
                    titleLabSize = 24)
    
    #Combine all plots into a grid (2 rows, 3 columns)
    combined_plot <- (p1 | p2 | p3) / 
      (p4 | p5 | p6)
    
    #Show the combined plot
    combined_plot
