dim_clust <- function(
    
    # Seurat 
    so, 
    assay = NULL,
    # PCA features
    blacklist_genes = NULL, 
    # Dim reduction
    reduction = NULL, # Default name_pca
    dims_pca  = 30,  # Default 20 - RunPCA dims
    dims_umap = 30,  # Default NULL - RunUMAP dims 
    dims_tsne = 20,  # Default 1:5
    min_dist  = 0.3,  # Default 0.3  - RunUMAP dmin.ist - controls how tightly the embedding
    # Cluster 
    dims_cluster = 30,  # Default 1:10 - FindNeighbors dims  
    cluster_res  = 0.8 # Default 0.8 - FindClusters resoluton - above(below) 1.0 to obtain larger(smaller) number of communities
) 
{
    
    # Set default assay 
    DefaultAssay(so) <- assay
    
    # Set prefix
    prefix    <- tolower(assay)
    name_pca  <- paste0(prefix, "_pca")
    name_umap <- paste0(prefix, "_umap")
    name_tsne <- paste0(prefix, "_tsne")
    name_nn   <- paste0(prefix, "_nn")
    name_snn  <- paste0(prefix, "_snn")
    name_nno  <- paste0(prefix, "_nno")
    
    # Set reduction 
    if(is.null(reduction)) {reduction <- name_pca} else {reduction <- reduction}

    # Subset variable features by blacklist genes
    if(is.null(blacklist_genes)) {features_pca <- NULL} else {features_pca <- VariableFeatures(so)[!VariableFeatures(so) %in% blacklist_genes]}

    DefaultAssay(so) <- assay
    so <- RunPCA(
        object         = so,
        assay          = assay,
        features       = features_pca,
        npcs           = dims_pca,
        reduction.name = name_pca
    )

    so <- FindNeighbors(
        object          = so,
        assay           = assay,
        dims            = 1:dims_cluster,
        k.param         = 20,
        reduction       = reduction,
        graph.name      = c(name_nn, name_snn)
    )

    so <- FindNeighbors(
        object          = so,
        assay           = assay,
        dims            = 1:dims_cluster,
        k.param         = 20,
        reduction       = reduction,
        return.neighbor = TRUE,
        graph.name      = name_nno
    )

    so <- FindClusters(
        object         = so,
        assay          = assay,
        resolution     = cluster_res,
        graph.name     = name_snn
    )

    so <- RunUMAP(
        object         = so,
        assay          = assay,
        n.neighbors    = 20,
        nn.name        = name_nno,
        reduction.name = paste0(name_umap, "_nno")
    )

    so <- RunUMAP(
        object         = so,
        assay          = assay,
        n.neighbors    = 20,
        dims           = 1:dims_umap,
        metric         = "euclidean",
        min.dist       = min_dist,
        reduction      = reduction,
        reduction.name = name_umap
    )

    so <- RunTSNE(
        object           = so,
        assay            = assay,
        reduction        = reduction,
        reduction.name   = name_tsne,
        dims             = 1:dims_tsne,
        check_duplicates = FALSE
    )

    return(so)

}
