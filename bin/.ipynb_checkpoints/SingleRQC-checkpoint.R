library_load <- suppressMessages(
    list(
        
        # Single R
        library(SingleCellExperiment),
        library(SingleR)
        
    )
)

#####################
### SingleRSeurat ###
#####################

setGeneric("SingleRSeurat", valueClass=c("Seurat"), def=function(so="Seurat", ref="SummarizedExperiment", cluster="character", ...) {standardGeneric("SingleRSeurat")})

setMethod("SingleRSeurat", signature=c(so="Seurat", ref="SummarizedExperiment", cluster="character"), definition=function(so, ref, cluster=NULL) {
    
    
    # Seurat object to SingleCellExperiment
    sce <- SummarizedExperiment(list(counts=so@assays$RNA@counts), colData = so@meta.data)
    
    if(is.null(cluster)) {
        
        message("Running SingleR for each cell")
        
        # Predict labels
        label_main <- SingleR::SingleR(test=sce, ref=ref, labels=ref$label.main, assay.type.test="counts", de.method="classic")
        label_fine <- SingleR::SingleR(test=sce, ref=ref, labels=ref$label.fine, assay.type.test="counts", de.method="classic")
        
        # Add labels to Seurat object
        label_main_meta <- as.data.frame(label_main) %>% dplyr::select(pruned.labels, tuning.scores.first) %>% dplyr::rename(main_labels_cluster=pruned.labels, main_delta_score_cluster=tuning.scores.first)
        label_fine_meta <- as.data.frame(label_fine) %>% dplyr::select(pruned.labels, tuning.scores.first) %>% dplyr::rename(fine_labels_cluster=pruned.labels, fine_delta_score_cluster=tuning.scores.first)
        
        so <- AddMetaData(so, cbind(label_main_meta, label_fine_meta))
        
    } else if((!is.null(cluster)) & (cluster %in% colnames(so@meta.data))) {
        
        message("Running SingleR in cluster mode")
        
        # Predict labels
        label_main <- SingleR::SingleR(test=sce, ref=ref, labels=ref$label.main, assay.type.test="counts", de.method="classic", clusters=sce$SCT_snn_res.0.8)
        label_fine <- SingleR::SingleR(test=sce, ref=ref, labels=ref$label.fine, assay.type.test="counts", de.method="classic", clusters=sce$SCT_snn_res.0.8)
        
        # Rename and combine SingleR label 
        label_main_meta <- as.data.frame(label_main) %>% dplyr::select(pruned.labels, tuning.scores.first) %>% dplyr::rename(main_labels_cluster=pruned.labels, main_delta_score_cluster=tuning.scores.first)
        label_fine_meta <- as.data.frame(label_fine) %>% dplyr::select(pruned.labels, tuning.scores.first) %>% dplyr::rename(fine_labels_cluster=pruned.labels, fine_delta_score_cluster=tuning.scores.first)
        
        label_main_meta[, cluster] <- rownames(label_main_meta)
        label_fine_meta[, cluster] <- rownames(label_fine_meta)
        
        label_main_meta <- dplyr::left_join(so@meta.data[, cluster, drop=FALSE], label_main_meta, by = cluster)
        label_fine_meta <- dplyr::left_join(so@meta.data[, cluster, drop=FALSE], label_fine_meta, by = cluster)
        
        # Add labels to Seurat object
        so <- AddMetaData(so, c(label_main_meta, label_fine_meta))
        
    } else if(!(is.null(cluster)) & !(cluster %in% colnames(so@meta.data))) {
        
        warning("Running SinglerR in cluster mode: Cluster variable not found in Seurat meta data")
        
    }
        
    return(so)
    
}
         )

###########################
### SingleRScoreHeatMap ###
###########################

setGeneric("SingleRScoreHeatMap", valueClass=c("gtree"), def=function(singler_obj="Seurat", color="character", main="character", ...) {standardGeneric("SingleRScoreHeatMap")})

setMethod("SingleRScoreHeatMap", signature=c(singler_obj="DFrame", color="character", main="character"), definition=function(singler_obj, color, main) {
    
    # SingleR score for matrix 
    singler_score <- singler_obj$scores
    rownames(singler_score) <- rownames(singler_obj)

    singler_score <- t(as.data.frame(singler_score))
    singler_score <- singler_score[match(names(color), rownames(singler_score)), ]
    
    # SingleR labels for top annotation 
    singler_labes <- factor(singler_obj$pruned.labels, levels=names(color))
    names(singler_labes) <- rownames(singler_obj)
    
    top_annotation <- HeatmapAnnotation(
        value=singler_labes, 
        col=list(value=color), 
        annotation_legend_param=list(value=list(title = "Label")),
        show_annotation_name=FALSE
    )
    
    # Colors
    color_break <- seq(min(singler_score), max(singler_score), 0.01)
    color_ramp <- viridis(length(color_break), option="magma")
    
    hm <- grid.grabExpr(
        draw(
            ComplexHeatmap::pheatmap(
                mat=singler_score,
                main=main, 
                name="Score", 
                fontsize_row=10,
                scale="none",
                cluster_rows=FALSE,
                cluster_cols=TRUE,
                cellwidth=0.08, 
                cellheight=10, 
                treeheight_row=30,
                treeheight_col=30,
                show_rownames=TRUE,
                show_colnames=FALSE,
                color=color_ramp,
                breaks=color_break, 
                top_annotation=top_annotation, 
                border_color=NA)
            )
        )
    
    return(hm)
    
}
         )