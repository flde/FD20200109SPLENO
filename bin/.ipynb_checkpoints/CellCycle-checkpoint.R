##############################
### CellCycle Class Object ###
##############################

setOldClass("igraph", igraph::make_empty_graph()) # Is S3 Class and needs to be reset for S4

setClass("CellCycle", slots=list(sample_name="character", so="Seurat", genes="data.frame", corr="list", graph="igraph"))

####################
### ToCorrMatrix ###
####################

setGeneric("CellCycle", valueClass="CellCycle", def=function(so="Seurat", genes="data.frame", ...) {standardGeneric("CellCycle")})

setMethod("CellCycle", signature=c(so="Seurat", genes="data.frame"), definition=function(sample_name="hello", so, genes) {
    
    cc <- new("CellCycle", sample_name=sample_name, so=so, genes=genes)
    
    return(cc)
    
    
}
         )

####################
### ToCorrMatrix ###
####################

setGeneric("ToCorrMatrix", valueClass="list", def=function(object) {standardGeneric("ToCorrMatrix")})

setMethod("ToCorrMatrix", signature="DFrame", definition=function(corr) {
    
    #' Takes output from correlatePairs and returns a correlation matrix for rho and fdr
    #' 
    #' @corr output DataFrame from correlatePairs
    
    library_load <- suppressMessages(library(reshape2))
    
    # Extend correlatePairs data frame to rho correlation matrix
    corr <- as.data.frame(corr)
    corr <- rbind(corr, data.frame(gene1=corr$gene2, gene2=corr$gene1, corr[, 3:5]))

    # Compute correlation matrix with rho and FDR values
    rho <- dcast(corr, gene1~gene2, value.var="rho")
    fdr <- dcast(corr, gene1~gene2, value.var="FDR")
    
    # Add gene rownames and remove gene name column 
    rownames(rho) <- rho$gene1
    rho <- rho[, -1]
    rownames(fdr) <- fdr$gene1
    fdr <- fdr[, -1]
    
    # Remove genes which are not correlated at all 
    rho <- rho[!apply(rho, 1, function(x) all(is.na(x))), !apply(rho, 2, function(x) all(is.na(x)))]
    fdr <- fdr[rownames(fdr) %in% rownames(rho), colnames(fdr) %in% colnames(rho)]

    # Set the diagonal of the correlation matrix for rho to 1 and for FDR to 0
    diag(rho) <- 1
    diag(fdr) <- 0

    return(list(rho=rho, fdr=fdr))
}
         )

##################
### CorrMatrix ###
##################

setGeneric("CorrMatrix", valueClass="CellCycle", def=function(object) {standardGeneric("CorrMatrix")})

setMethod("CorrMatrix", signature="CellCycle", definition=function(object) {
    
    # library 
    library_load <- suppressMessages(library(scran))
    
    # Get normalized count data from Seurat object 
    cnt <- GetAssayData(object@so, assay="RNA", slot="data")
    cnt <- cnt[rownames(cnt) %in% object@genes$gene, ]

    # CorrelatePairs
    corr <- scran::correlatePairs(cnt)

    # CorrelatePairs to list of correlation matrix for rho and fdr 
    corr <- ToCorrMatrix(corr)
              
    object@corr <- corr
    
    return(object)
}
         )

##################
### CorrClust ###
##################

setGeneric("CorrClust", valueClass="hclust", def=function(corr="list") {standardGeneric("CorrClust")})

setMethod("CorrClust", signature=c(corr="list"), definition=function(corr) {
    
    #' Compute hierarchical cluster from correlation matrix with dissimilarity euclidean 1-R distance matrix
    #' 
    #' @param corr List from corr with rho and FDR
    
    # Library 
    library_load <- suppressMessages(library(data.table))

    # Get rho correlation matrix 
    rho <- corr[["rho"]]
                                                                     
    # Compute euclidean distance of 1-R 
    corr_dist <- dist(1-rho, method="euclidean")
        
    # Compute hierarchical clustering 
    corr_clust <- hclust(corr_dist, method="complete")
    
    return(corr_clust)
}
         )
          
###################
### CorrHeatMap ###
###################
          
setGeneric("CorrHeatMap", valueClass="gTree", def=function(object="CellCycle", value_var="character") {standardGeneric("CorrHeatMap")})

setMethod("CorrHeatMap", signature=c(object="CellCycle", value_var="character"), definition=function(object, value_var="rho") {
    
    #' Compute heatmap based on correlation and hierarchical clustering 
    #' 
    #' @param corr List of correlation matrix for gene sets.
    #' @param corr_clust List of hierarchical cluster corresponding to corr.
    
    library_load <- suppressMessages(list(library(ComplexHeatmap), library(circlize)))

    # Compute dendrogram from hierarchical clustering 
    dend <- as.dendrogram(CorrClust(object@corr))

    # Get data matrix 
    mat <- object@corr[[value_var]]
        
    # Set data matrix for value_var
    if(value_var == "rho") {
        
        # Color 
        color_break <- seq(-1, 1, 0.01)
        color_ramp <- viridis(length(color_break), option="plasma")
        
        # Legend title 
        title <- "rho"
            
    }
        
    if(value_var == "fdr") {
        
        # -log10
        mat[mat == 0] <- .Machine$double.xmin
        mat <- -log10(mat)
            
        # Color
        color_break <- seq(0, -log10(.Machine$double.xmin), 0.01)
        color_ramp <- viridis(length(color_break), option="plasma")
        
        # Legend title 
        title <- "fdr"
            
    }           
        
    # Expression mean and frequency annotation plots
    cnt <- GetAssayData(object@so, assay="RNA", slot="counts")
    cnt <- cnt[rownames(cnt) %in% rownames(mat), ]
    cnt_freq <- 100 * rowSums(cnt >= 1) / ncol(cnt)
        
    cnt_norm <- GetAssayData(object@so, assay="RNA", slot="data")
    cnt_norm <- cnt_norm[rownames(cnt_norm) %in% rownames(mat), ]
    cnt_norm[cnt <= 1] <- NA
    cnt_norm <- apply(cnt_norm, 1, function(x) if(all(is.na(x))) {x=rep(0, length(x))} else{x}) %>% t()
    cnt_norm_means <- rowMeans(cnt_norm, na.rm=TRUE)

    cnt_freq_col=circlize::colorRamp2(seq(0, max(cnt_freq), length=100), rev(hcl.colors(100,"Purples")))
    cnt_norm_means_col=circlize::colorRamp2(seq(0, max(cnt_norm_means), length=100), rev(hcl.colors(100,"Blues")))
        
    df <- data.frame(ExpMean=cnt_norm_means, Freq=cnt_freq)
        
    right_annotation <- rowAnnotation(df=df, col=list(ExpMean=cnt_norm_means_col, Freq=cnt_freq_col), show_annotation_name=FALSE)
    
    # Cell cycle phase class annotation
    genes <- object@genes
    genes <- genes[genes$gene %in% rownames(mat), ]
    genes <- genes[match(rownames(mat), genes$gene), ]
    
    cc_phase_class <- factor(toupper(genes$cc_phase_class), levels=names(color$cc_phase_class))
    names(cc_phase_class) <- genes$gene
    cc_phase_class <- cc_phase_class[labels(dend)]
    
    top_annotation <- HeatmapAnnotation(
        value=cc_phase_class, 
        col=list(value=color$cc_phase_class), 
        annotation_legend_param=list(value=list(title = "Phase")), 
        show_annotation_name=FALSE, 
        show_legend=FALSE
    )
                      
    left_annotation <- rowAnnotation(
        value=cc_phase_class, 
        col=list(value=color$cc_phase_class), 
        annotation_legend_param=list(value=list(title = "Phase")), 
        show_annotation_name=FALSE
    )
                      
    # Plot heatmap 
    corr_hm <- grid.grabExpr(
        draw(
            ComplexHeatmap::pheatmap(
                mat=mat,
                main=object@sample_name,
                name=title, 
                fontsize_row=10,
                scale="none",
                cluster_rows=dend,
                cluster_cols=dend,
                cellwidth=1, 
                cellheight=1, 
                treeheight_row=20,
                treeheight_col=20,
                show_rownames=FALSE,
                show_colnames=FALSE,
                color=color_ramp,
                breaks=color_break, 
                right_annotation=right_annotation, 
                top_annotation=top_annotation, 
                left_annotation=left_annotation, 
                border_color =NA)
        )
    )
    
    return(corr_hm)
}
         )
    
#################
### CorrGraph ###
#################
          
setGeneric("CorrGraph", valueClass="CellCycle", def=function(object="CellCycle") {standardGeneric("CorrGraph")})

setMethod("CorrGraph", signature=c(object="CellCycle"), definition=function(object) {
    
    library_load <- suppressMessages(list(library(igraph)))
    
    rho <- object@corr[["rho"]]
    fdr <- object@corr[["fdr"]]

    rho <- as.matrix(rho)
    fdr <- as.matrix(fdr)

    # Get igraph object for rho and fdr matrix
    graph_rho <- graph.adjacency(rho, weighted=TRUE, mode="undirected", diag=FALSE)
    graph_fdr <- graph.adjacency(fdr, weighted=TRUE, mode="undirected", diag=FALSE)
        
    # Add fdr to graph_rho edge 
    E(graph_rho)$rho <- E(graph_rho)$weight
    E(graph_rho)$fdr <- E(graph_fdr)$weight
    
    # Color edges 
    E(graph_rho)[which(E(graph_rho)$weight < 0)]$color <- RColorBrewer::brewer.pal(11, "RdBu")[2]
    E(graph_rho)[which(E(graph_rho)$weight > 0)]$color <- RColorBrewer::brewer.pal(11, "RdBu")[10]
    E(graph_rho)[which(E(graph_rho)$weight == 0)]$color <- "gray"

    # Scale size of verdicts by frequency expression in single cells 
    cnt <- GetAssayData(object@so, assay="RNA", slot="counts")
    cnt <- cnt[rownames(cnt) %in% V(graph_rho)$name, ]
    cnt <- cnt[match(V(graph_rho)$name, rownames(cnt)), ]
    cnt_freq <- rowSums(cnt >= 1) / ncol(cnt)
    V(graph_rho)$cnt_freq <- cnt_freq
    V(graph_rho)$size <- 20*cnt_freq
    
    # Set color of verdicts by cc_phase_class
    genes <- object@genes[object@genes$gene %in% V(graph_rho)$name, ]
    genes <- genes[match(V(graph_rho)$name, genes$gene), ]
    V(graph_rho)$cc_phase_class <- toupper(genes$cc_phase_class)
    V(graph_rho)$color <- color$cc_phase_class[V(graph_rho)$cc_phase_class]
    
    # Update CellCycle object
    object@graph <- graph_rho
    
    return(object)
    
}
         )

####################
### VertexSelect ###
####################
    
setGeneric("ModuleSelect", valueClass=c(object="CellCycle"), def=function(object="CellCycle", module="character") {standardGeneric("ModuleSelect")})

setMethod("ModuleSelect", signature=c(object="CellCycle", module="character"), definition=function(object, module) {
    
    genes <- object@genes[toupper(object@genes$cc_phase_class) == module, ]
    corr <- lapply(object@corr, function(corr) {corr[rownames(corr) %in% genes$gene, colnames(corr) %in% genes$gene]})
    
    object_module <- CellCycle(sample_name=object@sample_name, so=object@so, genes=genes)
    object_module@corr <- corr

    return(object_module)
    
}
         )
    
######################
### CorrGraphPrune ###
######################
    
setGeneric("CorrGraphPrune", valueClass=c(object="CellCycle"), def=function(object="CellCycle", ...) {standardGeneric("CorrGraphPrune")})

setMethod("CorrGraphPrune", signature=c(object="CellCycle"), definition=function(object, fdr_max=0.01, rho_min=0, rho_direction=NULL, degree_min=0, cnt_freq_min=0) {
    
    # Get graph
    graph <- object@graph
    
    # Delte edges 
    graph <- delete_edges(graph, E(graph)[which(E(graph)$fdr > fdr_max)])
    graph <- delete_edges(graph, E(graph)[which(abs(E(graph)$rho) < rho_min)])
    if(!is.null(rho_direction)) {graph <- delete_edges(graph, E(graph)[which(sign(E(graph)$rho) != rho_direction)])}
    
    # Delte verdicts 
    graph <- delete_vertices(graph, V(graph)[which(V(graph)$cnt_freq < cnt_freq_min)])
    graph <- delete_vertices(graph, V(graph)[which(degree(graph) < degree_min)])
    
    # Update CellCycle object
    object@graph <- graph

    return(object)
}
         )
    
############################
### CorrGraphModulePrune ###
############################
    
setGeneric("CorrGraphModulePrune", valueClass=c(object="CellCycle"), def=function(object="CellCycle", ...) {standardGeneric("CorrGraphModulePrune")})

setMethod("CorrGraphModulePrune", signature=c(object="CellCycle"), definition=function(object, fdr_max=0.01, rho_min=0, rho_direction=NULL, degree_min=0) {

    modules <- unique(V(object@graph)$cc_phase_class)

    genes_prune <- list()
    for(module in modules) {

        object_module <- ModuleSelect(object, module=module) 
        object_module <- CorrGraph(object_module)
        object_module <- CorrGraphPrune(object_module, fdr_max=fdr_max, rho_min=rho_min, rho_direction=rho_direction, degree_min=degree_min)

        genes_prune[[module]] <- object_module@genes[object_module@genes$gene %in% V(object_module@graph)$name, ]

    }

    genes_prune <- do.call("rbind", genes_prune); rownames(genes_prune) <- NULL

    # Run CellCycle pipeline on pruned object 
    object_prune <- CellCycle(sample_name=object@sample_name, so=object@so, genes=genes_prune)
    object_prune <- CorrMatrix(object_prune)
    object_prune <- CorrGraph(object_prune)
    object_prune <- CorrGraphPrune(object_prune)
    
    return(object_prune)
    
}
         )
    
    
    

    
#####################
### CorrGraphGrob ###
#####################

setGeneric("CorrGraphGrob", valueClass=c(object="grob"), def=function(object="CellCycle", ...) {standardGeneric("CorrGraphGrob")})

setMethod("CorrGraphGrob", signature=c(object="CellCycle"), definition=function(object) {

    # as.grob is looking for object in the global environment
    object <<- object
    
    graph_grob <- as.grob(expression(plot(object@graph, layout=layout_on_sphere, vertex.label.color=NA, vertex.label=NA, main=object@sample_name)))
    legend_grob <- legendGrob(names(color$cc_phase_class), pch=19, gp=gpar(col=color$cc_phase_class, fill=color$cc_phase_class))
    
    grob <- packGrob(packGrob(frameGrob(), graph_grob), legend_grob, height=unit(1, "null"), side="left")


    return(grob)

}
         )
    
###############
### ccScore ###
###############

setGeneric("ccScore", valueClass=c(object="CellCycle"), def=function(object="CellCycle", ...) {standardGeneric("ccScore")})

setMethod("ccScore", signature=c(object="CellCycle"), definition=function(object) {
    
    # Get slots from CellCycle object 
    so <- object@so
    genes <- object@genes
    
    # Remove meta data from Seurat object which are cell cycle related to avoid doubling 
    col_remove <- colnames(so@meta.data)[grep("msG1S|msS|msG2|msG2M", colnames(so@meta.data))]
    if(length(col_remove) > 0) {
        
        so@meta.data <- so@meta.data[, -which(colnames(so@meta.data) %in% col_remove)]
    }
    
    # Convert genes to list by cc_phase_class    
    genes <- split(genes, f=genes$cc_phase_class)
    genes <- lapply(genes, function(x) x[, 1])
    
    # Score cell cycle 
    for (cc_phase_class in names(genes)) {
        
        so <- AddModuleScore(so, features=genes[cc_phase_class], assays="RNA", slot="data", ctrl=100, nbin=25, name=paste0("ms", toupper(cc_phase_class), "_RNA"))
        
    }

    # Compute cc_phase_class
    meta <- so@meta.data
    cc_meta <- meta[, (ncol(meta)-length(genes)+1):(ncol(meta))]
    
    # Needed for cc_phase_class_select
    cc_phase_class <- colnames(cc_meta)
    cc_phase_class <- gsub("ms|_RNA1","", cc_phase_class)

    cc_phase_class_select <- function(cc_score) {
        
        if(all(cc_score <= 0)) {
            cc_phase <- "G1"
        } else {
            cc_phase <- cc_phase_class[which.max(cc_score)]
        }

        return(cc_phase)

    }

    so$cc_phase_class <- apply(cc_meta, 1, cc_phase_class_select)
                    
    # Update CellCycle object
    object@so <- so
    
    return(object)
    
}
         )
    
######################
### ccScoreScatter ###
######################

setGeneric("ccScoreScatter", valueClass=c(object="ggplot"), def=function(object="CellCycle", ...) {standardGeneric("ccScoreScatter")})

setMethod("ccScoreScatter", signature=c(object="CellCycle"), function(object) {
    
    so <- object@so
    
    meta <- so@meta.data
    meta$cc_phase_class <- factor(meta$cc_phase_class, levels=names(color$cc_phase_class))
    
    cc_plot <- ggplot(meta, aes(x=msS_RNA1, y=msG2M_RNA1, color=cc_phase_class)) + 
        geom_point() + 
        geom_vline(aes(xintercept=0), color="black", linetype="longdash") +
        geom_hline(aes(yintercept=0), color="black", linetype="longdash") + 
        xlim(min(c(so$msS_RNA1, so$msG2M_RNA1)), max(c(so$msS_RNA1, so$msG2M_RNA1))) + 
        ylim(min(c(so$msS_RNA1, so$msG2M_RNA1)), max(c(so$msS_RNA1, so$msG2M_RNA1))) + 
        ggtitle(object@sample_name) + xlab("S-phase [Score]") + ylab("G2M [Score]") + 
        scale_color_manual(values=color$cc_phase_class) + 
        theme(aspect.ratio=1)
    
    return(cc_plot)
}
         )
    
######################
### ccScoreHeatMap ###
#####################      

setGeneric("ccScoreHeatMap", valueClass=c(object="gTree"), def=function(object="CellCycle", ...) {standardGeneric("ccScoreHeatMap")})

setMethod("ccScoreHeatMap", signature=c(object="CellCycle"), function(object) {
    
    so <- object@so
    
    # Get meta data from Seurat object
    meta <- so@meta.data
    
    # Get cell cycle meta information 
    cc_meta <- cbind(meta[, grep("msG1S|msS|msG2|msG2M", colnames(meta))], meta[, "cc_phase_class", drop=FALSE])
    cc_meta <- arrange(cc_meta, factor(cc_phase_class, levels=names(color$cc_phase_class)))

    # Sort cc_meta by cell cycle score 
    cc_meta <- split(cc_meta, f=cc_meta$cc_phase_class)
    cc_meta <- cc_meta[names(color$cc_phase_class)[names(color$cc_phase_class) %in% names(cc_meta)]]
    
    for(cc_phase_class in names(cc_meta)) {
        
        if(cc_phase_class == "G1") {
            
            cc_meta[[cc_phase_class]]$cc_score_min <- apply(cc_meta[[cc_phase_class]][, -ncol(cc_meta[[cc_phase_class]])], 1, min)
            cc_meta[[cc_phase_class]] <- arrange(cc_meta[[cc_phase_class]], cc_score_min)
            cc_meta[[cc_phase_class]]$cc_score_min <- NULL

        } else {
            
            cc_meta[[cc_phase_class]]$cc_score_max <- apply(cc_meta[[cc_phase_class]][, -ncol(cc_meta[[cc_phase_class]])], 1, max)
            cc_meta[[cc_phase_class]] <- arrange(cc_meta[[cc_phase_class]], cc_score_max)
            cc_meta[[cc_phase_class]]$cc_score_max <- NULL
        }

    }

    cc_meta <- do.call("rbind", cc_meta)
    
    # Extract cell cycle score 
    cc_phase_score <- cc_meta[, grep("msG1S|msS|msG2|msG2M", colnames(cc_meta))]
    
    cc_phase_score <- t(cc_phase_score)
    rownames(cc_phase_score) <- gsub("ms|_RNA1", "", rownames(cc_phase_score))
    cc_phase_score <- cc_phase_score[names(color$cc_phase_class)[names(color$cc_phase_class) %in% rownames(cc_phase_score)], ]

    # Extract cell cycle phase class 
    cc_phase_class <- cc_meta[, grep("cc_phase_class", colnames(cc_meta)), drop=FALSE]
    cc_phase_class$cc_phase_class <- factor(cc_phase_class$cc_phase_class, levels=names(color$cc_phase_class))

    # Top annotation 
    top_annotation <- HeatmapAnnotation(
        df=cc_phase_class, 
        col=list(cc_phase_class=color$cc_phase_class), 
        annotation_legend_param=list(cc_phase_class=list(title = "Phase")), 
        show_annotation_name=FALSE
    )

    # Colors 
    color_break <- seq(-(round(max(abs(cc_phase_score)), digits=1) + 0.1), (round(max(abs(cc_phase_score)), digits=1) + 0.1), 0.01)
    color_ramp <- rev(RColorBrewer::brewer.pal(12, "RdBu"))

    # Heat map 
    hm <- grid.grabExpr(
        draw(
            ComplexHeatmap::pheatmap(
                mat=cc_phase_score,
                main=object@sample_name, 
                name="Score", 
                fontsize_row=10,
                scale="none",
                cluster_rows=FALSE,
                cluster_cols=FALSE,
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

    