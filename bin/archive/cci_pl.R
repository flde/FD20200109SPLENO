#########################
### Heat map net slot ###
#########################
hm_net <- function(cellchat, slot, title, diff=NULL, celltype=NULL, cluster=FALSE) {
    
    if(is.null(diff)) {
        
        mat <- cellchat@net[[slot]]
        if(!is.null(celltype)) {mat <- mat[celltype, celltype]}
        mat_name <- "Count"
        mat_color <- viridis::magma(100)
        mat_breaks <- seq(0, max(mat), length.out=length(mat_color))
        mat_color_function <- circlize::colorRamp2(mat_breaks, mat_color) 
        
    } else {
        
        mat_1 <- cellchat@net[[1]][["count"]]
        mat_2 <- cellchat@net[[2]][["count"]]
        
        mat_1 <- mat_1[intersect(rownames(mat_1), rownames(mat_2)), intersect(rownames(mat_1), rownames(mat_2))]
        mat_2 <- mat_2[intersect(rownames(mat_1), rownames(mat_2)), intersect(rownames(mat_1), rownames(mat_2))]
        
        mat <- log2((mat_2+1)/(mat_1+1))
        if(!is.null(celltype)) {mat <- mat[celltype, celltype]}
        
        mat_name <- "log2FC"
        mat_color <- c(rev(RColorBrewer::brewer.pal(11,"RdBu"))[1:5], "#ffffff", rev(RColorBrewer::brewer.pal(11,"RdBu"))[7:11])
        mat_breaks <- seq(-max(abs(mat)), max(abs(mat)), length.out=length(mat_color))
        mat_color_function <- circlize::colorRamp2(mat_breaks, mat_color) 
        
    }
    
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        name=mat_name, 
        
        col=mat_color_function, 
        
        width=ncol(mat)*unit(6, "mm"), 
        height=ncol(mat)*unit(6, "mm"), 
        
        column_title=title, 
        column_title_gp=gpar(fontsize=18, fontface="bold"), 
        
        row_names_gp=gpar(fontsize=18, fontface="plain"), 
        column_names_gp=gpar(fontsize=18, fontface="plain"), 
        
        row_labels=rownames(mat), 
        column_labels=colnames(mat), 
        
        border=FALSE, 
        
        cluster_rows=cluster, 
        cluster_columns=cluster,
        
        show_row_names=TRUE,
        show_column_names=TRUE, 
        
        row_gap=unit(1, "mm"), 
        column_gap=unit(1, "mm"), 
        
        row_title_rot=90, 
    
        rect_gp=gpar(col="white", lwd=1), 
        
        heatmap_legend_param=list(title_gp=gpar(fontsize=18, fontface="plain"), labels_gp=gpar(fontsize=16), legend_height=unit(5*6, "mm"), grid_width=unit(6, "mm"))

    )

    return(hm)
}
#############################
### Ligand log2FC heatmap ###
#############################
dea_res_hm <- function(l_res, pval_thr=0.05, color_mat=NULL, row_title=NULL, cluster_rows=FALSE, cluster_columns=FALSE, show_column_dend=FALSE, show_column_names=TRUE, width=1.0, height=1.0, use_raster=FALSE, fontsize_select=1) {
    
    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Ensure unique ligand symbol in case of ligand complex 
    l_res[[1]] <- l_res[[1]] %>% dplyr::mutate(ligand.symbol=as.character(ligand.symbol)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub(ligand.symbol, "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("^_", "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("_$", "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("_", ", ", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand.symbol=ifelse(ligand.symbol!=ligand_complex & ligand_complex!="", paste0(ligand.symbol, " (", ligand_complex, ")"), ligand.symbol))

    l_res[[2]] <- l_res[[2]] %>% dplyr::mutate(ligand.symbol=as.character(ligand.symbol)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub(ligand.symbol, "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("^_", "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("_$", "", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand_complex=gsub("_", ", ", ligand_complex)) %>% 
        dplyr::rowwise() %>% dplyr::mutate(ligand.symbol=ifelse(ligand.symbol!=ligand_complex & ligand_complex!="", paste0(ligand.symbol, " (", ligand_complex, ")"), ligand.symbol))

    # Get ligand log2fc mat 
    l_log2fc_mat <- l_res[[1]] %>% dplyr::select(-ligand.symbol, -ligand_complex, -sender) %>% dplyr::distinct() %>% tibble::column_to_rownames("interaction_name")

    # Get row labels 
    row_labels <- l_res[[1]] %>% dplyr::select(interaction_name, ligand.symbol) %>% dplyr::filter(interaction_name %in% rownames(l_log2fc_mat)) %>% dplyr::distinct() %>% tibble::column_to_rownames("interaction_name")
    row_labels <- row_labels[rownames(l_log2fc_mat), , drop=FALSE]$ligand.symbol
    row_labels <- gsub("_", ", ", row_labels)

    # Column annotation     
    top_annotation <- HeatmapAnnotation(

        df=data.frame(celltype=rep(l_res[[1]]$sender[1], ncol(l_log2fc_mat))), 
        col=list(celltype=color$celltype_low), 
        simple_anno_size=unit(fontsize_scale*5, "mm"), 
        show_annotation_name=FALSE, 
        show_legend=FALSE, 
        border=TRUE

    )

    # Color mat 
    if(is.null(color_mat)) {
        
        
        color <- viridis::mako(length(breaks))
        
    } else {

        color <- c(color_mat[2], "#FFFFFF", color_mat[1])
        
    }

    breaks <- seq(-2, 2, length.out=3)
    color_function_mat <- colorRamp2(breaks, color)

    # Get ligand adj pval mat 
    l_adj_pval_mat <- l_res[[2]] %>% dplyr::select(-ligand.symbol, -ligand_complex, -sender) %>% dplyr::distinct() %>% tibble::column_to_rownames("interaction_name")
    l_adj_pval_mat <- l_adj_pval_mat[rownames(l_log2fc_mat), , drop=FALSE]
    
    # TF heatmap 
    hm <- Heatmap(

        matrix=l_log2fc_mat, 

        col=color_function_mat, 
        
        width=fontsize_scale*5*ncol(l_log2fc_mat)*unit(width, "mm"),
        height=fontsize_scale*5*nrow(l_log2fc_mat)*unit(height, "mm"), 

        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        column_title="Ligand", 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=cluster_rows, 
        row_labels=row_labels, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=NULL, 
        row_title=row_title , 
        row_gap=unit(0.5, "mm"),
        show_row_names=TRUE,
        row_names_side="left",

        cluster_columns=FALSE,
        clustering_distance_columns="pearson", 
        cluster_column_slices=FALSE, 
        show_column_dend=TRUE, 
        column_split=NULL,
        column_gap=unit(0.5, "mm"), 
        column_dend_height=unit(fontsize_scale*3, "mm"), 
        show_column_names=show_column_names, 

        top_annotation=top_annotation, 

        heatmap_legend_param=list(title="log2FC", at=c(min(breaks), 0, max(breaks)), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 

        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt")),

        use_raster=FALSE, raster_by_magick=TRUE, raster_resize_mat=mean, 

        cell_fun=function(j, i, x, y, w, h, fill) {
            
            if(l_adj_pval_mat[i, j] <= pval_thr) {
                
                grid.text("*", x, y-unit(fontsize_scale*1, "mm"), gp=gpar(fontsize=fontsize[1]))
            
            }
        
        }
            
    
    )

    return(hm)
    
}
    
#######################
### LR prob heatmap ###
#######################
lr_res_hm <- function(lr_res, pval_prob_thr=0.05, column_split_order=NULL, row_title=NULL, cluster_rows=FALSE, cluster_columns=FALSE, show_column_dend=FALSE, show_column_names=TRUE, width=1.0, height=1.0, use_raster=FALSE, fontsize_select=1) {

    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Get receptor mat 
    if("bottom_annotation" %in% colnames(lr_res)) {

        lr_pval_mat <- lr_res %>% dplyr::select(interaction_name, pval, sample_group, col_name, target, bottom_annotation) %>% 
            dplyr::mutate(pval=cut(pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), labels=c("***", "**", "*", "ns"))) %>% 
            dplyr::arrange(col_name) %>% 
            dplyr::mutate(col_group=paste0(col_name, ":", target, ":", bottom_annotation)) %>% dplyr::select(-sample_group, -col_name, -target, -bottom_annotation)
        
    } else {

        lr_pval_mat <- lr_res %>% dplyr::select(interaction_name, pval, sample_group, col_name, target) %>% 
            dplyr::mutate(pval=cut(pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), labels=c("***", "**", "*", "ns"))) %>% 
            dplyr::arrange(col_name) %>% 
            dplyr::mutate(col_group=paste0(col_name, ":", target)) %>% dplyr::select(-sample_group, -col_name, -target)
        
    }

    lr_pval_mat <- lr_pval_mat %>% pivot_wider(., names_from=col_group, values_from=pval) %>% tibble::column_to_rownames("interaction_name")
    
    # Get row labels 
    row_labels <- lr_res %>% dplyr::select(interaction_name, receptor_complex) %>% dplyr::filter(interaction_name %in% rownames(lr_pval_mat)) %>% dplyr::distinct() %>% tibble::column_to_rownames("interaction_name")
    row_labels <- row_labels[rownames(lr_pval_mat), , drop=FALSE]$receptor_complex
    row_labels <- gsub("_", ", ", row_labels)

    # Get column labels
    col_labels <- word(colnames(lr_pval_mat), 1, sep=":")
    
    # Column annotation
    top_annotation <- HeatmapAnnotation(

        df=data.frame(celltype=word(colnames(lr_pval_mat), 2, sep=":")), 
        col=list(celltype=color$celltype_low), 
        simple_anno_size=unit(fontsize_scale*5, "mm"), 
        show_annotation_name=FALSE, 
        show_legend=FALSE, 
        border=TRUE
        
    )

    # Bottom annotation
    if("bottom_annotation" %in% colnames(lr_res)) {

        bottom_annotation <- HeatmapAnnotation(

            df=data.frame(genotype=word(colnames(lr_pval_mat), 3, sep=":")), 
            col=list(genotype=c("(+/+)"="#000000", "(cre/+)"="#FFFFFF")), 
            simple_anno_size=unit(fontsize_scale*5, "mm"), 
            show_annotation_name=FALSE, 
            show_legend=FALSE, 
            border=TRUE
            
        )
        
    } else {

        bottom_annotation <- NULL
        
    }
    
    # Column split based on annotation
    col_split <- word(colnames(lr_pval_mat), 2, sep=":")
    if(!is.null(column_split_order)) {col_split <- factor(col_split, levels=column_split_order)}

    # Color mat 
    color_function_mat <- c("ns"="#FFFFFF", "*"="#ECE0FF", "**"="#B9AEF5", "***"="#595491")

    # TF heatmap 
    hm <- Heatmap(

        matrix=lr_pval_mat, 

        col=color_function_mat, 
        
        width=fontsize_scale*5*ncol(lr_pval_mat)*unit(width, "mm"),
        height=fontsize_scale*5*nrow(lr_pval_mat)*unit(height, "mm"), 

        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        column_title="Receptor", 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=cluster_rows, 
        row_labels=row_labels, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=NULL, 
        row_title=row_title , 
        row_gap=unit(0.5, "mm"),
        show_row_names=TRUE,
        row_names_side="right",

        cluster_columns=FALSE,
        column_label=col_labels,
        clustering_distance_columns="pearson", 
        cluster_column_slices=FALSE, 
        show_column_dend=TRUE, 
        column_split=col_split,
        column_gap=unit(1.0, "mm"), 
        column_dend_height=unit(fontsize_scale*3, "mm"), 
        show_column_names=show_column_names, 

        top_annotation=top_annotation,
        bottom_annotation=bottom_annotation, 

        heatmap_legend_param=list(title="p-value", at=c("ns", "*", "**", "***"), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm"), border = "black"), 

        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt")),

        use_raster=FALSE, raster_by_magick=TRUE, raster_resize_mat=mean
    
    )

    return(hm)
    
}