#########################
### Heat map net slot ###
#########################
hm_net <- function(cellchat, slot, title, mode="single", source=NULL, target=NULL, cluster=FALSE, breaks_max=NULL, width=1.2, height=1.2, fontsize_select=1) {

    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]
    
    if(mode=="single") {
        
        mat <- cellchat@net[[slot]]

        if(!is.null(source)) {mat <- mat[, source, drop=FALSE]}
        if(!is.null(target)) {mat <- mat[target, , drop=FALSE]}
        
        # Color mat 
        color <- c("#FFFFFF", rev(RColorBrewer::brewer.pal(11,"RdBu"))[11])
        if(is.null(breaks_max)) {breaks_max=max(mat)}
        breaks <- seq(0, breaks_max, length.out=2)
        mat_color_function <- circlize::colorRamp2(breaks, color)

        mat_name <- slot
        
    } else if(mode=="diff") {
        
        mat_1 <- cellchat@net[[1]][[slot]]
        mat_2 <- cellchat@net[[2]][[slot]]
        
        mat_1 <- mat_1[intersect(rownames(mat_1), rownames(mat_2)), intersect(rownames(mat_1), rownames(mat_2))]
        mat_2 <- mat_2[intersect(rownames(mat_1), rownames(mat_2)), intersect(rownames(mat_1), rownames(mat_2))]
        
        mat <- mat_2-mat_1
        
        if(!is.null(source)) {mat <- mat[, source, drop=FALSE]}
        if(!is.null(target)) {mat <- mat[target, , drop=FALSE]}
        
        mat_name <- "log2FC"
        mat_color <- c(rev(RColorBrewer::brewer.pal(11,"RdBu"))[1:5], "#ffffff", rev(RColorBrewer::brewer.pal(11,"RdBu"))[7:11])
        mat_breaks <- seq(-max(abs(mat)), max(abs(mat)), length.out=length(mat_color))
        mat_color_function <- circlize::colorRamp2(mat_breaks, mat_color) 
        
    }
    
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        name=mat_name, 
        
        col=mat_color_function, 
        
        width=fontsize_scale*5*ncol(mat)*unit(width, "mm"),
        height=fontsize_scale*5*nrow(mat)*unit(height, "mm"), 
        
        column_title=title, 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 
        
        row_names_gp=gpar(fontsize=fontsize[1], fontface="italic"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        row_labels=rownames(mat), 
        column_labels=colnames(mat), 
        
        cluster_rows=cluster, 
        cluster_columns=cluster,
        
        show_row_names=TRUE,
        show_column_names=TRUE, 
        
        row_gap=unit(1, "mm"), 
        column_gap=unit(1, "mm"), 
        
        row_title_rot=90, 

        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt")),
        
        heatmap_legend_param=list(title=mat_name, at=c(min(breaks), max(breaks)), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm"))

    )

    return(hm)
}

###################
### DEA heatmap ###
###################
dea_hm <- function(mat_1, mat_2, row_labels, row_splits, celltype, lr_type, dea_p_value_thr=0.01, width=1.2, height=1.2, use_raster=FALSE, fontsize_select=1) {
    
    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]
    
    # Column annotation     
    top_annotation <- HeatmapAnnotation(

        df=data.frame(celltype=rep(celltype, ncol(mat_1))), 
        col=list(celltype=color$celltype_low), 
        simple_anno_size=unit(fontsize_scale*5, "mm"), 
        show_annotation_name=FALSE, 
        show_legend=FALSE, 
        border=TRUE

    )

    # Color mat 
    color <- c(rev(RColorBrewer::brewer.pal(11,"RdBu"))[1], "#FFFFFF", rev(RColorBrewer::brewer.pal(11,"RdBu"))[11])
    breaks <- seq(-2, 2, length.out=3)
    color_function_mat <- circlize::colorRamp2(breaks, color)
    
    # TF heatmap 
    hm <- Heatmap(

        matrix=mat_1, 

        col=color_function_mat, 
        
        width=fontsize_scale*5*ncol(mat_1)*unit(width, "mm"),
        height=fontsize_scale*5*nrow(mat_1)*unit(height, "mm"), 

        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        column_title=lr_type, 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="italic"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=FALSE, 
        row_labels=row_labels, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=row_splits, 
        row_title=NA, 
        row_title_rot=0, 
        row_gap=unit(fontsize_scale*2, "mm"),
        show_row_names=TRUE,
        row_names_side=ifelse(lr_type=="Ligand", "left", "right"),

        cluster_columns=FALSE,
        clustering_distance_columns="pearson", 
        cluster_column_slices=FALSE, 
        show_column_dend=TRUE, 
        column_split=NULL,
        column_gap=unit(fontsize_scale*2, "mm"), 
        column_dend_height=unit(fontsize_scale*3, "mm"), 
        show_column_names=TRUE, 

        top_annotation=top_annotation, 

        heatmap_legend_param=list(title="log2FC", at=c(min(breaks), 0, max(breaks)), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 

        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt")),

        use_raster=FALSE, raster_by_magick=TRUE, raster_resize_mat=mean, 

        cell_fun=function(j, i, x, y, w, h, fill) {
            
            if(mat_2[i, j] <= dea_p_value_thr) {
                
                grid.text("*", x, y-unit(fontsize_scale*1, "mm"), gp=gpar(fontsize=fontsize[1]))
            
            }
        
        }
              
    )

    return(hm)
    
}

######################
### DEA LR heatmap ###
######################
dea_lr_hm <- function(lr_res, dea_res, source, target, dea_p_value_thr=0.01, dea_p_value_filter_thr=1, lr_de_check=FALSE, fontsize_select=1) {

    # Select DEA result for source 
    dea_res_l <- dea_res[[source]]
    dea_res_r <- dea_res[[target]]

    # Get DEA genes and LR    
    dea_genes_l <- lapply(dea_res_l, function(x) {x %>% dplyr::filter(p_val <= dea_p_value_filter_thr) %>% rownames()}) %>% do.call(c, .) %>% unique()
    dea_genes_r <- lapply(dea_res_r, function(x) {x %>% dplyr::filter(p_val <= dea_p_value_filter_thr) %>% rownames()}) %>% do.call(c, .) %>% unique()
    dea_genes <- c(dea_genes_l, dea_genes_r) %>% unique()
    dea_lr <- lr_res %>% dplyr::filter(ligand_symbol %in% dea_genes | receptor_symbol %in% dea_genes) %>% dplyr::pull(interaction_name)

    # Transform LR results 
    lr_res <- lr_res %>% dplyr::filter(interaction_name %in% dea_lr) %>% dplyr::select(pathway_name, ligand_symbol, receptor_symbol) %>% dplyr::distinct() %>% dplyr::mutate(interaction_id=paste0(pathway_name, "_", ligand_symbol, "_", receptor_symbol))

    # Prepare ligand log2FC and p-value matrix 
    mat_l_1 <- lapply(names(dea_res_l), function(i) {
    
        dea_res_l[[i]] %>% tibble::rownames_to_column("ligand_symbol") %>% dplyr::select(ligand_symbol, avg_log2FC) %>% dplyr::rename(!!i:=avg_log2FC)
        
    }
                   )
    
    mat_l_2 <- lapply(names(dea_res_l), function(i) {
    
        dea_res_l[[i]] %>% tibble::rownames_to_column("ligand_symbol") %>% dplyr::select(ligand_symbol, p_val) %>% dplyr::rename(!!i:=p_val)
        
    }
                   )

    # Prepare receptor log2FC and p-value matrix 
    mat_r_1 <- lapply(names(dea_res_r), function(i) {
    
        dea_res_r[[i]] %>% tibble::rownames_to_column("receptor_symbol") %>% dplyr::select(receptor_symbol, avg_log2FC) %>% dplyr::rename(!!i:=avg_log2FC)
        
    }
                   )
    
    mat_r_2 <- lapply(names(dea_res_r), function(i) {
    
        dea_res_r[[i]] %>% tibble::rownames_to_column("receptor_symbol") %>% dplyr::select(receptor_symbol, p_val) %>% dplyr::rename(!!i:=p_val)
        
    }
                   )

    # Join results 
    mat_l_1 <- purrr::reduce(mat_l_1, ~ full_join(.x, .y, by = "ligand_symbol"))
    mat_l_2 <- purrr::reduce(mat_l_2, ~ full_join(.x, .y, by = "ligand_symbol"))
    mat_r_1 <- purrr::reduce(mat_r_1, ~ full_join(.x, .y, by = "receptor_symbol"))
    mat_r_2 <- purrr::reduce(mat_r_2, ~ full_join(.x, .y, by = "receptor_symbol"))

    # Match results with LR matrix
    mat_l_1 <- left_join(lr_res, mat_l_1, by=join_by(ligand_symbol)) %>% dplyr::select(-pathway_name, -ligand_symbol, -receptor_symbol) %>% tibble::column_to_rownames("interaction_id")
    mat_l_2 <- left_join(lr_res, mat_l_2, by=join_by(ligand_symbol)) %>% dplyr::select(-pathway_name, -ligand_symbol, -receptor_symbol) %>% tibble::column_to_rownames("interaction_id")
    mat_r_1 <- left_join(lr_res, mat_r_1, by=join_by(receptor_symbol)) %>% dplyr::select(-pathway_name, -ligand_symbol, -receptor_symbol) %>% tibble::column_to_rownames("interaction_id")
    mat_r_2 <- left_join(lr_res, mat_r_2, by=join_by(receptor_symbol)) %>% dplyr::select(-pathway_name, -ligand_symbol, -receptor_symbol) %>% tibble::column_to_rownames("interaction_id")

    # Check if all data are in order
    lr_res <- lr_res %>% tibble::column_to_rownames("interaction_id")

    if (!all(rownames(lr_res)==rownames(mat_l_1)) & all(rownames(lr_res)==rownames(mat_l_2))) {message("Something wrong with the interaction matching")}
    if (!all(rownames(lr_res)==rownames(mat_r_1)) & all(rownames(lr_res)==rownames(mat_r_2))) {message("Something wrong with the interaction matching")}

    # Impute missing values 
    mat_l_1[is.na(mat_l_1)] <- 0
    mat_l_2[is.na(mat_l_2)] <- 1
    mat_r_1[is.na(mat_r_1)] <- 0
    mat_r_2[is.na(mat_r_2)] <- 1
    
    # LR DE Check and filter 
    if(lr_de_check) {

        lr_de_check <- ((rowSums(mat_l_2<=dea_p_value_thr)>0) | (rowSums(mat_r_2<=dea_p_value_thr)>0))
        
    } else {

        lr_de_check <- rep(TRUE, nrow(lr_res))
        
    }
    

    mat_l_1 <- mat_l_1[lr_de_check, , drop=FALSE]
    mat_l_2 <- mat_l_2[lr_de_check, , drop=FALSE]
    mat_r_1 <- mat_r_1[lr_de_check, , drop=FALSE]
    mat_r_2 <- mat_r_2[lr_de_check, , drop=FALSE]

    lr_res <- lr_res[lr_de_check, , drop=FALSE]

    # Get Heatmap components 
    row_split <- lr_res$pathway_name
    row_labels_l <- lr_res$ligand_symbol
    row_labels_r <- lr_res$receptor_symbol

    # Make Heatmaps 
    dea_l_hm <- dea_hm(mat_l_1, mat_l_2, row_labels_l, row_split, source, "Ligand", dea_p_value_thr=dea_p_value_thr, width=1.0, height=1.2, use_raster=FALSE, fontsize_select=fontsize_select)
    dea_r_hm <- dea_hm(mat_r_1, mat_r_2, row_labels_r, row_split, target, "Receptor", dea_p_value_thr=dea_p_value_thr, width=1.0, height=1.2, use_raster=FALSE, fontsize_select=fontsize_select)

    return(draw(dea_l_hm + dea_r_hm))
    
}