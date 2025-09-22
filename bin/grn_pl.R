auc_hm <- function(auc_res, celltype_level, celltype_color, celltype_filter=NULL, rank_thr=NULL, p_val_adj_thr=NULL, module=NULL, column_title_rot=0, fontsize_select=1) {

    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    if(!is.null(celltype_filter)) {auc_res <- auc_res[auc_res$celltype_low %in% celltype_filter, ]}
    if(!is.null(rank_thr)) {auc_res <- auc_res[auc_res$module %in% auc_res[auc_res$rank_1 <=rank_thr | auc_res$rank_2 >=rank_thr, ]$module, ]}
    if(!is.null(p_val_adj_thr)) {auc_res <- auc_res[auc_res$module %in% auc_res[auc_res$p_val_adj<=p_val_adj_thr,  ]$module, ]}
    if(!is.null(module)) {auc_res <- auc_res[auc_res$module %in% module, ]}
    
    # Re-set cell types in case some are missing from DEA
    celltype_level <- celltype_level[celltype_level %in% auc_res$celltype_low]

    # auc_res$avg_log2FC <- -sign(auc_res$avg_log2FC)*log10(auc_res$p_val_adj)
    
    # Get log2FC matrix 
    mat_1 <- auc_res %>% 
        dplyr::select(avg_log2FC, celltype_low, module) %>% 
        reshape2::dcast(., module ~ celltype_low, value.var="avg_log2FC") %>% 
        replace(is.na(.), 0)
    
    # Get padj matrix 
    mat_2 <- auc_res %>% 
        dplyr::select(p_val_adj, celltype_low, module) %>% 
        reshape2::dcast(., module ~ celltype_low, value.var="p_val_adj") %>% 
        replace(is.na(.), 1)
    
    # # Combine group_select with mat and preserve module names
    # mat_1 <- left_join(auc_res, mat_1, by="module")
    # mat_2 <- left_join(auc_res, mat_2, by="module")
    
    # # Store column split and column name 
    # column_split_1 <- mat_1$main_group
    # column_split_2 <- mat_2$main_group
    
    column_names_1 <- mat_1$module
    column_names_2 <- mat_2$module
    
    # Remove cell type and module names from matirx 
    mat_1 <- mat_1 %>% dplyr::select(-module) %>% as.matrix() %>% replace(is.na(.), 0)
    mat_2 <- mat_2 %>% dplyr::select(-module) %>% as.matrix() %>% replace(is.na(.), 1)
    
    # Transform matrix and set column names
    mat_1 <- t(mat_1)
    mat_2 <- t(mat_2)

    
    colnames(mat_1) <- column_names_1
    colnames(mat_2) <- column_names_2
    
    # Color function - expression z-score  
    color_ramp_mat <- c(color$sample_group["Bl6_NaCl_D6"], "#ffffff", color$sample_group["Bl6_CpG_D6"])
    mat_quant <- quantile(abs(mat_1), 0.99)
    mat_quant <- 0.5

    breaks_mat <- seq(-mat_quant, mat_quant, length.out=length(color_ramp_mat))
    color_function_mat <- circlize::colorRamp2(breaks_mat, color_ramp_mat) 
    
    # Top annotation - cell type
    row_annotation <- rowAnnotation(

        df=data.frame(

            celltype_annotation=celltype_level


        ), 

        annotation_legend_param=list(

            celltype_annotation=list(title="Celltype", title_gp=gpar(fontsize=fontsize[2], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), grid_width=unit(fontsize_scale*5.5, "mm"), grid_height=unit(fontsize_scale*5.5, "mm"))

        ), 

        col=list(

            celltype_annotation=celltype_color

        ), 

        gp=gpar(col="black", lwd=unit(0.6667, "pt")), 
        simple_anno_size=unit(fontsize_scale*6, "mm"), 
        border=FALSE,
        show_annotation_name=FALSE, 
        show_legend=FALSE

    )
    
    mat_1 <- mat_1[celltype_level, ]
    mat_2 <- mat_2[celltype_level, ]
    
    # Heat map 
    ComplexHeatmap::Heatmap(

        matrix=mat_1, 
        # name="log2(FC)"

        col=color_function_mat, 
        na_col="white", 

        width=unit(fontsize_scale*5.5*ncol(mat_1), "mm"), 
        height=unit(fontsize_scale*5.5*nrow(mat_1), "mm"), 

        # column_title=NULL, 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="italic"), 

        cluster_rows=FALSE, 
        show_row_dend=FALSE,    
        row_split=NULL,
        row_title_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        row_gap=unit(fontsize_scale*2, "mm"),
        show_row_names=TRUE,
        row_names_side="right",

        cluster_columns=TRUE,
        cluster_column_slices=FALSE, 
        column_title_rot=column_title_rot, 
        show_column_dend=FALSE, 
        # column_split=group_select$main_group,
        column_gap=unit(fontsize_scale*2, "mm"), 
        show_column_names=TRUE, 

        # top_annotation=NULL, 
        right_annotation=row_annotation, 

        rect_gp=gpar(col=NA, lwd=0, alpha=1), 

        heatmap_legend_param=list(title="log2FC", title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*5*3, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 

        border=TRUE, 
        border_gp=gpar(col="black", lwd=unit(0.6667, "pt")), 

        use_raster=FALSE, 
        
        cell_fun=function(j, i, x, y, w, h, fill) {
            
        if(mat_2[i, j] < p_val_adj_thr) {
            
            grid.text("*", x, y-unit(fontsize_scale*1, "mm"), gp=gpar(fontsize=fontsize[1]))
        
            }
        
        }

    ) 
    
    
}