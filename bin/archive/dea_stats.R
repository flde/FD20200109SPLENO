#########################
### Fetch DEA results ###
#########################
dea_res_fetch <- function(dea_res, row_order, p_adj_thr=0.05, column_title=NULL, name=NULL, grouping_var_query=NULL, grouping_var_names=NULL) {

    if(is.null(grouping_var_query)) {
        
        grouping_var_query <- c(
            
            "Tx-Baseline", 
            "D14-Baseline", 
            "D100-Baseline", 
            "GVHD-Baseline"
        
        )
    }

    if(is.null(grouping_var_names)) {
        
        grouping_var_names <- c(
            
            "Tx", 
            "D14", 
            "D100", 
            "GVHD"
        
        )
        
    }
    
    results_n <- lapply(dea_res, function(result) {
    
        # Count significant DEA results 
        dea_res_n <- lapply(result, function(x) {nrow(x[x$p_val_adj <= p_adj_thr, ])})
        
        grouping_var_failed <- grouping_var_query[!grouping_var_query %in% names(dea_res_n)]
    
        if(length(grouping_var_failed)!=0) {
    
            for(n in grouping_var_failed) {

                dea_res_n[[n]] <- NA

            }

        }

        # Order results 
        dea_res_n <- dea_res_n[grouping_var_query]
        
        # Rename results 
        names(dea_res_n) <- grouping_var_names
        
        dea_res_n <- do.call(rbind, dea_res_n)
        
        return(dea_res_n)

    }
                       )

    mat <- do.call(cbind, results_n)
    colnames(mat) <- names(dea_res)
    
    # Normalize and transform count matrix
    mat <- t(mat)
    
    # Order rows 
    return(mat)
    
}

#################################
### Plot heatmap DEA results ####
#################################
dea_res_hm <- function(mat, fontsize_select=1) {

    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]
    
    color_ramp_mat <- c("white", rev(RColorBrewer::brewer.pal(11,"RdBu"))[11])
    breaks_mat <- seq(0,  max(mat, na.rm=TRUE), length.out=length(color_ramp_mat))
    color_function_mat <- circlize::colorRamp2(breaks_mat, color_ramp_mat) 
    
    row_labels <- parse(text = rownames(mat))

    breaks <- round(max(breaks_mat))
    if (breaks %% 2 != 0) breaks <- breaks + 1
    
    hm <- Heatmap(
        
        mat,
    
        col=color_function_mat, 
        na_col="#d3d3d3", 
        
        width=fontsize_scale*5*ncol(mat)*unit(1, "mm"),
        height=fontsize_scale*5*nrow(mat)*unit(1, "mm"), 
    
        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=FALSE,
        show_row_names=TRUE,
        row_names_side="left",
        row_labels=row_labels, 
        
        cluster_columns=FALSE,
        show_column_names=TRUE, 
    
        show_heatmap_legend=TRUE, 

        heatmap_legend_param=list(title="log10(DEG)", at=c(0, breaks/2, breaks), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 
        
        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt"))

    )

    return(hm)
}