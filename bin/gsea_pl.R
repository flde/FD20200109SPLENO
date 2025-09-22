library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(
        
            library(dplyr)

        )
    )
)

#####################
### GSEA dot plot ###
#####################
gsea_pl <- function(gsea_res, adj_pval_thr=0.05, title=NULL, color=c(RColorBrewer::brewer.pal(8, "Set1")[1], RColorBrewer::brewer.pal(8, "Set1")[2]), color_names=c("Pos", "Neg"), size_range=5, pathway_suffix=NULL, top=20) {
    
    # Set GSEA data frame 
    gsea_res <- as.data.frame(gsea_res)
    gsea_res <- na.omit(gsea_res) 
    
    # Set color names 
    color <- setNames(color, color_names)
    
    # Fix pathway names
    if(!is.null(pathway_suffix)) {gsea_res$pathway <- gsub(pathway_suffix, "", gsea_res$pathway)}
    gsea_res$pathway <- gsub("_", " ", gsea_res$pathway)
    
    # Filter hits 
    gsea_up <- gsea_res[sign(gsea_res$NES)==+1, ]
    gsea_down <- gsea_res[sign(gsea_res$NES)==-1, ]
    
    gsea_up <- gsea_up[order(gsea_up$padj), ][1:top, ]
    gsea_down <- gsea_down[order(gsea_down$padj), ][1:top, ]
    
    gsea_res <- rbind(gsea_up, gsea_down)
    gsea_res <- na.omit(gsea_res)
    gsea_res <- distinct(gsea_res)
    
    # Add color 
    gsea_res$color <- ifelse(sign(gsea_res$ES)==1, names(color)[1], names(color)[2])
    gsea_res$color <- ifelse(gsea_res$padj<=adj_pval_thr, gsea_res$color, NA)

    gsea_res$sign_log_pval_values <- -log10(gsea_res$padj) * sign(gsea_res$ES)
    
    # Order hits 
    gsea_res <- gsea_res[rev(order(gsea_res$sign_log_pval_values)), ]
    gsea_res$pathway <- substr(gsea_res$pathway, 1, 50)
    gsea_res$pathway <- factor(gsea_res$pathway, levels=rev(gsea_res$pathway))

    x_max <- max(abs(gsea_res$sign_log_pval_values))
    if(x_max<abs(log10(adj_pval_thr))) {x_max <- abs(log10(adj_pval_thr))}
    x_max <- ceiling(ceiling(x_max))
    
    # Plot 
    plot <- ggplot(gsea_res, aes(x=sign_log_pval_values, y=pathway, color=color)) + 
        
        geom_vline(xintercept=-log10(adj_pval_thr), linetype="dashed") + 
        geom_vline(xintercept=log10(adj_pval_thr), linetype="dashed") +
    
        geom_point(aes(size=abs(NES))) +

        ggtitle(title) +
        xlab("Signed -log10 adj. p-value") + ylab("") + 
        scale_x_continuous(breaks=c(-x_max, 0, x_max), limits=c(-x_max, x_max)) +
        scale_color_manual(values=color, na.value="black", drop=FALSE) +
        scale_size(range=c(0, size_range)) + 
        guides(
            
            color=guide_legend(order=1, title="Agent", keywidth=0.75, keyheight=0.75, override.aes=list(size=4)), 
            size=guide_legend(order=2, title="Abs. (NES)", keywidth=0.75, keyheight=0.75)
            
        ) +
    
        theme(
            
            legend.position="bottom", 
            legend.justification="top", 
            axis.text.y=element_text(size=14, hjust=1, vjust=0.5, face="plain", margin=margin(t=0, r=2, b=0, l=0), color="black")
            
        ) 
    
    return(plot)
    
}

##############################
### GSEA heat map dot plot ###
##############################
gsea_hm <- function(mat_1, mat_2, col_label, col_split, pathway_suffix=NULL, row_split=NULL, row_genes=NULL, width=1.2, height=1.2, use_raster=FALSE, fontsize_select=1){

    # Fontsize select
    fontsize <- list(c(16,18), c(6,8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Check if mat matche
    if(!is.null(row_genes)){
        
        keep <- intersect(rownames(mat_2), row_genes)
        mat_1 <- mat_1[keep, , drop=FALSE]
        mat_2 <- mat_2[keep, , drop=FALSE]
        
    }

    # Row split 
    if(!is.null(row_split)) {row_split <- row_split[rownames(mat_1), ]$row_split_label}

    # Pathway suffix 
    if(!is.null(pathway_suffix)) {

        # Fix pathway names
        if(!is.null(pathway_suffix)) {rownames(mat_1) <- gsub(pathway_suffix, "", rownames(mat_1))}
        rownames(mat_1) <- gsub("_", " ", rownames(mat_1))
    
    }
    
    # Bottom annotation 
    bottom_annotation <- HeatmapAnnotation(
        
        df=data.frame(col_label=col_label),
        col=list(col_label=color$celltype_low),
        simple_anno_size=unit(fontsize_scale*5, "mm"),
        show_annotation_name=FALSE, show_legend=FALSE, border=TRUE
    
    )

    colnames(mat_1) <- col_label
    colnames(mat_2) <- col_label

    max_abs <- max(abs(mat_2), na.rm=TRUE)
    if(!is.finite(max_abs) || max_abs==0) max_abs <- 1
    
    nes_col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c(rev(RColorBrewer::brewer.pal(11,"RdBu"))[1], "#FFFFFF", rev(RColorBrewer::brewer.pal(11,"RdBu"))[11]))


    star_levels <- c("ns","*","**","***")
    p_bins <- cut(mat_1, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), labels=c("***","**","*","ns"), include.lowest=TRUE, right=TRUE)
    p_bins <- factor(p_bins, levels=star_levels, ordered=TRUE)

    size_map <- c(
        
        "ns"=fontsize_scale*2.5*unit(width, "mm")*0, 
        "*"=fontsize_scale*2.5*unit(width, "mm")*0.50, 
        "**"=fontsize_scale*2.5*unit(width, "mm")*0.75, 
        "***"=fontsize_scale*2.5*unit(width, "mm")*1.00
    )
    
    rad_mat <- matrix(size_map[as.character(p_bins)], nrow=nrow(mat_1), dimnames=dimnames(mat_1)) * 0.80

    base_mat <- matrix(0, nrow=nrow(mat_1), ncol=ncol(mat_1), dimnames=list(rownames(mat_1), colnames(mat_1)))

    hm <- Heatmap(
        
        matrix=base_mat,
        
        col=function(z) "white",
        
        width=fontsize_scale*5*ncol(base_mat)*unit(width, "mm"),
        height=fontsize_scale*5*nrow(base_mat)*unit(height, "mm"),
        
        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"),
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"),
        row_names_gp=gpar(fontsize=fontsize[1], fontface="italic"),
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"),
        
        cluster_rows=FALSE, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,
        row_split=row_split, 
        # row_title=NULL, 
        row_title_rot=0, 
        row_gap=unit(fontsize_scale*2, "mm"),
        show_row_names=TRUE, 
        row_names_side="right",
        
        cluster_columns=FALSE, 
        clustering_distance_columns="pearson",
        cluster_column_slices=FALSE, 
        show_column_dend=TRUE,
        column_split=col_split, 
        column_gap=unit(fontsize_scale*2, "mm"),
        column_dend_height=unit(fontsize_scale*3, "mm"), 
        show_column_names=TRUE,
        
        bottom_annotation=bottom_annotation,
        
        show_heatmap_legend=FALSE,
        
        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")),
        
        use_raster=use_raster, raster_by_magick=TRUE, raster_resize_mat=mean,
        
        layer_fun=function(j, i, x, y, w, h, ...) {
            
            r  <- rad_mat[cbind(i, j)]
            nes <- mat_2[cbind(i, j)]
            ok <- !is.na(r) & !is.na(nes) & r > 0
            
            if(any(ok)){
                
                x_ok <- x[ok]
                y_ok <- y[ok]
                r_ok <- r[ok]
                nes_ok <- nes[ok]
                fill_ok <- nes_col_fun(nes_ok)
                
                for(k in seq_along(r_ok)){
                    
                  grid.circle(x_ok[k], y_ok[k], r=unit(r_ok[k], "mm"), gp=gpar(fill=fill_ok[k], col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")))
                    
                }
            }
        }
    
    )

    draw(hm)
    
    nes_lgd <- Legend(
        
        title="NES", 
        col_fun=nes_col_fun,
        title_gp=gpar(fontsize=fontsize[1], fontface="plain"),
        labels_gp=gpar(fontsize=fontsize[1]),
        grid_width=unit(fontsize_scale*3, "mm"),
        legend_height=unit(fontsize_scale*15, "mm")
    )
    
    size_lgd <- Legend(
        
        title="Adj. p-value", 
        type="points",
        at=star_levels,
        legend_gp=gpar(fill="white", col="black"),
        size=unit(unname(size_map[star_levels]), "mm"),
        title_gp=gpar(fontsize=fontsize[1], fontface="plain"),
        labels_gp=gpar(fontsize=fontsize[1])
    
    )
    
    draw(hm, heatmap_legend_list=list(nes_lgd, size_lgd))

}