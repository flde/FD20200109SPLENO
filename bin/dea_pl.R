####################
### Volcano plot ###
####################
v_pl <- function(dea_res, log2_thr=1, adj_pval_thr=0.05, point_size=4, top_label=10, label=NULL, label_merge=TRUE, label_size=5, y_limit=NULL, title=NULL, color_pos=c("pos"="#0000ffff"), color_neg=c("neg"="#fd8008ff"), aspect_ratio=1) {
    
    # Set rownames to genes
    if("gene" %in% colnames(dea_res)) {rownames(dea_res) <- dea_res$gene}
    
    # Annotate entries significance by log2_thr and adj_pval_thr
    dea_res$p_val_adj <- ifelse(dea_res$p_val_adj == 0, min(dea_res$p_val_adj), dea_res$p_val_adj)
    dea_res$sig <- ifelse((abs(dea_res$avg_log2FC) >= log2_thr) & (dea_res$p_val_adj <= adj_pval_thr), "s", "ns")
    
    # Set color based on significance and direction of dea e.g. positive and negative 
    dea_res$color <- ifelse(dea_res$sig == "s" & dea_res$avg_log2FC > +log2_thr, names(color_pos), "n.s.")
    dea_res$color <- ifelse(dea_res$sig == "s" & dea_res$avg_log2FC < -log2_thr, names(color_neg), dea_res$color)

    dea_res$color <- factor(dea_res$color, levels=c(names(color_pos), names(color_neg), "n.s."))
    
    color <- c(color_pos, "#7f7f7f", color_neg)
    names(color) <- c(names(color_pos), "n.s.", names(color_neg))
    
    # Create labels based log2FC and p_val_adj
    dea_pos <- dea_res[dea_res$avg_log2FC > 0 & dea_res$sig == "s", ]
    dea_neg <- dea_res[dea_res$avg_log2FC < 0 & dea_res$sig == "s", ]

    pos_labels <- dea_pos[rev(order(dea_pos$avg_log2FC)), ] %>% rownames()
    neg_labels <- dea_neg[order(dea_neg$avg_log2FC), ] %>% rownames()
    
    pos_top_labels <- dea_pos[rev(order(dea_pos$avg_log2FC)), ][1:top_label, ] %>% rownames()
    neg_top_labels <- dea_neg[order(dea_neg$avg_log2FC), ][1:top_label, ] %>% rownames()
    
    # Set labels 
    if(is.null(label)) {

        label_select <- c(pos_top_labels, neg_top_labels)
        
        
    } else if (!is.null(label) & !label_merge) {

        label_select <- label[label %in% c(pos_labels, neg_labels)]
        
    } else if (!is.null(label) & label_merge) {

        label_select <- c(label[label %in% c(pos_labels, neg_labels)], c(pos_top_labels, neg_top_labels)) %>% unique()
        
    }

    dea_res$label <- ifelse(rownames(dea_res) %in% label_select, rownames(dea_res), NA)

    # Average expression 
    dea_res$avg_expression_ratio <- rowMeans(dea_res[, c("pct.1", "pct.2")])

    # Set colors
    color_values <- c(color_pos, color_neg, "n.s."="#7f7f7f")
    # dea$color <- factor(dea$color, levels=names(color_values))

    # Limit log2FC 
    if(!is.null(y_limit)) {dea_res <- dplyr::mutate(dea_res, avg_log2FC=ifelse(abs(avg_log2FC)>y_limit, sign(avg_log2FC)*y_limit, avg_log2FC))}
    if(is.null(y_limit)) {y_limit <- max(abs(dea_res$avg_log2FC))}

    # Plot
    volcano_plot <- ggplot(dea_res, aes(x=avg_expression_ratio, y=avg_log2FC, color=color, label=label), alpha=1) + 
    
        geom_point(data=dea_res[dea_res$color=="n.s.", ], size=point_size, shape=19) + 
        geom_point(data=dea_res[dea_res$color!="n.s.", ], size=point_size, shape=19) +
        geom_hline(aes(yintercept=+log2_thr), linetype="dotted", colour="black") +
        geom_hline(aes(yintercept=-log2_thr), linetype="dotted", colour="black") +
        ggrepel::geom_text_repel(segment.color="black", force=10, force_pull=1, max.overlaps=getOption("ggrepel.max.overlaps", default=100), size=label_size, alpha=1, guide="none", segment.size=0.1, color='black', min.segment.length=0.1, fontface="italic") + 
        ylim(-y_limit-1.0, y_limit+1.0) +  
        xlim(0, 1) + 
        ggtitle(title) + xlab("Expression [ratio]") + ylab("log2FC") + 
        scale_color_manual(values=c(color_pos, color_neg, "n.s."="#7f7f7f"), name="DEA") + 
    
        guides(
            
            color=guide_legend(order=1, title="Group", size=2, keywidth=0.75, keyheight=0.75), 
            alpha="none"
            
        ) + 
    
        theme(
            
            legend.position="right", 
            aspect.ratio=aspect_ratio
            
        )  +
    
        annotate("text", x=Inf, y=Inf, label=paste0("N[cells]~", dea_res$n_cells_1[1]), hjust=1.1, vjust=1.5, size=label_size, parse=TRUE) +
        annotate("text", x=Inf, y=-Inf, label=paste0("N[cells]~", dea_res$n_cells_2[1]), hjust=1.1, vjust=-0.5, size=label_size, parse=TRUE)
    
    return(volcano_plot)
    
}


#################################
### Plot heatmap DEA results ####
#################################
dea_res_hm <- function(mat, col_split=NULL, use_raster=FALSE, fontsize_select=1) {

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
    
        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        # column_title=col_label, 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=FALSE,
        show_row_names=TRUE,
        row_names_side="left",
        row_labels=row_labels, 
        
        cluster_columns=FALSE,
        column_split=col_split, 
        show_column_names=TRUE, 
    
        show_heatmap_legend=TRUE, 

        heatmap_legend_param=list(title="log10(DEG)", at=c(0, breaks/2, breaks), title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*15, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 
        
        border=TRUE, 
        rect_gp=gpar(col="black", lwd=unit(fontsize_scale*2*0.6667, "pt")), 
        border_gp=gpar(col="black", lwd=unit(fontsize_scale*3*0.6667, "pt")), 

        use_raster=use_raster, raster_by_magick=TRUE, raster_resize_mat=mean

    )

    return(hm)
}