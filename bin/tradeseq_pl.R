library_load <- suppressMessages(
    
    list(
        
        library(tradeSeq), 
        library(SingleCellExperiment), 

        # Data 
        library(dplyr), 
        library(tidyverse), 
        library(data.table), 

        library(ggplot2), 
        library(ggforce)
    )

)
###############################
### Plort smooth expression ###
###############################
plot_smooth <- function(fitgam, gene, n_points=50, point=FALSE, line=TRUE, condition_color=NULL, condition_line_type=FALSE, line_size=3, label_size=2, label_parse=FALSE, y_scale_full=FALSE) {

    cnt <- assays(fitgam)$counts[gene, , drop=FALSE]
    cnt <- as.data.frame(t(cnt))
    colnames(cnt) <- "exp"

    cnt$pseudotime <- colData(fitgam)$crv$pseudotime
    cnt$condition <- colData(fitgam)$tradeSeq$conditions
    
    cnt_smooth <- predictSmooth(fitgam, gene, nPoints=n_points, tidy=TRUE)

    if(y_scale_full) {

        cnt_smooth <- cnt_smooth[, c("yhat", "time", "condition")]
        colnames(cnt_smooth) <- c("exp", "pseudotime", "condition")

        max_exp <- max(log1p(cnt_smooth$exp), na.rm=TRUE)
        max_break <- floor(max_exp * 10) / 10
    
        if(max_exp<0.1) {
    
            max_exp <- 0.1
            max_break <- 0.1
            
        }

        breaks <- c(0.0, max_break)

        cnt_smooth <- cnt_smooth[cnt_smooth$condition %in% names(condition_color), ]

        
    } else {

        cnt_smooth <- cnt_smooth[cnt_smooth$condition %in% names(condition_color), ]
        cnt_smooth <- cnt_smooth[, c("yhat", "time", "condition")]
        colnames(cnt_smooth) <- c("exp", "pseudotime", "condition")
    
        max_exp <- max(log1p(cnt_smooth$exp), na.rm=TRUE)
        max_break <- floor(max_exp * 10) / 10
    
        if(max_exp<0.1) {
    
            max_exp <- 0.1
            max_break <- 0.1
            
        }
        
        breaks <- c(0.0, max_break)
        
    }
    
    if(isTRUE(condition_line_type)) {
        
        if(length(condition_color)==3) {condition_line_type <- setNames(c("solid", "solid", "dashed"), names(condition_color))}
        if(length(condition_color)==4) {condition_line_type <- setNames(c("solid", "solid", "dashed", "dashed"), names(condition_color))}
        if(length(condition_color)==6) {condition_line_type <- setNames(c("solid", "solid", "dashed", "solid", "solid", "dashed"), names(condition_color))}
    
    } else if (isFALSE(condition_line_type)) {

        condition_line_type <- setNames(rep("solid", length(condition_color)), names(condition_color))
        
    } else {

        condition_line_type <- condition_line_type

    }
    
    p <- ggplot(NULL, aes(x=pseudotime, y=log1p(exp), color=condition, linetype=condition)) +
        scale_color_manual(values=condition_color) + 
        scale_y_continuous(breaks=breaks, labels=format(breaks, nsmall=1), limits=c(0, max_exp), expand=c(0, 0.01)) + 
        scale_x_continuous(breaks=c(0.0, 0.5, 1.0), labels=c("0.0", "0.5", "1.0"), limits=c(0, 1), expand=c(0, 0)) + 
        scale_linetype_manual(values=condition_line_type) 
        # annotate("text", x=Inf, y=-Inf, label=gene, hjust=1.1, vjust=-0.5, size=label_size, parse=label_parse)
    
    if(point) {p <- p + geom_point(data=cnt, size=1)}
    if(line) {p <- p + geom_line(data=cnt_smooth, size=line_size, alpha=1)}
        
    
    return(p)
    
}

##############################
### Plort smooth prototype ###
##############################
plot_prototype <- function(fitgam, genes, condition_vec, n_points=50, line=TRUE, condition_color=NULL, condition_line_type=FALSE, line_size=3,  y_scale_full=FALSE) {

    cnt <- assays(fitgam)$counts[genes, , drop=FALSE]
    cnt <- as.data.frame(t(cnt))
    colnames(cnt) <- "exp"

    cnt$pseudotime <- colData(fitgam)$crv$pseudotime
    cnt$condition <- colData(fitgam)$tradeSeq$conditions
    
    cnt_smooth <- predictSmooth(fitgam, genes, nPoints=n_points, tidy=TRUE)
    
    # Get expression mat
    mat <- cnt_smooth %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")

    # Scale data 
    mat <- t(apply(mat, 1, scales::rescale))

    # Select plotting condition 
    mat <- mat[, grepl(paste0(paste0(condition_vec, "$"), collapse = "|"), colnames(mat))]

    # Compute prototype stats
    mat <- t(mat) %>% as.data.frame() %>% tibble::rownames_to_column("sample_group") %>% 
    
        tidyr::pivot_longer(
            
            cols=-sample_group,
            names_to="gene",
            values_to="exp"
            
        ) %>% 
    
        dplyr::mutate(pseudotime=str_split(sample_group, ":", simplify = TRUE)[, 1], condition=str_split(sample_group, ":", simplify = TRUE)[, 2]) %>% 
    
        group_by(condition, pseudotime) %>%
    
        dplyr::summarise(
            
            mean_exp=mean(exp, na.rm=TRUE),
            sd_exp=sd(exp, na.rm=TRUE),

            median_exp=median(exp, na.rm=TRUE),
            q25_exp=quantile(exp, 0.25, na.rm=TRUE),
            q75_exp=quantile(exp, 0.75, na.rm=TRUE), 
            .groups="drop"
            
        ) %>% 
    
        dplyr::ungroup() %>% 

        dplyr::mutate(

            ymin=mean_exp-sd_exp, 
            ymax=mean_exp+sd_exp
            
        ) %>% 

        dplyr::mutate(pseudotime=as.numeric(pseudotime))

    # Plot
    p <- ggplot(mat, aes(x=pseudotime, y=mean_exp, color=condition, fill=condition)) +
        geom_line(size=1.5/2.141959, alpha=1) + 
        geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.2, colour=NA) + 
        scale_y_continuous(breaks=c(0.0, 1.0), labels=c(0.0, 1.0), limits=c(0, 1), expand=c(0, 0)) + 
        scale_x_continuous(breaks=c(0.0, 0.5, 1.0), labels=c("0.0", "0.5", "1.0"), limits=c(0, 1), expand=c(0, 0)) + 
        scale_color_manual(values=color_mat) + 
        scale_fill_manual(values=color_mat) + 
        annotate("text", x=0.5, y=0, label=paste0("N=", length(genes)), hjust=0.5, vjust=-1, size=2) + 
        theme(legend.position="none")
            
    return(p)
    
}


####################################
### PT smooth expression heatmap ###
####################################
pt_hm <- function(fitgam, genes, condition_mat, condition_scale, n_points, color_mat, row_split=NULL, cluster_rows=FALSE, show_row_names=FALSE, cluster_columns=FALSE, show_column_dend=FALSE, bottom_annotation=TRUE, width=0.75, height=75, use_raster=TRUE, fontsize_select=1) {
    
    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Get mat
    if(is.null(condition_scale)) {condition_scale <- condition_mat}

    mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)
    mat <- mat %>% dplyr::filter(condition %in% condition_scale)
    mat <- mat %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")

    # Scale data 
    mat <- t(apply(mat, 1, scales::rescale))
    
    # Set color
    breaks_mat_exp <- seq(0, 1, length.out=3)
    color_mat_exp <- viridis::mako(length(breaks_mat_exp))
    color_function_mat_exp <- circlize::colorRamp2(breaks_mat_exp, color_mat_exp)

    # Select plotting condition 
    mat <- mat[, grepl(paste0(condition_mat, "$"), colnames(mat))] # Be careful that only works for single sample names

    # Order mat by row split 
    if(!is.null(row_split)) {mat <- mat[names(row_split), ]}

    # Column annotation 
    top_annotation <- HeatmapAnnotation(

        df=data.frame(condition=rep(names(color_mat), n_points)), 
        col=list(condition=color_mat), 
        simple_anno_size=unit(fontsize_scale*5, "mm"), 
        show_annotation_name=FALSE, 
        show_legend=FALSE

    )

    # Bottom annotation with DPT
    if(bottom_annotation) {

        # Get DPT mat
        dpt_mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)[1:n_points, ]$time

        # Set color
        breaks_dpi <- seq(0, 1, length.out=n_points)
    
        color_function_dpt <- viridis::rocket(length(breaks_dpi))
        names(color_function_dpt) <- breaks_dpi

        # DPT bottom annotation 
        bottom_annotation <- HeatmapAnnotation(
            
            DPT=dpt_mat, 
            col=list(DPT=color_function_dpt), 
            simple_anno_size=unit(fontsize_scale*5, "mm"), 
            show_annotation_name=FALSE, 
            show_legend=FALSE, 
            border=TRUE
            
        )
    
    } else {

        bottom_annotation=NULL
        
    }
    
    # DPT heatmap 
    hm <- Heatmap(
    
        matrix=mat, 
    
        col=color_function_mat_exp, 
        # na_col="white", 
    
        width=fontsize_scale*n_points*unit(width, "mm"),
        height=unit(fontsize_scale*height, "mm"),
    
        # row_title=NULL, 
        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 
    
        column_title="DPT exp.", 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 
    
        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
    
        cluster_rows=cluster_rows, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=row_split, 
        row_gap=unit(0.5, "mm"),
        row_names_side="left",
        show_row_names=show_row_names,
    
        cluster_columns=cluster_columns,
        cluster_column_slices=FALSE, 
        show_column_dend=show_column_dend, 
        column_split=NULL,
        column_gap=unit(0.5, "mm"), 
        show_column_names=FALSE, 

        top_annotation=top_annotation, 
    
        bottom_annotation=bottom_annotation, 
        right_annotation=NULL, 
    
        rect_gp=gpar(col=NA), 
    
        heatmap_legend_param=list(title="Scale Exp.", title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*5*3, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 
    
        border=TRUE, 
        border_gp=gpar(col="black", lwd=unit(0.6667, "pt")), 
    
        use_raster=use_raster, raster_by_magick=TRUE
    
    )

    return(hm)
    
}

####################################
### PT smooth expression heatmap ###
####################################
pt_tf_hm <- function(tf_tg_mat, row_split=NULL, row_split_order=NULL, row_title=NULL, cluster_rows=FALSE, cluster_columns=FALSE, show_column_dend=FALSE, show_column_names=TRUE, bottom_annotation=TRUE, width=0.75, height=75, use_raster=TRUE, fontsize_select=1) {
    
    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Set row split
    if(!is.null(row_split)) {tf_tg_mat <- tf_tg_mat[names(row_split), ]}

    # TF heatmap 
    hm <- Heatmap(

        matrix=tf_tg_mat, 

        col=structure(c("#FFFFFF", "#000038"), names=c(0, 1)), 
        
        width=fontsize_scale*5*ncol(tf_tg_mat)*unit(width, "mm"),
        height=unit(fontsize_scale*height, "mm"), 

        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        column_title="TF", 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 

        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        
        cluster_rows=cluster_rows, 
        # clustering_distance_rows="pearson", 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=row_split, 
        row_title=row_title , 
        row_gap=unit(0.5, "mm"),
        show_row_names=FALSE,
        row_names_side="right",

        cluster_columns=cluster_columns,
        # clustering_distance_columns="pearson", 
        cluster_column_slices=FALSE, 
        show_column_dend=TRUE, 
        column_split=NULL,
        column_gap=unit(0.5, "mm"), 
        column_dend_height = unit(fontsize_scale*3, "mm"), 
        show_column_names=show_column_names, 

        heatmap_legend_param=list(title="TF", title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*5*3, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 

        rect_gp=gpar(col=NA, lwd=1), 

        border=TRUE, 
        border_gp=gpar(col="black", lwd=unit(0.6667, "pt")),

        use_raster=use_raster, raster_by_magick=TRUE
    
    )

    return(hm)
    
}

##########################
### PT diff FC heatmap ###
##########################
pt_diff_hm <- function(fitgam, genes, condition_qry, condition_ref, n_points, color_mat, row_split=NULL, cluster_rows=FALSE, show_row_names=FALSE, cluster_columns=FALSE, show_column_dend=FALSE, bottom_annotation=TRUE, width=0.75, height=75, use_raster=TRUE, fontsize_select=1) {
    
    # Set font size 
    fontsize <- list(size_1=c(16, 18), size_2=c(6, 8))[[fontsize_select]]
    fontsize_scale <- c(1, 0.5)[[fontsize_select]]

    # Smoothed expression 
    mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)
    
    mat_qry <- mat %>% dplyr::filter(condition==condition_qry) %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")
    mat_ref <- mat %>% dplyr::filter(condition==condition_ref) %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")

    mat <- log2((mat_qry+0.1)/(mat_ref+0.1))

    # Set color
    breaks_mat_fc <- seq(-2, 2, length.out=3)
    color_mat_fc <- c(color_mat[[condition_ref]], "white", color_mat[[condition_qry]])
    color_function_mat_fc <- circlize::colorRamp2(breaks_mat_fc, color_mat_fc)

    # Order mat by row split 
    if(!is.null(row_split)) {mat <- mat[names(row_split), ]}

    # Column annotation 
    top_annotation <- HeatmapAnnotation(

        df=data.frame(condition=rep(condition_qry, n_points)), 
        col=list(condition=color_mat), 
        simple_anno_size=unit(fontsize_scale*5, "mm"), 
        show_annotation_name=FALSE, 
        show_legend=FALSE

    )

    # Bottom annotation with DPT
    if(bottom_annotation) {

        # Get DPT mat
        dpt_mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)[1:n_points, ]$time

        # Set color
        breaks_dpi <- seq(0, 1, length.out=n_points)
    
        color_function_dpt <- viridis::rocket(length(breaks_dpi))
        names(color_function_dpt) <- breaks_dpi

        # DPT bottom annotation 
        bottom_annotation <- HeatmapAnnotation(
            
            DPT=dpt_mat, 
            col=list(DPT=color_function_dpt), 
            simple_anno_size=unit(fontsize_scale*5, "mm"), 
            show_annotation_name=FALSE, 
            show_legend=FALSE, 
            border=TRUE
            
        )
    
    } else {

        bottom_annotation=NULL
        
    }
    
    # DPT heatmap 
    hm <- Heatmap(
    
        matrix=mat, 
    
        col=color_function_mat_fc, 
        # na_col="white", 
    
        width=fontsize_scale*n_points*unit(width, "mm"),
        height=unit(fontsize_scale*height, "mm"),
    
        # row_title=NULL, 
        row_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 
    
        column_title="DPT exp.", 
        column_title_gp=gpar(fontsize=fontsize[1], fontface="bold"), 
    
        row_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
        column_names_gp=gpar(fontsize=fontsize[1], fontface="plain"), 
    
        cluster_rows=cluster_rows, 
        cluster_row_slices=FALSE, 
        show_row_dend=FALSE,   
        row_split=row_split, 
        row_gap=unit(0.5, "mm"),
        row_names_side="left",
        show_row_names=show_row_names,
    
        cluster_columns=cluster_columns,
        cluster_column_slices=FALSE, 
        show_column_dend=show_column_dend, 
        column_split=NULL,
        column_gap=unit(0.5, "mm"), 
        show_column_names=FALSE, 

        top_annotation=top_annotation, 
    
        bottom_annotation=bottom_annotation, 
        right_annotation=NULL, 
    
        rect_gp=gpar(col=NA), 
    
        heatmap_legend_param=list(title="Scale Exp.", title_gp=gpar(fontsize=fontsize[1], fontface="plain"), labels_gp=gpar(fontsize=fontsize[1]), legend_height=unit(fontsize_scale*5*3, "mm"), grid_width=unit(fontsize_scale*3, "mm")), 
    
        border=TRUE, 
        border_gp=gpar(col="black", lwd=unit(0.6667, "pt")), 
    
        use_raster=use_raster, raster_by_magick=TRUE
    
    )

    return(hm)
    
}
































